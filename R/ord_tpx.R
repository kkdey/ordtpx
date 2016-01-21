#######  Undocumented "tpx" utility functions #########

## ** Only referenced from topics.R

## check counts (can be an object from tm, slam, or a simple co-occurance matrix)
CheckCounts <- function(fcounts){
  if(class(fcounts)[1] == "TermDocumentMatrix"){ fcounts <- t(fcounts) }
  if(is.null(dimnames(fcounts)[[1]])){ dimnames(fcounts)[[1]] <- paste("doc",1:nrow(fcounts)) }
  if(is.null(dimnames(fcounts)[[2]])){ dimnames(fcounts)[[2]] <- paste("wrd",1:ncol(fcounts)) }
  empty <- row_sums(fcounts) == 0
  if(sum(empty) != 0){
    fcounts <- fcounts[!empty,]
    cat(paste("Removed", sum(empty), "blank documents.\n")) }
  return(as.simple_triplet_matrix(fcounts))
}


## theta initialization

## ** main workhorse function.  Only Called by the above wrappers.
## topic estimation for a given number of topics (taken as ncol(theta))
tpxfit <- function(fcounts, X, param_set, del_beta, a_mu, b_mu, ztree_options, tol, verb,
                   admix, grp, tmax, wtol, qn)
{
  ## inputs and dimensions
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix") }
  mu_tree_set <- mu_tree_build_set(param_set);
  K <- length(param_set);
  levels <- length(mu_tree_set[[1]]);
  theta <- do.call(cbind, lapply(1:K, function(l) mu_tree_set[[l]][[levels]]/mu_tree_set[[l]][[1]]));
  n <- nrow(X)
  p <- ncol(X)
  m <- row_sums(X)

  ## recycle these in tpcweights to save time
  xvo <- X$v[order(X$i)]
  wrd <- X$j[order(X$i)]-1
  doc <- c(0,cumsum(as.double(table(factor(X$i, levels=c(1:nrow(X)))))))

  ## Initialize
  omega <- tpxweights(n=n, p=p, xvo=xvo, wrd=wrd, doc=doc, start=tpxOmegaStart(X,theta), theta=theta)

  ## tracking
  iter <- 0
  dif <- tol+1+qn
  update <- TRUE
  if(verb){
    cat("log posterior increase: " )
    digits <- max(1, -floor(log(tol, base=10))) }

  Y <- NULL # only used for qn > 0
  Q0 <- col_sums(X)/sum(X)
  L <- tpxlpost(fcounts, omega=omega, param_set=param_set, del_beta, a_mu, b_mu, ztree_options=1);
 # if(is.infinite(L)){ L <- sum( (log(Q0)*col_sums(X))[Q0>0] ) }

  ## Iterate towards MAP
  while( update  && iter < tmax ){

    ## sequential quadratic programming for conditional Y solution

    if(admix && wtol > 0){ Wfit <- tpxweights(n=nrow(X), p=ncol(X), xvo=xvo, wrd=wrd, doc=doc,
                                start=omega, theta=theta,  verb=0, nef=TRUE, wtol=wtol, tmax=20) }
    else{ Wfit <- omega }

    ## Construct the MRA of z-values given the current iterates of omega /theta

    z_tree <- z_tree_construct(fcounts, omega_iter = Wfit, theta_iter = t(theta), ztree_options = 1);

    ## Extract the beta and mu_0 parameters from the MRA tree

    param_set_fit <- param_extract_ztree(z_tree, del_beta, a_mu, b_mu);


    ## Build a MRA of mu-tree sets (set of clusters)

    mu_tree_set_fit <- mu_tree_build_set(param_set_fit);

    ## Extract the theta updates from the MRA tree

    levels <- length(mu_tree_set_fit[[1]]);
    theta_fit <- do.call(cbind, lapply(1:nclus, function(l) mu_tree_set_fit[[l]][[levels]]/mu_tree_set_fit[[l]][[1]]));
    move <- list(theta=theta_fit, omega=Wfit);

    ## joint parameter EM update
    ## move <- tpxEM(X=X, m=m, theta=theta, omega=Wfit, alpha=alpha, admix=admix, grp=grp)

    ## quasinewton-newton acceleration
    QNup <- tpxQN(move=move, fcounts=fcounts, Y=Y, X=X, del_beta=del_beta, a_mu=a_mu, b_mu=b_mu,
                  ztree_options=ztree_options, verb=verb, admix=admix, grp=grp, doqn=qn-dif)
    move <- QNup$move
    Y <- QNup$Y

    if(QNup$L < L){  # happens on bad Wfit, so fully reverse
      if(verb > 10){ cat("_reversing a step_") }
      ##move <- tpxEM(X=X, m=m, theta=theta, omega=omega, alpha=alpha, admix=admix, grp=grp)
      z_tree <- z_tree_construct(fcounts, omega_iter = omega, theta_iter = t(theta), ztree_options = 1);
      param_set_fit <- param_extract_ztree(z_tree, del_beta, a_mu, b_mu);
      mu_tree_set_fit <- mu_tree_build_set(param_set_fit);
      levels <- length(mu_tree_set_fit[[1]]);
      theta_fit <- do.call(cbind, lapply(1:nclus, function(l) mu_tree_set_fit[[l]][[levels]]/mu_tree_set_fit[[l]][[1]]));
      move <- list(theta=theta_fit, omega=omega);
      QNup$L <-  tpxlpost(fcounts, move$omega, param_set_fit, del_beta, a_mu, b_mu, ztree_options=1) }

    ## calculate dif
    dif <- (QNup$L-L)

    L <- QNup$L


    ## check convergence
    if(abs(dif) < tol){
      if(sum(abs(theta-move$theta)) < tol){ update = FALSE } }

    ## print
    if(verb>0 && (iter-1)%%ceiling(10/verb)==0 && iter>0){
      cat( paste( round(dif,digits), #" (", sum(abs(theta-move$theta)),")",
                 ", ", sep="") ) }

    ## heartbeat for long jobs
    if(((iter+1)%%1000)==0){
          cat(sprintf("p %d iter %d diff %g\n",
                nrow(theta), iter+1,round(dif))) }

    ## iterate
    iter <- iter+1
    theta <- move$theta;
    theta_tree_set <- lapply(1:K, function(k) mra_bottom_up(theta[,k]));
    param_set <- param_extract_mu_tree(theta_tree_set)
    omega <- move$omega

  }

  ## final log posterior
  L <- tpxlpost(fcounts, omega, param_set, del_beta, a_mu, b_mu, ztree_options=1);

  ## summary print
  if(verb>0){
    cat("done.")
    if(verb>1) { cat(paste(" (L = ", round(L,digits), ")", sep="")) }
    cat("\n")
  }

  out <- list(param_set=param_set, omega=omega, K=K, L=L, iter=iter)
  invisible(out) }


## ** called from topics.R (predict) and tpx.R
## Conditional solution for topic weights given theta
tpxweights <- function(n, p, xvo, wrd, doc, start, theta, verb=FALSE, nef=TRUE, wtol=10^{-5}, tmax=1000)
{
  K <- ncol(theta)
  start[start == 0] <- 0.1/K
  start <- start/rowSums(start)
  omega <- .C("Romega",
              n = as.integer(n),
              p = as.integer(p),
              K = as.integer(K),
              doc = as.integer(doc),
              wrd = as.integer(wrd),
              X = as.double(xvo),
              theta = as.double(theta),
              W = as.double(t(start)),
              nef = as.integer(nef),
              tol = as.double(wtol),
              tmax = as.integer(tmax),
              verb = as.integer(verb),
              PACKAGE="ordtpx")
  return(t(matrix(omega$W, nrow=ncol(theta), ncol=n))) }

## ** Called only in tpx.R


## Quasi Newton update for q>0
tpxQN <- function(move, fcounts, Y, X, del_beta, a_mu, b_mu, ztree_options, verb, admix, grp, doqn)
{
  ## always check likelihood
  theta_tree_set_in <- lapply(1:K, function(k) mra_bottom_up(move$theta[,k]));
  param_set_in <- param_extract_mu_tree(theta_tree_set_in)
  L <- tpxlpost(fcounts, move$omega, param_set_in, del_beta, a_mu, b_mu, ztree_options)

  if(doqn < 0){ return(list(move=move, L=L, Y=Y)) }

  ## update Y accounting
  Y <- cbind(Y, tpxToNEF(theta=move$theta, omega=move$omega))
  if(ncol(Y) < 3){ return(list(Y=Y, move=move, L=L)) }
  if(ncol(Y) > 3){ warning("mis-specification in quasi-newton update; please report this bug.") }

  ## Check quasinewton secant conditions and solve F(x) - x = 0.
  U <- as.matrix(Y[,2]-Y[,1])
  V <- as.matrix(Y[,3]-Y[,2])
  sUU <- sum(U^2)
  sVU <- sum(V*U)
  Ynew <- Y[,3] + V*(sVU/(sUU-sVU))
  qnup <- tpxFromNEF(Ynew, n=nrow(move$omega),
                     p=nrow(move$theta), K=ncol(move$theta))

  ## check for a likelihood improvement
  theta_tree_set_nup <- lapply(1:K, function(k) mra_bottom_up(qnup$theta[,k]));
  param_set_nup <- param_extract_mu_tree(theta_tree_set_nup)
  Lqnup <- try(tpxlpost(X=X, qnup$omega, param_set_nup, del_beta, a_mu, b_mu, ztree_options), silent=TRUE)


  if(inherits(Lqnup, "try-error")){
    if(verb>10){ cat("(QN: try error) ") }
    return(list(Y=Y[,-1], move=move, L=L)) }

  if(verb>10){ cat(paste("(QN diff ", round(Lqnup-L,3), ")\n", sep="")) }

  if(Lqnup < L){
    return(list(Y=Y[,-1], move=move, L=L)) }
  else{
    L <- Lqnup
    Y <- cbind(Y[,2],Ynew)
    return( list(Y=Y, move=qnup, L=L) )
  }
}


tpxOmegaStart <- function(X, theta)
  {
    if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }
    omega <- try(tcrossprod_simple_triplet_matrix(X, solve(t(theta)%*%theta)%*%t(theta)), silent=TRUE )
    if(inherits(omega,"try-error")){ return( matrix( 1/ncol(theta), nrow=nrow(X), ncol=ncol(theta) ) ) }
    omega[omega <= 0] <- .5
    return( normalize(omega, byrow=TRUE) )
  }


## fast computation of sparse P(X) for X>0
tpxQ <- function(theta, omega, doc, wrd){

  if(length(wrd)!=length(doc)){stop("index mis-match in tpxQ") }
  if(ncol(omega)!=ncol(theta)){stop("theta/omega mis-match in tpxQ") }

  out <- .C("RcalcQ",
            n = as.integer(nrow(omega)),
            p = as.integer(nrow(theta)),
            K = as.integer(ncol(theta)),
            doc = as.integer(doc-1),
            wrd = as.integer(wrd-1),
            N = as.integer(length(wrd)),
            omega = as.double(omega),
            theta = as.double(theta),
            q = double(length(wrd)),
            PACKAGE="ordtpx" )

  return( out$q ) }

## model and component likelihoods for mixture model
tpxMixQ <- function(X, omega, theta, grp=NULL, qhat=FALSE){
  if(is.null(grp)){ grp <- rep(1, nrow(X)) }
  K <- ncol(omega)
  n <- nrow(X)
  mixhat <- .C("RmixQ",
               n = as.integer(nrow(X)),
               p = as.integer(ncol(X)),
               K = as.integer(K),
               N = as.integer(length(X$v)),
               B = as.integer(nrow(omega)),
               cnt = as.double(X$v),
               doc = as.integer(X$i-1),
               wrd = as.integer(X$j-1),
               grp = as.integer(as.numeric(grp)-1),
               omega = as.double(omega),
               theta = as.double(theta),
               Q = double(K*n),
               PACKAGE="ordtpx")
  ## model and component likelihoods
  lQ <- matrix(mixhat$Q, ncol=K)
  lqlhd <- log(row_sums(exp(lQ)))
  lqlhd[is.infinite(lqlhd)] <- -600 # remove infs
  if(qhat){
    qhat <- exp(lQ-lqlhd)
    ## deal with numerical overload
    infq <- row_sums(qhat) < .999
    if(sum(infq)>0){
      qhat[infq,] <- 0
      qhat[n*(apply(matrix(lQ[infq,],ncol=K),1,which.max)-1) + (1:n)[infq]] <- 1 }
  }
  return(list(lQ=lQ, lqlhd=lqlhd, qhat=qhat)) }


## functions to move theta/omega to and from NEF.
tpxToNEF <- function(theta, omega){
  n <- nrow(omega)
  p <- nrow(theta)
  K <- ncol(omega)
  return(.C("RtoNEF",
            n=as.integer(n), p=as.integer(p), K=as.integer(K),
            Y=double((p-1)*K + n*(K-1)),
            theta=as.double(theta), tomega=as.double(t(omega)),
            PACKAGE="ordtpx")$Y)
}

## 'From' NEF representation back to probabilities
tpxFromNEF <- function(Y, n, p, K){
  bck <- .C("RfromNEF",
            n=as.integer(n), p=as.integer(p), K=as.integer(K),
            Y=as.double(Y), theta=double(K*p), tomega=double(K*n),
            PACKAGE="ordtpx")
  return(list(omega=t( matrix(bck$tomega, nrow=K) ), theta=matrix(bck$theta, ncol=K)))
}

