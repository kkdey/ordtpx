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
ord.tpxfit <- function(fcounts, X, param_set, del_beta, a_mu, b_mu, ztree_options, tol, verb,
                   admix, grp, tmax, wtol, qn, acc, adapt.method)
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
  omega <- ord.tpxweights(n=n, p=p, xvo=xvo, wrd=wrd, doc=doc,
                          start=ord.tpxOmegaStart(X,theta), theta=theta)

 # omega <- matrix(1/K, ncol=K, nrow=n)
  ## tracking
  iter <- 0
  dif <- tol+1+qn
  update <- TRUE
  if(verb){
    cat("log posterior increase: " )
    digits <- max(1, -floor(log(tol, base=10))) }

  Y <- NULL # only used for qn > 0
  Q0 <- col_sums(X)/sum(X)
  L <- ord.tpxlpost(fcounts, omega_iter = omega, theta_iter = theta,
                    del_beta, a_mu, b_mu, ztree_options=1,
                    adapt.method=adapt.method);
 # if(is.infinite(L)){ L <- sum( (log(Q0)*col_sums(X))[Q0>0] ) }

  ## Iterate towards MAP
  #tmax <- 5000

  while( update  && iter < tmax ){

    ## sequential quadratic programming for conditional Y solution

    move <- ord.tpxEM(X=X, m=m, theta=theta, omega=omega, method_admix = 1)


    if(admix && wtol > 0){ Wfit <- ord.tpxweights(n=nrow(X), p=ncol(X), xvo=xvo, wrd=wrd, doc=doc,
                                start=move$omega, theta=move$theta,  verb=0, nef=TRUE, wtol=wtol, tmax=20) }
    if(!admix | wtol <=0){ Wfit <- move$omega }


    if(!acc){
      #   z_tree <- z_tree_construct(fcounts, omega_iter = move$omega, theta_iter = t(move$theta), ztree_options = 1);
      #    param_set_fit <- param_extract_ztree(z_tree, del_beta, a_mu, b_mu);
      L_new <- ord.tpxlpost(fcounts, move$omega, move$theta,
                            del_beta, a_mu, b_mu, ztree_options=1,
                            adapt.method=adapt.method)
      QNup <- list("move"=move, "L"=L_new, "Y"=NULL)
      Y <- QNup$Y
    }
    ## joint parameter EM update
    ## move <- tpxEM(X=X, m=m, theta=theta, omega=Wfit, alpha=alpha, admix=admix, grp=grp)

    if(acc){
      ## quasinewton acceleration
      QNup <- ord.tpxQN(move=move, fcounts=fcounts, Y=Y, del_beta=del_beta, a_mu=a_mu, b_mu=b_mu,
                        ztree_options=ztree_options, adapt.method = adapt.method,
                        verb=verb, admix=admix, grp=grp, doqn=qn-dif)
      move <- QNup$move
      Y <- QNup$Y
    }


    if(adapt.method=="beta"){
    ## Construct the MRA of z-values given the current iterates of omega /theta

    z_tree <- z_tree_construct(fcounts, omega_iter = move$omega,
                               theta_iter = t(move$theta),
                               ztree_options = 1);

    ## Extract the beta and mu_0 parameters from the MRA tree

    param_set_fit <- param_extract_ztree(z_tree, del_beta, a_mu, b_mu);


    ## Build a MRA of mu-tree sets (set of clusters)

    mu_tree_set_fit <- mu_tree_build_set(param_set_fit);

    ## Extract the theta updates from the MRA tree

    levels <- length(mu_tree_set_fit[[1]]);
    theta_fit <- do.call(cbind, lapply(1:K,
                  function(l) mu_tree_set_fit[[l]][[levels]]/mu_tree_set_fit[[l]][[1]]));
    move <- list(theta=move$theta, omega=move$omega);
    }

    if(adapt.method=="smash"){
      row_total <- rowSums(fcounts);
      z_leaf_est <- round(sweep(move$theta, MARGIN=2, colSums(sweep(move$omega, MARGIN = 1,
                                                               row_total, "*")), "*"));
      #  plot(z_leaf_est[,1])
      #  plot(z_leaf_est[,2])
      z_leaf_smoothed <- do.call(cbind, lapply(1:dim(z_leaf_est)[2], function(k)
      {
        out <- suppressMessages(smashr::smash.poiss(z_leaf_est[,k]))
        return(out)
      }))
      #  plot(z_leaf_smoothed[,1])
      #  plot(z_leaf_smoothed[,2])
      theta_smoothed <- ordtpx::ord.normalizetpx(z_leaf_smoothed, byrow=FALSE)
      move <- list(theta=theta_smoothed, omega=move$omega)

      QNup$L <- ord.tpxlpost(fcounts, move$omega, move$theta,
                             del_beta, a_mu, b_mu, ztree_options=1,
                             adapt.method=adapt.method)
    }
    if(adapt.method=="bash"){
      row_total <- rowSums(fcounts);
      z_leaf_est <- round(sweep(move$theta, MARGIN=2, colSums(sweep(move$omega, MARGIN = 1, row_total, "*")), "*"));
      z_leaf_smoothed <- do.call(cbind, lapply(1:dim(z_leaf_est)[2], function(k)
      {
        out <- suppressMessages(binshrink(z_leaf_est[,k])$est)
        return(out)
      }))
      theta_smoothed <- ordtpx::ord.normalizetpx(z_leaf_smoothed, byrow=FALSE)
      move <- list(theta=theta_smoothed, omega=omega)
      QNup$L <-  ord.tpxlpost(fcounts, move$omega, move$theta,
                              del_beta, a_mu, b_mu, ztree_options=1,
                              adapt.method=adapt.method)

    }



#   plot(move$theta[,1], type="l")
#    plot(move$theta[,2], type="l")
#    barplot(t(move$omega),
#            col = 2:(K+1),
#            axisnames = F, space = 0, border = NA,
#            main=paste("No. of clusters=", K),
#            las=1, ylim=c(0,1), cex.axis=1.5, cex.main=1.4)


#  if(adapt.method!="bash"){
    if(QNup$L < L){  # happens on bad Wfit, so fully reverse
      if(verb > 10){ cat("_reversing a step_") }
 #     cat("We enter Qnup$L < L")
#      cat("\n")
      ##move <- tpxEM(X=X, m=m, theta=theta, omega=omega, alpha=alpha, admix=admix, grp=grp)
      if(adapt.method=="beta"){
        move <- ord.tpxEM(X=X, m=m, theta=theta, omega=omega)
        z_tree <- z_tree_construct(fcounts, omega_iter = move$omega, theta_iter = t(move$theta), ztree_options = 1);
        param_set_fit <- param_extract_ztree(z_tree, del_beta, a_mu, b_mu);
        mu_tree_set_fit <- mu_tree_build_set(param_set_fit);
        levels <- length(mu_tree_set_fit[[1]]);
        theta_fit <- do.call(cbind, lapply(1:K, function(l) mu_tree_set_fit[[l]][[levels]]/mu_tree_set_fit[[l]][[1]]));
        move <- list(theta=theta_fit, omega=omega);
        QNup$L <-  ord.tpxlpost(fcounts, move$omega, move$theta,
                                del_beta, a_mu, b_mu, ztree_options=1,
                                adapt.method=adapt.method)
      }
      if(adapt.method=="smash"){
        move <- ord.tpxEM(X=X, m=m, theta=theta, omega=omega)
        z_leaf_est <- round(sweep(move$theta, MARGIN=2, colSums(sweep(move$omega, MARGIN = 1, row_total, "*")), "*"));
        z_leaf_smoothed <- do.call(cbind, lapply(1:dim(z_leaf_est)[2], function(k)
        {
          out <- suppressMessages(smashr::smash.poiss(z_leaf_est[,k]))
          return(out)
        }))
        theta_smoothed <- ordtpx::ord.normalizetpx(z_leaf_smoothed, byrow=FALSE)
        move <- list(theta=theta_smoothed, omega=omega)
        QNup$L <-  ord.tpxlpost(fcounts, move$omega, move$theta,
                                del_beta, a_mu, b_mu, ztree_options=1,
                                adapt.method=adapt.method)
      }
      if(adapt.method=="bash"){
        move <- ord.tpxEM(X=X, m=m, theta=theta, omega=omega)
        z_leaf_est <- round(sweep(move$theta, MARGIN=2, colSums(sweep(move$omega, MARGIN = 1, row_total, "*")), "*"));
        z_leaf_smoothed <- do.call(cbind, lapply(1:dim(z_leaf_est)[2], function(k)
        {
          out <- suppressMessages(binshrink(z_leaf_est[,k])$est)
          return(out)
        }))
        theta_smoothed <- ordtpx::ord.normalizetpx(z_leaf_smoothed, byrow=FALSE)
        move <- list(theta=theta_smoothed, omega=omega)
        QNup$L <-  ord.tpxlpost(fcounts, move$omega, move$theta,
                                del_beta, a_mu, b_mu, ztree_options=1,
                                adapt.method=adapt.method)

      }

    }
#}


    ## calculate dif
    dif <- (QNup$L-L)

    L <- QNup$L


    ## check convergence
    if(abs(dif) < tol){
      if(sum(abs(theta-move$theta)) < tol){ update = FALSE } }

    ## print
    if(verb>0 && iter>0){
      cat( paste( round(abs(dif),digits), #" (", sum(abs(theta-move$theta)),")",
                 ", ", sep="") ) }

    ## heartbeat for long jobs
    if(((iter+1)%%1000)==0){
          cat(sprintf("p %d iter %d diff %g\n",
                nrow(theta), iter+1,round(dif))) }

    ## iterate
    iter <- iter+1
    theta <- move$theta;
    omega <- move$omega

  }

  ## final log posterior
  L <- ord.tpxlpost(fcounts, omega, theta, del_beta, a_mu, b_mu,
                    ztree_options=1, adapt.method=adapt.method);

  ## summary print
  if(verb>0){
    cat("done.")
    if(verb>1) { cat(paste(" (L = ", round(L,digits), ")", sep="")) }
    cat("\n")
  }

  out <- list(theta=theta, omega=omega, K=K, L=L, iter=iter)
  invisible(out) }


## ** called from topics.R (predict) and tpx.R
## Conditional solution for topic weights given theta
ord.tpxweights <- function(n, p, xvo, wrd, doc, start, theta, verb=FALSE, nef=TRUE, wtol=10^{-5}, tmax=1000)
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

ord.tpxEM <- function(X, m, theta, omega, method_admix=1){
  n <- nrow(X)
  p <- ncol(X)
  K <- ncol(theta)

  if(method_admix==1){
    lambda <- omega%*%t(theta);
    counts2 <- as.matrix(X);
    temp <- counts2/lambda;
    t_matrix <- (t(temp) %*% omega)*theta;
    w_matrix <- (temp %*% theta)*omega;

    theta <- normalizetpx(t_matrix, byrow=FALSE)
    omega <- normalizetpx(w_matrix+(1/(n*K)), byrow=TRUE)
    full_indices <- which(omega==1, arr.ind=T)
    full_indices_rows <- unique(full_indices[,1]);
    omega[full_indices_rows,] <- omega[full_indices_rows,] + (1/(n*K));
    omega <- normalizetpx(omega, byrow=TRUE)

  }

  if(method_admix==2){ Xhat <- (X$v/tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j))*(omega[X$i,]*theta[X$j,])
  Zhat <- .C("Rzhat", n=as.integer(n), p=as.integer(p), K=as.integer(K), N=as.integer(nrow(Xhat)),
             Xhat=as.double(Xhat), doc=as.integer(X$i-1), wrd=as.integer(X$j-1),
             zj = as.double(rep(0,K*p)), zi = as.double(rep(0,K*n)), PACKAGE="maptpx")
  theta <- normalizetpx(matrix(Zhat$zj, ncol=K), byrow=FALSE)
  omega <- normalizetpx(matrix(Zhat$zi+1/K, ncol=K)) }

  return(list(theta=theta, omega=omega))
  }


## Quasi Newton update for q>0
ord.tpxQN <- function(move, fcounts, Y, del_beta, a_mu, b_mu,
                      ztree_options, adapt.method,
                      verb, admix, grp, doqn)
{
  ## always check likelihood
  K <- ncol(move$theta);
  L <- ord.tpxlpost(fcounts, move$omega, move$theta, del_beta,
                    a_mu, b_mu, ztree_options, adapt.method=adapt.method)

  if(doqn < 0){ return(list(move=move, L=L, Y=Y)) }

  temp_omega <- move$omega;
  temp_theta <- move$theta;

  temp_omega[temp_omega >= 1 - 1e-14]=1 - 1e-14
  temp_omega[temp_omega <= 1e-14]=1e-14

  temp_theta[temp_theta >= 1 - 1e-14]=1 - 1e-14
  temp_theta[temp_theta <= 1e-14]=1e-14

  temp_omega <- ord.normalizetpx(temp_omega, byrow=TRUE)
  temp_theta <- ord.normalizetpx(temp_theta, byrow=FALSE)

  ## update Y accounting
  Y <- cbind(Y, tpxToNEF(theta=temp_theta, omega=temp_omega))
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
  Lqnup <- try(ord.tpxlpost(fcounts, qnup$omega, qnup$theta,
                            del_beta, a_mu, b_mu, ztree_options,
                            adapt.method=adapt.method), silent=TRUE)


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


ord.tpxOmegaStart <- function(X, theta)
  {
    if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }
    omega <- try(tcrossprod_simple_triplet_matrix(X, solve(t(theta)%*%theta)%*%t(theta)), silent=TRUE )
    if(inherits(omega,"try-error")){ return( matrix( 1/ncol(theta), nrow=nrow(X), ncol=ncol(theta) ) ) }
    omega[omega <= 0] <- .5
    return( ord.normalizetpx(omega, byrow=TRUE) )
  }


## fast computation of sparse P(X) for X>0
ord.tpxQ <- function(theta, omega, doc, wrd){

  if(length(wrd)!=length(doc)){stop("index mis-match in ord.tpxQ") }
  if(ncol(omega)!=ncol(theta)){stop("theta/omega mis-match in ord.tpxQ") }

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
ord.tpxMixQ <- function(X, omega, theta, grp=NULL, qhat=FALSE){
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
            PACKAGE="maptpx")
  return(list(omega=t( matrix(bck$tomega, nrow=K) ), theta=matrix(bck$theta, ncol=K)))
}

ord.tpxinit <- function(fcounts, X, K1, alpha, verb, param_set, del_beta, a_mu, b_mu,
                    ztree_options, tol, admix, grp, tmax, wtol, qn, acc){
  ## initheta can be matrix, or c(nK, tmax, tol, verb)
  ini_mu_tree_set <- mu_tree_build_set(param_set);
  levels <- length(ini_mu_tree_set[[1]]);
  initheta <- do.call(cbind, lapply(1:K1, function(l) ini_mu_tree_set[[l]][[levels]]/ini_mu_tree_set[[l]][[1]]));

  if(is.matrix(initheta)){
    if(ncol(initheta)!=K1){ stop("mis-match between initheta and K.") }
    if(prod(initheta>0) != 1){ stop("use probs > 0 for initheta.") }
    return(ord.normalizetpx(initheta, byrow=FALSE)) }

  if(is.matrix(alpha)){
    if(nrow(alpha)!=ncol(X) || ncol(alpha)!=K1){ stop("bad matrix alpha dimensions; check your K") }
    return(ord.normalizetpx(alpha, byrow=FALSE)) }

  if(is.null(initheta)){ ilength <- K1-1 }
  else{ ilength <- initheta[1] }
  if(ilength < 1){ ilength <- 1 }

  ## set number of initial steps
  if(length(initheta)>1){ tmax <- initheta[2] }else{ tmax <- 3 }
  ## set the tolerance
  if(length(initheta)>2){ tol <- initheta[3] }else{ tol <- 0.5 }
  ## print option
  if(length(initheta)>3){ verb <- initheta[4] }else{ verb <- 0 }


  if(verb){ cat("Building initial topics")
    if(verb > 1){ cat(" for K = ") }
    else{ cat("... ") } }


    ## Solve for map omega in NEF space
    fit <- ord.tpxfit(fcounts=fcounts, X=X, param_set=param_set, del_beta=del_beta, a_mu=a_mu, b_mu=b_mu,
                  ztree_options=ztree_options, tol=tol, verb=verb, admix=TRUE, grp=NULL, tmax=tmax, wtol=-1,
                  qn=-1, acc = acc);

  initheta <- fit$theta;
  return(initheta)
}
