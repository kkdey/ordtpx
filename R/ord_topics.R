##### Estimation for Topic Models ######
## intended main function; provides defaults and fits topic model for the user defined K
ord_topics <- function(counts,
                       K,
                       shape=NULL,
                       initopics=NULL,
                       tol=0.1,
                       ord=TRUE,
                       del_beta=NULL,
                       a_mu=NULL,
                       b_mu=NULL,
                       ztree_options=1,
                       verb=1,
                       reflect=TRUE,
                       tmax=10000,
                       tmax_start=10,
                       wtol=10^(-4),
                       qn=100,
                       grp=NULL,
                       admix=TRUE,
                       nonzero=FALSE,
                       dcut=-10,
                       acc=TRUE,
                       burn_trials=5,
                       init_method = c("mra", "taddy", "kmeans"),
                       adapt.method=c("beta", "smash", "bash"))
{
  ceil <- ceiling(log(dim(counts)[2])/log(2));
  if(log(dim(counts)[2])%%log(2)!=0) {
    cat(sprintf("number of features not a power of 2"));
    if(reflect){
      fcounts <- cbind(counts, counts[,dim(counts)[2]-(1:(2^{ceil}-dim(counts)[2]))]);
    }
    if(!reflect){
      fcounts <- cbind(counts, matrix(0, dim(counts)[1], 2^{ceil}-dim(counts)[2]));
    }}else{
    fcounts <- counts;
    }

    if(adapt.method=="smash" | adapt.method=="bash"){
      del_beta <- NULL
      a_mu <- NULL
      b_mu <- NULL
    }

  levels <- ceil+1;

  ## Convert the counts matrix to triplet matrix
  X <- CheckCounts(fcounts)
  p <- ncol(X)
  n <- nrow(X)

  if(verb > 0)
    cat(sprintf("\nEstimating on a %d document collection.\n", nrow(X)))

  ## check the prior parameters for theta
  if(prod(del_beta>0) != 1){ stop("use veta parameters > 0\n") }

  ## check the list of candidate K values
  if(K <=1){ stop(cat("use K values >= 2")) }
  cat(sprintf("\nFitting a ordered topic model with %d topics \n", K))

  ## Null model log probability
  sx <- sum(X)
  qnull <- col_sums(X)/sx
  null <- sum( X$v*log(qnull[X$j]) ) - 0.5*(n+p)*(log(sx) - log(2*pi))

  ## initialize

  loglik_list <- vector(mode = "list", length = burn_trials)
  init_theta_list <- vector(mode = "list", length = burn_trials)

  iter.max_val <- floor(seq(5,min(burn_trials*5, dim(fcounts)[1]), length.out=burn_trials))


  for(init.repeat in 1:burn_trials){
    if(init.repeat==1){
      inds <- 1: min(ceiling(nrow(X)*.05),100);
    }else{
      inds <- sample(1:dim(fcounts)[1],min(ceiling(nrow(X)*.05),100), replace=FALSE);
    }
    if(init_method=="mra"){
      theta_tree_start <- lapply(1:K, function(s) return(mra_tree_prior_theta(levels,del_beta)));
      param_set_start <- param_extract_mu_tree(theta_tree_start);
      initopics_theta <- ord.tpxinit(fcounts[inds,],
                                     X[inds,], K, shape,
                                     verb, param_set_start, del_beta, a_mu, b_mu, ztree_options, tol,
                                     admix, grp, tmax, wtol, qn, acc);
    }
    if(init_method=="taddy"){
      initopics_theta <-   tpxinit(X[inds,], initopics, K[1],
                                   shape, verb, nbundles=1, use_squarem=FALSE, init.adapt = TRUE)
    }
    if(init_method=="kmeans"){
      kmeans.init=kmeans(fcounts,K,nstart=5, iter.max=iter.max_val[init.repeat])
      phi=t(kmeans.init$centers)
      initopics_theta = apply(phi, 2, function(x) return(x/sum(x)))
    }

    initopics_theta[initopics_theta==0] <- 1e-14;
    initopics_theta[initopics_theta==0] <- 1- 1e-14;
    initopics_theta <- ord.normalizetpx(initopics_theta, byrow = FALSE)

    initopics_theta_tree_set <- lapply(1:K, function(k) mra_bottom_up(initopics_theta[,k]));
    initopics_param_set <- param_extract_mu_tree(initopics_theta_tree_set)

    suppressMessages(fit <- ord.tpxfit(fcounts=fcounts, X=X, param_set=initopics_param_set, del_beta=del_beta, a_mu=a_mu, b_mu=b_mu,
                      ztree_options=ztree_options, tol=tol, verb=verb, admix=admix, grp=grp, tmax=tmax_start, wtol=wtol,
                      qn=qn, acc=acc, adapt.method=adapt.method));
    loglik_list[[init.repeat]] <- fit$L
    init_theta_list[[init.repeat]] <- initopics_theta;
  }

  initopics_theta <- init_theta_list[[which.max(unlist(loglik_list))]];

  # if(init_method=="mra"){
  #    theta_tree_start <- lapply(1:K, function(s) return(mra_tree_prior_theta(levels,del_beta)));
  #    param_set_start <- param_extract_mu_tree(theta_tree_start);
  #    initopics_theta <- ord.tpxinit(fcounts[1:min(ceiling(nrow(X)*.05),100),], X[1:min(ceiling(nrow(X)*.05),100),], K, shape,
  #                            verb, param_set_start, del_beta, a_mu, b_mu, ztree_options, tol,
  #                            admix, grp, tmax, wtol, qn, acc);
  # }else{
  #   initopics_theta <-   tpxinit(X[1:min(ceiling(nrow(X)*.05),100),], initopics, K[1],
  #                        shape, verb, nbundles=1, use_squarem=FALSE, init.adapt = TRUE)
  # }

  ## initialize
  #initopics_theta_2 <- ord.tpxinit(X[1:min(ceiling(nrow(X)*.05),100),], initopics, K[1]+3, shape, verb)
  #initopics_theta_2 <- initopics[,sort(sample(1:(K[1]+2), K, replace=FALSE))];

  initopics_theta[initopics_theta==0] <- 1e-14;
  initopics_theta[initopics_theta==0] <- 1- 1e-14;
  initopics_theta <- ord.normalizetpx(initopics_theta, byrow = FALSE)

  initopics_theta_tree_set <- lapply(1:K, function(k) mra_bottom_up(initopics_theta[,k]));
  initopics_param_set <- param_extract_mu_tree(initopics_theta_tree_set)

  fit <- ord.tpxfit(fcounts=fcounts, X=X, param_set=initopics_param_set, del_beta=del_beta, a_mu=a_mu, b_mu=b_mu,
                ztree_options=ztree_options, tol=tol, verb=verb, admix=admix, grp=grp, tmax=tmax, wtol=wtol,
                qn=qn, acc=acc, adapt.method=adapt.method);

  ## clean up and out
  if(ord){ worder <- order(col_sums(fit$omega), decreasing=TRUE) } else{ worder <- 1:K }
  ## Main parameters
  theta=matrix(fit$theta[,worder], ncol=K, dimnames=list(phrase=dimnames(X)[[2]], topic=paste(1:K)) )
  omega=matrix(fit$omega[,worder], ncol=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  if(nrow(omega)==nrow(X)){ dimnames(omega)[[1]] <- dimnames(X)[[1]] }
  ## topic object
  out <- list(K=K, theta=theta, omega=omega, loglik=fit$L, X=X, null=null)
  class(out) <- "topics"
  invisible(out) }


