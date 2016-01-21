##### Estimation for Topic Models ######
## intended main function; provides defaults and fits topic model for the user defined K
ord_topics <- function(counts, K, shape=NULL, initopics=NULL, tol=0.1,
                  ord=TRUE, del_beta, a_mu, b_mu, verb=TRUE, ...)
  ## tpxselect defaults: tmax=10000, wtol=10^(-4), qn=100, grp=NULL, admix=TRUE, nonzero=FALSE, dcut=-10
{
  if(log(dim(counts)[2])%%log(2)!=0) stop("number of features not a power of 2")

  ## Convert the counts matrix to triplet matrix
  X <- CheckCounts(counts)
  p <- ncol(X)
  n <- nrow(X)

  if(verb)
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
  levels <- log(dim(counts)[2])/log(2)+1;
  theta_tree_start <- lapply(1:K, function(s) return(mra_tree_prior_theta(levels,del_beta)));

  param_set_start <- param_extract_mu_tree(theta_tree_start);

  fit <- tpxfit(counts=counts, X=X, param_set=param_set_start, del_beta=del_beta, a_mu=a_mu, b_mu=b_mu,
                tol=tol, verb=verb, admix=admix, grp=grp, tmax=tmax, wtol=wtol, qn=qn);


  #initopics <- tpxinit(X[1:min(ceiling(nrow(X)*.05),100),], initopics, K[1], shape, verb)

  ## either search for marginal MAP K and return bayes factors, or just fit
  ## tpx <- tpxSelect(X, K, bf, initopics, alpha=shape, tol, kill, verb, ...)
  ## K <- tpx$K

  ## clean up and out
  if(ord){ worder <- order(col_sums(fit$omega), decreasing=TRUE) } # order by decreasing usage
  else{ worder <- 1:K }
  ## Main parameters
  mu_tree_set <- mu_tree_build_set(fit$param_set);
  theta <- do.call(cbind, lapply(1:nclus, function(l) mu_tree_set[[l]][[levels]]/mu_tree_set[[l]][[1]]));
  theta=matrix(theta[,worder], ncol=K, dimnames=list(phrase=dimnames(X)[[2]], topic=paste(1:K)) )
  omega=matrix(fit$omega[,worder], ncol=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  if(nrow(omega)==nrow(X)){ dimnames(omega)[[1]] <- dimnames(X)[[1]] }

  ## topic object
  out <- list(K=K, theta=theta, omega=omega, param_set=fit$param_set, loglik=fit$L, X=X, null=null)
  class(out) <- "topics"
  invisible(out) }


