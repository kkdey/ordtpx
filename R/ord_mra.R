
## build a prior theta MRA tree

# del_beta is a vector of length (#Levels-1) where #Levels is the number of levels in the MRA tree

mra_tree_prior_theta <- function(S, del_beta)
{
  if(length(del_beta) !=(S-1)) stop("The length of flow proportions vector (del_beta) should be 1 less than the number of levels in MRA tree")
  theta <- vector(mode="list", length=S);
  theta[[1]] <- 1;
  for(s in 2:S){
    beta_vec <- rbeta(length(theta[[(s-1)]]), del_beta[(s-1)], del_beta[(s-1)]);
    theta[[s]] <- as.vector(rbind(beta_vec*theta[[(s-1)]], (1-beta_vec)*theta[[(s-1)]]));
  }
  return(theta)
}


## build a prior mu MRA tree

mra_tree_prior_mu <- function(S, del_beta, a_mu, b_mu)
{
  if(length(del_beta) !=(S-1)) stop("The length of flow proportions vector (del_beta) should match the number of levels in MRA tree")
  mu <- vector(mode="list", length=S);
  mu <- lapply(mra_tree_prior_theta(S,del_beta), "*", rgamma(1,a_mu,b_mu));
  return(mu)
}

## extract the theta tree from the mu MRA tree

extract_theta_tree_from_mu <- function(mu)
{
  if(!is.list(mu)) stop("Argument must be a list")
  theta <- lapply(mu,"/",mu[[1]]);
  return(theta)
}

## Construct the MRA tree bottom up given the leaf nodes values known

mra_bottom_up <- function(leaf_val)
{
  if(log(length(leaf_val))%%log(2)!=0) stop("size of vector not a power of 2- wavelet not constructed")
  wavefit <- wavethresh::wd(leaf_val,  filter.number=1,   family="DaubExPhase")
  S <- floor(log(length(leaf_val))/log(2));
  out <- lapply(0:S, function(s) wavethresh::accessC(wavefit,level=s));
  if(length(out)!=(S+1)) stop("size of the list does not match with number of levels")
  scaled_out <- lapply(1:(S+1), function(s) return(out[[s]]*(sqrt(2))^(S+1-s)))
  return(scaled_out)
}

## Construct the Z value MRA tree from the counts data and current iterates of omega and theta

z_tree_construct <- function(counts, omega_iter, theta_iter, ztree_options=c(1,2))
{
  if(ztree_options==2){
  row_total <- rowSums(counts);
  z_leaf_est <- round(sweep(theta_iter, MARGIN=1, colSums(sweep(omega_iter, MARGIN = 1, row_total, "*")), "*"));
  }
  if(ztree_options==1){
    z_leaf_est <- (t(omega_iter) %*% (counts/(omega_iter %*% theta_iter)))*theta_iter;
  }
  z_tree_out <- vector(mode="list",length=dim(omega_iter)[2])
  for(k in 1: dim(omega_iter)[2]){
  z_tree_out[[k]] <- mra_bottom_up(z_leaf_est[k,]);
  }
  return(z_tree_out)
}

## Extract the beta parameters and the root or top level  mu parameter in a hierarchy of lists
## from Z value MRA tree

param_extract_ztree <- function(z_tree_in, del_beta, a_mu, b_mu)
{
  if(!is.list(z_tree_in) | !is.list(z_tree_in[[1]])) stop("z_tree input must be a list of lists")
  K <- length(z_tree_in);
  beta_set <- vector(mode="list",length=K);
  param_set <- parallel::mclapply(1:K, function(k)
                {
                    intree <- z_tree_in[[k]];
                    S <- length(intree);
                    beta_out <- vector(mode="list",length=S-1);
                    for(s in 2:S){
                      beta_out[[(s-1)]] <- (intree[[s]][c(TRUE,FALSE)] + del_beta[(s-1)]-1)/(intree[[(s-1)]]+ 2*(del_beta[(s-1)]-1));
                    }
                    mu_out <- (intree[[1]]+a_mu - 1)/(b_mu+1);
                    out_list <- list("beta_tree"=beta_out,"mu"=mu_out);
                    return(out_list)
  }, mc.cores=parallel::detectCores())
  return(param_set)
}

## Build a set of mu trees across the topics, given the hierarchical loist structure of the beta values
## and the top level mu value.

mu_tree_build <- function(beta_tree, mu_top)
{
  S <- (length(beta_tree)+1);
  mu_tree <- vector(mode="list", length=S);
  mu_tree[[1]] <- mu_top;
  for(s in 2:S){
    mu_tree[[s]] <- as.vector(rbind(beta_tree[[(s-1)]]*mu_tree[[(s-1)]], (1-beta_tree[[(s-1)]])*mu_tree[[(s-1)]]));
  }
  return(mu_tree)
}


mu_tree_build_set <- function(param_set)
{
  K <- length(param_set)
  mu_tree_set <- lapply(1:K, function(k)
              {
                out <- mu_tree_build(param_set[[k]]$beta_tree,param_set[[k]]$mu);
                return(out)
  })
  return(mu_tree_set)
}

## Build a set of theta trees from mu trees

theta_tree_build_set <- function(mu_tree_set)
{
  K <- length(mu_tree_set)
  theta_tree_set <- lapply(1:K, function(k)
                                {
                                      out <- extract_theta_tree_from_mu(mu_tree_set[[k]])
                                      return(out)
  })
  return(theta_tree_set)
}

## Extract the parameters (the beta values and top mu values) for all clusters in a hierarchical list structure
## given the MRA tree structure of mu values for different clusters

param_extract_mu_tree <- function(mu_tree_set)
{
  if(!is.list(mu_tree_set) | !is.list(mu_tree_set[[1]])) stop("mu_tree input must be a list of lists")
  K <- length(mu_tree_set);
  beta_set <- vector(mode="list",length=K);
  param_set <- parallel::mclapply(1:K, function(k)
  {
    intree <- mu_tree_set[[k]];
    S <- length(intree);
    beta_out <- vector(mode="list",length=S-1);
    for(s in 2:S){
      beta_out[[(s-1)]] <- (intree[[s]][c(TRUE,FALSE)])/(intree[[(s-1)]]);
    }
    mu_out <- intree[[1]];
    out_list <- list("beta_tree"=beta_out,"mu"=mu_out);
    return(out_list)
  }, mc.cores=parallel::detectCores());
  return(param_set)
}



## Prior calculation given the beta values and the top mu value (parameters we can extract from param_extract_mu_tree
## or param_extract_ztree)

prior_calc_fn <- function(param_set_in, del_beta, a_mu, b_mu)
{
  nlevels <- length(del_beta);
  beta_set_across_classes <- vector(mode="list",length=nlevels)

  for(l in 1:(nlevels))
  {
    beta_set_across_classes[[l]] <- unlist(lapply(1:nclus, function(k) return(param_set[[k]]$beta_tree[[l]])))
  }

  mu_set_across_classes <- unlist(lapply(1:nclus, function(k) return(param_set[[k]]$mu)));

  prior_calc <- 0;
  for(l in 1:nlevels)
  {
    prior_calc <- prior_calc + sum(dbeta(beta_set_across_classes[[l]],del_beta[l],del_beta[l],log=TRUE));
  }

  prior_calc <- prior_calc + sum(dgamma(mu_set_across_classes,a_mu,b_mu,log=TRUE));
  return(prior_calc)
}

## Log likelihood calcultion given the set of Z values trees for different clusters and the parameter trees.

loglik_fn <- function(z_tree, param_set)
{
  fn_out <- sum(unlist(parallel::mclapply(1:nclus, function(k)
  {
    intree <- z_tree[[k]];
    S <- length(intree);
    loglik_out <- 0;
    for(s in 2:S){
      loglik_out <- loglik_out + sum(dbinom(round(intree[[s]][c(TRUE,FALSE)]), round(intree[[(s-1)]]), param_set[[k]]$beta_tree[[(s-1)]],log=TRUE));
    }
    loglik_out <- loglik_out + dpois(round(intree[[1]]),param_set[[k]]$mu, log=TRUE);
    return(loglik_out)
  }, mc.cores=parallel::detectCores())))
  return(fn_out)
}

## Posterior computation: similar to tpxlpost() function in maptpx package of Matt Taddy

tpxlpost <- function(counts, omega_iter, param_set, del_beta, a_mu, b_mu, ztree_options=c(1,2))
{
  z_tree <- z_tree_construct(counts, omega_iter = omega_iter, theta_iter = theta_iter, ztree_options = 1)
  loglik_value <- loglik_fn(z_tree, param_set = param_set);
  prior_calc <- prior_calc_fn(param_set, del_beta, a_mu, b_mu);
  prior_omega <- (1/dim(omega_iter)[2])*sum(log(omega_iter[omega_iter > 0]))
  posterior <- prior_calc + prior_omega + loglik_value;
  return(posterior)
}

