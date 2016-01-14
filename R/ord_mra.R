
mra_tree_prior_theta <- function(S, del_beta)
{
  if(length(del_beta) !=S) stop("The length of flow proportions vector (del_beta) should match the number of levels in MRA tree")
  theta <- vector(mode="list", length=S);
  theta[[1]] <- 1;
  for(s in 2:S){
    beta_vec <- rbeta(length(theta[[(s-1)]]), del_beta[s], del_beta[s]);
    theta[[s]] <- as.vector(rbind(beta_vec*theta[[(s-1)]], (1-beta_vec)*theta[[(s-1)]]));
  }
  return(theta)
}

mra_tree_prior_mu <- function(S, del_beta, a_mu, b_mu)
{
  if(length(del_beta) !=S) stop("The length of flow proportions vector (del_beta) should match the number of levels in MRA tree")
  mu <- vector(mode="list", length=S);
  mu <- rgamma(1,a_mu,b_mu) * mra_tree_prior_theta(S,del_beta);
  return(mu)
}

