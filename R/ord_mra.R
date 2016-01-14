
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
  mu <- lapply(mra_tree_prior_theta(S,del_beta), "*", rgamma(1,a_mu,b_mu));
  return(mu)
}

extract_theta_tree_from_mu <- function(mu)
{
  if(!is.list(mu)) stop("Argument must be a list")
  theta <- lapply(mu,"/",mu[[1]]);
  return(theta)
}

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
