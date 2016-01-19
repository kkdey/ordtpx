
## test code for ordtpx

del_beta <- c(2,3,4,2);
levels <- 4
theta_tree <- mra_tree_prior_theta(levels, del_beta);

a_mu <- 3;
b_mu <- 4;

mu_tree <- mra_tree_prior_mu(levels, del_beta, a_mu, b_mu);

nclus <- 5;
mu_tree_set <- lapply(1:nclus, function(s) return(mra_tree_prior_mu(levels,del_beta, a_mu, b_mu)));
