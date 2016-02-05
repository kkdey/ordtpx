
## test code for ordtpx

del_beta <- c(2,3,4);
levels <- 4
theta_tree <- mra_tree_prior_theta(levels, del_beta);

a_mu <- 3;
b_mu <- 4;

mu_tree <- mra_tree_prior_mu(levels, del_beta, a_mu, b_mu);

nclus <- 2;
mu_tree_set <- lapply(1:nclus, function(s) return(mra_tree_prior_mu(levels,del_beta, a_mu, b_mu)));
param_set <- param_extract_mu_tree(mu_tree_set)

prior_calc <- prior_calc_fn(param_set, del_beta, a_mu, b_mu)

theta_iter <- do.call(rbind, lapply(1:nclus, function(l) mu_tree_set[[l]][[levels]]/mu_tree_set[[l]][[1]]));
omega_iter <- rbind(c(0.5,0.5),c(0.8,0.2),c(0.1,0.9));

counts <- t(do.call(cbind,lapply(1:dim(omega_iter)[1], function(x) rmultinom(1,1000,prob=omega_iter[x,]%*%theta_iter))));

system.time(ord_topics <- ordtpx::ord_topics(counts,K=2, del_beta = c(2,9,5), a_mu=2, b_mu=3, ztree_options=2, tol=0.001));
system.time(map_topics <- maptpx::topics(counts, K=2, tol=0.001));


z_tree <- z_tree_construct(counts, omega_iter = omega_iter, theta_iter = theta_iter, ztree_options = 1)

loglik_value <- loglik_fn(z_tree, param_set);

posterior <- ord_tpxlpost(counts, omega_iter, param_set, del_beta, a_mu, b_mu, ztree_options = 1)
