---
title: 'Validating ordtpx: Ordered topic model'
author: "Kushal K Dey"
date: "April 11, 2016"
output: html_document
---

## Overview

In this script, we perform test runs to validate the *ordtpx* package. We take out chunks from the MRA analysis of the *ordtpx* code and then try to dry run or test run these code chunks.

```{r}
theta_tree <- mra_tree_prior_theta(4, c(1,20,200))
print(theta_tree)
plot(theta_tree[[4]], type="l")
```

Note that the $\Theta$ matrix is pretty smooth because the lower level $\beta$ values are generated from a very concentrated Beta distribution $Beta(200,200)$. If we want less smooth theta, one may take a less concentrated prior at the bottom levels.

```{r}
theta_tree <- mra_tree_prior_theta(4, c(1,10,20))
print(theta_tree)
plot(theta_tree[[4]], type="l")
```

We now build the mu tree over the theta tree we generate above.

```{r}
mu_tree <- mra_tree_prior_mu(4, c(1,20,200),2,2);
print(mu_tree)
plot(mu_tree[[4]], type="l")
```

We extract the theta tree from the mu tree we generated. The shape of the extracted theta would be similar to that of the mu (lower level)

```{r}
theta_tree_extract <- extract_theta_tree_from_mu(mu_tree)
plot(theta_tree_extract[[4]], type="l")
```


If we know the leaf node theta/ mu values, we can construct the tree by adding the adjacent values together.


```{r}
theta_tree <- mra_tree_prior_theta(4, c(1,20,200))
print(theta_tree)
plot(theta_tree[[4]], type="l")
```

```{r}
theta_tree_construct <- mra_bottom_up(theta_tree[[4]])
print(theta_tree)
print(theta_tree_construct) ## should be same as theta_tree
```


Now we generate counts data using the Multinomial model.

```{r}
nclus <- 2;
del_beta <- c(2,10,50, 100, 500, 1000, 2000);
levels <- 8

a_mu <- 30;
b_mu <- 80;

mu_tree_set <- lapply(1:nclus, function(s) return(mra_tree_prior_mu(levels,del_beta, a_mu, b_mu)));
param_set <- param_extract_mu_tree(mu_tree_set)

prior_calc <- prior_calc_fn(param_set, del_beta, a_mu, b_mu)

theta_sim <- do.call(rbind, lapply(1:nclus, function(l) mu_tree_set[[l]][[levels]]/mu_tree_set[[l]][[1]]));

n.out <- 200
omega_sim <- rbind( cbind( rep(1, n.out), rep(0, n.out)), 
                  cbind( rep(0, n.out), rep(1, n.out)),
                  cbind( seq(0.6, 0.4, length.out = n.out), 
                         1- seq(0.6, 0.4,length.out=n.out)) )
dim(omega_sim)
K <- dim(omega_sim)[2]

barplot(t(omega_sim),
        col = 2:(K+1),
        axisnames = F, space = 0, border = NA, 
        main=paste("No. of clusters=", K),
        las=1, ylim=c(0,1), cex.axis=1.5, cex.main=1.4)


counts <- t(do.call(cbind,lapply(1:dim(omega_sim)[1], function(x) rmultinom(1,1000,prob=omega_sim[x,]%*%theta_sim))));

```

We first apply **maptpx** package.

```{r}
topic_clus <- maptpx::topics(counts, K=2, tol=0.1)
K <- 2
barplot(t(topic_clus$omega),
        col = 2:(K+1),
        axisnames = F, space = 0, border = NA, 
        main=paste("No. of clusters=", K),
        las=1, ylim=c(0,1), cex.axis=1.5, cex.main=1.4)
```

We also plot the theta matrices 

```{r}
plot(topic_clus$theta[,1], type="l")
plot(topic_clus$theta[,2], type="l")
```

We now perform *ordtpx*.

```{r}
library(slam)
initopics=NULL
tol=0.1
ord=TRUE
ztree_options=1
verb=1
reflect=TRUE
tmax=10000
wtol=10^(-4)
qn=100
grp=NULL
admix=TRUE
nonzero=FALSE
dcut=-10

K <- 2

del_beta <- c(2,10,50, 100, 500, 1000, 2000);
levels <- 8

a_mu <- 30;
b_mu <- 80;

del_beta_1 <- c(2,2,4,4,6,6,7)
system.time(ord_topics <- ord_topics(counts, K=2, del_beta = del_beta_1, a_mu=3, b_mu=4, ztree_options=2, tol=0.001));

barplot(t(omega),
        col = 2:(K+1),
        axisnames = F, space = 0, border = NA, 
        main=paste("No. of clusters=", K),
        las=1, ylim=c(0,1), cex.axis=1.5, cex.main=1.4)
```

```{r}
plot(theta[,1], type="l")
plot(theta[,2], type="l")
```

