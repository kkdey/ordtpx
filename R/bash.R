
## Bash : Binomial adaptive shrinkage

postmean = function(m, betahat,sebetahat,v){
  UseMethod("postmean")
}

#' @export
postmean.default = function(m, x, flow_count_child){
  colSums(comppostprob(m, x, flow_count_child) * comp_postmean(m, x, flow_count_child))
}

comppostprob=function(m,x,s,v){
  UseMethod("comppostprob")
}

comppostprob.default = function(m,x,flow_count_child,FUN="+"){
  lpost = log_compdens_conv(m,x,flow_count_child,FUN="+") + log(m$pi) # lpost is k by n of log posterior prob (unnormalized)
  lpmax = apply(lpost,2,max) #dmax is of length n
  tmp = exp(t(lpost)-lpmax) #subtracting the max of the logs is just done for numerical stability
  tmp = tmp/rowSums(tmp)
  ismissing = is.na(x)
  tmp[ismissing,]=m$pi
  t(tmp)
}


matrix_dens = function(flow_count_child, alphavec){
  k = length(alphavec)
  flow.prop <- flow_count_child[c(TRUE,FALSE)]/ (flow_count_child[c(TRUE,FALSE)]+ flow_count_child[c(FALSE,TRUE)])
  n = length(flow.prop)
  f1 <- outer(flow_count_child[c(TRUE,FALSE)], alphavec, "+")
  f2 <- outer(flow_count_child[c(FALSE,TRUE)], alphavec, "+")
  ldens = dbeta(flow.prop,f1,f2,log=TRUE)
  maxldens = apply(ldens, 1, max)
  ldens = ldens - maxldens
  return(exp(ldens))
}




initpi = function(k,n,null.comp,randomstart){
  if(randomstart){
    pi = rgamma(k,1,1)
  } else {
    if(k<n){
      pi=rep(1,k)/n #default initialization strongly favours null; puts weight 1/n on everything except null
      pi[null.comp] = (n-k+1)/n #the motivation is data can quickly drive away from null, but tend to drive only slowly toward null.
    } else {
      pi=rep(1,k)/k
    }
  }
  pi=normalize(pi)
  return(pi)
}

setprior=function(prior,k,nullweight,null.comp){
  if(!is.numeric(prior)){
    if(prior=="nullbiased"){ # set up prior to favour "null"
      prior = rep(1,k)
      prior[null.comp] = nullweight #prior 10-1 in favour of null by default
    }else if(prior=="uniform"){
      prior = rep(1,k)
    } else if(prior=="unit"){
      prior = rep(1/k,k)
    }
  }
  if(length(prior)!=k | !is.numeric(prior)){
    stop("invalid prior specification")
  }
  return(prior)
}





autoselect.mixsd = function(betahat,sebetahat,mult){
  sebetahat=sebetahat[sebetahat!=0] #To avoid exact measure causing (usually by mistake)
  sigmaamin = min(sebetahat)/10 #so that the minimum is small compared with measurement precision
  if(all(betahat^2<=sebetahat^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  }else{
    sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2
  }
  if(mult==0){
    return(c(0,sigmaamax/2))
  }else{
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}

autoselect.mixalpha <- function(flow.prophat, z_level, mult=sqrt(2)){
  sebetahat <- sqrt(flow.prophat*(1-flow.prophat)/z_level)
  betahat_filtered <- flow.prophat[which(sebetahat!=0 & sebetahat!=Inf)]
  sebetahat_filtered <- sebetahat[which(sebetahat!=0 & sebetahat!=Inf)];
  autoselect_sd <- autoselect.mixsd(betahat = betahat_filtered,
                                     sebetahat=sebetahat_filtered,
                                     mult=mult);
  autoselect_alpha <- 0.5*((1/(4*autoselect_sd^2)) - 1);
  autoselect_alpha[autoselect_alpha > 10^7] = 10^7
  return(unique(autoselect_alpha))
}

mixalpha <- autoselect.mixalpha(flow.prophat = flow.prophat,
                                z_level = z_level,
                                mult=sqrt(2))
k = length(mixalpha)
g=beta.mix(pi,mixalpha)


############################### METHODS FOR betamix class ###########################

#' @title Constructor for normalmix class
#'
#' @description Creates an object of class normalmix (finite mixture of univariate normals)
#'
#' @details None
#'
#' @param pi vector of mixture proportions
#' @param mean vector of means
#' @param sd vector of standard deviations
#'
#' @return an object of class normalmix
#'
#' @export
#'
#' @examples normalmix(c(0.5,0.5),c(0,0),c(1,2))
#'


beta.mix = function(pi, alpha){
  structure(data.frame(pi, alpha),class="beta.mix")
}

#' @title comp_sd.normalmix
#' @description returns sds of the normal mixture
#' @param m a normal mixture distribution with k components
#' @return a vector of length k
#' @export
#'
comp_alpha.betamix = function(m){
  m$alpha
}

compdens.betamix = function(m,y,log=FALSE){
  k=ncomp(m)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(dbeta(d, m$alpha, m$alpha, log),nrow=k))
}

log_compdens_conv.betamix = function(m,x,flow_count_child, FUN="+"){
  if(length(s)==1){s=rep(s,length(x))}
  shape1 <- outer(m$alpha, flow_count_child[c(TRUE,FALSE)], FUN); ## posterior shape1 parameter
  shape2 <- outer(m$alpha, flow_count_child[c(FALSE, TRUE)], FUN); ## posterior shape2 parameter
  return(t(dbeta(x, shape1, shape2, log=TRUE)))
}


compdens_conv.betamix = function(m,x,flow_count_child, FUN="+"){
  if(length(s)==1){s=rep(s,length(x))}
  shape1 <- outer(flow_count_child[c(TRUE,FALSE)], m$alpha, FUN); ## posterior shape1 parameter
  shape2 <- outer(flow_count_child[c(FALSE, TRUE)], m$alpha, FUN); ## posterior shape2 parameter
  return(t(dbeta(x, shape1, shape2)))
}

comp_postmean.betamix = function(m, x, flow_count_child){
  num <- outer(flow_count_child[c(TRUE,FALSE)], m$alpha, FUN);
  den <- num + outer(flow_count_child[c(FALSE, TRUE)], m$alpha, FUN);
  tmp <- num/den;
  tmp[ismissing,]=m$mean
  t(tmp)
}

comp_postsd.betamix = function(m ,x, flow_count_child){
  num <- outer(flow_count_child[c(TRUE,FALSE)], m$alpha, FUN);
  den <- num + outer(flow_count_child[c(FALSE, TRUE)], m$alpha, FUN);
  tmp <- (num*den)/((num+den)^2 * (num + den +1))
}

comp_postmean2.betamix = function(m, x, flow_count_child){
  comp_postsd(m, x, flow_count_child)^2 + comp_postmean(m, x, flow_count_child)^2
}

dens_conv_mixlik = function(m,x,flow_count_child,pilik,FUN="+"){
  UseMethod("dens_conv_mixlik")
}
dens_conv_mixlik.default = function(m,x,flow_count_child, pilik,FUN="+"){
  l=dim(pilik)[2]
  colSums(rep(m$pi,l) * compdens_conv_mixlik(m,x,flow_count_child,FUN))
}


comppostprob_mixlik=function(m,x,flow_count_child,pilik){
  UseMethod("comppostprob_mixlik")
}

comppostprob_mixlik.default = function(m,x,flow_count_child,pilik){
  l=dim(pilik)[2]
  k=length(m$pi)
  tmp= (t(rep(m$pi,l) * compdens_conv_mixlik(m,x,flow_count_child,pilik))/dens_conv_mixlik(m,x,flow_count_child,pilik))
  group=rep(1:k,l)
  return(rowsum(t(tmp),group))
}


compdens_conv_mixlik = function(m, x, flow_count_child, pilik, FUN="+"){
  UseMethod("compdens_conv_mixlik")
}
compdens_conv_mixlik.default = function(m, x, flow_count_child, pilik, FUN="+"){
  dens=NULL
  for (i in 1:dim(pilik)[2]){
    dens=rbind(dens,pilik[,i]*compdens_conv(m,x,s[,i],v[i],FUN))
  }
  return(dens)
}

comppostprob_mixlik2=function(m,x,s,v,pilik){
  UseMethod("comppostprob_mixlik2")
}
#' @export
comppostprob_mixlik2.default = function(m,x,flow_count_child,pilik){
  l=dim(pilik)[2]
  k=length(m$pi)
  tmp= (t(rep(m$pi,l) * compdens_conv_mixlik(m,x,flow_count_child,pilik))/dens_conv_mixlik(m,x,flow_count_child,pilik))
  ismissing = (is.na(x) | apply(is.na(flow_count_child[c(TRUE,FALSE)]),1,sum))
  tmp[ismissing,]=rep(m$pi,l)/l;
  return(t(tmp))
}

loglik_conv_mixlik = function(m,flow_count_child,v,pilik,FUN="+"){
  UseMethod("loglik_conv_mixlik")
}

loglik_conv_mixlik.default = function(m,flow_count_child,pilik,FUN="+"){
  flow_prop <- flow_count_child[c(TRUE,FALSE)]/(flow_count_child[c(TRUE,FALSE)]+ flow_count_child[c(FALSE,TRUE)]);
  sum(log(dens_conv_mixlik(m,flow_prop,flow_count_child, pilik,FUN)))
}



negpenloglik = function(pi,matrix_lik,prior){return(-penloglik(pi,matrix_lik,prior))}

penloglik = function(pi, matrix_lik, prior){
  pi = normalize(pmax(0,pi))
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  loglik = sum(log(m.rowsum))
  subset = (prior != 1.0)
  priordens = sum((prior-1)[subset]*log(pi[subset]))
  return(loglik+priordens)
}

normalizer = function(x){return(x/sum(x))}

fixpoint = function(pi, matrix_lik, prior){
  pi = normalize(pmax(0,pi)) #avoid occasional problems with negative pis due to rounding
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  classprob = m/m.rowsum #an n by k matrix
  pinew = normalize(colSums(classprob) + prior - 1)
  return(pinew)
}



mixEM = function(matrix_lik,prior,pi_init=NULL,control=list()){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)

  k=dim(matrix_lik)[2]
  if(is.null(pi_init)){
    pi_init = rep(1/k,k)# Use as starting point for pi
  }
  res = squarem(par=pi_init,fixptfn=fixpoint, objfn=negpenloglik,matrix_lik=matrix_lik, prior=prior, control=controlinput)
  return(list(pihat = normalize(pmax(0,res$par)), B=res$value.objfn,
              niter = res$iter, converged=res$convergence))
}

mixIP = function(matrix_lik, prior, pi_init = NULL, control = list()){
  if(!require("REBayes",quietly=TRUE)){stop("mixIP requires installation of package REBayes")}
  n = nrow(matrix_lik)
  k = ncol(matrix_lik)
  #A = matrix_lik
  A = rbind(diag(length(prior)),matrix_lik) # add in observations corresponding to prior
  w = c(prior-1,rep(1,n))
  A = A[w!=0,]    #remove zero weight entries, as these otherwise cause errors
  w = w[w!=0]
  #w = rep(1,n+k)
  res = REBayes::KWDual(A, rep(1,k), ashr:::normalize(w), control=control)
  return(list(pihat = normalize(res$f), niter = NULL, converged=(res$status=="OPTIMAL")))
}

mixVBEM = function(matrix_lik, prior, pi_init = NULL,control=list()){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)

  k=ncol(matrix_lik)
  if(is.null(pi_init)){  pi_init = rep(1,k)  }# Use as starting point for pi
  res = squarem(par=pi_init,fixptfn=VBfixpoint, objfn=VBnegpenloglik,matrix_lik=matrix_lik, prior=prior, control=controlinput)

  return(list(pihat = res$par/sum(res$par), B=res$value.objfn, niter = res$iter, converged=res$convergence,post=res$par))
}


VBfixpoint = function(pipost, matrix_lik, prior){
  n=nrow(matrix_lik)
  k=ncol(matrix_lik)
  avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
  classprob = avgpipost*matrix_lik
  classprob = classprob/rowSums(classprob) # n by k matrix
  pipostnew = colSums(classprob) + prior
  return(pipostnew)
}

VBnegpenloglik=function(pipost,matrix_lik,prior){
  return(-VBpenloglik(pipost,matrix_lik,prior))
}

VBpenloglik = function(pipost, matrix_lik, prior){
  n=nrow(matrix_lik)
  k=ncol(matrix_lik)
  avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
  classprob = avgpipost*matrix_lik
  classprob = classprob/rowSums(classprob) # n by k matrix

  B= sum(classprob*log(avgpipost*matrix_lik),na.rm=TRUE) - diriKL(prior,pipost) - sum(classprob*log(classprob))
  return(B)
}


estimate_mixprop = function(counts,
                            g,
                            prior,
                            optmethod=c("mixEM","mixVBEM","mixIP"),
                            null.comp=1,
                            control=list()){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  optmethod=match.arg(optmethod)
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)

  controlinput$tol = min(0.1/n,1.e-7) # set convergence criteria to be more stringent for larger samples

  if(controlinput$trace==TRUE){tic()}

  counts_odd <- counts[c(TRUE,FALSE)]
  counts_even <- counts[c(FALSE, TRUE)]

  flow.prophat <- counts_odd/(counts_odd+counts_even)

  matrix_llik = t(log_compdens_conv.betamix(g, flow.prophat, counts)) #an n by k matrix
  matrix_llik = matrix_llik - apply(matrix_llik,1, max) #avoid numerical issues by subtracting max of each row
  matrix_lik = exp(matrix_llik)

  if(optmethod=="mixVBEM" || max(prior[-1])>1 || min(gradient(matrix_lik)+prior[1]-1,na.rm=TRUE)<0){
    fit=do.call(optmethod,args = list(matrix_lik= matrix_lik, prior=prior, pi_init=pi_init, control=controlinput))
  } else {
    fit = list(converged=TRUE,pihat=c(1,rep(0,k-1)))
  }

  ## check if IP method returns negative mixing proportions. If so, run EM.
  if (optmethod == "mixIP" & (min(fit$pihat) < -10 ^ -12)) {
    message("Interior point method returned negative mixing proportions.\n Switching to EM optimization.")
    optmethod <- "mixEM"
    fit = do.call(optmethod, args = list(matrix_lik = matrix_lik,
                                         prior = prior, pi_init = pi_init,
                                         control = controlinput))
  }

  if(!fit$converged){
    warning("Optimization failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.")
  }

  pi = fit$pihat
  converged = fit$converged

  loglik.final =  penloglik(pi,matrix_lik,1) #compute penloglik without penalty
  null.loglik = sum(log(matrix_lik[,null.comp]))
  g$pi=pi
  if(controlinput$trace==TRUE){toc()}

  return(list(loglik=loglik.final,null.loglik=null.loglik,
              matrix_lik=matrix_lik,converged=converged,g=g))

}


bash.workhorse = function(counts,
                          gridmult=sqrt(2),
                          randomstart=FALSE,
                          pointmass = TRUE,
                          optmethod = c("mixIP","mixEM","mixVBEM"),
                          nullweight=10,
                          prior=c("nullbiased","uniform","unit"),
                          g=NULL,
                          control=list()
){
  if(log(length(counts))%%log(2)!=0){
    stop("The length of the counts vector must be a power of 2")
  }

  control.default=list(K = 1, method=3, square=TRUE, step.min0=1,
                       step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07,
                       maxiter=5000, trace=FALSE)
  if(n>50000){control.default$trace=TRUE}
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)
  if(controlinput$maxiter==0){
    stop("option control$maxiter=0 deprecated; used fixg=TRUE instead")
  }

  if(n==0){
    stop("Error: all input values are missing")
  }

  counts_odd <- counts[c(TRUE,FALSE)]
  counts_even <- counts[c(FALSE, TRUE)]

  flow.prophat <- counts_odd/(counts_odd+counts_even)

  mixalpha <- autoselect.mixalpha(flow.prophat = flow.prophat,
                                  z_level = (counts_even+counts_odd),
                                  mult=gridmult)
  if(pointmass){ mixsd = c(10^10,mixsd) }

  null.comp <- which.max(alpha)
  k <- length(mixalpha)
  pi = initpi(k,n,null.comp,randomstart)
  prior = setprior(prior,k,nullweight,null.comp)


  g <- betamix(pi, mixalpha)
  pi.fit <- estimate_mixprop(counts, g , prior, null.comp=null.comp,
                          optmethod=optmethod, control=controlinput)

  flow.posteriormean <- postmean(pi.fit$g, flow.prophat, flow_count_child)
  counts_odd_shrunk <- flow.posteriormean*(counts_even+counts_odd);
  counts_even_shrunk <- (1-flow.posteriormean)*(counts_even+counts_odd)

  counts_shrunk <- c(rbind(counts_odd_shrunk, counts_even_shrunk));
  loglik_final <- pi.fit$loglik;
  ll <- list("counts"=counts_shrunk, "loglik"=loglik_final);
  return(ll)

}


