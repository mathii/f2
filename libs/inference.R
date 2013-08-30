## Functions used to infer the ages of chunks.
library(MASS)

########################################################################################################
## loglikelihood for age t, conditional on number of singletons and genetic length
## t: Age (coalescent units)
## Lg: Genetic length
## Ne: Effective pop size
## D: Number of doubletons
## pf: Function giving the probability that a recombination changes the tmrca
## error.params: (shape, rate) for the gamma overestimate in Lg
## S.params - list: S, theta, Lp, Ep: singletons, theta, phyical length, expected singletons
########################################################################################################

loglikelihood.age <-function(t, Lg, Ne, D, pf, error.params=NA, S.params=NA){
  inner <- function(x){return(NA)}
  t1 <- 0
  
  if(all(is.na(error.params))){
    t1 <- (dgamma(Lg, shape=1+1/D, rate=4*Ne*t*pf(t), log=TRUE))
  } else{
    inner <- function(x){
      dgamma(x, shape=1+1/D, rate=4*Ne*t*pf(t))*dgamma(Lg-x,shape=error.params[1],rate=error.params[2])
    }
    ## This is to ensure that we are not integrating over large zero regions
    upper=max(0,Lg-1e-6)                #the 1e-6 is because when shape<0 the gamma dist is not defined at x=0
    new.upper=upper
    while(inner(new.upper)==0 & upper>Lg/16){
      upper=new.upper
      new.upper=upper/2
    }
    t1 <- log(integrate(inner, lower=0, upper=upper, stop.on.error=FALSE)$value)
  }

  t2 <- 0
  if(!all(is.na(S.params))){            #If we supplied singleton information, use that. 
    ## t2 <- dpois(S.params$S, lambda=S.params$theta*S.params$Lp*t + 2*S.params$Ep, log=TRUE) 
    ## Sum of poisson and negative binomial...
    S <- S.params$S
    s.vals <- 0:S
    nbprob <- 1/(1+S.params$Ep)         #Ep is theta*Lp/2/n, here. 
    probs <- dpois(s.vals, lambda=S.params$theta*S.params$Lp*t)*dnbinom(S-s.vals, size=2, prob=nbprob)
    t2 <- log(sum(probs))
  }

  ## bit hacky until we figure out how to deal with this underflow...
  if(t1+t2< -1e4){                       #This should be small enough!
    return(-1e4-t)
  }
  return(t1+t2)
}

########################################################################################################
## Compute the approximate probability that a recombination when the tmrca is tau
## changes the number of lineages when there are n lineages
########################################################################################################

pn <- function(tau, n){
  inner <- function(u, t, n){
    return(exp(-(2+wedding.cake(u, n-2))*(u-t)))
  }
  inner.vec <- Vectorize(inner, vectorize.args="u")

  outer <- function(t, tau, n){
    return( integrate(inner.vec, lower=t, upper=tau, n=n, t=t)$value )
  }
  outer.vec <- Vectorize(outer, vectorize.args="t")

  value <- integrate(outer.vec, lower=0, upper=tau, tau=tau, n=n)$value/tau
  return(1-value)
}

########################################################################################################
## Return approximated function that a recombination when the tmrca is tau
## changes the number of lineages when there are n lineages
########################################################################################################

make.pnfn <- function(n, t.grid=c((1:10)/10, 2:10)){
  res=rep(0, length(t.grid))
  for( i in 1:length(t.grid)){
    res[i] <- pn(t.grid[i], n)
  }
  return(approxfun(t.grid, res, rule=2))
}

########################################################################################################
## Compute the likelihood matrix of n observations at g grid points
## Lgs: Vector of genetic lengths
## Ds: Vector of numbers of doubletons
## Ne: Effective pop size
## pf: function returning probability that a recombination doesn't change the tmrca
## logt.grid: grid of values to compute likelihood at
## error.params: (shape, rate) for the gamma overestimate in Lg
## S.params - data.frame: S, theta, Lp, Ep: singletons, theta, phyical length, expected singletons
## verbose: report progress
########################################################################################################

compute.ll.matrix <- function(Lgs, Ds, Ne, pf, logt.grid=(0:60)/10, error.params=NA, S.params=NA, verbose=FALSE){
  ll.mat <- matrix(0, nrow=length(Lgs), ncol=length(logt.grid))
  
  for(i in 1:length(Lgs)){
    for(j in 1:length(logt.grid)){
      if(verbose){cat(paste("\r", i, j))}
      s.p <- NA
      if(!all(is.na(S.params))){s.p <- S.params[i,]}
      ll.mat[i,j] <- loglikelihood.age(t=(10^(logt.grid[j]))/2/Ne, Lg=Lgs[i], D=Ds[i], Ne=Ne, pf=pf,  error.params=error.params, S.params=s.p)
    }
  }
  if(verbose){cat("\n")}
  return(ll.mat)
}

########################################################################################################
## Estimate the density function, approximating as a discrete distribution. 
## Lgs: Vector of genetic lengths
## Ds: Vector of numbers of doubletons
## Ne: Effective pop size
## pf: function returning probability that a recombination doesn't change the tmrca
## logt.grid: grid of values to discretise at
## error.params: (shape, rate) for the gamma overestimate in Lg
## S.params - data.frame: S, theta, Lp, Ep: singletons, theta, phyical length, expected singletons
## verbose: show mcmc estimates
## prior: prior density
## ll.mat: use this matrix if precalculated
## n.sims, burn.in, thin: mcmc parameters
## alpha: precision of prior
########################################################################################################

estimate.t.density.mcmc <- function(Lgs, Ds, Ne, pf, logt.grid=(0:60)/10, verbose=FALSE, prior=function(x){return(1+0*x)}, ll.mat=NA, error.params=NA, S.params=NA, n.sims=1000, burn.in=100, thin=10, alpha=1){
  k <- length(logt.grid)
  
  if(all(is.na(ll.mat))){
    ll.mat <- compute.ll.matrix( Lgs, Ds, Ne, pf, logt.grid=logt.grid, error.params=error.params, S.params=S.params)
  }
  f.mat <- exp(ll.mat)
  for(i in 1:NROW(f.mat)){
    if(all(f.mat[i,]==0)){
      f.mat[i,] <- 1
    }
  }  

  prior.probs <- prior(logt.grid)
  prior.probs <- prior.probs/sum(prior.probs)

  all.counts <- mcmc.sample(f.mat, prior.probs, n.sims, burn.in, thin, alpha)

  all.probs <- all.counts/sum(all.counts)
  d.fun <- approxfun(logt.grid, all.probs)
  K <- integrate(d.fun, lower=min(logt.grid), upper=max(logt.grid))$value
  d.fun <- approxfun(logt.grid, all.probs/K)
  
  return(d.fun)
}

########################################################################################################
## Wrapper for the c function which does the mcmc sampling. If the compiled c sampler is available
## use that, otherwise fall through to the R sampler. The c sampler is defined in inference.c, which
## should be compiled with "R CMD SHLIB inference.c". The R version is of course *much* slower and
## also unsupported. Really only for small examples, or debugging. 
## 
## f.mat: matrix of likelihoods
## prior.probs: prior probabilities for each bin
## n_sims, burn_in, thin, alpha: mcmc sampling parameters
########################################################################################################

mcmc.sample <- function(f.mat, prior.probs, n.sims, burn.in, thin, alpha){
  if(is.loaded("mcmc_density_sampler")){
    dims <- as.integer(dim(f.mat))
    f.mat <- as.double(c(f.mat))
    prior <- as.double(prior.probs)
    mcmc.params <- as.integer(c(n.sims, burn.in, thin, alpha))
    seed <- as.integer(1234)
    all.counts <- as.integer(0*prior.probs)
    res <- .C("mcmc_density_sampler", f.mat, dims, prior, mcmc.params, seed, all.counts)
    all.counts <- res[[6]]
  } else {
    cat("Warning: Using unsupported R sampler. Results may be unreliable\n")
    warning("Using unsupported R sampler. Results may be unreliable")
    n <- NROW(f.mat)
    k <- NCOL(f.mat)
    
    all.counts <- 0*logt.grid                  #record all counts.
    probs <- prior.probs
    plot(logt.grid, probs, col="blue", type="l", ylim=c(0,2*max(probs)))
    total.sims <- n.sims+burn.in
    for(sim in 1:total.sims){
      counts <- 0*logt.grid
      for(i in 1:n){
        pr=probs*f.mat[i,]
        ti <- sample(1:k, size=1, prob=pr)
        counts[ti] <- counts[ti]+1
      }
      probs <- counts/(alpha+n)+alpha*prior.probs/(alpha+n)
      lines(logt.grid, probs, col=rgb(sim/total.sims, 0, 1-sim/total.sims, 10*thin/n.sims))
      if(sim>burn.in & !((sim-burn.in)%%thin)){all.counts <- all.counts+counts}
    }
  }
  return(all.counts)
}


########################################################################################################
## A_n(t), the expected number of lineges remaining in the coalescent at time t
## given that there are n lineages at time 0
########################################################################################################

wedding.cake <- function(t, n){
  if(n==0){return(0)}
  
  i <- 1:n
  terms <- exp(-t*i*(i-1)/2)*sapply(i, wedding.cake.ratio, n=n)
  return(sum(terms))
}

########################################################################################################
## Helper function for wedding.cake
## Note, if you need to call this a lot, tabulate it. 
########################################################################################################

wedding.cake.ratio <- function(i,n){
  return(exp(log(2*i-1)+sum(log((n-i+1):n))-sum(log(n:(n+i-1)))))
}

########################################################################################################
## just a weak prior for the log of the coalescent time
########################################################################################################

norm.2.p <- function(x){return(dnorm(x, mean=2))}
