## Functions for the v2 analysis.

############################################################
##
## Derivative of the positive poission log likelihood with
## respect to mu
##
############################################################

d.pospoisson <- function(mu, y, x){
    N <- length(x)
    if(length(y)!=N){stop("x and y must be the same length")}
    if(length(mu)!=1){stop("mu must be length 1")}
    return(sum(y/mu - x - x*exp(-mu*x)/(1-exp(-mu*x))))
}

############################################################
##
## Derivative of the poission log likelihood with
## respect to mu
##
############################################################

d.poisson <- function(mu, y, x){
    N <- length(x)
    if(length(y)!=N){stop("x and y must be the same length")}
    if(length(mu)!=1){stop("mu must be length 1")}
    return(sum(y/mu-x))
}

############################################################
##
## Compute likelihhod profile of a bunch of haps with
## respect to mutation rate
##
############################################################

likelihood.profile.mu <- function(haps, Ne, mu.grid, error.params=NA,
                                  verbose=TRUE, tol=1/2/Ne,
                                  hap.detect.rate=1, shape=1.5,
                                  v2=TRUE, theta.estimates){

    med.t.hat <- 0*mu.grid
    ll <- matrix(NA, nrow=NROW(haps), ncol=length(mu.grid))

    likelihood.func <- loglikelihood.age
    if(v2){
        likelihood.func <- loglikelihood.age.v2
    }

    for(i in 1:length(mu.grid)){
        mu=mu.grid[i]
        cat("mu=", mu, "\n")

        S.params <- haps[,c("f1", "hap.len")]
        names(S.params) <- c("S", "Lp")
        S.params$theta <- 4*Ne*mu
        S.params$Ep <- S.params$Lp*(theta.estimates[haps$ID1]+theta.estimates[haps$ID2])

        t.hats <- MLE.from.haps(haps, Ne, S.params=S.params, error.params=error.params, verbose=verbose, tol=tol, mu=hap.detect.rate, ignoreErrors=TRUE, v2=v2, shape=shape)
        med.t.hat[i] <- median(t.hats, na.rm=TRUE)
        for(j in 1:NROW(haps)){
            e <- tryCatch({
                ## Remembe that the t.hats are in generations, but the loglikelihhod is in coalescent time
                ll[j,i] <- likelihood.func(t.hats[j]/2/Ne, haps$map.len[j], Ne,haps$f2[j], pf=function(x){1},error.params=error.params, S.params=S.params[j,], mu=hap.detect.rate, shape=shape)}
                        , error=function(e){e})
        }
    }

    exclude <- apply(is.na(ll), 1, any)
    cat(paste0("Excluded ", sum(exclude), " haplotypes"))
    ll <- ll[!exclude,]
    ll <- colSums(ll)
    
    return(list(x=mu.grid, y=ll, median.t.hat=med.t.hat))
}
