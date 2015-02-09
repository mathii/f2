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
