## Analyse results from 1kg for a single chromosome. Should be fairly easy to modify to run
## on general data. 

######################################################################################################

library(RColorBrewer)
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==7){
  code.dir <- args[1]
  res.dir <- args[2]
  Ne <- as.numeric(args[3])
  nseq <- as.numeric(args[4])
  mu <- as.numeric(args[5])
  max.log <- as.numeric(args[6])
  bins <- as.numeric(args[7])
  plots <- TRUE
} else{
  stop("Need to specify 7 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
source(paste(code.dir, "/analysis/1kgsetup.R", sep=""))
haps <- read.table(paste(res.dir, "/f2_haplotypes.txt.gz", sep=""), as.is=TRUE, header=TRUE)
error.params <- scan(paste(res.dir, "error_params.txt", sep="/"), quiet=TRUE)
theta.estimates <- scan(paste(res.dir, "theta_estimates.txt", sep="/"), quiet=TRUE)

######################################################################################################

#if(!("p.fun" %in% ls())){
#  cat("Making pn function\n")
#  p.fun <- make.pnfn(nseq*2)              #probably need to speed this up.
#  ## p.fun <- make.pnfn(250)              #probably need to speed this up. 
#}

p.fun<-function(t){return(1)}

haps$hap.len <- haps$hap.right-haps$hap.left
haps <- haps[haps$map.len>0,]

## Singleton parameters
S.params <- haps[,c("f1", "hap.len")]
names(S.params) <- c("S", "Lp")
S.params$theta <- 4*Ne*mu
S.params$Ep <- S.params$Lp*(theta.estimates[haps$ID1]+theta.estimates[haps$ID2])

## Compute likelihood
cat("Computing likelihood\n")
logt.grid <- seq(0, max.log, length.out=bins)
ll.mat <- compute.ll.matrix( haps$map.len, haps$f2, Ne, p.fun, logt.grid=logt.grid, error.params=error.params, S.params=S.params, verbose=FALSE)

save.image(paste(res.dir, "/ll_environment.Rdata", sep=""))

## estimate densities, by population.
populations <- sort(unique(pop.map))
npop <- length(populations)             #14
densities <- rep(list(list()),npop)

ID1.pop <- pop.map[haps$ID1]
ID2.pop <- pop.map[haps$ID2]
for(i in 1:(npop)){
  for(j in i:npop){
    include <- (ID1.pop==populations[i]&ID2.pop==populations[j])|(ID1.pop==populations[j]&ID2.pop==populations[i])
    if(sum(include)>0){
      alpha <- round(0.05*sum(include))
      dens <- estimate.t.density.mcmc(0 ,0, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mat[include,])
      densities[[i]][[j]] <- densities[[j]][[i]] <- dens
    }else{
      densities[[i]][[j]] <- densities[[j]][[i]] <- function(x){return(0)}
    }
  }
}

## plots. One plot of all within-group densities, and one of all densities in total.
density.summary.plots(densities, populations, pop.cols, res.dir, xlim=c(1,4), ylim=c(0,3) )
