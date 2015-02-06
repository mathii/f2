## Analyse results from data for a single chromosome. Should be fairly easy to modify to run
## on general data. 

######################################################################################################

library(RColorBrewer)
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==8){
  code.dir <- args[1]
  res.dir <- args[2]
  Ne <- as.numeric(args[3])
  nseq <- as.numeric(args[4])
  mu <- as.numeric(args[5])
  max.log <- as.numeric(args[6])
  bins <- as.numeric(args[7])
  setup.file <- args[8]
  plots <- TRUE
} else{
  stop("Need to specify 8 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
source(setup.file)
haps <- read.table(paste(res.dir, "/f2_haplotypes.txt.gz", sep=""), as.is=TRUE, header=TRUE)
error.params <- scan(paste(res.dir, "error_params.txt", sep="/"), quiet=TRUE)
theta.estimates <- scan(paste(res.dir, "theta_estimates.txt", sep="/"), quiet=TRUE)

######################################################################################################

p.fun<-function(t){return(1)}

## Shouldn't need this, but for backwards compatability.
if("ibd.len" %in% names(haps)){names(haps)[which(names(haps)=="ibd.len")] <- "hap.len"}
haps <- haps[haps$map.len>0,]

## Singleton parameters
S.params <- haps[,c("f1", "hap.len")]
names(S.params) <- c("S", "Lp")
S.params$theta <- 4*Ne*mu
## These theta estimates are theta/2n
S.params$Ep <- S.params$Lp*(theta.estimates[haps$ID1]+theta.estimates[haps$ID2])
t.hats <- MLE.from.haps(haps, Ne, S.params=S.params,  error.params=error.params, verbose=TRUE)

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
    if(sum(include)>1){
      dens <- density(log10(t.hats[include]))
      densities[[i]][[j]] <- densities[[j]][[i]] <- approxfun(dens, rule=2)
    }else{
      densities[[i]][[j]] <- densities[[j]][[i]] <- function(x){return(0*x)}
    }
  }
}

## plots. One plot of all within-group densities, and one of all densities in total.
density.summary.plots(densities, populations, pop.cols, res.dir, max.log=max.log, xlim=c(1,5), ylim=c(0,2) )
