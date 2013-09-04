## Estimate the overestimate in length of the genetic length by sampling haplotypes at random. 

######################################################################################################

library(RColorBrewer)                   
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==9){
  code.dir <- args[1]
  hap.dir <- args[2]
  res.dir <- args[3]
  sample.file <- args[4]
  map.file <- args[5]
  pairs <- as.numeric(args[6])
  each <- as.numeric(args[7])
  direction <- args[8]                  #either "one.way" or "two.way"
  nbp <- as.numeric(args[9])                        #length of sequence. 
} else{
  stop("Need to specify 9 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
samples <- scan(sample.file, quiet=TRUE, what="")

######################################################################################################

## Length overestimate parameters 
error.params <- fit.gamma.to.error(paste(hap.dir, "by_sample", sep="/"), samples , map.file, pairs=pairs, each=each, verbose=TRUE, direction=direction)
write.table(error.params, paste(res.dir, "error_params.txt", sep="/"), sep="\t", col.names=FALSE, row.names=FALSE)

cat("Estimating theta\n")
## Estimate theta per sample
N <- length(samples)
singletons <- count.singletons.from.positions(cbind(1:N, 1:N), rep(0, length(samples)), rep(nbp, length(samples)), paste(hap.dir, "pos.idx.f1.gz", sep="/") )
thetas <- singletons/2/nbp              #expected number of singletons per base, per chromosome. 
write.table(thetas, paste(res.dir, "theta_estimates.txt", sep="/"), sep="\t", col.names=FALSE, row.names=FALSE)
