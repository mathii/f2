## Estimate the overestimate in length of the genetic length by sampling haplotypes at random. 

######################################################################################################

library(RColorBrewer)                   
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==10){
  code.dir <- args[1]
  hap.dir <- args[2]
  res.dir <- args[3]
  sample.file <- args[4]
  map.file <- args[5]
  pairs <- as.numeric(args[6])
  each <- as.numeric(args[7])
  direction <- args[8]                  #either "one.way" or "two.way"
  nbp <- as.numeric(args[9])                        #length of sequence.
  setup.file <- args[10]
} else{
  stop("Need to specify 10 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
source(setup.file)
samples <- scan(sample.file, quiet=TRUE, what="")

######################################################################################################

## Length overestimate parameters
pops <- unique(pop.map)
npop <- length(unique(pop.map))
shapes <- matrix(0, npop, npop)
rates <- matrix(0, npop, npop)

for(i in 1:npop){
  for(j in i:npop){
    cat(paste0(pops[i], "-", pops[j]))
    samples1 <- samples[pop.map==pops[i]]
    samples2 <- samples[pop.map==pops[j]]

    error.params <- fit.gamma.to.error(paste(hap.dir, "by_sample", sep="/"), map.file, samples1, samples2, pairs=pairs, each=each, verbose=TRUE, direction=direction)
    shapes[i,j] <- shapes[j,i] <- error.params[1]
    rates[i,j] <- rates[j,i] <- error.params[2]
  }
}
    write.table(shapes, paste(res.dir, "error_params_bypop_shapes.txt", sep="/"), sep="\t", col.names=FALSE, row.names=FALSE)
    write.table(rates, paste(res.dir, "error_params_bypop_rates.txt", sep="/"), sep="\t", col.names=FALSE, row.names=FALSE)

## cat("\nEstimating theta\n")
## ## Estimate theta per sample
## N <- length(samples)
## singletons <- count.singletons.from.positions(cbind(1:N, 1:N), rep(0, length(samples)), rep(nbp, length(samples)), paste(hap.dir, "pos.idx.f1.gz", sep="/") )
## thetas <- singletons/2/nbp              #expected number of singletons per base, per chromosome. 
## write.table(thetas, paste(res.dir, "theta_estimates.txt", sep="/"), sep="\t", col.names=FALSE, row.names=FALSE)
