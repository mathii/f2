## Estimate the overestimate in length of the genetic length by sampling haplotypes at random. 

######################################################################################################

library(RColorBrewer)                   
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==8){
  code.dir <- args[1]
  hap.dir <- args[2]
  res.dir <- args[3]
  sample.file <- args[4]
  map.file <- args[5]
  pairs <- as.numeric(args[6])
  each <- as.numeric(args[7])
  direction <- args[8]                  #either "one.way" or "two.way"
} else{
  stop("Need to specify 8 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
samples <- scan(sample.file, quiet=TRUE, what="")

######################################################################################################

error.params <- fit.gamma.to.error(paste(hap.dir, "by_sample", sep="/"), samples , map.file, pairs=pairs, each=each, verbose=TRUE, direction=direction)
write.table(error.params, paste(res.dir, "error_params.txt", sep="/"), sep="\t", col.names=FALSE, row.names=FALSE)
