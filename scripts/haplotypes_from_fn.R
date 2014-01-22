## find the haplotypes using the fn variants.
## args should be 1) the directory with the haplotypes in (f1, f2 and by sample) 2) The output file and 3) The map file/rate
## 4) Singleton power, 5) mutation rate
## 6) which n.
##
## Right now we are just using the lengths

args <- commandArgs(TRUE)

if(length(args)==7){
  code.dir <- args[1]
  hap.root <- args[2]
  res.dir <- args[3]
  map.file <- args[4]
  verbose <- as.numeric(args[5])
  n <- as.numeric(args[6])
  Ne <- as.numeric(args[7])
  
} else{
  stop("Need to specify 7 arguments")
}

source(paste(code.dir, "/libs/include.R", sep=""))

fn.file <- paste0(hap.root, "/pos.idx.f", n, ".gz")

pos.file <- paste(hap.root, "/by_sample/pos.gz", sep="")
by.sample.gt.root <- paste(hap.root, "/by_sample", sep="")
samples <- scan(paste(hap.root, "/samples.txt", sep=""), quiet=TRUE, what="")
error.params <- scan(paste(res.dir, "error_params.txt", sep="/"), quiet=TRUE)

sim.pop.map <- rep("ALL", length(samples))
names(sim.pop.map) <- samples

haps <- find.haplotypes.from.fn(fn.file, pos.file, by.sample.gt.root, sim.pop.map, "ALL", map.file, verbose=verbose)

## Just remove the ones with zero genetic length.
cat(paste0("Removed ", sum(haps$map.len==0), " haplotypes with length 0\n"))
haps <- haps[haps$map.len>0,]

## from the coalescent_fns.R file - need to compute for higher n
qns <- list(function(t){1}, q2, q3, q4, q5)
## estimate ages
t.hats <- MLE.from.haps(haps, Ne, S.params=NA,  error.params=error.params, verbose=TRUE, p.fun=qns[[n]])

haps <- cbind(haps, t.hat=t.hats)

write.table(haps, paste0(res.dir, "/f", n, "_results.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE)
