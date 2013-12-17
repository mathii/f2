## find the haplotypes using the f2 variants.
## args should be 1) the directory with the haplotypes in (f1, f2 and by sample) 2) The output file and 3) The map file/rate
## 4) Singleton power, 5) mutation rate

args <- commandArgs(TRUE)

if(length(args)==5){
  code.dir <- args[1]
  hap.root <- args[2]
  out.file <- args[3]
  map.file <- args[4]
  verbose <- as.numeric(args[5])
} else{
  stop("Need to specify 4 arguments")
}

source(paste(code.dir, "/libs/include.R", sep=""))

f1.file <- paste(hap.root, "/pos.idx.f1.gz", sep="")
f2.file <- paste(hap.root, "/pos.idx.f2.gz", sep="")
pos.file <- paste(hap.root, "/by_sample/pos.gz", sep="")
by.sample.gt.root <- paste(hap.root, "/by_sample", sep="")
samples <- scan(paste(hap.root, "/samples.txt", sep=""), quiet=TRUE, what="")
sim.pop.map <- rep("ALL", length(samples))
names(sim.pop.map) <- samples

haps <- find.haplotypes.from.f2(f1.file, f2.file, pos.file, by.sample.gt.root, sim.pop.map, "ALL", map.file, verbose=verbose)

## Just remove the ones with zero length.
cat(paste0("Removed ", sum(haps$map.len>0), " haplotypes with length 0\n"))
haps <- haps[haps$map.len>0,]

write.table(haps, out.file, row.names=FALSE, col.names=TRUE, quote=FALSE)
