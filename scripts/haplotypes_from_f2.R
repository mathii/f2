## find the haplotypes using the f2 variants.
## args should be 1) the directory with the haplotypes in (f1, f2 and by sample) 2) The output file and 3) The map file/rate
## 4) Singleton power, 5) mutation rate

args <- commandArgs(TRUE)

if(length(args)==4){
  code.dir <- args[1]
  hap.root <- args[2]
  out.file <- args[3]
  map.file <- args[4]
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

all.results <- find.haplotypes.from.f2(f1.file, f2.file, pos.file, by.sample.gt.root, sim.pop.map, "ALL", map.file)

write.table(all.results, out.file, row.names=FALSE, col.names=TRUE, quote=FALSE)
