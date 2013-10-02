## Plot the median by chromosome. See if they vary...
code.dir <- "~/f2/code/"
path.to.results <- "~/f2/results/"
source(paste0(code.dir, "analysis/1kgsetup.R"))

within <- matrix(0,14,22)

for(chr in 1:22){
  data <- read.table(paste0(path.to.results, "chr", chr, "/results/q50.txt"))
  for(i in 1:14){
    within[i,chr] <- data[i,i]
  }
}
rownames(within) <- rownames(data)

plot(within[1,], col=pop.cols[rownames(within)[1]], pch=16, bty="n", xlab="Chromosome", ylab="median", ylim=c(0,300), xaxt="n" )
axis(1, at=1:22)
for(i in 2:14){
  points(within[i,], col=pop.cols[rownames(within)[i]], pch=16)
}

alls <- read.table(paste0(path.to.results, "all/q50.txt"))
for(i in 1:14){
  abline(h=alls[i,i], col=pop.cols[rownames(alls)[i]])        
}

## formal test
med=c(within)
chr=as.factor(rep(1:22, each=14))
pop=as.factor(rep(rownames(within), 22))
