library(RColorBrewer)
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==3){
  chr.res.dir <- args[1]
  res.dir <- args[2]
  code.dir <- args[3]
  plots <- TRUE
} else{
  stop("Need to specify 3 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
chrs <- 1:22

######################################################################################################

ll.mat <- list()
subenv <- new.env()
logt.grid=NULL

## load saved results
cat("Loading data\n")
i=1
for(chr in chrs){
  cat(paste("\r", chr))
  load(paste(chr.res.dir, "/chr", chr, "/results/ll_environment.Rdata", sep=""), envir=subenv)
  ll.mat[[i]] <- subenv$ll.mat

  if(i==1){
    logt.grid <- subenv$logt.grid
  }else if (!all(subenv$logt.grid==logt.grid)){
    stop("grids do not match")
  }
  i=i+1
}
cat("\n")

ll.mat <- do.call("rbind", ll.mat)

## the following is cnp'd from run_1kg_analysis_chr.R. Put it in a function. 

## estimate densities, by population.
populations <- sort(unique(pop.map))
npop <- length(populations)             #14
densities <- rep(list(list()),npop)

ID1.pop <- pop.map[haps$ID1]
ID2.pop <- pop.map[haps$ID2]
for(i in 1:(npop)){
  for(j in i:npop){
    include <- (ID1.pop==populations[i]&ID2.pop==populations[j])|(ID1.pop==populations[j]&ID2.pop==populations[i])
    alpha <- round(0.05*sum(include))
    dens <- estimate.t.density.mcmc(0 ,0, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mat[include,])
    densities[[i]][[j]] <- dens
    densities[[j]][[i]] <- dens    
  }
}

## plots. One plot of all within-group densities, and one of all densities in total.
within.list <- list()
for(i in 1:npop){
  within.list[[i]] <- densities[[i]][[i]]
}
names(within.list) <- populations
if(plots){pdf(paste(res.dir, "/within.pdf", sep=""))}else{dev.new()}
plot.densities(within.list, logt.grid, cols=pop.cols[populations], xlim=c(1,5), ylim=c(0,2.5), main="Within")
if(plots){dev.off()}

for(i in 1:npop){
  between.list <- list()
  for(j in 1:npop){
    between.list[[j]] <- densities[[i]][[j]]
  }
  names(between.list) <- populations

  if(plots){pdf(paste(res.dir, "/between_", populations[i], ".pdf", sep=""))}else{dev.new()}
  plot.densities(between.list, logt.grid, cols=pop.cols[populations], xlim=c(1,5), ylim=c(0,2.5), main=populations[i])
  if(plots){dev.off()}
}

medians <- matrix(0, npop, npop)
for(i in 1:npop){
  for(j in 1:npop){
    medians[i,j] <- 10^quantile.density(densities[[i]][[j]], 0.5)
  }
}
colnames(medians) <- rownames(medians) <- populations
write.table(medians, paste(res.dir, "/medians.txt", sep=""), row.names=TRUE, col.names=TRUE, sep="\t")
