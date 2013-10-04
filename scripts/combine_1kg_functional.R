## Script to combine the analsis of different functional classes.
## Just by functional class, and within vs between. 

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
source(paste(code.dir, "/analysis/1kgsetup.R", sep=""))
chrs <- c(1:22)

######################################################################################################

subenv <- new.env()
logt.grid=NULL

ll.mats.by.class <- list()
haps.by.class <- list()

classes <- c("lof", "coding", "noncoding")

for( cls in classes){

  haps <- list()
  ll.mat <- list()

  ## load saved results
  cat("Loading data\n")
  i=1
  for(chr in chrs){
    cat(paste("\r", chr))
    load(paste0(chr.res.dir, "/chr", chr, "/", cls, "/results/ll_environment.Rdata"), envir=subenv)
    ll.mat[[i]] <- subenv$ll.mat

    haps[[i]] <- subenv$haps  
    if(i==1){
      logt.grid <- subenv$logt.grid
    }else if (!all(subenv$logt.grid==logt.grid)){
      stop("grids do not match")
    }
    i=i+1
  }
  cat("\n")

  ll.mats.by.class[[cls]] <- do.call("rbind", ll.mat)
  haps.by.class[[cls]] <- do.call("rbind", haps)
}


## the following is cnp'd from run_1kg_analysis_chr.R.
## estimate densities, by class for within/between.
densities <- rep(list(list()),3)

i=1
for(cls in classes){
  haps <- haps.by.class[[cls]]
  ll.mat <- ll.mats.by.class[[cls]]
  
  ID1.pop <- pop.map[haps$ID1]
  ID2.pop <- pop.map[haps$ID2]

  within <- ID1.pop==ID2.pop
  between <- !within

  alpha.w <- round(0.05*sum(within))
  alpha.b <- round(0.05*sum(between))
  
  dens.w <- estimate.t.density.mcmc(0 ,0, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha.w,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mat[within,])
  dens.b <- estimate.t.density.mcmc(0 ,0, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha.b,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mat[between,])
  
  densities[[i]][["within"]] <- dens.w
  densities[[i]][["between"]] <- dens.b
  i=i+1
}

q5 <- q50 <- q95 <- matrix(0,nrow=2,ncol=3)

for( i in 1:3){
  for(j in 1:2){
    a <- quantile.density(densities[[i]][[j]], 0.05)
    b <- quantile.density(densities[[i]][[j]], 0.5)
    c <- quantile.density(densities[[i]][[j]], 0.95)
    q5[j,i] <- a
    q50[j,i] <- b
    q95[j,i] <- c    
  }
}

## Just plot the medians
pdf(paste0(res.dir, "/by_func.pdf"))
xs <- c(1,2,3,4,5,6)
ys <- c(q50[1,], q50[2,])
plot(xs, ys, pch=16, col=c("blue", "blue", "blue", "red", "red", "red"), bty="n", xaxt="n", ylab=expression(Age (Log[10]~generations)), ylim=c(0,4), xlab="")
segments(1:3,q5[1,],1:3,q95[1,], col="blue")
segments(4:6,q5[2,],4:6,q95[2,], col="red")
axis(1, at=c(1:6), labels=c("LOF", "coding", "noncoding", "LOF", "coding", "noncoding"))
legend("topleft", c("Within", "Between"), col=c("blue", "red"), pch=16, bty="n")
abline(v=3.5, lty=3)
dev.off()
