## Script to combine the analsis of different functional classes.
## Just by functional class, and within vs between. 

library(RColorBrewer)
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==4){
  chr.res.dir <- args[1]
  res.dir <- args[2]
  code.dir <- args[3]
  setup.file <- args[4]
  plots <- FALSE
} else{
  stop("Need to specify 4 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
source(setup.file)
chrs <- c(1:22)

######################################################################################################

subenv <- new.env()
logt.grid=NULL

ll.mats.by.class <- list()
haps.by.class <- list()

classes <- c("lof", "coding", "noncoding")
ns <- rep(0,2*length(classes))

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
## densities <- list()

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
  
  dens.w <- estimate.t.density.mcmc(0 ,0, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha.w,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mat[within,], return.density.only=FALSE)
  dens.b <- estimate.t.density.mcmc(0 ,0, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha.b,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mat[between,], return.density.only=FALSE)

  densities[[i]] <- dens.w$density
  densities[[i+length(classes)]] <- dens.b$density

  ns[i] <- sum(within)
  ns[i+length(classes)] <- sum(between)

  i=i+1
}

col=c(rep("#377EBA", 3), rep("#E41A1C", 3))
border=col
fill=paste0(col, "80")
x.pos=c(1.25,2.25,3.25,4.75,5.75,6.75)

viola.plot(densities, x.pos=x.pos, eps=2e-2, col=col, border=border, fill=fill, labels=rep(c("LOF", "Coding", "Noncoding"),2), xlim=c(1,7), ylab=expression(Age~(Log[10]~generations)) )
abline(v=4, lty=3)
mtext(paste0("(",format(ns, big.mark=",", trim=TRUE),")"), 1, at=x.pos, line=1)
mtext(c("Within population", "Between populations"), 3, at=c(2.25, 5.75), line=0)

