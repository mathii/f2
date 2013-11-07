## compare msmsc gene flow estimates to the between population
## f2 density estimates for CEU and CHB populations.

sim.type="CEU_CHB" 

msmc <- read.table("~/f2_age/1000g/CHB_CEU_4combined_0,1,4,5_msmc.final.txt", header=T, as.is=T)

######################################################################################################

## source(paste("~/f2_age/code/libs/include.R", sep=""))
dyn.load("~/f2_age/code/libs/inference.so")
code.dir="~/f2_age/code/"
chr.res.dir="~/f2_age/1000g/results/"
source(paste("~/f2_age/code/analysis/1kgsetup.R", sep=""))
chrs <- c(1:22)

######################################################################################################

if(! "d.b" %in% ls() ){
  haps <- list()
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
    
    haps[[i]] <- subenv$haps
    
    if(NROW(ll.mat[[i]])!=NROW(haps[[i]])){stop(paste0("Error in chr", chr))}
  
    if(i==1){
      logt.grid <- subenv$logt.grid
    }else if (!all(subenv$logt.grid==logt.grid)){
      stop("grids do not match")
    }
  i=i+1
  }
  cat("\n")
  
  ll.mat <- do.call("rbind", ll.mat)
  haps <- do.call("rbind", haps)
  
  ## the following is cnp'd from run_1kg_analysis_chr.R.
  ## estimate densities, by population.
  populations <- sort(unique(pop.map))
  npop <- length(populations)             #14
  densities <- rep(list(list()),npop)

  ID1.pop <- pop.map[haps$ID1]
  ID2.pop <- pop.map[haps$ID2]
  
  include <- (ID1.pop=="CEU"&ID2.pop=="CHB")|(ID2.pop=="CEU"&ID1.pop=="CHB")
  
  alpha <- round(0.05*sum(include))
  d.b <- estimate.t.density.mcmc(0 ,0, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mat[include,])
}

par(mar=c(5.1,4.1,2.1,4.1))
plot(msmc$left_time_boundary/1.25e-8, 2*msmc$lambda_01/(msmc$lambda_00+msmc$lambda_11), type="s", log="x", ylim=c(0,1), xlim=c(100,20000), bty="n", xlab="Generations", ylab="MSMC gene flow estimate", col="#E41A1C", yaxt="n", main="CEU_CHB")
  axis(2, col="#E41A1C")
  scale <- max(d.b(logt.grid))
lines(10^logt.grid, d.b(logt.grid)/scale, col="#377EBA")
axis(4, col="#377EBA", labels=FALSE)
mtext(expression(f[2]~age~density), 4, line=2)

## chose other features here
## abline(v=224, lty=3)
## abline(v=560*2, lty=3)
## abline(v=1120, lty=3)
## abline(v=560, lty=3)
