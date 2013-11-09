## compare msmsc gene flow estimates to the between population
## f2 density estimates for CEU and CHB populations.

sim.type="CEU_CHB" 

msmc <- read.table("~/f2_age/1000g/CHB_CEU_4combined_0,1,4,5_msmc.final.txt", header=T, as.is=T)

######################################################################################################

## source(paste("~/f2_age/code/libs/include.R", sep=""))
dyn.load("~/f2_age/code/libs/inference.so")
code.dir="~/f2_age/code/"
chr.res.dir="~/f2_age/1000g/results/"
source(paste("~/f2_age/code/analysis/1kg_setup.R", sep=""))
source(paste0(code.dir, "/libs/include.R"))
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
plot(0.03*msmc$left_time_boundary/1.25e-8, 2*msmc$lambda_01/(msmc$lambda_00+msmc$lambda_11), type="s", log="x", ylim=c(0,1), bty="n", xlab="Time (kya)", ylab="MSMC gene flow estimate", col="#E41A1C", yaxt="n", lwd=2, xaxt="n")
axis(2, col="#E41A1C")
axis(1, at=c(seq(5,10,1), seq(20,100,10), seq(200,500,100)), labels=FALSE)
mtext(c(5,10,20,50,100,200,500), 1,at= c(5,10,20,50,100,200,500), line=1)
scale <- max(d.b(logt.grid))
lines(0.03*10^logt.grid, d.b(logt.grid)/scale, col="#377EBA", lwd=2)
axis(4, col="#377EBA", labels=FALSE, lwd=2, at=c(0,1))
mtext(expression(f[2]~age~density), 4, line=2)
abline(v=40, lty=3, col="#E41A1C")
abline(v=11, lty=3, col="#377EBA")

