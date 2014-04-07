## compare msmsc gene flow estimates to the between population
## f2 density estimates for CEU and CHB populations.

sim.type="CEU_CHB" 

msmc <- read.table("~/f2_age/1000g/CHB_CEU_4combined_0,1,4,5_msmc.final.txt", header=T, as.is=T)
msmc$flow <- 2*msmc$lambda_01/(msmc$lambda_00+msmc$lambda_11)

######################################################################################################

## source(paste("~/f2_age/code/libs/include.R", sep=""))
dyn.load("~/f2_age/code/libs/inference.so")
code.dir="~/f2_age/code/"
chr.res.dir="~/f2_age/1000g/results/"
source(paste("~/f2_age/code/analysis/1kg_setup.R", sep=""))
source(paste0(code.dir, "/libs/include.R"))
chrs <- c(1:22)

######################################################################################################

if(TRUE ){
  t.hats <- list()
  subenv <- new.env()
  logt.grid=NULL

  ## load saved results
  cat("Loading data\n")
  i=1
  if(TRUE){                            #To load from chr
    for(chr in chrs){
      cat(paste("\r", chr))
      load(paste(chr.res.dir, "/chr", chr, "/results/ll_environment.Rdata", sep=""), envir=subenv)

      haps <- subenv$haps

      ID1.pop <- pop.map[haps$ID1]
      ID2.pop <- pop.map[haps$ID2]
      include <- (ID1.pop=="CEU"&ID2.pop=="CHB")|(ID2.pop=="CEU"&ID1.pop=="CHB")
      
      t.hats[[i]]<-MLE.from.haps(haps[include,], subenv$Ne, S.params=subenv$S.params[include,], error.params=subenv$error.params, verbose=TRUE)
    
      i=i+1
    }
    cat("\n")
    
    t.hats <- do.call("c", t.hats)
    t.hat.dens <- density(log10(t.hats))
  } else{                               #Load from new code
     load(paste(chr.res.dir, "/all/all_results.Rdata", sep=""))
     ID1.pop <- pop.map[haps$ID1]
     ID2.pop <- pop.map[haps$ID2]
     include <- (ID1.pop=="CEU"&ID2.pop=="CHB")|(ID2.pop=="CEU"&ID1.pop=="CHB")
     t.hats <- t.hats[include]
     haps <- haps[include,]
     t.hat.dens <- density(log10(t.hats))
  }
}

par(mar=c(5.1,4.1,2.1,4.1))
plot(0.03*msmc$left_time_boundary/1.25e-8, 2*msmc$lambda_01/(msmc$lambda_00+msmc$lambda_11), type="s", log="x", ylim=c(0,1), bty="n", xlab="Time (kya)", ylab="MSMC gene flow estimate", col="#E41A1C", yaxt="n", lwd=2, xaxt="n", xlim=c(1,500))
axis(2, col="#E41A1C")
axis(1, at=c(seq(1,10,1), seq(20,100,10), seq(200,500,100)), labels=FALSE)
mtext(c(1,2,5,10,20,50,100,200,500), 1,at= c(1,2,5,10,20,50,100,200,500), line=1)
scale <- max(t.hat.dens$y)
lines(0.03*10^t.hat.dens$x, t.hat.dens$y/scale, col="#377EBA", lwd=2)
axis(4, col="#377EBA", labels=FALSE, lwd=2, at=c(0,1))
mtext(expression(f[2]~age~density), 4, line=2)
abline(v=0.03*msmc$left_time_boundary[max(which(msmc$flow<0.5))+1]/1.25e-8, lty=3, col="#E41A1C")
abline(v=median(0.03*t.hats), lty=3, col="#377EBA")

