## compare msmsc gene flow estimates to the between population
## f2 density estimates for CEU and CHB populations.

sim.type="CHB_MXLNAT" 

msmc <- read.table("~/f2_age/1000g/CHB_MXLNAT_4combined_0,1,4,5_msmc.final.txt", header=T, as.is=T)
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
hom.tracts <- list()
for(sample in names(pop.map[pop.map=="MXL"])){
  data <- read.table(paste0("/Users/mathii/f2_age/1000g/MXL/", sample, ".bed"), header=FALSE, as.is=TRUE)
  hom.tracts[[sample]] <- data[data[,4]==6,]
}

subenv <- new.env()
t.hats <- list()

## load saved results
cat("Loading data\n")
i=1
for(chr in chrs){
  cat(paste("\r", chr))
  load(paste(chr.res.dir, "/chr", chr, "/results/ll_environment.Rdata", sep=""), envir=subenv)
  

  haps <- subenv$haps
  haps$chr <- chr                         #duh

  ID1.pop <- pop.map[haps$ID1]
  ID2.pop <- pop.map[haps$ID2]
  include <- (ID1.pop=="CHB"&ID2.pop=="MXL")|(ID2.pop=="CHB"&ID1.pop=="MXL")
  
  which.MXL <- ifelse(ID1.pop=="MXL", haps$ID1, haps$ID2)
  for(i in 1:NROW(haps)){
    if(!include[i]){next}
    cat(paste0("\r", i))
    this.MXL <- names(pop.map)[which.MXL[i]]
    my.tracts <- hom.tracts[[this.MXL]]
    
    n.tracts <- sum((haps$chr[i]==my.tracts[,1]) & (haps$hap.left[i]>my.tracts[,2]) & (haps$hap.right[i]<my.tracts[,3]))

    if(!n.tracts){include[i] <- FALSE}
  }

  t.hats[[i]]<-MLE.from.haps(haps[include,], subenv$Ne, S.params=subenv$S.params[include,], error.params=subenv$error.params, verbose=TRUE)

  i=i+1
}
cat("\n")
  
t.hats <- do.call("c", t.hats)
t.hat.dens <- density(log10(t.hats))
par(mar=c(5.1,4.1,2.1,4.1))

plot(0.03*msmc$left_time_boundary/1.25e-8, msmc$flow, type="s", log="x", ylim=c(0,1), bty="n", xlab="Time (kya)", ylab="MSMC gene flow estimate", col="#E41A1C", yaxt="n", lwd=2, xaxt="n", xlim=c(1,500))
axis(2, col="#E41A1C")
axis(1, at=c(seq(1,10,1), seq(20,100,10), seq(200,500,100)), labels=FALSE)
mtext(c(1,2,5,10,20,50,100,200,500), 1,at= c(1,2,5,10,20,50,100,200,500), line=1)
scale <- max(t.hat.dens$y)
lines(0.03*10^t.hat.dens$x, t.hat.dens$y/scale, col="#377EBA", lwd=2)
axis(4, col="#377EBA", labels=FALSE, lwd=2, at=c(0,1))
mtext(expression(f[2]~age~density), 4, line=2)
abline(v=0.03*msmc$left_time_boundary[max(which(msmc$flow<0.5))+1]/1.25e-8, lty=3, col="#E41A1C")
abline(v=median(0.03*t.hats), lty=3, col="#377EBA")

