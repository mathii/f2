## plot the distibution of the MLE for different values of the parameter K.

######################################################################################################

library(RColorBrewer)
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==3){
  chr.res.dir <- args[1]
  res.dir <- args[2]
  code.dir <- args[3]
  plots <- TRUE
  publication.plots <- FALSE
  include.densities <- FALSE
} else{
  stop("Need to specify 3 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
chrs <- 1:22
Ks <- c(1,1.5,2.0)

######################################################################################################

all.t.hats <- list()
all.matched.haps <- list()

subenv <- new.env()
logt.grid=NULL

## load saved results
cat("Loading data\n")
i=1
for(chr in chrs){
  cat(paste("\r", chr))
  
  load(paste(chr.res.dir, "/chr", chr, "/results/ll_and_density_environment.Rdata", sep=""), envir=subenv)

  this.t.hats <- matrix(NA, nrow=NROW(subenv$matched), ncol=length(Ks))
  for(j in 1:length(Ks)){
      this.t.hats[,j]<-MLE.from.haps(subenv$matched, subenv$Ne, S.params=subenv$S.params, error.params=subenv$error.params, shape=Ks[j], verbose=TRUE)
  }
  all.t.hats[[i]] <- this.t.hats
  all.matched.haps[[i]] <- subenv$matched

  i=i+1
}
cat("\n")

t.hats <- do.call("rbind", all.t.hats)
matched <- do.call("rbind", all.matched.haps)

rm(all.t.hats)
rm(all.matched.haps)

pdf("~/Dropbox/f2_figures/SFigure10.pdf")
cols <- brewer.pal(length(Ks), "Set1")
qqplot(matched$Age, t.hats[,1], col=cols[1], type="l", lwd=2, bty="n", log="xy", xlab="True age (generations)", ylab="Estimated age (generations)")
for(i in 1:length(Ks)){
    qq <- qqplot(matched$Age, t.hats[,i], plot.it=FALSE)
    lines(qq, col=cols[i], lwd=2)
}
legend("topleft", paste0("K=", Ks), lwd=2, bty="n", col=cols)
abline(0,1, lty=2, col="black")
dev.off()
