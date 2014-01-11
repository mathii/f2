## Make plots for additional scenarios - used for supplementary figure 5.

############################################################################################
## ancient split


load("~/f2_age/simulations/ancient_split/chr20/results/ll_and_density_environment.Rdata")
dyn.load("~/f2_age/code/libs/inference.so")
code.dir<-"~/f2_age/code/"
source("~/f2_age/code/libs/include.R")

within <- (matched$ID1<51 & matched$ID2<51)|(matched$ID1>50 & matched$ID2>50)
between <- !within

t.hats.b<-MLE.from.haps(matched[between,], Ne, S.params=S.params[between,], error.params=error.params, verbose=TRUE)
t.hat.b.dens <- density(log10(t.hats.b))
t.hats.w<-MLE.from.haps(matched[within,], Ne, S.params=S.params[within,], error.params=error.params, verbose=TRUE)
t.hat.w.dens <- density(log10(t.hats.w))

t.true.w.dens <- density(log10(matched$Age[within]))
t.true.b.dens <- density(log10(matched$Age[between]))

grid=10^seq(0,5,length.out=1000)

pdf("~/Desktop/ancient_split_density.pdf")
plot(10^t.hat.w.dens$x, t.hat.w.dens$y, col="#377EBA", type="l", bty="n",log="x", xlab="Age (generations)", ylab="Density", ylim=c(0,4), lwd=2, xlim=c(10,11000))
lines(10^t.hat.b.dens$x, t.hat.b.dens$y, col="#E41A1C", lwd=2)

lines(10^t.true.b.dens$x, t.true.b.dens$y, col="#E41A1C", lty=2, lwd=2)
lines(10^t.true.w.dens$x, t.true.w.dens$y, col="#377EBA", lty=2, lwd=2)

abline(v=1120, lty=3, lwd=2)
abline(v=median(t.hats.b), lty=3, col="#E41A1C", lwd=2)
legend("topleft", c("Within", "Between", "Split time"), col=c("#377EBA", "#E41A1C", "black"), lty=c(1,1,3), bty="n", lwd=2)

dev.off()

pdf("~/Desktop/ancient_split_compare.pdf")
plot.mle.and.density(matched$Age, t.hats[,6], dd, main="", xlim=c(0.5,4), ylim=c(0.5,4), col="#000000", alpha="80")
abline(v=log10(560), lty=3)
dev.off()

############################################################################################
## bottleneck split

load("~/f2_age/simulations/bottleneck/chr20/results/ll_and_density_environment.Rdata")
dyn.load("~/f2_age/code/libs/inference.so")
code.dir<-"~/f2_age/code/"
source("~/f2_age/code/libs/include.R")

t.hats<-MLE.from.haps(matched, Ne, S.params=S.params, error.params=error.params, verbose=TRUE)
t.hat.dens <- density(log10(t.hats))
t.true.dens <- density(log10(matched$Age), adjust=0.25)


grid=10^seq(0,5,length.out=1000)

pdf("~/Desktop/bottleneck_density.pdf")
plot(10^t.hat.dens$x, t.hat.dens$y, col="#377EBA", type="l", bty="n",log="x", xlab="Age (generations)", ylab="Density", lwd=2, xlim=c(1,10000), ylim=c(0,3))
lines(10^t.true.dens$x, t.true.dens$y, col="#377EBA", lty=2, lwd=2) 

abline(v=14, lty=3, lwd=2)
abline(v=16.8, lty=3, lwd=2)

dev.off()
