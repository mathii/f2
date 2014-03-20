## Compute confidence intervals for a whole-genome simulation.

######################################################################################################

library(RColorBrewer)
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==3){
  chr.res.dir <- args[1]
  res.dir <- args[2]
  code.dir <- args[3]
} else{
  stop("Need to specify 3 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
chrs <- 1:22

######################################################################################################

t.hats <- list()
t.true <- list()
cis <- list()
subenv <- new.env()

## load saved results
cat("Loading data\n")
i=1
for(chr in chrs){
  cat(paste("\r", chr))

  load(paste(chr.res.dir, "/chr", chr, "/results/ll_and_density_environment.Rdata", sep=""), envir=subenv)
  t.hats[[i]] <- subenv$t.hats

  n.haps <- NROW(subenv$matched)
  this.ci <- matrix(0, nrow=n.haps, ncol=2)
  t.hats[[i]] <- MLE.from.haps(subenv$matched, subenv$Ne,S.params=subenv$S.params,  error.params=subenv$error.params, verbose=TRUE, tol=1e-10)
  for(j in 1:n.haps){
      cat(paste0("\r",chr, ":", j, "/", n.haps))
      this.ci[j,] <- confidence.interval( Lg=subenv$matched$map.len[j],  Ne=subenv$Ne, pf=subenv$p.fun, D=subenv$matched$f2[j], error.params=subenv$error.params, S.params=subenv$S.params[j,], alpha=0.95, max.search=10)
  }
  cis[[i]] <- this.ci
  t.true[[i]] <- subenv$matched$Age
  i=i+1
}

t.hats <- do.call("c", t.hats)
t.true <- do.call("c", t.true)
cis <- do.call("rbind", cis)

save.image(paste0(res.dir, "/all_ci.Rdata"))

plot( t.true , cis[,2], pch=".", log="xy", col="red", xlim=c(1,50000), ylim=c(1,50000))
abline(0,1,lty=2)
points( t.true, cis[,1], pch=".",  col="blue")
points( t.true, t.hats, pch=".",  col="grey")
ql <- qqplot(t.true, cis[,1], plot.it=FALSE)
lines(ql, col="blue")
qu <- qqplot( t.true, cis[,2], plot.it=FALSE)
lines(qu, col="red")
qt <- qqplot( t.true, t.hats, plot.it=FALSE)
lines(qt, col="grey")
