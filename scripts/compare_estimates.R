## compute and compare the age estimates. Need to generate the matched haplotypes
## using compare_haplotypes.R first, saved in res.dir/matched_haps.txt.gz
## Script compares results using genetic length with and without singletons, 
## and with and without corrections. Geneitc length corrections us a gamma
## distribution, which was estimated from simulated human data. 

######################################################################################################

library(RColorBrewer)
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==7){
  code.dir <- args[1]
  res.dir <- args[2]
  Ne <- as.numeric(args[3])
  nseq <- as.numeric(args[4])/2
  mu <- as.numeric(args[5])
  max.log <- as.numeric(args[6])
  bins <- as.numeric(args[7])
  plots <- TRUE
  verb <- FALSE
  publication.plots <- TRUE
} else{
  stop("Need to specify 7 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
error.params <- scan(paste(res.dir, "error_params.txt", sep="/"), quiet=TRUE)
theta.estimates <- scan(paste(res.dir, "theta_estimates.txt", sep="/"), quiet=TRUE)
######################################################################################################

p.fun <- function(t){return(1)}
logt.grid <- seq(0, max.log, length.out=bins)

######################################################################################################

matched <- NA
if(!file.exists( paste(res.dir, "/matched_haps.txt.gz", sep=""))){
  cat("Salvaging likelihood\n")
  matched <- read.table(paste(res.dir, "/f2_haplotypes.txt.gz", sep=""),  as.is=TRUE, header=TRUE)
  matched$hap.len <- matched$hap.right-matched$hap.left
  matched <- matched[matched$ID1!=matched$ID2,]
  matched <- matched[matched$map.len>0,]      
  ll.mats <- list()
  S.params <- matched[,c("f1", "hap.len")]
  names(S.params) <- c("S", "Lp")
  S.params$theta <- 4*Ne*mu
  S.params$Ep <- S.params$Lp*(theta.estimates[matched$ID1]+theta.estimates[matched$ID2])
  ll.mats[[6]] <- compute.ll.matrix( matched$map.len, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=error.params, S.params=S.params, verbose=verb)
  save.image(paste(res.dir, "/ll_and_density_environment.Rdata", sep=""))
}
matched <- read.table(paste(res.dir, "/matched_haps.txt.gz", sep=""), as.is=TRUE, header=TRUE)

######################################################################################################


## Set up singleton params
t.hats <- matrix(0, nrow=NROW(matched), ncol=6)
S.params <- matched[,c("f1", "hap.len")]
names(S.params) <- c("S", "Lp")
S.params$theta <- 4*Ne*mu
S.params$Ep <- S.params$Lp*(theta.estimates[matched$ID1]+theta.estimates[matched$ID2])

## Real true S and Lp, for comparison. 
S.params.true <- matched[,c("true.f1", "true.len")]
names(S.params.true) <- c("S", "Lp")
S.params.true$theta <- 4*Ne*mu
S.params.true$Ep <-  S.params.true$Lp*(theta.estimates[matched$ID1]+theta.estimates[matched$ID2])

cat("Calculating MLE\n")
for(j in 1:NROW(matched)){
  cat(paste("\r", j))
  ## Using true values
  max.search <- (10^max.log)/2/Ne
  t.hats[j,1] <- 2*Ne*optimize(loglikelihood.age, interval=c(0,max.search), maximum=TRUE,  Lg=matched$true.map[j]/100,  Ne=Ne, pf=p.fun, D=matched$f2[j], error.params=NA)$maximum
  t.hats[j,4] <- 2*Ne*optimize(loglikelihood.age, interval=c(0,max.search), maximum=TRUE,  Lg=matched$true.map[j]/100,  Ne=Ne, pf=p.fun, D=matched$f2[j], error.params=NA, S.params=S.params.true[j,])$maximum

  ## Using observed values
  t.hats[j,2] <- 2*Ne*optimize(loglikelihood.age, interval=c(0,max.search), maximum=TRUE,  Lg=matched$map.len[j],  Ne=Ne, pf=p.fun,  D=matched$f2[j], error.params=NA)$maximum
  t.hats[j,5] <- 2*Ne*optimize(loglikelihood.age, interval=c(0,max.search), maximum=TRUE,  Lg=matched$map.len[j],  Ne=Ne, pf=p.fun,  D=matched$f2[j], error.params=NA, S.params=S.params[j,])$maximum

  ## Obesrved values with corrected likelihood
  max.search <- 100*t.hats[j,2]/2/Ne
  if(log10(t.hats[j,2])<0){max.search <- max.search/50} #ugh
  t.hats[j,3] <- 2*Ne*optimize(loglikelihood.age, interval=c(0,max.search), maximum=TRUE,  Lg=matched$map.len[j],  Ne=Ne, pf=p.fun, D=matched$f2[j], error.params=error.params)$maximum
  t.hats[j,6] <- 2*Ne*optimize(loglikelihood.age, interval=c(0,max.search), maximum=TRUE,  Lg=matched$map.len[j],  Ne=Ne, pf=p.fun, D=matched$f2[j], error.params=error.params, S.params=S.params[j,])$maximum
}
cat("\rCalculating ll matrices\n")

norm.2.p <- function(x){return(dnorm(x, mean=2))}

denss <- rep(list(NA),6)
ll.mats <- list()
ll.mats[[1]] <- compute.ll.matrix( matched$true.map/100, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=NA, S.params=NA, verbose=verb)
ll.mats[[2]] <- compute.ll.matrix( matched$map.len, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=NA, S.params=NA, verbose=verb)
ll.mats[[3]] <- compute.ll.matrix( matched$map.len, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=error.params, S.params=NA, verbose=verb)
ll.mats[[4]] <- compute.ll.matrix( matched$true.map/100, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=NA, S.params=S.params.true, verbose=verb)
ll.mats[[5]] <- compute.ll.matrix( matched$map.len, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=NA, S.params=S.params, verbose=verb)
ll.mats[[6]] <- compute.ll.matrix( matched$map.len, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=error.params, S.params=S.params, verbose=verb)

cat("Sampling Densities\n")
alpha <- round(0.05*NROW(matched))
for(i in 1:6){
  ## Lg and D are dummy entries... the ll matrix is the only information we use... 
  denss[[i]] <- estimate.t.density.mcmc(0*matched$map.len ,0*matched$f2, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mats[[i]])
}

save.image(paste(res.dir, "/ll_and_density_environment.Rdata", sep=""))

labels=c("True (Lg)", "Observed (Lg)", "Corrected (Lg)", "True (Lg+S)", "Observed (Lg+S)", "Corrected (Lg+S)")
if(plots){pdf(paste(res.dir, "/compare_estimates.pdf", sep=""), height=12, width=18)}else{dev.new()}
par(mfrow=c(2,3))
for(i in 1:6){
  plot.mle.and.density(matched$Age, t.hats[,i], denss[[i]], main=labels[i], xlim=c(1,10^max.log), ylim=c(1,10^max.log))
}
if(plots){dev.off()}

if(publication.plots){
    for(i in 1:6){
        png(paste0(res.dir, "/estimate.", labels[i], ".png"), height=600, width=600)
        plot.mle.and.density(matched$Age, t.hats[,i], NA, xlim=c(1,10^max.log), ylim=c(1,10^max.log), cex=2, alpha="40", col="#101010", line.col="red")
        dev.off()
    }
}


if(plots){pdf(paste(res.dir, "/coverage.pdf", sep=""), height=6, width=6)}else{dev.new()}
plot.coverage.curves(matched$Age, t.hats, labels=labels, lty=rep(1,6), lwd=1, cols=c("black", brewer.pal(5, "Set1")))
if(plots){dev.off()}

## Now do confidence intervals.
cis <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
lower <- matrix(0, nrow=NROW(matched), ncol=length(cis))
upper <- matrix(0, nrow=NROW(matched), ncol=length(cis))
for(i in 1:length(cis)){
  for(j in 1:NROW(matched)){
    cat(paste0("\r",i, ":", j))
    ci <- confidence.interval( Lg=matched$map.len[j],  Ne=Ne, pf=p.fun, D=matched$f2[j], error.params=error.params, S.params=S.params[j,], alpha=cis[i], max.search=10)
    lower[j,i] <- ci[1]
    upper[j,i] <- ci[2]
  }
}

gt.lower <- 0*lower
lt.upper <- 0*upper

for(i in 1:length(cis)){
  gt.lower[,i] <- matched$Age>lower[,i]
  lt.upper[,i] <- matched$Age<upper[,i]

}
probs <- colMeans(gt.lower&lt.upper)

if(plots){pdf(paste(res.dir, "/CI_coverage.pdf", sep=""), height=6, width=6)}else{dev.new()}
plot(c(0,1), c(0,1), lwd=2, lty=2, type="l", col="grey80", bty="n", xlab=expression(alpha), ylab=expression(Empirical~alpha~CI))
lines(cis, probs, type="b", col="#e41a1c",lwd=2, pch=16)
if(plots){dev.off()}
