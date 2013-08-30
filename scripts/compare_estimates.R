## compute and compare the age estimates. Need to generate the matched haplotypes
## using compare_haplotypes.R first, saved in res.dir/matched_haps.txt.gz
## Script compares results using genetic length with and without singletons, 
## and with and without corrections. Geneitc length corrections us a gamma(1.2, 375M)
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
} else{
  stop("Need to specify 7 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
matched <- read.table(paste(res.dir, "/matched_haps.txt.gz", sep=""), as.is=TRUE, header=TRUE)
error.params <- scan(paste(res.dir, "error_params.txt", sep="/"), quiet=TRUE)

######################################################################################################

cat("Calculating pn(t)\n")
if(!("p.fun" %in% ls())){
  p.fun <- make.pnfn(2*nseq)
}
                    
## Set up singleton params
t.hats <- matrix(0, nrow=NROW(matched), ncol=6)
S.params <- matched[,c("f1", "hap.len")]
names(S.params) <- c("S", "Lp")
S.params$theta <- 4*Ne*mu
S.params$Ep <- S.params$Lp*4*Ne*mu/2/nseq

## Real true S and Lp, for comparison. 
S.params.true <- matched[,c("true.f1", "true.len")]
names(S.params.true) <- c("S", "Lp")
S.params.true$theta <- 4*Ne*mu
S.params.true$Ep <- S.params.true$Lp*4*Ne*mu/2/nseq

cat("Calculating MLE\n")
for(j in 1:NROW(matched)){
  ## cat(paste("\r", j))
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

denss <- list()
ll.mats <- list()
logt.grid <- seq(0, max.log, length.out=bins)
ll.mats[[1]] <- compute.ll.matrix( matched$true.map/100, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=NA, S.params=NA, verbose=FALSE)
ll.mats[[2]] <- compute.ll.matrix( matched$map.len, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=NA, S.params=NA, verbose=FALSE)
ll.mats[[3]] <- compute.ll.matrix( matched$map.len, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=error.params, S.params=NA, verbose=FALSE)
ll.mats[[4]] <- compute.ll.matrix( matched$true.map/100, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=NA, S.params=S.params.true, verbose=FALSE)
ll.mats[[5]] <- compute.ll.matrix( matched$map.len, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=NA, S.params=S.params, verbose=FALSE)
ll.mats[[6]] <- compute.ll.matrix( matched$map.len, matched$f2, Ne, p.fun, logt.grid=logt.grid, error.params=error.params, S.params=S.params, verbose=FALSE)

cat("Sampling Densities\n")
par(mfrow=c(2,3))
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
  plot.mle.and.density(matched$Age, t.hats[,i], denss[[i]], main=labels[i], xlim=c(0,max.log), ylim=c(0,max.log))
}
if(plots){dev.off()}

if(plots){pdf(paste(res.dir, "/coverage.pdf", sep=""), height=6, width=6)}else{dev.new()}
plot.coverage.curves(matched$Age, t.hats, labels=labels, lty=rep(1,6), lwd=1, cols=c("black", brewer.pal(5, "Set1")))
if(plots){dev.off()}
