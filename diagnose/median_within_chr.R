## plot the median in bins along one chromosome.
code.dir <- "~/f2_age/code/"
path.to.results <- "~/f2_age/1000g/results/"
source(paste0(code.dir, "analysis/1kgsetup.R"))

## Which CHR?
chrs <- c(8:9)
bin.w <- 10e6

start <- 0
end <- 150e6

bins.l <- seq(start, end, bin.w)
bins.r <- c(bins.l[2:length(bins.l)], Inf)
medians <- matrix(NA, nrow=length(chrs), ncol=length(bins.l))
  
for(k in 1:length(chrs)){
  chr <- chrs[k]
  load(paste0(path.to.results, "/chr", chr, "/results/ll_environment.Rdata"))
  ## dyn.load(paste0(code.dir, "/libs/inference.so"))
  
  for( i in 1:length(bins.l)){
    cat(paste0("Iteration", i, "\n"))  
    include <- haps$hap.left>=bins.l[i] & haps$hap.right<bins.r[i]
    if(sum(include)>100){
      alpha <- round(0.05*sum(include))
      dens <- estimate.t.density.mcmc(0 ,0, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=1000, thin=10, ll.mat=ll.mat[include,])
      medians[k,i] <- quantile.density(dens, 0.5)
    }
  }  
}

cols <- rainbow(length(chrs))
plot(bins.l, medians[1,], type="s", xlab="Position", ylab="Median", bty="n", col=cols[1], ylim=c(2.1,2.4), xlim=c(0,150e6))
if(length(chrs)>1){
  for(j in 2:length(chrs)){
    lines(bins.l, medians[j,], type="s", col=cols[j])
  }
}
legend("topright", paste0("Chromosome ", chrs), col=cols, lty=1, bty="n")
