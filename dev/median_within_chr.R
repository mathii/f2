## plot the median in bins along one chromosome.
code.dir <- "~/f2_age/code/"
path.to.results <- "~/f2_age/1000g/results/"
source(paste0(code.dir, "analysis/1kg_setup.R"))
source(paste0(code.dir, "libs/include.R"))

## Which CHR?
chrs <- c(9)
bin.w <- 1e6

## Which pops?
pops.include <- list("CEU"=c("CEU", "GBR", "FIN", "IBS", "TSI"), "ASN"=c("CHB", "CHS", "JPT"), "AFR"=c("YRI", "LWK")) 

start <- 0
end <- 150e6

bins.l <- seq(start, end, bin.w)
bins.r <- c(bins.l[2:length(bins.l)], Inf)
medians <- matrix(NA, nrow=length(pops.include), ncol=length(bins.l))



for(k in 1:length(chrs)){
  chr <- chrs[k]
  load(paste0(path.to.results, "/chr", chr, "/results/ll_environment.Rdata"))
  ## dyn.load(paste0(code.dir, "/libs/inference.so"))
  these.t.hats <-  MLE.from.haps(haps, Ne, S.params,  error.params, verbose=TRUE, tol=1e-10)
  for( i in 1:length(bins.l)){
      include <- haps$hap.left>=bins.l[i] & haps$hap.right<bins.r[i]
      for(j in 1:length(pops.include)){          
          this.include <- pop.map[haps$ID1] %in% pops.include[[j]] & pop.map[haps$ID2] %in% pops.include[[j]]
      medians[j,i] <- median(these.t.hats[include & this.include])
      }
  }  
}

## cols <- rainbow(length(chrs))
cols <- c("blue", "darkgreen", "brown")
plot(bins.l, medians[1,], log="y", type="s", xlab="Position", ylab="Median", bty="n", col=cols[1], ylim=c(100,4000))
if(length(pops.include)>1){
  for(j in 2:length(pops.include)){
    lines(bins.l, medians[j,], type="s", col=cols[j])
  }
}
legend("topleft", names(pops.include), col=cols, lty=1, bty="n")
