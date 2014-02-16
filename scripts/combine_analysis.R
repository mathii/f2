library(RColorBrewer)
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==4){
  chr.res.dir <- args[1]
  res.dir <- args[2]
  code.dir <- args[3]
  setup.file <- args[4]
  plots <- TRUE
} else{
  stop("Need to specify 4 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
source(setup.file)
chrs <- c(1:22)

######################################################################################################

t.hats <- list()
haps <- list()
subenv <- new.env()
logt.grid=NULL

## load saved results
cat("Loading data\n")
i=1
for(chr in chrs){
  cat(paste("\r", chr))
  load(paste(chr.res.dir, "/chr", chr, "/results/ll_environment.Rdata", sep=""), envir=subenv)

  if(is.null(subenv$t.hats)){   #backwards compatibility. 
    t.hats[[i]] <- MLE.from.haps(subenv$haps, subenv$Ne,S.params=subenv$S.params,  error.params=subenv$error.params, verbose=TRUE, tol=1e-10)
  } else{
    t.hats[[i]] <- subenv$t.hats
  }

  subenv$haps$chr <- chr
  haps[[i]] <- subenv$haps
  i=i+1
}

cat("\n")
rm(subenv)
t.hats <- do.call("c", t.hats)
haps <- do.call("rbind", haps)

## the following is cnp'd from run_1kg_analysis_chr.R.
## estimate densities, by population.
populations <- sort(unique(pop.map))
npop <- length(populations)             #14
densities <- rep(list(list()),npop)

ID1.pop <- pop.map[haps$ID1]
ID2.pop <- pop.map[haps$ID2]
for(i in 1:(npop)){
  for(j in i:npop){
    include <- (ID1.pop==populations[i]&ID2.pop==populations[j])|(ID1.pop==populations[j]&ID2.pop==populations[i])
    dens <- density(log10(t.hats[include]))
    densities[[i]][[j]] <- densities[[j]][[i]] <- approxfun(dens, rule=2)
  }
}

l.o <- c("ASW", "LWK", "YRI", "CLM", "MXL", "PUR",  "CHB", "CHS", "JPT", "CEU", "FIN", "GBR", "IBS", "TSI")
legend.order=order(match(populations, l.o))

save.image(paste0(res.dir, "/all_results.RData"))

## plots. One plot of all within-group densities, and one of all densities in total.
density.summary.plots(densities, populations, pop.cols, res.dir, xlim=c(1,5), ylim=c(0,1.2), legend.order=legend.order )
haplotype.count.summary( ID1.pop, ID2.pop, populations, res.dir, legend.order=legend.order)


