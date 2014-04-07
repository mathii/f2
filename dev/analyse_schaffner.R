## Load and analyse the schaffner simulation.

load("~/f2_age/simulations/schaffner/chr10/results/ll_and_density_environment.Rdata")
code.dir <- "~/f2_age/code/"
source("~/f2_age/code/libs/include.R")

res.dir="~/f2_age/simulations/schaffner/chr10/results"

t.hats <- MLE.from.haps(matched, Ne, S.params=S.params,  error.params=error.params, verbose=TRUE, tol=1e-10)

pop.map <- c(rep("AFR",50), rep("EUR",50), rep("ASN",50), rep("AFM",50))
names(pop.map) <- c(paste0("EUR", 1:50), paste0("EUR", 1:50),paste0("ASN", 1:50),paste0("AFM", 1:50))
pop.cols <- c("brown", "blue", "green", "purple")
names(pop.cols) <- c("AFR", "EUR", "ASN", "AFM")
    # 1 africans
    # 2 europeans
    # 3 asians
    # 4 african-americans

populations <- sort(unique(pop.map))
npop <- length(populations)             #14
densities <- rep(list(list()),npop)

ID1.pop <- pop.map[matched$ID1]
ID2.pop <- pop.map[matched$ID2]
for(i in 1:(npop)){
  for(j in i:npop){
    include <- (ID1.pop==populations[i]&ID2.pop==populations[j])|(ID1.pop==populations[j]&ID2.pop==populations[i])
    dens <- density(log10(t.hats[include]))
    densities[[i]][[j]] <- densities[[j]][[i]] <- approxfun(dens, rule=2)
  }
}


legend.order=names(pop.cols)
## plots. One plot of all within-group densities, and one of all densities in total.
density.summary.plots(densities, populations, pop.cols, res.dir, xlim=c(1,5), ylim=c(0,1.2), legend.order=legend.order )
