library(RColorBrewer)
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

min.f2<-1 # Minimum number of f2 in haplotype
if(length(args)==4||length(args)==5){
  chr.res.dir <- args[1]
  res.dir <- args[2]
  code.dir <- args[3]
  setup.file <- args[4]
  plots <- TRUE
  if(length(args)==5){
    min.f2<-as.numeric(args[5])
  }
} else{
  stop("Need to specify 4 or 5 arguments")
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
  ## Legacy - some of my old datests have extra columns.
  haps[[i]] <- haps[[i]][,!(names(haps[[i]]) %in% c("ID.from", "ID.to" ))]
  i=i+1
}

cat("\n")
rm(subenv)
t.hats <- do.call("c", t.hats)
haps <- do.call("rbind", haps)

#Cut out everything with less than min.f2 variants
b4<-NROW(haps)
haps <- haps[haps$f2>=min.f2,]
a4<-NROW(haps)
cat(paste0("Removed ",  b4-a4, "haplotypes with too few f2 variants\n"))

## the following is cnp'd from run_1kg_analysis_chr.R.
## estimate densities, by population.
populations <- sort(unique(pop.map))
npop <- length(populations)             #14
densities <- rep(list(list()),npop)
q50.direct <- matrix(0, nrow=npop, ncol=npop)

ID1.pop <- pop.map[haps$ID1]
ID2.pop <- pop.map[haps$ID2]
bw <- c()
for(i in 1:(npop)){
    for(j in i:npop){
        include <- (ID1.pop==populations[i]&ID2.pop==populations[j])|(ID1.pop==populations[j]&ID2.pop==populations[i])
        if(sum(include)>1){
            dens <- density(log10(t.hats[include]))
            densities[[i]][[j]] <- densities[[j]][[i]] <- approxfun(dens, rule=2)
            q50.direct[i,j] <- q50.direct[j,i] <- median(t.hats[include])
        }else{
            densities[[i]][[j]] <- densities[[j]][[i]] <- function(x){return(0*x)}
            q50.direct[i,j] <- q50.direct[j,i] <- 999999
        }
    }
}

l.o <- c("ESN", "GWD", "LWK", "MSL", "YRI", "ACB", "ASW", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU")
legend.order=order(match(populations, l.o))

save.image(paste0(res.dir, "/all_results.RData"))

## plots. One plot of all within-group densities, and one of all densities in total.
density.summary.plots(densities, populations, pop.cols, res.dir, xlim=c(1,5), ylim=c(0,1.2), legend.order=legend.order, legend.cex=0.5 )
haplotype.count.summary( ID1.pop, ID2.pop, populations, res.dir, pop.counts=table(pop.map)[populations], legend.order=legend.order)

rownames(q50.direct) <- colnames(q50.direct) <- populations

write.table(10^(q50.direct[legend.order,legend.order]), paste0(res.dir, "/q50_direct.txt"), row.names=TRUE, col.names=TRUE, sep="\t")




