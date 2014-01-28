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

haps <- list(list(), list(), list(), list())
logt.grid=NULL

## load saved results
cat("Loading data\n")
j=1
for(f in 2:5){
    i=1
    for(chr in chrs){
        cat(paste("\r", chr))
        haps[[j]][[i]] <- read.table(paste(chr.res.dir, "/chr", chr, "/results/f", f, "_results.txt.gz", sep=""), as.is=TRUE, header=TRUE)
        i=i+1
    }
    haps[[j]] <- do.call("rbind", haps[[j]])
    j=j+1
}

cat("\n")

## the following is cnp'd from run_1kg_analysis_chr.R.
## estimate densities, by population.
populations <- sort(unique(pop.map))
npop <- length(populations)             #14
densities <- rep(list(list()),npop)

k=1
for(f in 2:5){
    cat(paste0(f, "\n"))
    ID1.pop <- pop.map[haps[[k]]$ID1]
    ID2.pop <- pop.map[haps[[k]]$ID2]

    pops <- matrix("", nrow=NROW(haps[[k]]), ncol=f)
    ID.labs <- paste0("ID", 1:f)
    ID.all <-haps[[k]][,ID.labs]

    for(foo in 1:NROW(ID.all)){
        for(bar in 1:NCOL(ID.all)){
            pops[foo,bar] <- pop.map[ID.all[foo,bar]]
        }
    }
    
    for(i in 1:(npop)){
        for(j in i:npop){
            include <- rep(FALSE, NROW(pops))
            for(pp in 1:NROW(pops)){
                include[pp] <- ( populations[i] %in% pops[pp,])& (populations[j] %in% pops[pp,])
            }
            dens <- density(log10(haps[[k]]$t.hat[include]))
            densities[[i]][[j]] <- densities[[j]][[i]] <- approxfun(dens, rule=2)
        }
    }
    
    l.o <- c("ASW", "LWK", "YRI", "CLM", "MXL", "PUR",  "CHB", "CHS", "JPT", "CEU", "FIN", "GBR", "IBS", "TSI")
    legend.order=order(match(populations, l.o))
    
## plots. One plot of all within-group densities, and one of all densities in total.
    density.summary.plots(densities, populations, pop.cols, res.dir, xlim=c(1,5), ylim=c(0,1.2), legend.order=legend.order, prefix=paste0("f", f, "_") )
    haplotype.count.summary( ID1.pop, ID2.pop, populations, res.dir, legend.order=legend.order, prefix=paste0("f", f, "_"))
    k=k+1
}

## Now do it by pair, but all f on same graph. 
