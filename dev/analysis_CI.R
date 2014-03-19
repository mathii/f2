## compute confidence intervals for each t.hat.

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
cis <- list()
haps <- list()
subenv <- new.env()
logt.grid=NULL

## load saved results
cat("Loading data\n")
i=1
for(chr in chrs){
  cat(paste("\r", chr))
  load(paste(chr.res.dir, "/chr", chr, "/results/ll_environment.Rdata", sep=""), envir=subenv)

  t.hats[[i]] <- MLE.from.haps(subenv$haps, subenv$Ne,S.params=subenv$S.params,  error.params=subenv$error.params, verbose=TRUE, tol=1e-10)
  
  ci[[i]] <- confidence.interval(subenv$haps$map.len, subenv$Ne, function(x){1}, subenv$haps$f2, subenv$error.params, subenv$haps$S.params, 0.95, max.search=10)

  subenv$haps$chr <- chr
  haps[[i]] <- subenv$haps
  ## Legacy - some of my old datests have extra columns.
  haps[[i]] <- haps[[i]][,!(names(haps[[i]]) %in% c("ID.from", "ID.to" ))]
  i=i+1
}

cat("\n")
rm(subenv)
t.hats <- do.call("c", t.hats.k1)
cis <- do.call("rbind", cis)
haps <- do.call("rbind", haps)

