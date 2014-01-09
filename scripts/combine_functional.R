## Script to combine the analsis of different functional classes.
## Just by functional class, and within vs between. 

library(RColorBrewer)
args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==4){
  chr.res.dir <- args[1]
  res.dir <- args[2]
  code.dir <- args[3]
  setup.file <- args[4]
  plots <- FALSE
} else{
  stop("Need to specify 4 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
source(setup.file)
chrs <- c(1:22)

######################################################################################################

subenv <- new.env()
logt.grid=NULL

ll.mats.by.class <- list()
haps.by.class <- list()
t.hats.by.class <- list()

classes <- c("lof", "coding", "noncoding", "intergenic")
ns <- rep(0,2*length(classes))

for( cls in classes){

  haps <- list()
  t.hats <- list()
  
  ## load saved results
  cat("Loading data\n")
  i=1
  for(chr in chrs){
    cat(paste("\r", chr))
    load(paste0(chr.res.dir, "/chr", chr, "/", cls, "/results/ll_environment.Rdata"), envir=subenv)

    haps[[i]] <- subenv$haps  

  ## if(is.null(subenv$t.hats)){   #backwards compatibility. 
      t.hats[[i]] <- MLE.from.haps(subenv$haps, subenv$Ne,S.params=subenv$S.params,  error.params=subenv$error.params, verbose=TRUE, tol=1e-12)
  ## } else{
  ##     t.hats[[i]] <- subenv$t.hats
  ## }
    
    i=i+1
  }
  cat("\n")

  haps.by.class[[cls]] <- do.call("rbind", haps)
  t.hats.by.class[[cls]] <- do.call("c", t.hats)
}


## the following is cnp'd from run_1kg_analysis_chr.R.
## estimate densities, by class for within/between.
densities <- list()
values <- list()

i=1
for(cls in classes){
  haps <- haps.by.class[[cls]]
  
  ID1.pop <- pop.map[haps$ID1]
  ID2.pop <- pop.map[haps$ID2]

  within <- ID1.pop==ID2.pop
  between <- !within

  dw <- density(log10(t.hats.by.class[[cls]][within]))
  db <- density(log10(t.hats.by.class[[cls]][between]))
  
  densities[[i]] <- approxfun(dw, rule=2)
  densities[[i+length(classes)]] <- approxfun(db, rule=2)

  values[[i]] <- log10(t.hats.by.class[[cls]][within])
  values[[i+length(classes)]] <- log10(t.hats.by.class[[cls]][between])
  
  ns[i] <- sum(within)
  ns[i+length(classes)] <- sum(between)

  i=i+1
}

save.image(paste0(res.dir, "/functional_results.RData"))

pdf(paste0(res.dir, "/functional_distribution.pdf"))
col=c(rep("#377EBA", 4), rep("#E41A1C", 4))
border=col
fill=paste0(col, "80")
x.pos=c(1.25,2.25,3.25,4.25, 5.75,6.75,7.75, 8.75)
viola.plot(densities, x.pos=x.pos, eps=2e-2, col=col, border=border, fill=fill, labels=rep(c("LOF", "Coding", "Noncoding", "Intergenic"),2), xlim=c(1,9), ylab=expression(Age~(Log[10]~generations)), scale=0.25 )
abline(v=5, lty=3)
mtext(paste0("(",format(ns, big.mark=",", trim=TRUE),")"), 1, at=x.pos, line=1)
mtext(c("Within population", "Between populations"), 3, at=c(2.75, 7.25), line=0)
dev.off()
