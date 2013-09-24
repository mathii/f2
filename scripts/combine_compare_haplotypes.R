## Combine the estimates from different chromosomes, plot overall power

args <- commandArgs(TRUE)
set.seed(12345)

######################################################################################################

if(length(args)==4){
  chr.res.dir <- args[1]
  res.dir <- args[2]
  code.dir <- args[3]
  dbl.pwr <- as.numeric(args[4])
  plots <- TRUE
} else{
  stop("Need to specify 4 arguments")
}

######################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
chrs <- 1:22
Ne <- 14000                             #variable
######################################################################################################

all.true.haps <- list()
all.matched.haps <- list()

## load saved results
i=1
for(chr in chrs){
  cat(paste("\r", chr))
  this.true.haps <- read.table(paste(chr.res.dir, "/chr", chr, "/haplotypes/NN_haplotypes.txt.gz", sep=""), as.is=TRUE, header=TRUE)
  map <- read.table(paste(chr.res.dir, "/chr", chr, "/map/cut.map.txt", sep=""), as.is=TRUE, header=TRUE)
  sim.mapfn <- approxfun(map[,2]-min(map[,2]), map[,4], rule=2 )
  this.true.haps$true.map <- sim.mapfn(this.true.haps$End)-sim.mapfn(this.true.haps$Start)
  this.true.haps <- this.true.haps[this.true.haps$ID1!=this.true.haps$ID2,]
  
  all.true.haps[[i]] <- this.true.haps
  this.matched.haps <- read.table(paste(chr.res.dir, "/chr", chr, "/results/matched_haps.txt.gz", sep=""), as.is=TRUE, header=TRUE)
  this.matched.haps$CHR <- chr
  all.matched.haps[[i]] <- this.matched.haps
  i=i+1
}
cat("\n")

true.haps <- do.call("rbind", all.true.haps)
matched <- do.call("rbind", all.matched.haps)
rm(all.true.haps)
rm(all.matched.haps)

pwr.matched <- matched[rbinom(NROW(matched), size=matched$f2, prob=dbl.pwr)>0,]

## Plot lengths
elements.Lg=list("All"=log10(true.haps$true.map/100), "Detected"=log10(matched$true.map/100), "66% power"=log10(pwr.matched$true.map/100) )
elements.Lp=list("All"=log10(true.haps$End-true.haps$Start), "Detected"=log10(matched$End-matched$Start), "66% power"=log10(pwr.matched$End-pwr.matched$Start) )
elements.t=list("All"=log10(true.haps$Age*4*Ne), "Detected"=log10(matched$Age), "66% power"=log10(pwr.matched$Age) )
blues=c("#377EBA40", "#377EBA80", "#377EBAB0")
border.blues=c("#377EBA", "#377EBA", "#377EBA")
reds=c("#E41A1C40", "#E41A1C80", "#E41A1CB0")
border.reds=c("#E41A1C", "#E41A1C", "#E41A1C")
if(plots){pdf(paste(res.dir, "/distribution_genetic_length.pdf", sep=""))}else{dev.new()}
overlapping.density.plot(elements.Lg, cols=reds, borders=border.reds, xlab=expression(Log[10]~"(Length (M))"), main="Genetic length")
if(plots){dev.off()}
if(plots){pdf(paste(res.dir, "/distribution_physical_length.pdf", sep=""))}else{dev.new()}
overlapping.density.plot(elements.Lp, cols=blues, borders=border.blues, xlab=expression(Log[10]~"(Length (b))"), main="Physical length")
if(plots){dev.off()}

## and power
x <- log10(matched$End - matched$Start)
y <- log10(matched$true.map/100)
z <- log10(matched$Age)
Lg.lim <- c(-5,0)
Lp.lim <- predict(lm(x~y), data.frame(y=Lg.lim))
## t.lim <- predict(lm(z~y), data.frame(y=Lg.lim))
t.lim=c(0,4)

if(plots){pdf(paste(res.dir, "/power.pdf", sep=""), width=6, height=6)}else{dev.new()}
par(list(mfrow=c(3,1), mar=c(3.8,4.1,4.1,2.1)))
overlapping.power.plot(elements.Lg, cols=c( "#377EBA", "#377EBA80"), xlab=expression(Log[10]~"(Length (M))"), xlim=Lg.lim, xaxt="n", yaxt="n", ylab="Power")
axis(1, at=seq(-6,1,1), line=-0.3)
axis(2, at=c(0,0.5,1))

## Power as a function of genetic length
par(list(mar=c(3.8,4.1,2.8,2.1)))
overlapping.power.plot(elements.Lp, cols=c( "#E41A1C", "#E41A1C80"), xlab=expression(Log[10]~"(Length (b))"), xlim=Lp.lim, xaxt="n", yaxt="n", ylab="Power")
axis(1, at=seq(0,10,0.5), line=-0.3)
axis(2, at=c(0,0.5,1))

## Power as a function of age. 
par(list(mar=c(5.1,4.1,2.8,2.1)))
overlapping.power.plot(elements.t, cols=c( "#4DAF4A", "#4DAF4A80"), xlab=expression(Log[10]~"(Age (generations))"), xlim=t.lim, xaxt="n", yaxt="n", ylab="Power", legend.pos="topright")
## axis(1, at=seq(-6,1,1), line=-0.3)
axis(1, at=seq(-1,5, 0.5), line=-0.3)
axis(2, at=c(0,0.5,1))

if(plots){dev.off()}



