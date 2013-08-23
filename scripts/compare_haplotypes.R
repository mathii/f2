## Compare the detected f2 haplotypes with the truth. Save the results for
## later analysis.

args <- commandArgs(TRUE)
set.seed(12345)

if(length(args)==7){
  code.dir <- args[1]
  hap.dir <- args[2]
  res.dir <- args[3]
  Ne <- as.numeric(args[4])
  dbl.pwr <- as.numeric(args[5])
  map.file <- args[6]
  nseq <- as.numeric(args[7])/2
  plots <- TRUE
} else{
  stop("Need to specify 7 arguments")
}

################################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
f2.haps <- read.table(paste(res.dir, "/f2_haplotypes.txt.gz", sep=""),  as.is=TRUE, header=TRUE)
map <- read.table(map.file, as.is=TRUE, header=TRUE)
sim.mapfn <- approxfun(map[,2]-min(map[,2]), map[,4] )
true.haps <- read.table(paste(hap.dir,"/NN_haplotypes.txt.gz", sep=""), as.is=TRUE, header=TRUE)

################################################################################################################

true.haps<-true.haps[order(true.haps$ID1, true.haps$ID2),]
true.haps$true.map <- sim.mapfn(true.haps$End)-sim.mapfn(true.haps$Start)
## We find a true hap if we fully cover it.
true.haps$hap.left <- NA
true.haps$hap.right <- NA
true.haps$f1 <- NA
true.haps$f2 <- NA
true.haps$t.hat <- NA
true.haps$map.len <- NA
true.haps$f2.hap.id <- NA
true.haps$Age <- true.haps$Age*Ne*4     #Is this a factor of 2.. ms says it is 4Ne?? 

## a bit of a hack - to make sure that we catch the ones at the end.
l.min <- min(f2.haps$hap.left)
r.max <- max(f2.haps$hap.right)
f2.haps$hap.left[f2.haps$hap.left==l.min] <- min(true.haps$Start)
f2.haps$hap.right[f2.haps$hap.right==r.max] <- max(true.haps$End)

for( i in 1:NROW(true.haps)){
      select <- f2.haps[f2.haps$ID1==true.haps[i,"ID1"]&f2.haps$ID2==true.haps[i,"ID2"]&f2.haps$pos>=true.haps[i,"Start"]&f2.haps$pos<=true.haps[i,"End"],]

  if(NROW(select)==1){
    true.haps[i,c("hap.left", "hap.right", "f1", "f2", "t.hat", "map.len", "f2.hap.id")] <- c(select[c("hap.left", "hap.right", "f1", "f2", "t.hat", "map.len" )], as.numeric(rownames(select)[1]))
  }
  if(NROW(select)>1){
    cat("this should never happen - overlapping f2 haplotypes\n")
    true.haps[i,c("hap.left", "hap.right", "f1", "f2", "t.hat", "map.len", "f2.hap.id")] <- c(select[1,c("hap.left", "hap.right", "f1", "f2", "t.hat",  "map.len")], as.numeric(rownames(select)[1]))
  }

}

## augment....
matched <- true.haps[!is.na(true.haps$f2.hap.id),]
matched$hap.len <- matched$hap.right-matched$hap.left
matched <- matched[matched$ID1!=matched$ID2,] #remove haplotypes you matched with yourself
matched$true.len<-matched$End-matched$Start
matched$true.f1 <- count.singletons.from.positions(matched[,c("ID1", "ID2")], matched$Start, matched$End, paste(hap.dir, "/pos.idx.f1.gz", sep=""))

write.table(matched, paste(res.dir, "/matched_haps.txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

## what would happen if we didn't see all doubletons
pwr.matched <- matched[rbinom(NROW(matched), size=matched$f2, prob=dbl.pwr)>0,]

all.map.hist <- hist(log10(true.haps$true.map), plot=FALSE, breaks=40)
match.map.hist <- hist(log10(matched$true.map), plot=FALSE, breaks=all.map.hist$breaks)
pwr.match.map.hist <- hist(log10(pwr.matched$true.map), plot=FALSE, breaks=all.map.hist$breaks)

all.hist <- hist(log10(true.haps$End-true.haps$Start), plot=FALSE, breaks=40)
match.hist <- hist(log10(matched$End-matched$Start), plot=FALSE, breaks=all.hist$breaks)
pwr.match.hist <- hist(log10(pwr.matched$End-pwr.matched$Start), plot=FALSE, breaks=all.hist$breaks)

if(plots){pdf(paste(res.dir, "/length_distribution.pdf", sep=""))}else{dev.new()}
## compute histogram to show power.
this.power <- length(unique(matched$f2.hap.id))/NROW(true.haps)
this.pwr.power <- length(unique(pwr.matched$f2.hap.id))/NROW(true.haps)
plot(all.hist$mids, all.hist$density, col="#377EBA", type="l", bty="n", xlab=expression(Log[10]~"(Length (b))"), ylab="Density", xaxt="n")
axis(1,at=c(1:8))
polygon(c(1,all.hist$mids,8), c(0,all.hist$density,0), col="#377EBA40", border="#377EBA" )
polygon(c(1,match.hist$mids,8), c(0,match.hist$density*this.power,0), col="#377EBA80", border="#377EBA")
polygon(c(1,pwr.match.hist$mids,8), c(0,pwr.match.hist$density*this.pwr.power,0), col="#377EBAB0", border="#377EBA")
if(plots){dev.off()}

## and power
x <- log10(matched$End - matched$Start)
y <- log10(matched$true.map)
lims <- predict(lm(x~y), data.frame(y=c(-2.5,1)))

if(plots){pdf(paste(res.dir, "/power.pdf", sep=""))}else{dev.new()}
par(list(mfrow=c(2,1), mar=c(3.8,4.1,4.1,2.1)))
plot(all.hist$mids ,match.hist$count/all.hist$count, col="#377EBA80", lwd=2, bty="n", type="l", pch=16, xlab=expression(Log[10]~"(Length (b))"), ylab="Power", xlim=lims, xaxt="n", yaxt="n")
lines(all.hist$mids ,pwr.match.hist$count/all.hist$count, col="#377EBAB0", lwd=2, bty="n", type="l", pch=16)
axis(1, at=seq(2,7,0.5), line=-0.3)
axis(2, at=c(0,0.5,1))
## Power as a function of genetic length

par(list(mar=c(5.1,4.1,2.8,2.1)))
plot(all.map.hist$mids ,match.map.hist$count/all.map.hist$count, col="#E41A1C80", lwd=2, bty="n", type="l", pch=16, xlab=expression(Log[10]~"(Length (cM))"), ylab="Power", xlim=c(-2.5,1), xaxt="n", yaxt="n")
lines(all.map.hist$mids ,pwr.match.map.hist$count/all.map.hist$count, col="#E41A1CB0", lwd=2, bty="n", type="l", pch=16)
axis(1, at=seq(-3,1.5,0.5), line=-0.3)
axis(2, at=c(0,0.5,1))
dev.off()

