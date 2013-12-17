## Compare the detected f2 haplotypes with the truth. Save the results for
## later analysis.

args <- commandArgs(TRUE)
set.seed(12345)

if(length(args)==8){
  code.dir <- args[1]
  hap.dir <- args[2]
  res.dir <- args[3]
  Ne <- as.numeric(args[4])
  dbl.pwr <- as.numeric(args[5])
  map.file <- args[6]
  nseq <- as.numeric(args[7])/2
  simulator <- args[8]
  plots <- TRUE
} else{
  stop("Need to specify 8 arguments")
}

################################################################################################################

source(paste(code.dir, "/libs/include.R", sep=""))
f2.haps <- read.table(paste(res.dir, "/f2_haplotypes.txt.gz", sep=""),  as.is=TRUE, header=TRUE)
map <- read.table(map.file, as.is=TRUE, header=TRUE)
sim.mapfn <- approxfun(map[,2]-min(map[,2]), map[,4], rule=2 )
true.haps <- read.table(paste(hap.dir,"/NN_haplotypes.txt.gz", sep=""), as.is=TRUE, header=TRUE)

################################################################################################################

true.haps<-true.haps[order(true.haps$ID1, true.haps$ID2),]
true.haps$true.map <- sim.mapfn(true.haps$End)-sim.mapfn(true.haps$Start)
## We find a true hap if we fully cover it.
true.haps$hap.left <- NA
true.haps$hap.right <- NA
true.haps$f1 <- NA
true.haps$f2 <- NA
true.haps$map.len <- NA
true.haps$f2.hap.id <- NA
if(simulator=="macs"){
  true.haps$Age <- true.haps$Age*Ne*4     #Is this a factor of 2.. ms says it is 4Ne?}
} else if(simulator=="fastsimcoal"){
} else{
  stop("I only know what to do with macs or fastsimcoal output")
}
true.haps <- true.haps[true.haps$ID1!=true.haps$ID2,]

## a bit of a hack - to make sure that we catch the ones at the end.
l.min <- min(f2.haps$hap.left)
r.max <- max(f2.haps$hap.right)
f2.haps$hap.left[f2.haps$hap.left==l.min] <- min(true.haps$Start)
f2.haps$hap.right[f2.haps$hap.right==r.max] <- max(true.haps$End)

for( i in 1:NROW(true.haps)){
      select <- f2.haps[f2.haps$ID1==true.haps[i,"ID1"]&f2.haps$ID2==true.haps[i,"ID2"]&f2.haps$pos>=true.haps[i,"Start"]&f2.haps$pos<=true.haps[i,"End"],]

  if(NROW(select)==1){
    true.haps[i,c("hap.left", "hap.right", "f1", "f2", "map.len", "f2.hap.id")] <- c(select[c("hap.left", "hap.right", "f1", "f2", "map.len" )], as.numeric(rownames(select)[1]))
  }
  if(NROW(select)>1){
    cat("this should never happen - overlapping f2 haplotypes\n")
    true.haps[i,c("hap.left", "hap.right", "f1", "f2", "map.len", "f2.hap.id")] <- c(select[1,c("hap.left", "hap.right", "f1", "f2",  "map.len")], as.numeric(rownames(select)[1]))
  }

}

## augment....
matched <- true.haps[!is.na(true.haps$f2.hap.id),]
matched$hap.len <- matched$hap.right-matched$hap.left
matched$true.len<-matched$End-matched$Start
matched$true.f1 <- count.singletons.from.positions(matched[,c("ID1", "ID2")], matched$Start, matched$End, paste(hap.dir, "/pos.idx.f1.gz", sep=""))
## filter...
matched <- matched[matched$ID1!=matched$ID2,] #remove haplotypes you matched with yourself
matched <- matched[matched$true.map>0,]       #remove things with 0 length (shouldn't really happen)

write.table(matched, paste(res.dir, "/matched_haps.txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

## what would happen if we didn't see all doubletons
pwr.matched <- matched[rbinom(NROW(matched), size=matched$f2, prob=dbl.pwr)>0,]


## Plot lengths
pwr.lab=paste(round(100*dbl.pwr), "% power", sep="")
elements.Lg=list("All"=log10(true.haps$true.map/100), "Detected"=log10(matched$true.map/100), "Power"=log10(pwr.matched$true.map/100) )
elements.Lp=list("All"=log10(true.haps$End-true.haps$Start), "Detected"=log10(matched$End-matched$Start), "Power"=log10(pwr.matched$End-pwr.matched$Start) )
names(elements.Lg)[[3]] <- pwr.lab
names(elements.Lp)[[3]] <- pwr.lab
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
Lg.lim <- c(-5,0)
Lp.lim <- predict(lm(x~y), data.frame(y=Lg.lim))

if(plots){pdf(paste(res.dir, "/power.pdf", sep=""))}else{dev.new()}
par(list(mfrow=c(2,1), mar=c(3.8,4.1,4.1,2.1)))
overlapping.power.plot(elements.Lp, cols=c( "#377EBA", "#377EBA80"), xlab=expression(Log[10]~"(Length (b))"), xlim=Lp.lim, xaxt="n", yaxt="n", ylab="Power")
axis(1, at=seq(2,7,0.5), line=-0.3)
axis(2, at=c(0,0.5,1))
## Power as a function of genetic length
par(list(mar=c(5.1,4.1,2.8,2.1)))
overlapping.power.plot(elements.Lg, cols=c( "#E41A1C", "#E41A1C80"), xlab=expression(Log[10]~"(Length (M))"), xlim=Lg.lim, xaxt="n", yaxt="n", ylab="Power")
axis(1, at=seq(-6,1,1), line=-0.3)
axis(2, at=c(0,0.5,1))
if(plots){dev.off()}

