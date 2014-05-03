## Just load the colour scheme and population assignments for the 1000 Genomes data
## The colours are the ones used in the 1000 genomes paper. 

## Actually, if you are making a new file for your own data, the only things
## you have to define are "pop.map", which maps sample names to populations
## and "pop.cols" which maps populations to colours. 

## eur<-c("GBR", "FIN", "IBS", "CEU", "TSI"); 
## asn<-c("CHS","CHB", "JPT"); 
## afr<-c("YRI", "LWK", "ASW"); 
## amr<-c("PUR", "CLM", "MXL"); 
## pop.list<-c(eur, asn, afr, amr); 
## pops.as.list<-list(eur, asn, afr, amr); 
## cont.list<-c("EUR", "ASN", "AFR", "AMR");
## names(pops.as.list) <- cont.list
## cols<-c(68, 638, 62, 26, 73, 254, 259, 258, 513, 488, 567, 498, 585, 503); 
## cols.cont<-c(592, 258, 649, 552);
## cols<-colours()[cols];
## names(cols) <- pop.list
## cols.cont<-colours()[cols.cont];
## names(cols.cont) <- cont.list

## cont.cols <- cols.cont

## pop.to.cont <- rep(NA, length(pop.list))
## names(pop.to.cont) <- pop.list
## for(i in 1:length(pop.to.cont)){
##   for(j in 1:length(cont.list)){
##     if(pop.list[i] %in% pops.as.list[[j]]){pop.to.cont[i] <- cont.list[j]}
##   }
## }

## cont.map <- pop.to.cont[pop.map]
## names(cont.map) <- names(pop.map)

cols <- c("goldenrod", "peru", "navajowhite3", "lightgoldenrod3", "palegoldenrod", "brown3", "orangered", "sienna1", "brown", "tomato3", "orange", "green3", "greenyellow", "green", "green4", "limegreen", "blue", "turquoise3", "cyan", "cornflowerblue", "darkblue", "darkmagenta", "darkviolet", "maroon", "deeppink2", "magenta")
pops <- c("ESN", "GWD", "LWK", "MSL", "YRI", "ACB", "ASW", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU")

panel <- read.table(paste(code.dir, "/analysis/1kgP3_panel.txt", sep=""), as.is=TRUE, header=FALSE, sep="\t")
pop.map <- panel[,2]
names(pop.map) <- panel[,1]
pop.cols <- cols
names(pop.cols) <- pops
