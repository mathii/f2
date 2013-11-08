## Just load the colour scheme and population assignments for the 1000 Genomes data
## The colours are the ones used in the 1000 genomes paper. 

## Actually, if you are making a new file for your own data, the only things
## you have to define are "pop.map", which maps sample names to populations
## and "pop.cols" which maps populations to colours. 

eur<-c("GBR", "FIN", "IBS", "CEU", "TSI"); 
asn<-c("CHS","CHB", "JPT"); 
afr<-c("YRI", "LWK", "ASW"); 
amr<-c("PUR", "CLM", "MXL"); 
pop.list<-c(eur, asn, afr, amr); 
pops.as.list<-list(eur, asn, afr, amr); 
cont.list<-c("EUR", "ASN", "AFR", "AMR");
names(pops.as.list) <- cont.list
cols<-c(68, 638, 62, 26, 73, 254, 259, 258, 513, 488, 567, 498, 585, 503); 
cols.cont<-c(592, 258, 649, 552);
cols<-colours()[cols];
names(cols) <- pop.list
cols.cont<-colours()[cols.cont];
names(cols.cont) <- cont.list

panel <- read.table(paste(code.dir, "/analysis/1kg_panel.txt", sep=""), as.is=TRUE, header=FALSE, sep="\t")
pop.map <- panel[,2]
names(pop.map) <- panel[,1]
pop.cols <- cols
cont.cols <- cols.cont

pop.to.cont <- rep(NA, length(pop.list))
names(pop.to.cont) <- pop.list
for(i in 1:length(pop.to.cont)){
  for(j in 1:length(cont.list)){
    if(pop.list[i] %in% pops.as.list[[j]]){pop.to.cont[i] <- cont.list[j]}
  }
}

cont.map <- pop.to.cont[pop.map]
names(cont.map) <- names(pop.map)

