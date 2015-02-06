## Just load the colour scheme and population assignments for the 1000 Genomes data
## The colours are the ones used in the 1000 genomes paper. 

## Actually, if you are making a new file for your own data, the only things
## you have to define are "pop.map", which maps sample names to populations
## and "pop.cols" which maps populations to colours. 

cols <- c("goldenrod", "peru", "navajowhite3", "lightgoldenrod3", "palegoldenrod", "brown3", "orangered", "sienna1", "brown", "tomato3", "orange", "green3", "greenyellow", "green", "green4", "limegreen", "blue", "turquoise3", "cyan", "cornflowerblue", "darkblue", "darkmagenta", "darkviolet", "maroon", "deeppink2", "magenta", "black", "grey50", "darkslategray1", "darkslategrey", "mediumaquamarine", "sandybrown")
pops <- c("ESN", "GWD", "LWK", "MSL", "YRI", "ACB", "ASW", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU", "Neanderthal", "Denisovan", "LBK", "Loschbour", "Ust_Ishim", "Clovis")

panel <- read.table(paste(code.dir, "/analysis/1kgP3_ancient_panel.txt", sep=""), as.is=TRUE, header=FALSE, sep="\t")
pop.map <- panel[,2]
names(pop.map) <- panel[,1]
pop.cols <- cols
names(pop.cols) <- pops
