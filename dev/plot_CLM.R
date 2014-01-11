## Run the first bit of /scripts/combine_analysis first. 

from <- which(populations=="GBR")
pops <- c( "CLM", "MXL", "PUR", "CEU", "IBS", "TSI")

pdf("~/Desktop/GBRvsCLM.pdf")
grid <- seq(1, 4, length.out=100)
d <- densities[[from]][[which(populations==pops[1])]]
plot(10^grid, d(grid), log="x", xlim=c(10,1000), ylim=c(0,1.2), lwd=2, type="l", col=pop.cols[pops[1]], bty="n", xlab="Age (generations)", ylab="Density")

for(i in 2:length(pops)){
    d <- densities[[from]][[which(populations==pops[i])]]
    lines(10^(grid), d(grid), lwd=2, col=pop.cols[pops[i]])
}

legend("topleft", pops, lwd=2, col=pop.cols[pops], bty="n")

dev.off()
