## compare the msmc results the other way.
library(RColorBrewer)

sim.types=c("ancient_split", "ancient_split_migration", "ancient_split_long_migration" )
mu <- 1.2e-8
cols <- brewer.pal(3, "Set1")
scale=1.8
for(outer.i in 1:length(sim.types)){
  sim.type=sim.types[outer.i]
  load(paste0("~/f2_age/simulations/", sim.type, "/chr20/results/ll_and_density_environment.Rdata"))
  dyn.load("~/f2_age/code/libs/inference.so")
  msmc <- read.table(paste0("~/f2_age/simulations/", sim.type, "/chr20/results/msmc_results_all.final.txt"), header=T, as.is=T)
     
  within <- (matched$ID1<51 & matched$ID2<51)|(matched$ID1>50 & matched$ID2>50)
  between <- !within

  alpha.b <- round(0.05*sum(between))

  d.b <- estimate.t.density.mcmc(0*matched$map.len[between,] ,0*matched$f2[between,], Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha.b ,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mats[[6]][between,])

  par(mar=c(5.1,4.1,2.1,4.1))
  if(outer.i==1){
    plot(msmc$left_time_boundary/2/mu, 2*msmc$lambda_01/(msmc$lambda_00+msmc$lambda_11), type="s", log="x", xlim=c(100,2000), ylim=c(0,1), bty="n", xlab="Generations", ylab="MSMC gene flow estimate", col=cols[outer.i], yaxt="n", main=sim.type, lty=1)
    axis(2)
    axis(4, lty=2, labels=FALSE, tick=FALSE)
    mtext(expression(f[2]~age~density), 4, line=2)
  } else{
    lines(msmc$left_time_boundary/2/mu, 2*msmc$lambda_01/(msmc$lambda_00+msmc$lambda_11), type="s", col=cols[outer.i], lty=1)

  }
    
  ## lines(10^logt.grid, d.b(logt.grid)/scale, lty=1, col=cols[outer.i])
  max.pt <- 10^logt.grid[d.b(logt.grid)==max(d.b(logt.grid))]
  abline(v=max.pt, col=cols[outer.i], lty=2)
  
## chose other features here
}

lines(c(224,1120), c(-0.02, -0.02), col=cols[3], type="b", pch=16, lwd=2)
lines(c(560,1120), c(-0.0, -0.0), col=cols[2], type="b", pch=16, lwd=2)
lines(c(1120,1120), c(0.02, 0.02), col=cols[1], type="b", pch=16, lwd=2)

## abline(v=1120, col="black", lty=3)
## abline(v=1110, col=cols[1], lty=3)
## abline(v=560, lty=3, col=cols[2])
## abline(v=224, lty=3, col=cols[3])

