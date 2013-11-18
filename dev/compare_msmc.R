## compare msmsc gene flow estimates to the between population
## f2 density estimates.

## sim.type="ancient_split_long_migration" 
## sim.type="ancient_split_migration" 
## sim.type="ancient_split" 
## sim.type="third_party" 
## sim.type="ancient_split_migration_growth" 
## sim.type="ancient_split_migration_bottleneck" 
sim.type="ancient_split_migration_200" 

chr="20"

load(paste0("~/f2_age/simulations/", sim.type, "/chr", chr, "/results/ll_and_density_environment.Rdata"))
dyn.load("~/f2_age/code/libs/inference.so")
msmc <- read.table(paste0("~/f2_age/simulations/", sim.type, "/chr", chr, "/results/msmc_results_all.final.txt"), header=T, as.is=T)

within <- (matched$ID1<51 & matched$ID2<51)|(matched$ID1>50 & matched$ID2>50)
within <- (matched$ID1<101 & matched$ID2<101)|(matched$ID1>100 & matched$ID2>100)
between <- !within

alpha.w <- round(0.05*sum(within))
alpha.b <- round(0.05*sum(between))
alpha <- round(0.05*NROW(matched))

d.b <- estimate.t.density.mcmc(0*matched$map.len[between,] ,0*matched$f2[between,], Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mats[[6]][between,])

par(mar=c(5.1,4.1,2.1,4.1))
plot(msmc$left_time_boundary/1.2e-8, 2*msmc$lambda_01/(msmc$lambda_00+msmc$lambda_11), type="s", log="x", ylim=c(0,1), xlim=c(200,2000), bty="n", xlab="Generations", ylab="MSMC gene flow estimate", col="#E41A1C", yaxt="n", main=sim.type)
axis(2, col="#E41A1C")
scale <- max(d.b(logt.grid))
lines(10^logt.grid, d.b(logt.grid)/scale, col="#377EBA")
axis(4, col="#377EBA", labels=FALSE)
mtext(expression(f[2]~age~density), 4, line=2)

## chose other features here
## abline(v=224, lty=3)
## abline(v=560*2, lty=3)
abline(v=1120, lty=3)
abline(v=560, lty=3)
