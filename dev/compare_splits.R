## make plots for some additional scenarios.

################################################################################
## Clean splits

load("~/f2_age/simulations/recent_split/chr20/results/ll_and_density_environment.Rdata")
dyn.load("~/f2_age/code/libs/inference.so")

within <- (matched$ID1<51 & matched$ID2<51)|(matched$ID1>50 & matched$ID2>50)
between <- !within

alpha.w <- round(0.05*sum(within))
alpha.b <- round(0.05*sum(between))
alpha <- round(0.05*NROW(matched))

  ## Lg and D are dummy entries... the ll matrix is the only information we use... 
d.w <- estimate.t.density.mcmc(0*matched$map.len[within,] ,0*matched$f2[within,], Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mats[[6]][within,])

d.b <- estimate.t.density.mcmc(0*matched$map.len[between,] ,0*matched$f2[within,], Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mats[[6]][between,])

dd <- estimate.t.density.mcmc(0*matched$map.len ,0*matched$f2, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mats[[6]])

pdf("~/f2_age/talk/figures/clean_split_density.pdf")
plot(logt.grid, d.w(logt.grid), col="#377EBA", type="l", bty="n", xlim=c(0,4), ylim=c(0,2), xlab=expression(Age~(log[10]~generations)), ylab="Density")
lines(logt.grid, d.b(logt.grid), col="#E41A1C")
abline(v=log10(560), lty=3)

legend("topleft", c("Within", "Between", "Split time"), col=c("#377EBA", "#E41A1C", "black"), lty=c(1,1,3), bty="n")
dev.off()

pdf("~/f2_age/talk/figures/clean_split_compare.pdf")
plot.mle.and.density(matched$Age, t.hats[,6], dd, main="", xlim=c(0.5,4), ylim=c(0.5,4), col="#000000", alpha="80")
abline(v=log10(560), lty=3)
dev.off()

################################################################################
## Dirty splits

load("~/f2_age/simulations/recent_split_migration/chr20/results/ll_and_density_environment.Rdata")
dyn.load("~/f2_age/code/libs/inference.so")

within <- (matched$ID1<51 & matched$ID2<51)|(matched$ID1>50 & matched$ID2>50)
between <- !within

alpha.w <- round(0.05*sum(within))
alpha.b <- round(0.05*sum(between))
alpha <- round(0.05*NROW(matched))

  ## Lg and D are dummy entries... the ll matrix is the only information we use... 
dm.w <- estimate.t.density.mcmc(0*matched$map.len[within,] ,0*matched$f2[within,], Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mats[[6]][within,])

dm.b <- estimate.t.density.mcmc(0*matched$map.len[between,] ,0*matched$f2[within,], Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mats[[6]][between,])

dmd <- estimate.t.density.mcmc(0*matched$map.len ,0*matched$f2, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mats[[6]])

pdf("~/f2_age/talk/figures/dirty_split_density.pdf")
plot(logt.grid, d.w(logt.grid), col="#377EBA", type="l", bty="n", xlim=c(0,4), ylim=c(0,2), xlab=expression(Age~(log[10]~generations)), ylab="Density")
lines(logt.grid, d.b(logt.grid), col="#E41A1C")
lines(logt.grid, dm.w(logt.grid), col="#377EBA", lty=2)
lines(logt.grid, dm.b(logt.grid), col="#E41A1C", lty=2)
legend("topleft", c("Within", "Between", "Split time"), col=c("#377EBA", "#E41A1C", "black"), lty=c(1,1,3), bty="n")

abline(v=log10(400), lty=3)
dev.off()

################################################################################
## Bottlenecks

load("~/f2_age/simulations/bottleneck/chr20/results/ll_and_density_environment.Rdata")
dd <- estimate.t.density.mcmc(0*matched$map.len ,0*matched$f2, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mats[[6]])
pdf("~/f2_age/talk/figures/bottleneck_density.pdf")
plot(logt.grid, dd(logt.grid), col="#377EBA", type="l", bty="n", xlim=c(0,4), ylim=c(0,2), xlab=expression(Age~(log[10]~generations)), ylab="Density")
legend("topleft", c("Density",  "Bottleneck"), col=c("#377EBA",  "black"), lty=c(1,3), bty="n")
abline(v=log10(14), lty=3)
abline(v=log10(16.8), lty=3)
dev.off()

pdf("~/f2_age/talk/figures/bottleneck_compare.pdf")
plot.mle.and.density(matched$Age, t.hats[,6], dd, main="", xlim=c(0.5,4), ylim=c(0.5,4), col="#000000", alpha="80")
abline(v=log10(14), lty=3)
abline(v=log10(16.8), lty=3)
dev.off()

################################################################################
## Expanding

load("~/f2_age/simulations/expanding/chr20/results/ll_and_density_environment.Rdata")
dd <- estimate.t.density.mcmc(0*matched$map.len ,0*matched$f2, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mats[[6]])

pdf("~/f2_age/talk/figures/expanding_compare.pdf")
plot.mle.and.density(matched$Age, t.hats[,6], dd, main="", xlim=c(0.5,4), ylim=c(0.5,4), col="#000000", alpha="80")
dev.off()

################################################################################
## Wrong Ne

load("~/f2_age/simulations/Ne_10/chr20/results/ll_and_density_environment.Rdata")
dd <- estimate.t.density.mcmc(0*matched$map.len ,0*matched$f2, Ne, p.fun, verbose=FALSE, logt.grid=logt.grid, prior=norm.2.p, alpha=alpha,error.params=NA, n.sims=10000, thin=100, ll.mat=ll.mats[[6]])

matched$Age <- matched$Age*10

pdf("~/f2_age/talk/figures/Ne_10_compare.pdf")
plot.mle.and.density(matched$Age, t.hats[,6], dd, main="", xlim=c(0,6), ylim=c(0,6), col="#000000", alpha="20")
dev.off()
