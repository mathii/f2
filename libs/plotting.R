## plot plot plot

############################################################################################################################
## Plot the coverage - how many of our estimates are within a given range of the truth
## true.age: vector of length n with true ages in generations
## est.ages: matrix of nxk with the k lines we want plotted
## labels: names of the lines
## cols, lty: colors and types of the lines.
############################################################################################################################

plot.coverage.curves <- function(true.age, t.hats, labels, cols, lty, lwd=1, grid=seq(0,1,length.out=100), ... ){
  n <- NROW(t.hats)
  k <- NCOL(t.hats)

  errors <- t.hats*0
  
  cdfs <- list()
  for(i in 1:k){
    cdfs[[i]] <- ecdf(abs(log10(t.hats[,i]/true.age)))
  }

  plot(grid, cdfs[[1]](grid), col=cols[1], lty=lty[1], lwd=lwd, bty="n", xlim=c(0,1), ylim=c(0,1), type="l")

  if(k>1){
    for(i in 2:k){
      lines(grid, cdfs[[i]](grid), col=cols[i], lty=lty[i], lwd=lwd)
    }
  }

  legend("topleft", labels, bty="n", lwd=lwd, lty=lty, col=cols)
}

############################################################################################################################
## Add a qq line to compare the estimate to the truth
## dens: density function, asa function of log(age)
## true.age: vector of length n with true ages in generations
## qs - quantile limits
############################################################################################################################

add.density.to.plot <- function(dens, true.age, qs=c(0.01,0.99), col="#E41A1C", xmax=4){
  xs=seq(0,xmax, length.out=500)
  cumdens <- approxfun(xs, cumsum(dens(xs))/max(cumsum(dens(xs))), rule=2)
  qnts=cumdens(xs)
  interior <- qnts>=qs[1] & qnts<=qs[2]
  qnts=qnts[interior]
  true.qnts=quantile(log10(true.age), qnts)
  lines(true.qnts, xs[interior], col=col)
}

############################################################################################################################
## Scatterplot of the mles, plus a density estimate
## t.hat: vector of length n with estimated ages in generations
## dens: density function, asa function of log(age)
## true.age: vector of length n with true ages in generations
############################################################################################################################

plot.mle.and.density <- function(true.age, t.hat, dens, alpha="15", cols="#377EBA", ...){
  plot(log10(true.age), log10(t.hat), col=paste(cols, alpha, sep=""), bty="n", xlab="True age", ylab="Estimated age", pch=16, ...)
  qq.mle <- qqplot(log10(true.age), log10(t.hat), plot.it=FALSE)
  lines(qq.mle, col="#4DAF4A")
  abline(0,1,col="black", lty=2)
  add.density.to.plot(dens, true.age)
}

############################################################################################################################
## Density plot, with different densities overlapping in different colours.
## Scales the densities so that they are proportional to the first
## item of elements. 
## elements: list of vectors of objects, names taken from name of list
## cols: colours to plot
## borders: border colours
############################################################################################################################

overlapping.density.plot <- function(elements, cols, borders, main="", xlab="", legend.pos="topleft", ...){
  n <- length(elements)
  denss <- list()
  for(i in 1:length(elements)){
    denss[[i]] <- density(elements[[i]])
  }

  plot(denss[[1]], col=borders[1], type="l", bty="n", main=main, xlab=xlab, ...)
  for(i in 1:n){
    scale <- length(elements[[i]])/length(elements[[1]])
    polygon(denss[[i]]$x, scale*denss[[i]]$y, col=cols[i], border=borders[i] )
  }

  if(!all(is.null(names(elements)))){
    legend(legend.pos, names(elements), fill=cols, border=borders, bty="n")
  }
}

############################################################################################################################
## Simillar to the overlapping density plot, but plots the scale factor
## as a function of the element value. elements has to be longer than 2, then
## elements: list of vectors of objects, names taken from name of list
## cols: colours to plot
## borders: border colours
############################################################################################################################

overlapping.power.plot <- function(elements, cols, main="", xlab="", legend.pos="topleft", lwd=2, lty=1, ...){
  n <- length(elements)
  denss <- list()
  for(i in 1:length(elements)){
    denss[[i]] <- density(elements[[i]])
  }

  for(i in 2:n){
    scale <- length(elements[[i]])/length(elements[[1]])
    dfn <- approxfun(denss[[i]]$x, scale*denss[[i]]$y, rule=2)
    power <- pmin(denss[[1]]$y, dfn(denss[[i]]$x)) #in case we go over 1 when the sample is small
    if(i==2){
      plot(denss[[1]]$x, power/denss[[1]]$y, col=cols[1], main=main, xlab=xlab, bty="n", lwd=lwd, lty=lty, type="l", ...)
    }else{
      lines(denss[[1]]$x, power/denss[[1]]$y, col=cols[i-1], main=main, xlab=xlab, bty="n", lwd=lwd, lty=lty)
    }
  }
  
  if(!all(is.null(names(elements)))){
    legend(legend.pos, names(elements)[2:n], col=cols, lwd=lwd, lty=lty, bty="n")
  }
}

############################################################################################################################
## Just plot a list of density estimates, where the densities are just functions
## densities: A list of densities, with names
## grid: Interpolate at these points
## cols: Colors to plot densities
## lwd, lty: standard graphical parameters
############################################################################################################################

plot.densities <- function(densities, grid, cols, lwd=1, lty=1, ...){
  names <- names(densities)
  n <- length(names)
  plot(grid, densities[[1]](grid), col=cols[1], lwd=lwd, lty=lty, bty="n", type="l", ...)
  if(n>1){
    for(i in 2:n){
      lines(grid, densities[[i]](grid), col=cols[i], lwd=lwd, lty=lty)
    }
  }
  legend("topright", names(densities), col=cols, lwd=lwd, lty=lty, bty="n")
}

############################################################################################################################
## Find quantile of density
## density: a function, which is a density
## quantile: the quantile to find
############################################################################################################################

quantile.density <- function(dens, quantile, lower=0, upper=6){
  ff <- function(x){
    return((integrate(dens, lower=lower, upper=x, stop.on.error=FALSE)$value-quantile)^2)
  }

  return(optimize(ff, interval=c(lower,upper))$minimum) 
}

############################################################################################################################
## Make summary plots of densities and report quantiles. 
## densities: a list of lists where the [[i]][[j]] entry is the density of i-j chunks
## populations: names of populations
## pop.cols: population colours, named vector
## res.dir: directory to output results. if NULL then open windows. 
############################################################################################################################

density.summary.plots <- function(densities, populations, pop.cols, res.dir=NULL, ...){
  plots <- !is.null(res.dir)
  npop <- length(populations)
  
  within.list <- list()
  for(i in 1:npop){
    within.list[[i]] <- densities[[i]][[i]]
  }
  names(within.list) <- populations
  if(plots){pdf(paste(res.dir, "/within.pdf", sep=""))}else{dev.new()}
  plot.densities(within.list, logt.grid, cols=pop.cols[populations], main="Within", ...)
  if(plots){dev.off()}
  
  for(i in 1:npop){
    between.list <- list()
  for(j in 1:npop){
    between.list[[j]] <- densities[[i]][[j]]
  }
    names(between.list) <- populations

    if(plots){pdf(paste(res.dir, "/between_", populations[i], ".pdf", sep=""))}else{dev.new()}
    plot.densities(between.list, logt.grid, cols=pop.cols[populations], main=populations[i], ...)
    if(plots){dev.off()}
  }

  if(plots){                            #save tables of quantiles. 
    qs <- c(0.05, 0.5, 0.95)
    for(q in qs){
      res <- matrix(0, npop, npop)
      for(i in 1:npop){
        for(j in 1:npop){
          res[i,j] <- 10^quantile.density(densities[[i]][[j]], q)
        }
      }
      
      colnames(res) <- rownames(res) <- populations
      write.table(res, paste(res.dir, "/q", round(q*100), ".txt", sep=""), row.names=TRUE, col.names=TRUE, sep="\t")
    }
  }
}
