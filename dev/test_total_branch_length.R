## test our computed expressions for the ratio of total branch length to 2t

source("~/f2_age/code/libs/coalescent_fns.R")
source("~/rare-var/coalescent.r")
library(scales)

#########################################
#Function to add branch lengths to tree
#########################################

annotate.tree<-function(tree) {
	colnames.old<-colnames(tree);
	tree<-cbind(tree, array(0, c(nrow(tree), 2)));
	colnames(tree)<-c(colnames.old, "Branch.length", "Parent");
	which.add<-ncol(tree)-1;
	which.add.2<-which.add+1;
	n<-(nrow(tree)+1)/2;
	for (j in (n+1):(2*n-1)) {
		tree[tree[j,3],which.add]<-tree[j,2]-tree[tree[j,3],2];
		tree[tree[j,4],which.add]<-tree[j,2]-tree[tree[j,4],2];
		tree[tree[j,3],which.add.2]<-j;
		tree[tree[j,4],which.add.2]<-j;
	}
	return(tree);
}

qns <- list(function(t){1}, q2, q3, q4, q5)
range <- c(0,4)
bins <- 100


grid <- seq(range[1], range[2], length.out=bins)


n.sims <- 10000

first <- TRUE
for(n in 2:5){
  t <- rep(0,n.sims)
  r <- rep(0,n.sims)
  for(i in 1:n.sims){
    tree <- annotate.tree(simulate.coalescent(sample=n, plot.tree=FALSE)$tree)
    t[i]<- rev(tree[,2])[1]
    r[i] <- sum(tree[,6])
  }

  ## From https://stat.ethz.ch/pipermail/r-help/2004-July/054984.html
  bpts  <- pretty(t, n=bins)
  INDEX <- cut(t, bpts, include.lowest = TRUE)
  bmeans <- tapply(r/2/t, INDEX, mean)
  bmids <- 0.5*(bpts[1:(length(bpts)-1)]+bpts[2:length(bpts)])
  
  if(first){
    plot(t, r/2/t, pch=".", col=alpha(n-1,min(1,1000/n.sims)), bty="n", xlim=range, ylim=c(1,2))
    first <- FALSE
  }else{
    points(t, r/2/t, pch=".", col=alpha(n-1,min(1,1000/n.sims)))
  }
  lines(bmids, bmeans, col=n-1, lty=1)  
  lines(grid, qns[[n]](grid), col=n-1, lty=3)
}
