PCA <- function(x,perc)
{
  p <- ncol(x);R <- cor(x);myEig <- eigen(R);sum.eigen <- sum(myEig$values)
  csum <- cumsum(myEig$values)/sum.eigen;index <- seq(1,p); d <- min(index[csum > perc])
  standardize <- function(x) {(x - mean(x))/sd(x)}; X <- apply(x, MARGIN=2, FUN=standardize)
  X.t <- X %*% myEig$vectors; x.s <- X.t[,1:d]
  list(x.s=x.s,d=d,load.matrix=myEig$vectors[,1:d])
}

PCA.f <- function(x,d)
{
  p <- ncol(x);R <- cor(x);myEig <- eigen(R);sum.eigen <- sum(myEig$values)
  csum <- cumsum(myEig$values)/sum.eigen;index <- seq(1,p); 
  standardize <- function(x) {(x - mean(x))/sd(x)}; X <- apply(x, MARGIN=2, FUN=standardize)
  X.t <- X %*% myEig$vectors; x.s <- X.t[,1:d]
  list(x.s=x.s)
}

###PCA using percentage#####
pca.group <- function(i,x,group,perc)
{
  p <- length(group);n <- nrow(x); pos_e <- cumsum(group);pos_s <- pos_e -group + 1; pe <- pos_e[i];ps <- pos_s[i]

  if ( group[i] != 1 ){
    xs <- x[,ps:pe]; xs <- xs[,colSums(xs)!=0];
    pca <- PCA(xs,perc);x.s <- pca$x.s; d <- pca$d
  }
  else if (group[i] == 1){
    x.s <- x[,ps:pe]; d <- 1
  }
  list(x.s=x.s,d=d)
}

###PCA using number predictors###
pca.group.f <- function(i,x,group,d)
{
  p <- length(group);n <- nrow(x); pos_e <- cumsum(group);pos_s <- pos_e -group + 1; pe <- pos_e[i];ps <- pos_s[i]

  if ( group[i] != 1 ){
    xs <- x[,ps:pe];xs <- xs[,colSums(xs)!=0];pca <- PCA.f(xs,d);x.s <- pca$x.s
  }
  else if (group[i] == 1){
    x.s <- x[,ps:pe]
  }
  list(x.s=x.s)
}

#source("/home/jin/CCT/Jin_matlab/Compound_bridge/PCA.r")