library(mvtnorm)

generate <- function(beta,k,n,covar,group,mag,quan)
{
 beta0 <- 0.5; p <- sum(group);  n.group <- length(group);

 x <- rmvnorm(n = n, mean = rep(0, p), covar) 
 
 index =seq(1,p)
 fun <- function(j,z)
 {
   t = z[,j]; n = dim(z)[1]; d = numeric(n);d[t <= -0.6744898] =0
   d[t>-0.6744898 & t<0.6744898]=1; d[t>=0.6744898]=2
   list (d=d)
 }
 x <- matrix(unlist(lapply(index,fun,x)),nrow=n);

 res <- rnorm(n = n, 0 , k); Ti <- beta0 + x%*%beta + res;
 
 Ce <- runif(n,min=0,max=quantile(exp(Ti),c(quan)))
 y <- pmin(exp(Ti), Ce); d <- as.numeric(exp(Ti) <= Ce); xs=x[order(y),]
 d=d[order(y)]; y=sort(y)

 w <- numeric(n); w[1]=d[1]/n
 for ( i in 2:n )
 {
   tmp = 1
   for ( j in 1: (i-1) )
     tmp = tmp*((n-j)/(n-j+1))^d[j]
     
   w[i]=d[i]/(n-i+1)*tmp
 }

 list(y=y,x=xs,d=d,w=w)
}
