dyn.load("gmcp.so")

GMCP_quant <- function(x,y,group,lambda,gamma)
{
  n <- nrow(x); p <- ncol(x); p.group <- length(group);
  
  end <- cumsum(group); start <- end -group + 1

  fun <- function(i,x)
  {
    sigma = t(x[,start[i]:end[i]])%*%x[,start[i]:end[i]]/n
    R = chol(sigma);x.t = x[,start[i]:end[i]]%*%solve(R)
    return(list(R=R,x.t=x.t))
  }
  result <- lapply(seq(1,p.group),fun,x) 

  fun <- function(j,result)
  {
    result[[j]]$x.t
  }
  x.t <- matrix(unlist(lapply(seq(1,p.group),fun,result)),nrow=n)
  
  beta <- numeric(p)
  param <- c(n, p.group, p,500,max(group))
  epsilon <- 1E-10

  fit <- .C("GMCP", y=as.double(y),x=as.double(t(x.t)),G=as.integer(group),param=as.integer(param),
            lambda=as.double(lambda),gamma=as.double(gamma),epsilon=as.double(epsilon),beta=as.double(beta))

  b <- fit$beta; beta <- numeric(p); 
  for ( j in 1: p.group)
  {
     beta[start[j]:end[j]] <- solve(result[[j]]$R)%*%b[start[j]:end[j]]
  }

  list(beta=beta,b=b,diff=fit$diff,x.t=x.t)
}
