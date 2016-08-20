dyn.load("glasso_weight.so")
com_bridge <- function(x,y,group.inner,group.outer,lambda,gamma)
{
  n <- nrow(x); p <- ncol(x); p.inner <- length(group.inner); p.outer <- length(group.outer)
  
  end <- cumsum(group.inner); start <- end -group.inner + 1

  fun <- function(i,x)
  {
    sigma = t(x[,start[i]:end[i]])%*%x[,start[i]:end[i]]/n
    R = chol(sigma);x.t = x[,start[i]:end[i]]%*%solve(R)
    return(list(R=R,x.t=x.t))
  }
  result <- lapply(seq(1,p.inner),fun,x) 

  fun <- function(j,result)
  {
    result[[j]]$x.t
  }
  x.t <- matrix(unlist(lapply(seq(1,p.inner),fun,result)),nrow=n)

  theta <- numeric(p.inner); tau_n <-  (lambda/((1-gamma)/gamma)^(gamma-1))^(1/(1-gamma))
  
  beta <- rep(0.01,p); 
  h <- numeric(1)
  param.path <- NULL; 
  repeat
  {
    theta.fun <- function(j,gamma,tau_n,group.outer,beta,start,end)
    {
     c_j <- group.outer[j]^(1-gamma); end.outer <- cumsum(group.outer); 
     start.outer <- end.outer -group.outer + 1; index <- c(start.outer[j]:end.outer[j])
     norm <- function(k,beta,start,end)
     {
       norm_2_k <- (end[k]-start[k]+1)*sum(beta[start[k]:end[k]]^2)
     }
     sum_norm <- sum(unlist(lapply(index,norm,beta,start,end))); theta_j <- c_j*((1-gamma)/gamma/tau_n)^gamma*sum_norm^gamma
     list(theta_j=theta_j)
    }
    index.outer <- seq(1,p.outer)
    theta <- unlist(lapply(index.outer,theta.fun,gamma,tau_n,group.outer,beta,start,end))

    weight <- function(j,group.outer,theta,gamma)
    {
     c_j <-  group.outer[j]^(1-gamma);  
     if ( theta[j] < 1e-5 )  theta[j] <- 1e-5
        
     w_j <- theta[j]^(1-1/gamma)*c_j^(1/gamma)
       
     if ( w_j >= 10000) w_j <- 10000
     list(w_j=w_j)
    }
    w.tmp <- unlist(lapply(index.outer,weight,group.outer,theta,gamma))*group.outer
    w <- rep(w.tmp,group.outer)

    param <- c(n,length(group.inner),p,200); epsilon <- 1E-05;
    fit <- .C("Glasso_weight", y=as.double(y),x=as.double(t(x.t)),G=as.integer(group.inner),param=as.integer(param),
            lambda=as.double(w),epsilon=as.double(epsilon),beta=as.double(beta),diff=as.double(2))
    param.path <- cbind(param.path,fit$beta);
    diff <- max(abs(beta-fit$beta)); beta <- fit$beta; h <- h + 1;
    if ( h > 100 | diff < 1E-5 )
    {
      break
    }
  }

  b <- beta; beta <- numeric(p); 
  for ( j in 1: p.inner)
  {
     beta[start[j]:end[j]] <- solve(result[[j]]$R)%*%b[start[j]:end[j]]
  }

  list(beta=beta,b=b,diff=diff,h=h,D=fit$diff,x.t=x.t,param.path=param.path)
}