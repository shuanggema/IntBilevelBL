source("GMCP.r")

df.fun <- function(betaf,group,x,y) {  
  pos_e = cumsum(group); pos_s = pos_e -group + 1; p = length(group); df = 0
  for ( i in 1:length(group) ){
    if ( betaf[pos_s[i]]!=0 ) {
      if ( group[i] != 1 ) {
        fit = lm(y~x[,pos_s[i]:pos_e[i]])
        ls.norm = sum(((fit$coefficients)[-1])^2)^0.5
        est.norm= sum((betaf[pos_s[i]:pos_e[i]])^2)^0.5
        df = df + 1 + (group[i]-1)*est.norm/ls.norm
      }
      else if ( group[i] == 1)
        df = df + 1
    }
  }
  list(df=df)
}

lammax <- function(x,y,group)
{
  p <- length(group)
  n <- dim(x)[1]
  totalp <- ncol(x)
  
  #calculate lambda max
  end = cumsum(group)
  start = end -group + 1
  inner <- rep(0,p)
  for ( j in 1:p )
  {
    if ( group[j] != 1 )
       inner[j] <- sum(colSums(x[,start[j]:end[j]]*as.vector(y)/n)^2)^0.5/group[j]^0.5
    else if ( group[j] == 1 )
       inner[j] <- abs(sum(x[,start[j]:end[j]]*as.vector(y)/n))
  }
  lambda.max = max(inner)
  list(lambda.max = lambda.max)
}

BIC.gMCP <- function (x,y,n.t,group,epsilon,n.step,gamma,g){ 
         
  n <- nrow(x); p <- ncol(x); p.group <- length(group); 
  pe <- cumsum(n.t); ps <- pe -n.t + 1; m <- length(n.t);
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

  lambda.max = lammax(x.t,y,group)$lambda.max
  lambda.min <- lambda.max*epsilon
  ss <- (log(lambda.max)-log(lambda.min))/(n.step-1)
  lambda.seq <- numeric(n.step)
  for ( i in 1:n.step ){
     lambda.seq[i] <- exp(log(lambda.max)-ss*(i-1))
  }

  lh <-numeric(n.step); penalty <- numeric(n.step); mm <- numeric(n.step); bic <- numeric(n.step); df <- numeric(n.step)

  for ( i in 1: n.step ){
     lambda <- lambda.seq[i]; fm <- GMCP_quant(x,y,group,lambda,gamma)
     beta.t <- fm$beta; y.h <- x%*%as.vector(beta.t); lh[i]<- n*log(sum((y- y.h)^2)/n)
     mm[i] <- length(beta.t[abs(beta.t)!=0]); tmp <- df.fun(beta.t,group,x,y); df[i] <- tmp$df
     penalty[i] <- df[i]*(log(n)+2*g*log(p))    #+2*log(d) sum(abs(betam))#
     bic[i] <- lh[i]+penalty[i]
   }
   fm <- GMCP_quant(x,y,group,min(lambda.seq[bic==min(bic)]),gamma)
   beta.opt <- fm$beta;
   list(bic=bic,lambda.seq=lambda.seq,lh=lh,df=df,penalty=penalty,mm=mm,beta.opt=beta.opt)
}

BIC.gMCP.ind <- function (x,y,group,epsilon,n.step,gamma,g){ 
         
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

  lambda.max = lammax(x.t,y,group)$lambda.max
  lambda.min <- lambda.max*epsilon
  ss <- (log(lambda.max)-log(lambda.min))/(n.step-1)
  lambda.seq <- numeric(n.step)
  for ( i in 1:n.step ){
     lambda.seq[i] <- exp(log(lambda.max)-ss*(i-1))
  }

  lh <-numeric(n.step); penalty <- numeric(n.step); mm <- numeric(n.step); bic <- numeric(n.step); df <- numeric(n.step)

  for ( i in 1: n.step ){
     lambda <- lambda.seq[i]; fm <- GMCP_quant(x,y,group,lambda,gamma)
     beta.t <- fm$beta; y.h <- x%*%as.vector(beta.t); lh[i]<- n*log(sum((y- y.h)^2)/n)
     mm[i] <- length(beta.t[abs(beta.t)!=0]); tmp <- df.fun(beta.t,group,x,y); df[i] <- tmp$df
     penalty[i] <- df[i]*(log(n)+2*g*log(p))    #+2*log(d) sum(abs(betam))#
     bic[i] <- lh[i]+penalty[i]
   }
   fm <- GMCP_quant(x,y,group,min(lambda.seq[bic==min(bic)]),gamma)
   beta.opt <- fm$beta;
   list(bic=bic,lambda.seq=lambda.seq,lh=lh,df=df,penalty=penalty,mm=mm,beta.opt=beta.opt)
}

BIC.comp.grlasso <- function (x,y,group,epsilon,n.step,gamma,g,n.t){ 
         
  n <- nrow(x); p <- ncol(x); p.group <- length(group); 
  end <- cumsum(group); start <- end -group + 1
  m <- length(n.t)
  
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

  lambda.max = lammax(x.t,y,group)$lambda.max
  lambda.min <- lambda.max*epsilon
  ss <- (log(lambda.max)-log(lambda.min))/(n.step-1)
  lambda.seq <- numeric(n.step)
  for ( i in 1:n.step ){
     lambda.seq[i] <- exp(log(lambda.max)-ss*(i-1))
  }

  lh <-numeric(n.step); penalty <- numeric(n.step); mm <- numeric(n.step); bic <- numeric(n.step); df <- numeric(n.step)

  for ( i in 1: n.step ){
     lambda <- lambda.seq[i]; fm <- GMCP_quant(x,y,group,lambda,gamma)
     beta.t <- fm$beta; y.h <- x%*%as.vector(beta.t); lh[i]<- n*log(sum((y- y.h)^2)/n)
     mm[i] <- length(beta.t[abs(beta.t)!=0]); tmp <- df.fun(beta.t,group,x,y); df[i] <- tmp$df
     penalty[i] <- df[i]*(sum(log(n.t))/m+2*g*log(p/m))    #+2*log(d) sum(abs(betam))#
     bic[i] <- lh[i]+penalty[i]
   }
   fm <- GMCP_quant(x,y,group,min(lambda.seq[bic==min(bic)]),gamma)
   beta.opt <- fm$beta;
   list(bic=bic,lambda.seq=lambda.seq,lh=lh,df=df,penalty=penalty,mm=mm,beta.opt=beta.opt)
}
