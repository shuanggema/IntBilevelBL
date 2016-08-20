source("CallC_com_bridge.r")

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

BIC.comp.bridge <- function (X,Y,n.t,group.inner,group.outer,epsilon,n.step,gamma,lambda.max,g){ 
         
  n <- nrow(X); p <- ncol(X); p.inner <- length(group.inner); p.outer <- length(group.outer);
  pe <- cumsum(n.t); ps <- pe -n.t + 1; m <- length(n.t);

  lambda.min <- lambda.max*epsilon
  ss <- (log(lambda.max)-log(lambda.min))/(n.step-1)
  lambda.seq <- numeric(n.step)
  for ( i in 1:n.step ){
     lambda.seq[i] <- exp(log(lambda.max)-ss*(i-1))
  }

  lh <-numeric(n.step); penalty <- numeric(n.step); mm <- numeric(n.step); bic <- numeric(n.step); df <- numeric(n.step)

  for ( i in 1: n.step ){
     lambda <- lambda.seq[i]; fm <- com_bridge(X,Y,group.inner,group.outer,lambda,gamma);
     beta.t <- fm$beta; Y.h <- X%*%as.vector(beta.t); lh[i]<- n*log(sum((Y- Y.h)^2)/n)
     mm[i] <- length(beta.t[abs(beta.t)!=0]); tmp <- df.fun(beta.t,group.inner,X,Y); df[i] <- tmp$df
     penalty[i] <- df[i]*(log(n)+2*g*log(p))    #+2*log(d) sum(abs(betam))#
     bic[i] <- lh[i]+penalty[i]
   }
   fm <- com_bridge(X,Y,group.inner,group.outer,min(lambda.seq[bic==min(bic)]),gamma);
   beta.opt <- fm$beta;
   list(bic=bic,lambda.seq=lambda.seq,lh=lh,df=df,penalty=penalty,mm=mm,beta.opt=beta.opt)
}

  
