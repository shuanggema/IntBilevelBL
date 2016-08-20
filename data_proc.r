
data.proc_1 <- function(k, n, quan, g.size, n.group, mag,beta1,beta2,beta3,beta4){
  p <- g.size*n.group; group <- rep(g.size,n.group); pos_e <- cumsum(group); 
  pos_s <- pos_e -group + 1; covar <- matrix(rep(0,p*p),nrow=p)

    for ( i in 1:p)
     for ( j in 1:p)
      covar[i,j] <- mag^(abs(i-j))

    for ( j in 1:n.group)
     for ( m in pos_s[j]:pos_e[j])
      for ( l in pos_s[j]:pos_e[j])
        covar[m,l]<- mag^(abs(m-l)) 

  for (i in 1:p)
    covar[i,i]=1
    
  data1 <- generate(beta1,k,n,covar,g.size,n.group,mag,quan);
  data2 <- generate(beta2,k,n,covar,g.size,n.group,mag,quan);
  data3 <- generate(beta3,k,n,covar,g.size,n.group,mag,quan);
  data4 <- generate(beta4,k,n,covar,g.size,n.group,mag,quan);
  x1=data1$x;y1=data1$y;d1=data1$d;w1=data1$w;
  x2=data2$x;y2=data2$y;d2=data2$d;w2=data2$w;
  x3=data3$x;y3=data3$y;d3=data3$d;w3=data3$w;
  x4=data4$x;y4=data4$y;d4=data4$d;w4=data4$w;
  
  list(x1=x1,x2=x2,x3=x3,x4=x4,y1=y1,y2=y2,y3=y3,y4=y4,d1=d1,d2=d2,d3=d3,d4=d4,w1=w1,w2=w2,w3=w3,w4=w4,covar=covar)
}

data.proc <- function(k, n, quan, g.size, n.group, mag,beta1,beta2,beta3,beta4){
  p <- g.size*n.group; group <- rep(g.size,n.group); pos_e <- cumsum(group); 
  pos_s <- pos_e -group + 1; covar <- matrix(rep(0,p*p),nrow=p)

  if ( mag == 1 ){
   for ( i in 1:p)
    for ( j in 1:p){
     if ( abs(i-j) == 1)
      covar[i,j]=0.2
     else if ( abs(i-j)==2)
      covar[i,j]=0.1
    }
  }
  else if ( mag == 2 ){
   for ( i in 1:p)
    for ( j in 1:p){
     if ( abs(i-j) == 1)
      covar[i,j]=0.5
     else if ( abs(i-j)==2)
      covar[i,j]=0.25
    }
  }
  else if ( mag == 3 ){
   for ( i in 1:p)
    for ( j in 1:p){
     if ( abs(i-j) == 1)
      covar[i,j]=0.6
     else if ( abs(i-j)==2)
      covar[i,j]=0.33
    }
  }

  else if ( mag !=1 & mag !=2 & mag !=3 ){
    for ( i in 1:p)
     for ( j in 1:p)
      covar[i,j] <- 0.2^(abs(i-j))

    for ( j in 1:n.group)
     for ( m in pos_s[j]:pos_e[j])
      for ( l in pos_s[j]:pos_e[j])
        covar[m,l]<- mag^(abs(m-l)) 
  }
  for (i in 1:p)
    covar[i,i]=1
    
  data1 <- generate(beta1,k,n,covar,g.size,n.group,mag,quan);
  data2 <- generate(beta2,k,n,covar,g.size,n.group,mag,quan);
  data3 <- generate(beta3,k,n,covar,g.size,n.group,mag,quan);
  data4 <- generate(beta4,k,n,covar,g.size,n.group,mag,quan);
  x1=data1$x;y1=data1$y;d1=data1$d;w1=data1$w;
  x2=data2$x;y2=data2$y;d2=data2$d;w2=data2$w;
  x3=data3$x;y3=data3$y;d3=data3$d;w3=data3$w;
  x4=data4$x;y4=data4$y;d4=data4$d;w4=data4$w;
  
  list(x1=x1,x2=x2,x3=x3,x4=x4,y1=y1,y2=y2,y3=y3,y4=y4,d1=d1,d2=d2,d3=d3,d4=d4,w1=w1,w2=w2,w3=w3,w4=w4,covar=covar)
}

data.proc.3 <- function(k, n, quan, g.size, n.group, mag,beta1,beta2,beta3){
  p <- g.size*n.group; group <- rep(g.size,n.group); pos_e <- cumsum(group); 
  pos_s <- pos_e -group + 1; covar <- matrix(rep(0,p*p),nrow=p)

  if ( mag == 1 ){
    for ( i in 1:p)
     for ( j in 1:p)
      covar[i,j] <- 0.2^(abs(i-j))

    for ( j in 1:n.group)
     for ( m in pos_s[j]:pos_e[j])
      for ( l in pos_s[j]:pos_e[j]){
       if ( abs(m-l) == 1)
         covar[m,l]=0.2
       else if ( abs(m-l)==2)
         covar[m,l]=0.1
       else 
         covar[m,l]=0
      }
  }
  else if ( mag == 2 ){
    for ( i in 1:p)
     for ( j in 1:p)
      covar[i,j] <- 0.2^(abs(i-j))

    for ( j in 1:n.group)
     for ( m in pos_s[j]:pos_e[j])
      for ( l in pos_s[j]:pos_e[j]){
       if ( abs(m-l) == 1)
         covar[m,l]=0.5
       else if ( abs(m-l)==2)
         covar[m,l]=0.25
       else 
         covar[m,l]=0
      }
  }
  else if ( mag == 3 ){
    for ( i in 1:p)
     for ( j in 1:p)
      covar[i,j] <- 0.2^(abs(i-j))

    for ( j in 1:n.group)
     for ( m in pos_s[j]:pos_e[j])
      for ( l in pos_s[j]:pos_e[j]){
       if ( abs(m-l) == 1)
         covar[m,l]=0.6
       else if ( abs(m-l)==2)
         covar[m,l]=0.33
       else 
         covar[m,l]=0
      }
  }

  else if ( mag !=1 & mag !=2 & mag !=3 ){
    for ( i in 1:p)
     for ( j in 1:p)
      covar[i,j] <- 0.2^(abs(i-j))

    for ( j in 1:n.group)
     for ( m in pos_s[j]:pos_e[j])
      for ( l in pos_s[j]:pos_e[j])
        covar[m,l]<- mag^(abs(m-l)) 
  }
  for (i in 1:p)
    covar[i,i]=1
    
  data1 <- generate(beta1,k,n,covar,g.size,n.group,mag,quan);
  data2 <- generate(beta2,k,n,covar,g.size,n.group,mag,quan);
  data3 <- generate(beta3,k,n,covar,g.size,n.group,mag,quan);
  x1=data1$x;y1=data1$y;d1=data1$d;w1=data1$w;
  x2=data2$x;y2=data2$y;d2=data2$d;w2=data2$w;
  x3=data3$x;y3=data3$y;d3=data3$d;w3=data3$w;

  list(x1=x1,x2=x2,x3=x3,y1=y1,y2=y2,y3=y3,d1=d1,d2=d2,d3=d3,w1=w1,w2=w2,w3=w3,covar=covar)
}

data.proc.4 <- function(k, n, quan, group, mag,beta1,beta2,beta3){
  p <- sum(group); pos_e <- cumsum(group);  n.group <- length(group); 
  pos_s <- pos_e -group + 1; covar <- matrix(rep(0,p*p),nrow=p)

  if ( mag == 1 ){
    for ( i in 1:p)
     for ( j in 1:p)
      covar[i,j] <- 0.2^(abs(i-j))

    for ( j in 1:n.group)
     for ( m in pos_s[j]:pos_e[j])
      for ( l in pos_s[j]:pos_e[j]){
       if ( abs(m-l) == 1)
         covar[m,l]=0.2
       else if ( abs(m-l)==2)
         covar[m,l]=0.1
       else 
         covar[m,l]=0
      }
  }
  else if ( mag == 2 ){
    for ( i in 1:p)
     for ( j in 1:p)
      covar[i,j] <- 0.2^(abs(i-j))

    for ( j in 1:n.group)
     for ( m in pos_s[j]:pos_e[j])
      for ( l in pos_s[j]:pos_e[j]){
       if ( abs(m-l) == 1)
         covar[m,l]=0.5
       else if ( abs(m-l)==2)
         covar[m,l]=0.25
       else 
         covar[m,l]=0
      }
  }
  else if ( mag == 3 ){
    for ( i in 1:p)
     for ( j in 1:p)
      covar[i,j] <- 0.2^(abs(i-j))

    for ( j in 1:n.group)
     for ( m in pos_s[j]:pos_e[j])
      for ( l in pos_s[j]:pos_e[j]){
       if ( abs(m-l) == 1)
         covar[m,l]=0.6
       else if ( abs(m-l)==2)
         covar[m,l]=0.33
       else 
         covar[m,l]=0
      }
  }

  else if ( mag !=1 & mag !=2 & mag !=3 ){
    for ( i in 1:p)
     for ( j in 1:p)
      covar[i,j] <- 0.2^(abs(i-j))

    for ( j in 1:n.group)
     for ( m in pos_s[j]:pos_e[j])
      for ( l in pos_s[j]:pos_e[j])
        covar[m,l]<- mag^(abs(m-l)) 
  }
  for (i in 1:p)
    covar[i,i]=1

  data1 <- generate(beta1,k,n,covar,group,mag,quan);
  data2 <- generate(beta2,k,n,covar,group,mag,quan);
  data3 <- generate(beta3,k,n,covar,group,mag,quan);
  x1=data1$x;y1=data1$y;d1=data1$d;w1=data1$w;
  x2=data2$x;y2=data2$y;d2=data2$d;w2=data2$w;
  x3=data3$x;y3=data3$y;d3=data3$d;w3=data3$w;

  list(x1=x1,x2=x2,x3=x3,y1=y1,y2=y2,y3=y3,d1=d1,d2=d2,d3=d3,w1=w1,w2=w2,w3=w3,covar=covar)
}

process_data <- function(X,Y,D){
  p <- ncol(X); n <- nrow(X);index <- seq(1,p)

  Xs<-X[order(Y),];D<-D[order(Y)];Y<-sort(Y)

  weight <- numeric(n);weight[1]=D[1]/n;
  for ( i in 2:n )
  {
     tmp = 1
     for ( j in 1: (i-1) )
      tmp = tmp*((n-j)/(n-j+1))^D[j]

     weight[i]=D[i]/(n-i+1)*tmp
  }

  list(X=Xs,Y=Y,weight=weight)
}

extract.x <- function(i,tmp){
  x <- tmp[[i]]$x.s
}

extract.d <- function(i,tmp){
  x <- tmp[[i]]$d
}

PCA.x <- function(x,perc,group){
  index <- seq(1,length(group));
  tmp <- lapply(index,pca.group,x,group,perc); n <- nrow(x); 
  x.new <- matrix(unlist(lapply(index,extract.x,tmp)),nrow=n);
  group.new <- unlist(lapply(index,extract.d,tmp));
  list(x.new=x.new,group.new=group.new)
}

PCA.x.f <- function(x,d,group){
  index <- seq(1,length(group));
  tmp <- lapply(index,pca.group.f,x,group,d); n <- nrow(x); 
  x.new <- matrix(unlist(lapply(index,extract.x,tmp)),nrow=n);
  group.new <- rep(d,length(group))
  list(x.new=x.new,group.new=group.new)
}

standard <- function (x) {
   p <- ncol(x)
   n <- nrow(x)
   x.mean <- matrix(rep(apply(x,2,mean),n),n,p,byrow=T)
   x.std <- t(t((x-x.mean))/(apply(x,2,sd)*(n-1)^0.5)*n^0.5)
}

weighted.standard.x <- function (x,w) {
   p <- ncol(x)
   n <- nrow(x)
   x.wmean <- matrix(rep(apply(x,2,weighted.mean,w),n),n,p,byrow=T)
   x.std <- (x-x.wmean)*w^0.5#/(apply(x,2,sd)*(n-1)^0.5)*n^0.5)
}
weighted.standard.y <- function (y,w) {
   n <- length(y)
   y.wmean <- weighted.mean(y,w)
   y.std <- (y-y.wmean)*w^0.5
}

transform <- function(x,y,w)
{
  x.sd=weighted.standard.x(x,w); 
  y.sd=weighted.standard.y(y,w);  
  list(x.sd=x.sd,y.sd=y.sd)
}

reweight <- function(x1,x2,x3,x4,y1,y2,y3,y4,w1,w2,w3,w4){
  w1.s <- w1[w1!=0];x1.s <- x1[w1!=0,];y1.s=y1[w1!=0];n1=length(w1.s);
  w2.s <- w2[w2!=0];x2.s <- x2[w2!=0,];y2.s=y2[w2!=0];n2=length(w2.s);
  w3.s <- w3[w3!=0];x3.s <- x3[w3!=0,];y3.s=y3[w3!=0];n3=length(w3.s);
  w4.s <- w4[w4!=0];x4.s <- x4[w4!=0,];y4.s=y4[w4!=0];n4=length(w4.s);

  data11<-transform(x1.s,log(y1.s),w1.s);
  data22<-transform(x2.s,log(y2.s),w2.s);
  data33<-transform(x3.s,log(y3.s),w3.s);
  data44<-transform(x4.s,log(y4.s),w4.s);

  x1.std<-standard(data11$x.sd);x2.std<-standard(data22$x.sd);x3.std<-standard(data33$x.sd);x4.std<-standard(data44$x.sd);
  y1.new<-data11$y.sd;y2.new<-data22$y.sd;
  y3.new<-data33$y.sd;y4.new<-data44$y.sd;

  y <- c(y1.new,y2.new,y3.new,y4.new);
  n.t<- c(n1,n2,n3,n4); 
  list(x1.sd=x1.std,x2.sd=x2.std,x3.sd=x3.std,x4.sd=x4.std,y=y,n.t=n.t)
}

reweight2 <- function(x1,x2,y1,y2,w1,w2){
  w1.s <- w1[w1!=0];x1.s <- x1[w1!=0,];y1.s=y1[w1!=0];n1=length(w1.s);
  w2.s <- w2[w2!=0];x2.s <- x2[w2!=0,];y2.s=y2[w2!=0];n2=length(w2.s);

  data11<-transform(x1.s,log(y1.s),w1.s);
  data22<-transform(x2.s,log(y2.s),w2.s);

  x1.std<-standard(data11$x.sd);x2.std<-standard(data22$x.sd);
  y1.new<-data11$y.sd;y2.new<-data22$y.sd;

  y <- c(y1.new,y2.new);
  n.t<- c(n1,n2); 
  list(x1.sd=x1.std,x2.sd=x2.std,y=y,n.t=n.t)
}

reweight3 <- function(x1,x2,x3,y1,y2,y3,w1,w2,w3){
  w1.s <- w1[w1!=0];x1.s <- x1[w1!=0,];y1.s=y1[w1!=0];n1=length(w1.s);
  w2.s <- w2[w2!=0];x2.s <- x2[w2!=0,];y2.s=y2[w2!=0];n2=length(w2.s);
  w3.s <- w3[w3!=0];x3.s <- x3[w3!=0,];y3.s=y3[w3!=0];n3=length(w3.s);
  y.s <- c(y1.s,y2.s,y3.s);
  
  data11<-transform(x1.s,log(y1.s),w1.s);
  data22<-transform(x2.s,log(y2.s),w2.s);
  data33<-transform(x3.s,log(y3.s),w3.s);

  x1.std<-standard(data11$x.sd);x2.std<-standard(data22$x.sd);x3.std<-standard(data33$x.sd);
  y1.new<-data11$y.sd;y2.new<-data22$y.sd;y3.new<-data33$y.sd;

  y <- c(y1.new,y2.new,y3.new);
  n.t<- c(n1,n2,n3); 
  list(y.s=y.s,x1.sd=x1.std,x2.sd=x2.std,x3.sd=x3.std,y=y,n.t=n.t)
}


process <- function(x1.sd,x2.sd,x3.sd,x4.sd,n.t,group1,group2,group3,group4)
{
  ntotal<- sum(n.t); n.group <- length(group1); pe<-cumsum(n.t);ps<-pe-n.t+1; 
  group.new <- as.vector(rbind(group1,group2,group3,group4)); pos_e_n <- cumsum(group.new); pos_s_n <- pos_e_n -group.new + 1;
  p.t <- c(ncol(x1.sd),ncol(x2.sd),ncol(x3.sd),ncol(x4.sd));
  
  x <- matrix(rep(0,ntotal*sum(p.t)),nrow=ntotal)
 
  for( j in 1:n.group )
    for ( k in 1:4 ){
      thiscommand<- paste("x.tmp <- x", k, ".sd",sep="");
      eval(parse(text=thiscommand));
      thiscommand<- paste("group.tmp <- group", k,sep="");
      eval(parse(text=thiscommand));
      pe.tmp <- cumsum(group.tmp); ps.tmp <- pe.tmp -group.tmp + 1;
      ind <- (j-1)*4+k;
      x[ps[k]:pe[k],pos_s_n[ind]:pos_e_n[ind]] <- x.tmp[,ps.tmp[j]:pe.tmp[j]];
    }

  x0=matrix(rep(0,ntotal*4),nrow=ntotal)
  for ( k in 1:4 ){
     start <- ps[k];end <- pe[k]
     x0[start:end,k] <- rep(1,end-start+1)
  }
  list(x0=x0,x=x,ps=ps,pe=pe,group.inner=group.new)
}

process2 <- function(x1.sd,x2.sd,n.t,group1,group2){
  ntotal<- sum(n.t); n.group <- length(group1); pe<-cumsum(n.t);ps<-pe-n.t+1; 
  group.new <- as.vector(rbind(group1,group2)); pos_e_n <- cumsum(group.new); pos_s_n <- pos_e_n -group.new + 1;
  p.t <- c(ncol(x1.sd),ncol(x2.sd));
  
  x <- matrix(rep(0,ntotal*sum(p.t)),nrow=ntotal)
 
  for( j in 1:n.group )
    for ( k in 1:2 ){
      thiscommand<- paste("x.tmp <- x", k, ".sd",sep="");
      eval(parse(text=thiscommand));
      thiscommand<- paste("group.tmp <- group", k,sep="");
      eval(parse(text=thiscommand));
      pe.tmp <- cumsum(group.tmp); ps.tmp <- pe.tmp -group.tmp + 1;
      ind <- (j-1)*2+k;
      x[ps[k]:pe[k],pos_s_n[ind]:pos_e_n[ind]] <- x.tmp[,ps.tmp[j]:pe.tmp[j]];
    }

  x0=matrix(rep(0,ntotal*2),nrow=ntotal)
  for ( k in 1:2 ){
     start <- ps[k];end <- pe[k]
     x0[start:end,k] <- rep(1,end-start+1)
  }
  list(x0=x0,x=x,ps=ps,pe=pe,group.inner=group.new)
}

process3 <- function(x1.sd,x2.sd,x3.sd,n.t,group1,group2,group3){
  ntotal<- sum(n.t); n.group <- length(group1); pe<-cumsum(n.t);ps<-pe-n.t+1; 
  group.new <- as.vector(rbind(group1,group2,group3)); pos_e_n <- cumsum(group.new); pos_s_n <- pos_e_n -group.new + 1;
  p.t <- c(ncol(x1.sd),ncol(x2.sd),ncol(x3.sd));
  
  x <- matrix(rep(0,ntotal*sum(p.t)),nrow=ntotal)
 
  for( j in 1:n.group )
    for ( k in 1:3 ){
      thiscommand<- paste("x.tmp <- x", k, ".sd",sep="");
      eval(parse(text=thiscommand));
      thiscommand<- paste("group.tmp <- group", k,sep="");
      eval(parse(text=thiscommand));
      pe.tmp <- cumsum(group.tmp); ps.tmp <- pe.tmp -group.tmp + 1;
      ind <- (j-1)*3+k;
      x[ps[k]:pe[k],pos_s_n[ind]:pos_e_n[ind]] <- x.tmp[,ps.tmp[j]:pe.tmp[j]];
    }

  x0=matrix(rep(0,ntotal*3),nrow=ntotal)
  for ( k in 1:3 ){
     start <- ps[k];end <- pe[k]
     x0[start:end,k] <- rep(1,end-start+1)
  }
  list(x0=x0,x=x,ps=ps,pe=pe,group.inner=group.new)
}

sp <- function(X,Y,group.inner,group.outer,lambda.seq,gamma.seq){
  n.lambda <- length(lambda.seq); n.gamma <- length(gamma.seq); BETA <- vector("list",n.gamma); 
  for ( k in 1: n.gamma)
    BETA[[k]] <- matrix(rep(0,n.lambda*ncol(X)),nrow=n.lambda)
    for ( j in 1: n.lambda){
       fm <- com_bridge(X,Y,group.inner,group.outer,lambda.seq[j],gamma.seq[k]);
       BETA[[k]][j,] <-  fm$beta;
    }
  list(BETA=BETA)
}

n.nonzero <- function(BETA,group){
  K <-  length(BETA); J <- nrow(BETA[[1]]);
  count <- matrix(rep(0,K*J),nrow=K)
  for ( k in 1: K)
    for ( j in 1: J){
       tmp <- BETA[[k]][j,];
       cou <- identify.group(tmp,group)$ident.g;
       if ( is.null(cou) == T )
          count[k,j] <- 0
       else if ( is.null(cou) != T ) 
          count[k,j] <- length(cou)
    }
  list(count=count)
}

#identify.group <- function(beta,group){
#   p.g <- length(group); pe<-cumsum(group);ps<-pe-group+1;
#   ident.g <- NULL;
#   for ( j in 1: p.g){
#     if (beta[ps[j]]!=0) ident.g <- c(ident.g,j)
#   }
#   list(ident.g=ident.g)
#}

identify.group <- function(beta,group){
   p.g <- length(group); pe<-cumsum(group);ps<-pe-group+1;
   ident.g <- NULL;
   for ( j in 1: p.g){
     tmp = beta[ps[j]:pe[j]];
     if ( sum(tmp^2)!=0) ident.g <- c(ident.g,j)
   }
   list(ident.g=ident.g)
}

n.nonzero.ind <- function(BETA){
  K <-  length(BETA); J <- nrow(BETA[[1]]);
  count <- matrix(rep(0,K*J),nrow=K)
  for ( k in 1: K)
    for ( j in 1: J){
       tmp <- BETA[[k]][j,];
       count[k,j] <- length(tmp[tmp!=0])
    }
  list(count=count)
}

##########correlation#################
cor_x <- function(x,group){
  pos_e <- cumsum(group); pos_s <- pos_e -group + 1; 
  p.group <- length(group); c.a <- numeric(p.group);
  for ( j in 1 : p.group){
    c.x <- cor(x[,pos_s[j]:pos_e[j]]); p.tmp <- pos_e[j]-pos_s[j]+1; sum.c <- numeric(1);
    for ( k in 1: (p.tmp-1))
       sum.c <- sum.c + c.x[k+1,k]
    c.a[j] <- sum.c / (p.tmp -1 )
  }
  list(rho = mean(c.a))
}

######identify outer group#############
identify.outer <- function(sel.group,group.outer){
   pe <- cumsum(group.outer); ps <- pe -group.outer + 1; p.outer <-length(group.outer); ind <- numeric(0);
   n.s <- length(sel.group)
   if ( n.s == 0 ) s.g.outer <- 0
   else if ( n.s >0 ) {
     for ( j in 1:n.s){
        k <- 1;
        repeat{ 
          mark <- 0     
          if ( sel.group[j] <= pe[k] & sel.group[j] >= ps[k] ) mark <- 1;
          if ( mark == 1 ) ind <- c(ind,k);
          k <- k + 1;
            if ( k > p.outer ){
              break
            }
        }  
     }
   }
   list(ind = ind)
}

#####degree of freedom###############
df.fun <- function(betaf,group.inner,group.outer,x,y)
{  
  pos_e = cumsum(group.inner); pos_s = pos_e -group.inner + 1;
  pe <- cumsum(group.outer); ps = pe -group.outer + 1;
  p.inner = length(group.inner); p.outer <- length(group.outer); df = 0;
  sel.group <- identify.group(betaf,group.inner)$ident.g;
  s.g.outer <- unique(identify.outer(sel.group,group.outer)$ind);
  if ( length(s.g.outer) == 0 ) df = 0
  else if ( length(s.g.outer) > 0 ){
    for ( j in 1: length(s.g.outer) ){
      index <- c(pos_s[ps[s.g.outer[j]]]:pos_e[pe[s.g.outer[j]]]); p.j <- pe[s.g.outer[j]]-ps[s.g.outer[j]]+1
      fit <- lm(y~x[,index]); 
      ls.norm <- sum(((fit$coefficients)[-1])^2)^0.5
      est.norm <- sum((betaf[index])^2)^0.5
      df <- df + 1 + (p.j-1)*est.norm/ls.norm
    }
  }
  list(df=df)
}

df.eval <- function(BETA,group.inner,group.outer,x,y){
  K <-  length(BETA); J <- nrow(BETA[[1]]);
  count <- matrix(rep(0,K*J),nrow=K)
  for ( k in 1: K)
    for ( j in 1: J){
       tmp <- BETA[[k]][j,];
       count[k,j] <- df.fun(tmp,group.inner,group.outer,x,y)$df
    }
  list(count=count)
}

