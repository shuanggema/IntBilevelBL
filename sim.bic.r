match.group <- function(b,matched,n.group)
{
  if (matched == 1|matched==4){
    #50% not matched
    a <- c(1:6,7,11,15,16,20,24)
  }
  else if (matched ==2|matched==5 ){
    #25% not matched
    a <- c(1:6,7,11,15,16,17,18)
  }
  else if (matched ==3|matched==6){
    #0 % not matched
    a <- c(1:6,13:18)
  }
  matched=match(a,b,nomatch=-1)
  tp <- length(matched[matched!=-1])
  fn <- length(a) -tp
  fp <- length(b) - tp
  tn <- n.group - tp -fn -fp
  list(tp=tp,fn=fn,fp=fp,tn=tn)
}

sim <- function(num,n.group,g.size,matched=1,quan,mag,k,g,perc,replic,epsilon,n.step,lambda.max,gamma.seq,example)
{
  tp <- numeric();tn <- numeric();fp <- numeric();fn <- numeric();sse.seq <- numeric();
  lambda.opt.seq <- numeric();lambda.min.seq <- numeric();lambda.max.seq <- numeric();
  group <- rep(g.size,n.group);
  
  for ( i in 1:replic)
  {
    if ( matched == 1){
    #50% not matched  : value =1
    beta1 <- c(rep(0.15,g.size),rep(0.1,g.size),rep(0.15,g.size),rep(0,2*g.size),rep(0.1,g.size),rep(0,g.size*(n.group-6)))
    beta2 <- c(rep(0.15,g.size),rep(0.1,g.size),rep(0,g.size),rep(0.15,g.size),rep(0,2*g.size),rep(0.1,g.size),rep(0,g.size*(n.group-7)))
    beta3 <- c(rep(0.1,g.size),rep(0.15,g.size),rep(0,2*g.size),rep(0.15,g.size),rep(0,2*g.size),rep(0.1,g.size),rep(0,g.size*(n.group-8)))
    }
    else if ( matched == 2){
    #25% not matched  : value =2
    beta1 <- c(rep(0.15,g.size),rep(0.1,g.size),rep(0.15,g.size),rep(0,2*g.size),rep(0.1,g.size),rep(0,g.size*(n.group-6)))
    beta2 <- c(rep(0.15,g.size),rep(0.1,g.size),rep(0,g.size),rep(0.15,g.size),rep(0,g.size),rep(0.1,g.size),rep(0,g.size*(n.group-6)))
    beta3 <- c(rep(0.1,g.size),rep(0.15,g.size),rep(0,2*g.size),rep(0.15,g.size),rep(0.1,g.size),rep(0,g.size*(n.group-6)))
    }
    else if ( matched ==3){
    #0 % not matched  : value =3
    beta1 <- c(rep(0.15,g.size),rep(0.1,g.size),rep(0,2*g.size),rep(0.15,g.size),rep(0.1,g.size),rep(0,g.size*(n.group-6)))
    beta2 <- c(rep(0.15,g.size),rep(0.1,g.size),rep(0,2*g.size),rep(0.15,g.size),rep(0.1,g.size),rep(0,g.size*(n.group-6)))
    beta3 <- c(rep(0.1,g.size),rep(0.15,g.size),rep(0,2*g.size),rep(0.15,g.size),rep(0.1,g.size),rep(0,g.size*(n.group-6)))
    }
    else if ( matched == 4){
    #50% not matched  : value =4
    beta1 <- c(rep(0.15,g.size-1),0,0,rep(0.1,g.size-1),rep(0.15,g.size-1),0,rep(0,2*g.size),0,rep(0.1,g.size-1),rep(0,g.size*(n.group-6)))
    beta2 <- c(rep(0.15,g.size-1),0,0,rep(0.1,g.size-1),rep(0,g.size),rep(0.15,g.size-1),0,rep(0,2*g.size),0,rep(0.1,g.size-1),rep(0,g.size*(n.group-7)))
    beta3 <- c(0,rep(0.1,g.size-1),rep(0.15,g.size-1),0,rep(0,2*g.size),rep(0.15,g.size-1),0,rep(0,2*g.size),0,rep(0.1,g.size-1),rep(0,g.size*(n.group-8)))
    }
    else if ( matched == 5){
    #25% not matched  : value =5
    beta1 <- c(rep(0.15,g.size-1),0,0,rep(0.1,g.size-1),rep(0.15,g.size-1),0,rep(0,2*g.size),0,rep(0.1,g.size-1),rep(0,g.size*(n.group-6)))
    beta2 <- c(rep(0.15,g.size-1),0,0,rep(0.1,g.size-1),rep(0,g.size),rep(0.15,g.size-1),0,rep(0,g.size),0,rep(0.1,g.size-1),rep(0,g.size*(n.group-6)))
    beta3 <- c(0,rep(0.1,g.size-1),rep(0.15,g.size-1),0,rep(0,2*g.size),rep(0.15,g.size-1),0,0,rep(0.1,g.size-1),rep(0,g.size*(n.group-6)))
    }
    else if ( matched ==6){
    #0 % not matched  : value =6; not all beta in a cluster are non-zero
    beta1 <- c(rep(0.15,g.size-1),0,0,rep(0.1,g.size-1),rep(0,2*g.size),rep(0.15,g.size-1),0,0,rep(0.1,g.size-1),rep(0,g.size*(n.group-6)))
    beta2 <- c(rep(0.15,g.size-1),0,0,rep(0.1,g.size-1),rep(0,2*g.size),rep(0.15,g.size-1),0,0,rep(0.1,g.size-1),rep(0,g.size*(n.group-6)))
    beta3 <- c(0,rep(0.1,g.size-1),rep(0.15,g.size-1),0,rep(0,2*g.size),rep(0.15,g.size-1),0,0,rep(0.1,g.size-1),rep(0,g.size*(n.group-6)))
    }
    
    source("/space/jliu20/Compound_bridge/data_proc.r");
    data.p <- data.proc.3(k, num, quan, g.size, n.group, mag,beta1,beta2,beta3)
    x1 <- data.p$x1; y1 <- data.p$y1; d1 <- data.p$d1; w1 <- data.p$w1;
    x2 <- data.p$x2; y2 <- data.p$y2; d2 <- data.p$d2; w2 <- data.p$w2;
    x3 <- data.p$x3; y3 <- data.p$y3; d3 <- data.p$d3; w3 <- data.p$w3;
    #rho <- cor_x(rbind(x1,x2,x3,x4),group)$rho
    
    source("/space/jliu20/Compound_bridge/PCA.r"); 
    #perc <- 0.6;
    t1<- PCA.x(x1,perc,group); t2<- PCA.x(x2,perc,group); t3<- PCA.x(x3,perc,group);  
    #t1<- PCA.x.f(x1,d,group); t2<- PCA.x.f(x2,d,group); t3<- PCA.x.f(x3,d,group);  
    x1<- t1$x.new ; x2 <- t2$x.new; x3 <- t3$x.new;
    group1 <- t1$group.new; group2 <- t2$group.new; group3 <- t3$group.new; 

    data.r <- reweight3(x1,x2,x3,y1,y2,y3,w1,w2,w3)
    x1.sd <- data.r$x1.sd; x2.sd <- data.r$x2.sd; x3.sd <- data.r$x3.sd; n.t <- data.r$n.t; y <- data.r$y;
    data.f <- process3(x1.sd,x2.sd,x3.sd,n.t,group1,group2,group3)
    
    x <- data.f$x; x0 <- data.f$x0; ps <- data.f$ps; pe <- data.f$pe; group.inner <- data.f$group.inner
    totalp <- ncol(x); m <- 3;
    group.outer <- rep(m,length(group.inner)/m);
    
    #fm <- com_bridge(x,y,group.inner,group.outer,0.0011,gamma.seq);
    #fm$beta[fm$beta!=0]
    #sel.group <- identify.group(fm$beta,group.inner)$ident.g;

    source("/space/jliu20/Compound_bridge/BIC.r")
    BIC.t <- BIC.comp.bridge(x,y,n.t,group.inner,group.outer,epsilon,n.step,gamma.seq,lambda.max,g)
    bic <- BIC.t$bic; lambda.seq <- BIC.t$lambda.seq; lambda.opt <- min(lambda.seq[bic==min(bic)])
    
    #source("/space/jliu20/Compound_bridge/cv_mult.r")
    #n.fold <- 5; #epsilon <- 0.05; n.step <- 40; lambda.max <- 0.0016;gamma.seq <- 0.9
    #cv <- cv.optim(cl,x,y,n.t,group.inner,group.outer,n.fold,epsilon,n.step,gamma.seq,lambda.max)
    #PE <- cv$PE; lambda.seq <- cv$lambda.seq; lambda.opt <- lambda.seq[PE==min(PE)]#opt <- min.opt(PE);

    #ss <- sp(x,y,group.inner,group.outer,lambda.seq,gamma.seq); n <- nrow(x);p <- ncol(x); 
    #ki <- log(p)/log(n); g <- 1-1/2/ki; 
    #mag.thresh <- mag; #mag.thresh[mag.thresh <0.35] <- 0;
    #cri <- df.eval(ss$BETA,group.inner,group.outer,x,y)$count*mag.thresh*(log(n)/n+2*g*log(p)/n)+log(PE/n);  
    #lambda.opt <- lambda.seq[cri==min(cri)] 

    lambda.opt.seq <- c(lambda.opt.seq,lambda.opt)
    lambda.min.seq <- c(lambda.min.seq,min(lambda.seq))
    lambda.max.seq <- c(lambda.max.seq,max(lambda.seq))

    fm <- com_bridge(x,y,group.inner,group.outer,lambda.opt,gamma.seq[1]);
    #index<-seq(1,ncol(x));index[fm$beta!=0];fm$beta[fm$beta!=0];
    sel.group <- identify.group(fm$beta,group.inner)$ident.g;
    if ( is.null(sel.group) == T) sel.group <- NA; 
    sse <- sum((y-x%*%fm$beta)^2)
    sse.seq <- c(sse.seq, sse)  
                                                                   
    outfile1 <- gsub("( )", "", paste("/space/jliu20/Compound_bridge/Sim/group_selected_",num,"_",k,"_",perc,"_",mag,"_",matched,"_",gamma.seq,"_",example,".txt"))
    write.table(matrix(sel.group,nrow=1),outfile1,sep=" ",quote=F, row.names=F,col.names=F,append=T);
    
    a <- match.group(sel.group,matched,n.group*3)
    tp = c(tp,a$tp);tn = c(tn,a$tn);fp = c(fp,a$fp);fn = c(fn,a$fn)
  }
  
  outfile2 <- gsub("( )", "", paste("/space/jliu20/Compound_bridge/Sim/sse_surv_",num,"_",k,"_",perc,"_",mag,"_",matched,"_",gamma.seq,"_",example,".txt"))
  write.table(sse.seq,outfile2,sep=",",quote=F, row.names=F,col.names=F,append=T);

  list(tp=tp,fp=fp,tn=tn,fn=fn,sse.seq=sse.seq,lambda.opt.seq=lambda.opt.seq,lambda.min.seq=lambda.min.seq,lambda.max.seq=lambda.max.seq)
}



#library(snow)
#cl <- makeCluster(5)


#source("/space/jliu20/Compound_bridge/sim.bic.r")
#source("/space/jliu20/Compound_bridge/generate.surv.categorical.r")
#sim <- function(num,n.group,g.size,matched=1,quan,mag,k,g,d,replic,epsilon,n.step,lambda.max,gamma.seq,example)
#num=100;n.group=200;g.size=5;matched=2;quan=0.9;mag=3;k=0.1;g=0;perc=0.9;epsilon=0.01;n.step=40;gamma.seq=0.9
#lambda.max=0.0020

#matched =2 25%
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0020,gamma.seq=0.9,"AR08")    #node01
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0015,gamma.seq=0.9,"AR05")    #node01
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0012,gamma.seq=0.9,"AR02")    #node01

#s=sim(100,200,5,matched=2,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"AR08")   #node02
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00020,gamma.seq=0.7,"AR05")   #node02
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00018,gamma.seq=0.7,"AR02")   #node02

#s=sim(100,200,5,matched=2,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00007,gamma.seq=0.5,"AR08")   #node03
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"AR05")   #node03   
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00003,gamma.seq=0.5,"AR02")   #node03

#s=sim(100,200,5,matched=2,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00120,gamma.seq=0.9,"Band")     #node04
#s=sim(100,200,5,matched=2,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00025,gamma.seq=0.7,"Band")     #node04   
#s=sim(100,200,5,matched=2,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.000035,gamma.seq=0.5,"Band")    #node04   

#s=sim(100,200,5,matched=2,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00160,gamma.seq=0.9,"Band")    #node06
#s=sim(100,200,5,matched=2,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"Band")    #node06
#s=sim(100,200,5,matched=2,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"Band")    #node06

#s=sim(100,200,5,matched=2,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00130,gamma.seq=0.9,"Band")    #node07
#s=sim(100,200,5,matched=2,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00017,gamma.seq=0.7,"Band")    #node07
#s=sim(100,200,5,matched=2,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00004,gamma.seq=0.5,"Band")    #node07

#matched =3 0%
#source("/space/jliu20/Compound_bridge/sim.bic.r")
#source("/space/jliu20/Compound_bridge/generate.surv.categorical.r")
#sum(s$lambda.opt.seq>s$lambda.min.seq)
#sum(s$lambda.opt.seq<s$lambda.max.seq)
#mean(s$fp);sd(s$fp);
#mean(s$tp);sd(s$tp);
#s=sim(100,200,5,matched=3,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0015,gamma.seq=0.9,"AR08")    #node08
#s=sim(100,200,5,matched=3,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0015,gamma.seq=0.9,"AR05")    #node08
#s=sim(100,200,5,matched=3,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0012,gamma.seq=0.9,"AR02")    #node08

#s=sim(100,200,5,matched=3,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"AR08")   #node09
#s=sim(100,200,5,matched=3,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"AR05")   #node09
#s=sim(100,200,5,matched=3,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00020,gamma.seq=0.7,"AR02")   #node09

#s=sim(100,200,5,matched=3,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"AR08")   #node10
#s=sim(100,200,5,matched=3,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"AR05")   #node10   
#s=sim(100,200,5,matched=3,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00004,gamma.seq=0.5,"AR02")   #node10

#s=sim(100,200,5,matched=3,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00150,gamma.seq=0.9,"Band")     #node11
#s=sim(100,200,5,matched=3,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00025,gamma.seq=0.7,"Band")     #node11   
#s=sim(100,200,5,matched=3,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00006,gamma.seq=0.5,"Band")     #node11   

#s=sim(100,200,5,matched=3,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00160,gamma.seq=0.9,"Band")    #node12
#s=sim(100,200,5,matched=3,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00025,gamma.seq=0.7,"Band")    #node12
#s=sim(100,200,5,matched=3,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"Band")    #node12

#s=sim(100,200,5,matched=3,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00140,gamma.seq=0.9,"Band")    #node13
#s=sim(100,200,5,matched=3,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"Band")    #node13
#s=sim(100,200,5,matched=3,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"Band")    #node13

#matched =4 50%; beta not all nonzero in an important cluster.
#source("/space/jliu20/Compound_bridge/sim.bic.r")
#source("/space/jliu20/Compound_bridge/generate.surv.categorical.r")
#s=sim(100,200,5,matched=4,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0015,gamma.seq=0.9,"AR08")    #node01
#s=sim(100,200,5,matched=4,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0015,gamma.seq=0.9,"AR05")    #node01
#s=sim(100,200,5,matched=4,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0012,gamma.seq=0.9,"AR02")    #node01

#s=sim(100,200,5,matched=4,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"AR08")   #node02
#s=sim(100,200,5,matched=4,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"AR05")   #node02
#s=sim(100,200,5,matched=4,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00020,gamma.seq=0.7,"AR02")   #node02

#s=sim(100,200,5,matched=4,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"AR08")   #node03
#s=sim(100,200,5,matched=4,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"AR05")   #node03   
#s=sim(100,200,5,matched=4,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00004,gamma.seq=0.5,"AR02")   #node03

#s=sim(100,200,5,matched=4,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00150,gamma.seq=0.9,"Band")     #node04
#s=sim(100,200,5,matched=4,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00025,gamma.seq=0.7,"Band")     #node04   
#s=sim(100,200,5,matched=4,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00006,gamma.seq=0.5,"Band")     #node04   

#s=sim(100,200,5,matched=4,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00160,gamma.seq=0.9,"Band")     #node06
#s=sim(100,200,5,matched=4,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00025,gamma.seq=0.7,"Band")     #node06
#s=sim(100,200,5,matched=4,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"Band")     #node06

#s=sim(100,200,5,matched=4,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00140,gamma.seq=0.9,"Band")     #node07
#s=sim(100,200,5,matched=4,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"Band")     #node07
#s=sim(100,200,5,matched=4,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"Band")     #node07

#matched =5 25%; beta not all nonzero in an important cluster.
#source("/space/jliu20/Compound_bridge/sim.bic.r")
#source("/space/jliu20/Compound_bridge/generate.surv.categorical.r")
#s=sim(100,200,5,matched=5,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0015,gamma.seq=0.9,"AR08")    #node08
#s=sim(100,200,5,matched=5,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0015,gamma.seq=0.9,"AR05")    #node08
#s=sim(100,200,5,matched=5,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0012,gamma.seq=0.9,"AR02")    #node08

#s=sim(100,200,5,matched=5,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"AR08")   #node09
#s=sim(100,200,5,matched=5,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"AR05")   #node09
#s=sim(100,200,5,matched=5,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00020,gamma.seq=0.7,"AR02")   #node09

#s=sim(100,200,5,matched=5,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"AR08")   #node10
#s=sim(100,200,5,matched=5,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"AR05")   #node10   
#s=sim(100,200,5,matched=5,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00004,gamma.seq=0.5,"AR02")   #node10

#s=sim(100,200,5,matched=5,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00150,gamma.seq=0.9,"Band")     #node11
#s=sim(100,200,5,matched=5,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00025,gamma.seq=0.7,"Band")     #node11   
#s=sim(100,200,5,matched=5,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00006,gamma.seq=0.5,"Band")     #node11   

#s=sim(100,200,5,matched=5,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00160,gamma.seq=0.9,"Band")     #node12
#s=sim(100,200,5,matched=5,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00025,gamma.seq=0.7,"Band")     #node12
#s=sim(100,200,5,matched=5,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"Band")     #node12

#s=sim(100,200,5,matched=5,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00140,gamma.seq=0.9,"Band")     #node13
#s=sim(100,200,5,matched=5,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"Band")     #node13
#s=sim(100,200,5,matched=5,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"Band")     #node13

#matched =6 0%; beta not all nonzero in an important cluster.
#source("/space/jliu20/Compound_bridge/sim.bic.r")
#source("/space/jliu20/Compound_bridge/generate.surv.categorical.r")
#s=sim(100,200,5,matched=6,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0015,gamma.seq=0.9,"AR08")    #node14
#s=sim(100,200,5,matched=6,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0015,gamma.seq=0.9,"AR05")    #node14
#s=sim(100,200,5,matched=6,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.0012,gamma.seq=0.9,"AR02")    #node14

#s=sim(100,200,5,matched=6,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"AR08")   #node15
#s=sim(100,200,5,matched=6,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"AR05")   #node15
#s=sim(100,200,5,matched=6,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00020,gamma.seq=0.7,"AR02")   #node15

#s=sim(100,200,5,matched=6,quan=0.9,mag=0.8,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"AR08")   #node16
#s=sim(100,200,5,matched=6,quan=0.9,mag=0.5,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"AR05")   #node16   
#s=sim(100,200,5,matched=6,quan=0.9,mag=0.2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00004,gamma.seq=0.5,"AR02")   #node16

#s=sim(100,200,5,matched=6,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00150,gamma.seq=0.9,"Band")     #node17
#s=sim(100,200,5,matched=6,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00025,gamma.seq=0.7,"Band")     #node17   
#s=sim(100,200,5,matched=6,quan=0.9,mag=1,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00006,gamma.seq=0.5,"Band")     #node17   

#s=sim(100,200,5,matched=6,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00160,gamma.seq=0.9,"Band")     #node18
#s=sim(100,200,5,matched=6,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00025,gamma.seq=0.7,"Band")     #node18
#s=sim(100,200,5,matched=6,quan=0.9,mag=2,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"Band")     #node18

#s=sim(100,200,5,matched=6,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00140,gamma.seq=0.9,"Band")     #node19
#s=sim(100,200,5,matched=6,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00030,gamma.seq=0.7,"Band")     #node19
#s=sim(100,200,5,matched=6,quan=0.9,mag=3,k=0.1,g=0,perc=0.9,5,epsilon=0.01,40,0.00005,gamma.seq=0.5,"Band")     #node19

#s$lambda.opt.seq> s$lambda.min.seq
###not use
#source("/space/jliu20/Compound_bridge/sim.bic.r")
#source("/space/jliu20/Compound_bridge/generate.surv.categorical.r")
#sim <- function(num,n.group,g.size,matched=1,quan,mag,k,g,d,replic,epsilon,n.step,lambda.max,gamma.seq,example)
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.8,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.0018,gamma.seq=0.9,"AR08")    #node01
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.5,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.0014,gamma.seq=0.9,"AR05")    #node01
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.2,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.0010,gamma.seq=0.9,"AR02")    #node01

#s=sim(100,200,5,matched=2,quan=0.9,mag=0.8,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00035,gamma.seq=0.7,"AR08")   #node02
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.5,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00019,gamma.seq=0.7,"AR05")   #node02
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.2,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00014,gamma.seq=0.7,"AR02")   #node02

#s=sim(100,200,5,matched=2,quan=0.9,mag=0.8,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00007,gamma.seq=0.5,"AR08")   #node03
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.5,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00004,gamma.seq=0.5,"AR05")   #node03
#s=sim(100,200,5,matched=2,quan=0.9,mag=0.2,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00003,gamma.seq=0.5,"AR02")   #node03

#s=sim(100,200,5,matched=2,quan=0.9,mag=1,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00120,gamma.seq=0.9,"Band")     #node04
#s=sim(100,200,5,matched=2,quan=0.9,mag=1,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00018,gamma.seq=0.7,"Band")     #node04
#s=sim(100,200,5,matched=2,quan=0.9,mag=1,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.000035,gamma.seq=0.5,"Band")    #node04

#s=sim(100,200,5,matched=2,quan=0.9,mag=2,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00160,gamma.seq=0.9,"Band")     #node08
#s=sim(100,200,5,matched=2,quan=0.9,mag=2,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00024,gamma.seq=0.7,"Band")     #node08
#s=sim(100,200,5,matched=2,quan=0.9,mag=2,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00004,gamma.seq=0.5,"Band")     #node08
                                      
#s=sim(100,200,5,matched=2,quan=0.9,mag=3,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00110,gamma.seq=0.9,"Band")     #node09
#s=sim(100,200,5,matched=2,quan=0.9,mag=3,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.00015,gamma.seq=0.7,"Band")     #node09
#s=sim(100,200,5,matched=2,quan=0.9,mag=3,k=0.1,g=0,d=3,5,epsilon=0.01,40,0.000027,gamma.seq=0.5,"Band")    #node09
