#50% not matched
match.s1 <- function(b,matched,n.group)
{
  if ( matched == 1|matched==4){
    #50% not matched
    a=c(1,4,7,10)
  }
  else if ( matched == 2|matched==5){
    #25% not matched
    a=c(1,4,7,10)
  }
  else if ( matched == 3|matched==6){
    # 0% not matched
    a=c(1,4,7,10)
  } 
  matched=match(a,b,nomatch=-1)
  tp <- length(matched[matched!=-1])
  fn <- length(a) -tp
  fp <- length(b) - tp
  tn <- n.group - tp -fn -fp
  list(tp=tp,fn=fn,fp=fp,tn=tn)
}

match.s2 <- function(b,matched,n.group)
{ 
  if ( matched == 1|matched==4){
    #50% not matched
    a=c(1,4,8,11)
  }
  else if ( matched == 2|matched==5){
    #25% not matched
    a=c(1,4,8,10)
  }
  else if ( matched == 3|matched==6){
    # 0% not matched
    a=c(1,4,7,10)
  } 
  matched=match(a,b,nomatch=-1)
  tp <- length(matched[matched!=-1])
  fn <- length(a) -tp
  fp <- length(b) - tp
  tn <- n.group - tp -fn -fp
  list(tp=tp,fn=fn,fp=fp,tn=tn)
}

match.s3 <- function(b,matched,n.group)
{  
  if ( matched == 1|matched==4){
    #50% not matched
    a=c(1,4,9,12)
  }
  else if ( matched == 2|matched==5){
    #25% not matched
    a=c(1,4,9,10)
  }
  else if ( matched == 3|matched==6){
    # 0% not matched
    a=c(1,4,7,10)
  } 
  matched=match(a,b,nomatch=-1)
  tp <- length(matched[matched!=-1])
  fn <- length(a) -tp
  fp <- length(b) - tp
  tn <- n.group - tp -fn -fp
  list(tp=tp,fn=fn,fp=fp,tn=tn)
}

sim <- function(num,matched=1,quan,g,mag,k,perc,replic,epsilon,n.step,gamma,example)
{
  tp <- numeric();tn <- numeric();fp <- numeric();fn <- numeric();sse.seq <- numeric();
  lambda.opt.seq1 <- numeric();lambda.min.seq1 <- numeric();lambda.max.seq1 <- numeric();
  lambda.opt.seq2 <- numeric();lambda.min.seq2 <- numeric();lambda.max.seq2 <- numeric();
  lambda.opt.seq3 <- numeric();lambda.min.seq3 <- numeric();lambda.max.seq3 <- numeric();
  group <- c(20,20,10,4,4,4,rep(5,3),6,6,6,rep(20,5),rep(10,10),rep(5,141)); n.group <- length(group);
  
  for ( i in 1:replic)
  {
    if ( matched == 1){
    #50% not matched  : value =1
    beta1 <- c(rep(0.1,20),rep(0,20),rep(0,10),rep(0.2,4),rep(0,4),rep(0,4),rep(0.1,5),rep(0,5),rep(0,5),rep(0.15,6),rep(0,6),rep(0,6),rep(0,905))
    beta2 <- c(rep(0.1,20),rep(0,20),rep(0,10),rep(0.2,4),rep(0,4),rep(0,4),rep(0,5),rep(0.1,5),rep(0,5),rep(0,6),rep(0.15,6),rep(0,6),rep(0,905))
    beta3 <- c(rep(0.1,20),rep(0,20),rep(0,10),rep(0.2,4),rep(0,4),rep(0,4),rep(0,5),rep(0,5),rep(0.1,5),rep(0,6),rep(0,6),rep(0.15,6),rep(0,905))
    }
    else if ( matched == 2){
    #25% not matched  : value =2
    beta1 <- c(rep(0.1,20),rep(0,20),rep(0,10),rep(0.2,4),rep(0,4),rep(0,4),rep(0.1,5),rep(0,5),rep(0,5),rep(0.15,6),rep(0,6),rep(0,6),rep(0,905))
    beta2 <- c(rep(0.1,20),rep(0,20),rep(0,10),rep(0.2,4),rep(0,4),rep(0,4),rep(0,5),rep(0.1,5),rep(0,5),rep(0.15,6),rep(0,6),rep(0,6),rep(0,905))
    beta3 <- c(rep(0.1,20),rep(0,20),rep(0,10),rep(0.2,4),rep(0,4),rep(0,4),rep(0,5),rep(0,5),rep(0.1,5),rep(0.15,6),rep(0,6),rep(0,6),rep(0,905))
    }
    else if ( matched ==3){
    #0 % not matched  : value =3
    beta1 <- c(rep(0.1,20),rep(0,20),rep(0,10),rep(0.2,4),rep(0,4),rep(0,4),rep(0.1,5),rep(0,5),rep(0,5),rep(0.15,6),rep(0,6),rep(0,6),rep(0,905))
    beta2 <- c(rep(0.1,20),rep(0,20),rep(0,10),rep(0.2,4),rep(0,4),rep(0,4),rep(0.1,5),rep(0,5),rep(0,5),rep(0.15,6),rep(0,6),rep(0,6),rep(0,905))
    beta3 <- c(rep(0.1,20),rep(0,20),rep(0,10),rep(0.2,4),rep(0,4),rep(0,4),rep(0.1,5),rep(0,5),rep(0,5),rep(0.15,6),rep(0,6),rep(0,6),rep(0,905))
    }
    else if ( matched == 4){
    #50% not matched  : value =4
    beta1 <- c(rep(0.1,15),rep(0,5),  rep(0,20),rep(0,10),rep(0.2,3),0, rep(0,4),rep(0,4),rep(0.1,3),0,0,rep(0,5),      rep(0,5),        rep(0.15,5),rep(0,7),rep(0,6),rep(0,905))
    beta2 <- c(rep(0,5),rep(0.1,15),  rep(0,20),rep(0,10),0,rep(0.2,3), rep(0,4),rep(0,4),rep(0,5),      0,0,rep(0.1,3),rep(0,5),        rep(0,6),rep(0.15,5),rep(0,7),rep(0,905))
    beta3 <- c(rep(0,5),rep(0.1,15),  rep(0,20),rep(0,10),0.2,0,0.2,0.2,rep(0,4),rep(0,4),rep(0,5),      rep(0,5),      0,rep(0.1,3),0,  rep(0,6),rep(0,7),rep(0.15,5),rep(0,905))
    }
    else if ( matched == 5){
    #25% not matched  : value =5
    beta1 <- c(rep(0.1,15),rep(0,5),  rep(0,20),rep(0,10),  rep(0.2,3),0,rep(0,4),rep(0,4), rep(0.1,3),0,0,rep(0,5),rep(0,5), rep(0.15,5),rep(0,7),rep(0,6),rep(0,905))
    beta2 <- c(rep(0,5),rep(0.1,15),  rep(0,20),rep(0,10),  0,rep(0.2,3),rep(0,4),rep(0,4), rep(0,5),0,0,rep(0.1,3),rep(0,5), rep(0.15,5),rep(0,7),rep(0,6),rep(0,905))
    beta3 <- c(rep(0,5),rep(0.1,15),  rep(0,20),rep(0,10),  0.2,0,0.2,0.2,rep(0,4),rep(0,4),rep(0,5),rep(0,5),0,rep(0.1,3),0, rep(0.15,5),rep(0,7),rep(0,6),rep(0,905))
    }
    else if ( matched ==6){
    #0 % not matched  : value =6; not all beta in a cluster are non-zero
    beta1 <- c(rep(0.1,15),rep(0,5),  rep(0,20),rep(0,10), rep(0.2,3),0,rep(0,4),rep(0,4),  rep(0.1,3),0,0,rep(0,5),rep(0,5),rep(0.15,5),rep(0,7),rep(0,6),rep(0,905))
    beta2 <- c(rep(0,5),rep(0.1,15),  rep(0,20),rep(0,10), 0,rep(0.2,3),rep(0,4),rep(0,4),  0,0,rep(0.1,3),rep(0,5),rep(0,5),rep(0.15,5),rep(0,7),rep(0,6),rep(0,905))
    beta3 <- c(rep(0,5),rep(0.1,15),  rep(0,20),rep(0,10), 0.2,0,0.2,0.2,rep(0,4),rep(0,4), 0,rep(0.1,3),0,rep(0,5),rep(0,5),rep(0.15,5),rep(0,7),rep(0,6),rep(0,905))    
    }

    source("/space/jliu20/Compound_bridge/data_proc.r");
    data.p <- data.proc.4(k, num, quan, group, mag,beta1,beta2,beta3)
    x1 <- data.p$x1; y1 <- data.p$y1; d1 <- data.p$d1; w1 <- data.p$w1;
    x2 <- data.p$x2; y2 <- data.p$y2; d2 <- data.p$d2; w2 <- data.p$w2;
    x3 <- data.p$x3; y3 <- data.p$y3; d3 <- data.p$d3; w3 <- data.p$w3;
    #rho <- cor_x(rbind(x1,x2,x3,x4),group)$rho
    
    source("/space/jliu20/Compound_bridge/PCA.r"); 
    t1<- PCA.x(x1,perc,group); t2<- PCA.x(x2,perc,group); t3<- PCA.x(x3,perc,group);  
    #t1<- PCA.x.f(x1,d,group); t2<- PCA.x.f(x2,d,group); t3<- PCA.x.f(x3,d,group);  
    x1<- t1$x.new ; x2 <- t2$x.new; x3 <- t3$x.new;
    group1 <- t1$group.new; group2 <- t2$group.new; group3 <- t3$group.new; 

    data.r <- reweight3(x1,x2,x3,y1,y2,y3,w1,w2,w3)
    x1.sd <- data.r$x1.sd; x2.sd <- data.r$x2.sd; x3.sd <- data.r$x3.sd; n.t <- data.r$n.t; y <- data.r$y;
    pe <- cumsum(n.t); ps <- pe -n.t + 1;
    y1.sd <- y[ps[1]:pe[1]]; y2.sd <- y[ps[2]:pe[2]]; y3.sd <- y[ps[3]:pe[3]]
    
    #data.f <- process3(x1.sd,x2.sd,x3.sd,n.t,group1,group2,group3)
    #x <- data.f$x; x0 <- data.f$x0; ps <- data.f$ps; pe <- data.f$pe; group.inner <- data.f$group.inner
    #totalp <- ncol(x); m <- 3;
    #group.outer <- rep(m,length(group.inner)/m);
    
    
    #source("/space/jliu20/Compound_bridge/cv_mult_gMCP.r")
    #n.fold <- 5; 
    #cv <- cv.optim(cl,x,y,n.t,group.inner,n.fold,epsilon,n.step,gamma)
    #PE <- cv$PE; lambda.seq <- cv$lambda.seq; lambda.opt <- lambda.seq[PE==min(PE)]
    source("/space/jliu20/Compound_bridge/BIC.gMCP.r")
    BIC.1t <-  BIC.gMCP.ind(x1.sd,y1.sd,group1,epsilon,n.step,gamma,g)
    bic1 <- BIC.1t$bic; lambda.seq1 <- BIC.1t$lambda.seq; lambda.opt1 <- min(lambda.seq1[bic1==min(bic1)]);b1.opt <- BIC.1t$beta.opt

    BIC.2t <-  BIC.gMCP.ind(x2.sd,y2.sd,group2,epsilon,n.step,gamma,g)
    bic2 <- BIC.2t$bic; lambda.seq2 <- BIC.2t$lambda.seq; lambda.opt2 <- min(lambda.seq2[bic2==min(bic2)]);b2.opt <- BIC.2t$beta.opt
    
    BIC.3t <-  BIC.gMCP.ind(x3.sd,y3.sd,group3,epsilon,n.step,gamma,g)
    bic3 <- BIC.3t$bic; lambda.seq3 <- BIC.3t$lambda.seq; lambda.opt3 <- min(lambda.seq3[bic3==min(bic3)]);b3.opt <- BIC.3t$beta.opt    

    lambda.opt.seq1 <- c(lambda.opt.seq1,lambda.opt1);lambda.min.seq1 <- c(lambda.min.seq1,min(lambda.seq1));lambda.max.seq1 <- c(lambda.max.seq1,max(lambda.seq1))
    lambda.opt.seq2 <- c(lambda.opt.seq2,lambda.opt2);lambda.min.seq2 <- c(lambda.min.seq2,min(lambda.seq2));lambda.max.seq2 <- c(lambda.max.seq2,max(lambda.seq2))
    lambda.opt.seq3 <- c(lambda.opt.seq3,lambda.opt3);lambda.min.seq3 <- c(lambda.min.seq3,min(lambda.seq3));lambda.max.seq3 <- c(lambda.max.seq3,max(lambda.seq3))
    #fm <- GMCP_quant(x,y,group.inner,lambda.opt,gamma)
    sel1.group <- identify.group(b1.opt,group1)$ident.g;
    sel2.group <- identify.group(b2.opt,group2)$ident.g;
    sel3.group <- identify.group(b3.opt,group3)$ident.g;
    
    if ( is.null(sel1.group) == T) sel1.group <- NA; 
    if ( is.null(sel2.group) == T) sel2.group <- NA; 
    if ( is.null(sel3.group) == T) sel3.group <- NA; 
    
    sse <- sum((y1.sd-x1.sd%*%b1.opt)^2)+sum((y2.sd-x2.sd%*%b2.opt)^2) + sum((y3.sd-x3.sd%*%b3.opt)^2)
    sse.seq <- c(sse.seq, sse)  
                                                                   
    outfile1 <- gsub("( )", "", paste("/space/jliu20/Compound_bridge/Sim/group_selected_gMCP_sc2_",num,"_",k,"_",perc,"_",mag,"_",matched,"_",gamma,"_",example,".txt"))
    write.table(matrix(sel1.group,nrow=1),outfile1,sep=" ",quote=F, row.names=F,col.names=F,append=T);
    write.table(matrix(sel2.group,nrow=1),outfile1,sep=" ",quote=F, row.names=F,col.names=F,append=T);
    write.table(matrix(sel3.group,nrow=1),outfile1,sep=" ",quote=F, row.names=F,col.names=F,append=T);

    a1 <- match.s1(sel1.group,matched,n.group);  a2 <- match.s2(sel2.group,matched,n.group);    a3 <- match.s3(sel3.group,matched,n.group)
    tp = c(tp,a1$tp+a2$tp+a3$tp);tn = c(tn,a1$tn+a2$tn+a3$tn);fp = c(fp,a1$fp+a2$fp+a3$fp);fn = c(fn,a1$fn+a2$fn+a3$fn)
  }
  
  outfile2 <- gsub("( )", "", paste("/space/jliu20/Compound_bridge/Sim/sse_surv_gMCP_sc2_",num,"_",k,"_",perc,"_",mag,"_",matched,"_",gamma,"_",example,".txt"))
  write.table(sse.seq,outfile2,sep=",",quote=F, row.names=F,col.names=F,append=T);

  list(tp=tp,fp=fp,tn=tn,fn=fn,sse.seq=sse.seq,lambda.opt.seq1=lambda.opt.seq1,lambda.min.seq1=lambda.min.seq1,lambda.max.seq1=lambda.max.seq1,
  lambda.opt.seq2=lambda.opt.seq2,lambda.min.seq2=lambda.min.seq2,lambda.max.seq2=lambda.max.seq2,
  lambda.opt.seq3=lambda.opt.seq3,lambda.min.seq3=lambda.min.seq3,lambda.max.seq3=lambda.max.seq3)
}


#library(snow)
#cl <- makeCluster(5)


#source("sim_gMCP.bic_uneven_sc2.r")
#source("generate.surv.categorical_uneven.r")
# matched == 1 50%
#s <- sim(100,matched=1,quan=0.9,g=0,mag=0.2,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR02") #node01
#s <- sim(100,matched=1,quan=0.9,g=0,mag=0.5,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR05") #node01
#s <- sim(100,matched=1,quan=0.9,g=0,mag=0.8,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR08") #node01

#s <- sim(100,matched=1,quan=0.9,g=0,mag=1,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node02
#s <- sim(100,matched=1,quan=0.9,g=0,mag=2,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node02
#s <- sim(100,matched=1,quan=0.9,g=0,mag=3,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node02

# matched == 2 25%
#s <- sim(100,matched=2,quan=0.9,g=0,mag=0.2,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR02") #node03
#s <- sim(100,matched=2,quan=0.9,g=0,mag=0.5,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR05") #node03
#s <- sim(100,matched=2,quan=0.9,g=0,mag=0.8,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR08") #node03

#s <- sim(100,matched=2,quan=0.9,g=0,mag=1,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node04
#s <- sim(100,matched=2,quan=0.9,g=0,mag=2,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node04
#s <- sim(100,matched=2,quan=0.9,g=0,mag=3,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node04

# matched == 3 0%
#s <- sim(100,matched=3,quan=0.9,g=0,mag=0.2,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR02") #node06
#s <- sim(100,matched=3,quan=0.9,g=0,mag=0.5,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR05") #node06
#s <- sim(100,matched=3,quan=0.9,g=0,mag=0.8,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR08") #node06

#s <- sim(100,matched=3,quan=0.9,g=0,mag=1,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node07
#s <- sim(100,matched=3,quan=0.9,g=0,mag=2,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node07
#s <- sim(100,matched=3,quan=0.9,g=0,mag=3,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node07

source("sim_gMCP.bic_uneven_sc2.r")
source("generate.surv.categorical_uneven.r")
#matched=4; group Lasso
#s <- sim(100,matched=4,quan=0.9,g=0,mag=0.2,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR02") #node20
#s <- sim(100,matched=4,quan=0.9,g=0,mag=0.5,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR05") #node20
#s <- sim(100,matched=4,quan=0.9,g=0,mag=0.8,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR08") #node20

#s <- sim(100,matched=4,quan=0.9,g=0,mag=1,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node20
#s <- sim(100,matched=4,quan=0.9,g=0,mag=2,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node21
#s <- sim(100,matched=4,quan=0.9,g=0,mag=3,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node21

#matched=5; group Lasso
#s <- sim(100,matched=5,quan=0.9,g=0,mag=0.2,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR02") #node21
#s <- sim(100,matched=5,quan=0.9,g=0,mag=0.5,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR05") #node21
#s <- sim(100,matched=5,quan=0.9,g=0,mag=0.8,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR08") #node01

#s <- sim(100,matched=5,quan=0.9,g=0,mag=1,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node02
#s <- sim(100,matched=5,quan=0.9,g=0,mag=2,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node03
#s <- sim(100,matched=5,quan=0.9,g=0,mag=3,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node04

#matched=6; group Lasso
#s <- sim(100,matched=6,quan=0.9,g=0,mag=0.2,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR02") #node06
#s <- sim(100,matched=6,quan=0.9,g=0,mag=0.5,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR05") #node07
#s <- sim(100,matched=6,quan=0.9,g=0,mag=0.8,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"AR08") #node08

#s <- sim(100,matched=6,quan=0.9,g=0,mag=1,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node09
#s <- sim(100,matched=6,quan=0.9,g=0,mag=2,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node10
#s <- sim(100,matched=6,quan=0.9,g=0,mag=3,0.1,0.9,5,epsilon=0.01,n.step=40,gamma=10000,"Band") #node11

#sum(s$lambda.opt.seq1 > s$lambda.min.seq1)
#sum(s$lambda.opt.seq2 > s$lambda.min.seq2)
#sum(s$lambda.opt.seq3 > s$lambda.min.seq3)
#mean(s$tp);sd(s$tp);mean(s$fp);sd(s$fp);

#source("sim_gMCP.bic_uneven_sc2.r")
#source("generate.surv.categorical.r")

#sim.outer(num,n.group,g.size,quan,mag,k,d,replic,epsilon,n.step,gamma,example)
#s <- sim.outer(100,100,10,quan=0.9,mag=0.2,k=0.1,3,5,epsilon=0.05,n.step=40,gamma=10000,"AR02") #node01
#s <- sim.outer(100,100,10,quan=0.9,mag=0.5,k=0.1,3,5,epsilon=0.02,n.step=40,gamma=10000,"AR05") #node02
#s <- sim.outer(100,100,10,quan=0.9,mag=0.8,k=0.1,3,5,epsilon=0.02,n.step=40,gamma=10000,"AR08") #node03

#s <- sim.outer(100,100,10,quan=0.9,mag=1,k=0.1,3,5,epsilon=0.05,n.step=40,gamma=10000,"Band") #node04
#s <- sim.outer(100,100,10,quan=0.9,mag=2,k=0.1,3,5,epsilon=0.04,n.step=40,gamma=10000,"Band") #node06
#s <- sim.outer(100,100,10,quan=0.9,mag=3,k=0.1,3,5,epsilon=0.02,n.step=40,gamma=10000,"Band") #node07

#s <- sim.outer(100,100,10,quan=0.9,mag=0.2,k=0.1,3,5,epsilon=0.05,n.step=40,gamma=6,"AR02") #node08
#s <- sim.outer(100,100,10,quan=0.9,mag=0.5,k=0.1,3,5,epsilon=0.02,n.step=40,gamma=6,"AR05") #node09
#s <- sim.outer(100,100,10,quan=0.9,mag=0.8,k=0.1,3,5,epsilon=0.02,n.step=40,gamma=6,"AR08") #node10

#s <- sim.outer(100,100,10,quan=0.9,mag=1,k=0.1,3,5,epsilon=0.02,n.step=40,gamma=6,"Band") #node11
#s <- sim.outer(100,100,10,quan=0.9,mag=2,k=0.1,3,5,epsilon=0.02,n.step=40,gamma=6,"Band") #node12
#s <- sim.outer(100,100,10,quan=0.9,mag=3,k=0.1,3,5,epsilon=0.02,n.step=40,gamma=6,"Band") #node13
