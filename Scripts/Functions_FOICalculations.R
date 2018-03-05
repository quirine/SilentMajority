

FOIsam.full = function(b, m, g, a, n, Ca, Cm, Cs, Xa, Xm, Xs, recov=1/MaxTime){
  (b * m * a ) *  ( (sum((recov*Ca) * Xa) + sum((recov*Cm) * Xm) + sum((recov*Cs) * Xs)) / (ga + (sum((recov*Ca) * Xa) + sum((recov*Cm) * Xm) + sum((recov*Cs) * Xs)))) * exp(-g*n) 
}

FOIsm.full = function(b, m, g, a, n, Cm, Cs, Xm, Xs, recov=1/MaxTime){
  (b * m * a ) *  ( (sum((recov*Cm) * Xm) + sum((recov*Cs) * Xs)) / (ga + (sum((recov*Cm) * Xm) + sum((recov*Cs) * Xs)))) * exp(-g*n) 
}

FOIs.full = function(b, m, g, a, n, Cm, Cs, Xm, Xs, recov=1/MaxTime){
  (b * m * a ) *  ( ( sum((recov*Cs) * Xs)) / (ga + (sum((recov*Cm) * Xm) + sum((recov*Cs) * Xs)))) * exp(-g*n) 
}

FOIsm.FOI.sam.full = function(Ca, Cm, Cs, Xa, Xm, Xs) {
  b =1 ; m = 1;  g = 1; a = 4; n = 1/7;  # these values cancel out in the devision
  FOIsm.full(b, m, g, a, n, Cm, Cs, Xm, Xs) / FOIsam.full(b, m, g, a, n, Ca, Cm, Cs, Xa, Xm, Xs)  
}

FOIs.FOI.sam.full = function(Ca, Cm, Cs, Xa, Xm, Xs) {
  b =1 ; m = 1;  g = 1; a = 4; n = 1/7;  # these values cancel out in the devision
  FOIs.full(b, m, g, a, n, Cm, Cs, Xm, Xs) / FOIsam.full(b, m, g, a, n, Ca, Cm, Cs, Xa, Xm, Xs)  
}

FOIsm.FOI.sam.sweep  = function(Ca, Cm, Cs, Xa, Xm, Xs, bites) {
  b =1 ; m = 1;  g = 1; a = 4; n = 1/7; # these values cancel out in the devision
  FOIsm.full(b, m, g, a, n, Cm, Cs, Xm, Xs) / FOIsam.full(b, m, g, a, n, Ca, Cm, Cs, Xa, Xm, Xs)  
}

Disease.prop.calculator = function(infec.hist, Prim.dis, Sec.dis, Tert.dis){
  Inf.prim<-numeric(); 
  for (ii in 1:dim(infec.hist)[1]) Inf.prim <- rbind(Inf.prim, infec.hist[ii,1] * Prim.dis)
  Inf.sec<-numeric(); 
  for (ii in 1:dim(infec.hist)[1]) Inf.sec <- rbind(Inf.sec, infec.hist[ii,2] * Sec.dis)
  Inf.tert<-numeric(); 
  for (ii in 1:dim(infec.hist)[1]) Inf.tert <- rbind(Inf.tert, infec.hist[ii,3] * Tert.dis)
  Inf.quart<-numeric(); 
  for (ii in 1:dim(infec.hist)[1]) Inf.quart <- rbind(Inf.quart, infec.hist[ii,4] * Tert.dis)
  
  Xa <- Xm <- Xs <- numeric()
  Xa <- cbind(Inf.prim[,1], Inf.sec[,1], Inf.tert[,1], Inf.quart[,1]) 
  Xm <- cbind(Inf.prim[,2], Inf.sec[,2], Inf.tert[,2], Inf.quart[,2]) 
  Xs <- cbind(Inf.prim[,3], Inf.sec[,3], Inf.tert[,3], Inf.quart[,3]) 
  return(data.frame(Xa,Xm,Xs))
}

Sympto.FOI.calculator = function(prev,X,Ca,Cm,Cs,infs,foi.vec,severe = FALSE){
  Prop.FOI = numeric()
  # browser()
  for (ii in 1:length(foi.vec)){
    Xa <- X[,1:4]; Xm <- X[,5:8]; Xs <- X[,9:12]; 
    x.tot <- prev[ii]  
    if (severe==FALSE){
      Prop.FOI <- c(Prop.FOI, FOIsm.FOI.sam.full(Ca, Cm, Cs, Xa[ii,1:infs]*x.tot, Xm[ii,1:infs]*x.tot, Xs[ii,1:infs]*x.tot))
    }
    else if (severe == TRUE){
      Prop.FOI <- c(Prop.FOI, FOIs.FOI.sam.full(Ca, Cm, Cs, Xa[ii,1:infs]*x.tot, Xm[ii,1:infs]*x.tot, Xs[ii,1:infs]*x.tot))
    }
  }
  return(Prop.FOI)
}

Sympto.FOI.sampler = function(prev,X,infs,foi.vec,Ca.all,Cm.all,Cs.all,cores=2,severe = FALSE){
  print(severe)
  Prop.FOI.all <- numeric()
  cl = registerDoParallel(cores=cores)
  Prop.FOI.all = foreach(ii = 1:dim(Ca.all)[1], .combine=cbind) %dopar% {
    Ca <- Ca.all[ii,]; Cm <- Cm.all[ii,]; Cs <- Cs.all[ii,]; 
    Sympto.FOI.calculator(prev,X,Ca,Cm,Cs,infs,foi.vec,severe)
  } 
  
return(Prop.FOI.all)
}



