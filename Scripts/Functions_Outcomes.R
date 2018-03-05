# functions to compile outcomes for final figures and paper

Derive.outcomes.X <- function(X, Prop.infectivity.prior.to.IIP.1, Prop.infectivity.prior.to.IIP.2, 
                            Prop.infectivity.prior.to.IIP.3, foi.vec, inf = 2)
{
  Xa <- X[,1:4]; Xm <- X[,5:8]; Xs <- X[,9:12]; 
  if (inf == 2){
    Xa.prim <- Xa[,1]; Xm.prim <- Xm[,1];  Xs.prim <- Xs[,1];
    Xa.sec <- Xa[,2]; Xm.sec <- Xm[,2];  Xs.sec <- Xs[,2];
    Xa <- rowSums(Xa[,1:2]); Xm <- rowSums(Xm[,1:2]);  Xs <- rowSums(Xs[,1:2]);
    Xa.norm <- Xa / (Xa + Xm + Xs);  Xm.norm <- Xm / (Xa + Xm + Xs);  Xs.norm <- Xs / (Xa + Xm + Xs);
    Xa.norm.prim <- Xa.prim / (Xa + Xm + Xs);  Xm.norm.prim <- Xm.prim / (Xa + Xm + Xs);  Xs.norm.prim <- Xs.prim/ (Xa + Xm + Xs);
    Xa.norm.sec <- Xa.sec / (Xa + Xm + Xs);  Xm.norm.sec <- Xm.sec / (Xa + Xm + Xs);  Xs.norm.sec <- Xs.sec / (Xa + Xm + Xs);
    Xinap.prim <- inapp.ratio.rev[1] - Xa.norm.prim
    Xinap.sec <- inapp.ratio.rev[2] - Xa.norm.sec
    Xinap <- (Xinap.prim * ( (Xa.prim + Xm.prim + Xs.prim )/(Xa + Xm + Xs) ) + Xinap.sec * ((Xa.sec + Xm.sec + Xs.sec ) / (Xa + Xm + Xs)) ) 
  }
  if (inf == 4){  
    Xa.prim <- Xa[,1]; Xm.prim <- Xm[,1];  Xs.prim <- Xs[,1];
    Xa.sec <- Xa[,2]; Xm.sec <- Xm[,2];  Xs.sec <- Xs[,2];
    Xa <- rowSums(X[,1:4]); Xm <- rowSums(X[,5:8]);  Xs <- rowSums(X[,9:12]);
    Xa.post <- rowSums(X[,3:4]); Xm.post <- rowSums(X[,7:8]);  Xs.post <- rowSums(X[,11:12]);
    Xa.norm <- Xa / (Xa + Xm + Xs );  Xm.norm <- Xm / (Xa + Xm + Xs );  Xs.norm <- Xs / (Xa + Xm + Xs );
    Xa.norm.prim <- Xa.prim / (Xa + Xm + Xs);  Xm.norm.prim <- Xm.prim / (Xa + Xm + Xs); Xs.norm.prim <- Xs.prim/ (Xa + Xm + Xs);
    Xa.norm.sec <- Xa.sec / (Xa + Xm + Xs);  Xm.norm.sec <- Xm.sec / (Xa + Xm + Xs);  Xs.norm.sec <- Xs.sec / (Xa + Xm + Xs);
    Xa.norm.post <- Xa.post / (Xa + Xm + Xs);  Xm.norm.post <- Xm.post / (Xa + Xm + Xs); Xs.norm.post <- Xs.post / (Xa + Xm + Xs);
    Xinap.prim <- inapp.ratio.rev[1] - Xa.norm.prim
    Xinap.sec <- inapp.ratio.rev[2] - Xa.norm.sec
    Xinap.post <- inapp.ratio.rev[3] - Xa.norm.post
    Xinap <- (Xinap.prim * ( (Xa.prim + Xm.prim + Xs.prim )/(Xa + Xm + Xs) ) + Xinap.sec * ( (Xa.sec + Xm.sec + Xs.sec )/(Xa + Xm + Xs) ) + Xinap.post * ( (Xa.post + Xm.post + Xs.post )/(Xa + Xm + Xs) )) 
  }
  Xsym <- Xm.norm + Xs.norm
  Xap <- 1 - Xa.norm - Xinap
  Prop.sym.inap <- Xinap / Xsym
  Prop.sym.ap <- Xap / Xsym
  Prop.sym.prim <- ( X[,5] + X[,9] ) / (Xm + Xs) 
  Prop.severe.prim <- X[,9]/Xs
  if (inf == 2){
    Prop.sym.sec <- ( X[,6] + X[,10] ) / (Xm + Xs) 
    Prop.severe.sec <- X[,10]/Xs
  }
  if (inf == 4){
    Prop.sym.sec <- ( rowSums(X[,6:8]) + rowSums(X[,10:12]) ) / (Xm + Xs) # and post-sec
    Prop.severe.sec <- ( rowSums(X[,10:12]) ) / (Xs) # and post-sec
  }
  
  Prop.sym.prim.presym <- Prop.sym.prim * median(Prop.infectivity.prior.to.IIP.1)
  Prop.sym.sec.presym <- Prop.sym.sec * median(Prop.infectivity.prior.to.IIP.2)  
  Prop.sym.presym <- Prop.sym.prim.presym + Prop.sym.sec.presym 
  
  Prop.severe.prim.presym <- Prop.severe.prim * median(Prop.infectivity.prior.to.IIP.3)  
  Prop.severe.sec.presym <- Prop.severe.sec * median(Prop.infectivity.prior.to.IIP.3)  
  Prop.severe.presym <- Prop.severe.prim.presym + Prop.severe.sec.presym 
  
  Prop.sym.prim.presym.low <- Prop.sym.prim * HDIofMCMC(Prop.infectivity.prior.to.IIP.1,0.95)[1] 
  Prop.sym.prim.presym.up <- Prop.sym.prim * HDIofMCMC(Prop.infectivity.prior.to.IIP.1,0.95)[2]
  Prop.sym.sec.presym.low <- Prop.sym.sec * HDIofMCMC(Prop.infectivity.prior.to.IIP.2,0.95)[1]
  Prop.sym.sec.presym.up <- Prop.sym.sec * HDIofMCMC(Prop.infectivity.prior.to.IIP.2,0.95)[2]
  Prop.sym.presym.low <- Prop.sym.prim.presym.low + Prop.sym.sec.presym.low
  Prop.sym.presym.up <- Prop.sym.prim.presym.up + Prop.sym.sec.presym.up
  
  Prop.sym.postsym <- 1 - Prop.sym.presym  
  Prop.severe.postsym <- 1 - Prop.severe.presym
  Prop.sym.ap.presym <- Prop.sym.presym * Prop.sym.ap
  Prop.sym.ap.postsym <- Prop.sym.postsym * Prop.sym.ap
  Prop.sym.inap.presym <- Prop.sym.presym * Prop.sym.inap
  Prop.sym.inap.postsym <- Prop.sym.postsym * Prop.sym.inap
  Outcomes <- data.frame(Xinap, Xsym, Xap, Xm, Xs, Prop.sym.inap, Prop.sym.ap, Prop.sym.prim, Prop.sym.sec,
                         Prop.sym.prim.presym, Prop.sym.sec.presym, Prop.sym.presym, Prop.sym.presym.low, Prop.sym.presym.up,
                         Prop.sym.postsym, Prop.sym.ap.presym, Prop.sym.ap.postsym, Prop.sym.inap.presym, Prop.sym.inap.postsym, Prop.severe.postsym )
  return(Outcomes)
}

Derive.outcomes.FOI <- function(Prop.FOI, XProps, inap = FALSE, Prop.FOI.DHF = FALSE) {
  Prop.FOI.sym.all <- Prop.FOI[,sort.list(Prop.FOI[NROW(Prop.FOI),])]
  Prop.FOI.asym.all <- 1 - Prop.FOI.sym.all
  if (inap == TRUE){
    Prop.FOI.asym.all <- (1 - Prop.FOI.sym.all) * Prim.dis[1]/inapp.ratio.rev[1]   
  }
  Prop.FOI.inap.all <- Prop.FOI.sym.all * XProps$Prop.sym.inap
  Prop.FOI.inap.asym.all <- Prop.FOI.asym.all + Prop.FOI.inap.all 
  Prop.FOI.presym.all <- Prop.FOI.sym.all * XProps$Prop.sym.ap.presym
  Prop.FOI.postsym.all <- Prop.FOI.sym.all * XProps$Prop.sym.ap.postsym
  Prop.FOI.inap.presym.all <- Prop.FOI.sym.all * XProps$Prop.sym.inap.presym
  Prop.FOI.inap.postsym.all <- Prop.FOI.sym.all * XProps$Prop.sym.inap.postsym
  idx <- which.min(abs(Prop.FOI[100,] - median(Prop.FOI[100,])))
  Prop.FOI.sym <- Prop.FOI[,idx]
  Prop.FOI.asym <- 1 - Prop.FOI.sym
  Prop.FOI.inap <- Prop.FOI.sym * XProps$Prop.sym.inap
  Prop.FOI.inap.asym <- Prop.FOI.asym + Prop.FOI.inap
  Prop.FOI.presym <- Prop.FOI.sym * XProps$Prop.sym.ap.presym
  Prop.FOI.postsym <- Prop.FOI.sym * XProps$Prop.sym.ap.postsym
  Prop.FOI.inap.presym <- Prop.FOI.sym * XProps$Prop.sym.inap.presym
  Prop.FOI.inap.postsym <- Prop.FOI.sym * XProps$Prop.sym.inap.postsym
  if (Prop.FOI.DHF == FALSE){
    Prop.FOI.postsym.detected <- Prop.FOI.postsym * detection.rate
    Prop.FOI.postsym.detected.low <- Prop.FOI.postsym * detection.rate.CI[1]
    Prop.FOI.postsym.detected.high <- Prop.FOI.postsym * detection.rate.CI[2]
  }
  else {
    idx <- which.min(abs(Prop.FOI.DHF[100,] - median(Prop.FOI.DHF[100,])))
    Prop.FOI.severe <- Prop.FOI.DHF[,idx]
    Prop.FOI.postsym <- (Prop.FOI.sym - Prop.FOI.severe) * XProps$Prop.sym.ap.postsym
    Prop.FOI.severe.postsym <- Prop.FOI.severe * XProps$Prop.severe.postsym  
    prop.severe <- XProps$Xs / (XProps$Xs + XProps$Xm)
    mild.detection.rate <- detection.rate - prop.severe
    mild.detection.rate.low <- detection.rate.CI[1] - prop.severe
    mild.detection.rate.high <- detection.rate.CI[2] - prop.severe
    Prop.FOI.postsym.detected <- Prop.FOI.postsym * mild.detection.rate + Prop.FOI.severe.postsym
    Prop.FOI.postsym.detected.low <- Prop.FOI.postsym * mild.detection.rate.low + Prop.FOI.severe.postsym
    Prop.FOI.postsym.detected.high <- Prop.FOI.postsym * mild.detection.rate.high + Prop.FOI.severe.postsym
  }
  
  Outcomes <- data.frame(Prop.FOI.sym, Prop.FOI.asym, Prop.FOI.inap, Prop.FOI.inap.asym,
                        Prop.FOI.presym, Prop.FOI.postsym, Prop.FOI.inap.presym, Prop.FOI.inap.postsym,Prop.FOI.postsym.detected,Prop.FOI.postsym.detected.low,Prop.FOI.postsym.detected.high)
}

Derive.outcomes.FOI.CIs <- function(Prop.FOI, XProps,Prop.FOI.DHF = FALSE) {
  # browser()
  Prop.FOI.sym.all <- Prop.FOI[,sort.list(Prop.FOI[NROW(Prop.FOI),])]
  Prop.FOI.asym.all <- 1 - Prop.FOI.sym.all
  Prop.FOI.inap.all <- Prop.FOI.sym.all * XProps$Prop.sym.inap
  Prop.FOI.inap.asym.all <- Prop.FOI.asym.all + Prop.FOI.inap.all 
  Prop.FOI.presym.all <- Prop.FOI.sym.all * XProps$Prop.sym.ap.presym
  Prop.FOI.postsym.all <- Prop.FOI.sym.all * XProps$Prop.sym.ap.postsym
  Prop.FOI.inap.presym.all <- Prop.FOI.sym.all * XProps$Prop.sym.inap.presym
  Prop.FOI.inap.postsym.all <- Prop.FOI.sym.all * XProps$Prop.sym.inap.postsym
  Prop.FOI.presym.of.ap.all <- Prop.FOI.presym.all / (Prop.FOI.presym.all + Prop.FOI.postsym.all)
  if (Prop.FOI.DHF == FALSE) {
  Prop.FOI.postsym.detected.all <- Prop.FOI.postsym.all * detection.rate
  Prop.FOI.postsym.detected.low.all <- Prop.FOI.postsym.all * detection.rate.CI[1]
  Prop.FOI.postsym.detected.high.all <- Prop.FOI.postsym.all * detection.rate.CI[2]
  }
  else {
    Prop.FOI.severe.all <- Prop.FOI.DHF[,sort.list(Prop.FOI[NROW(Prop.FOI.DHF),])]
    Prop.FOI.postsym.all <- (Prop.FOI.sym.all - Prop.FOI.severe.all) * XProps$Prop.sym.ap.postsym
    Prop.FOI.severe.postsym.all <- Prop.FOI.severe.all * XProps$Prop.severe.postsym  
    prop.severe <- XProps$Xs / (XProps$Xs + XProps$Xm)
    mild.detection.rate <- detection.rate - prop.severe
    mild.detection.rate.low <- detection.rate.CI[1] - prop.severe
    mild.detection.rate.high <- detection.rate.CI[2] - prop.severe
    Prop.FOI.postsym.detected.all <- Prop.FOI.postsym.all * mild.detection.rate + Prop.FOI.severe.postsym.all
    Prop.FOI.postsym.detected.low.all <- Prop.FOI.postsym.all * mild.detection.rate.low + Prop.FOI.severe.postsym.all
    Prop.FOI.postsym.detected.high.all <- Prop.FOI.postsym.all * mild.detection.rate.high + Prop.FOI.severe.postsym.all
  }
  # CIs for foi = 0.1
  Prop.FOI.sym.CI <- HDIofMCMC(Prop.FOI.sym.all[100,], 0.95)
  Prop.FOI.asym.CI <- HDIofMCMC(Prop.FOI.asym.all[100,], 0.95)
  Prop.FOI.inap.CI <- HDIofMCMC(Prop.FOI.inap.all[100,], 0.95)
  Prop.FOI.inap.asym.CI <- HDIofMCMC(Prop.FOI.inap.asym.all[100,], 0.95)
  Prop.FOI.presym.CI <- HDIofMCMC(Prop.FOI.presym.all[100,], 0.95)
  Prop.FOI.postsym.CI <- HDIofMCMC(Prop.FOI.postsym.all[100,], 0.95)
  Prop.FOI.inap.presym.CI <- HDIofMCMC(Prop.FOI.inap.presym.all[100,], 0.95)
  Prop.FOI.inap.postsym.CI <- HDIofMCMC(Prop.FOI.inap.postsym.all[100,], 0.95)
  Prop.FOI.presym.of.ap.all.CI <- HDIofMCMC(Prop.FOI.presym.of.ap.all[100,], 0.95)
  Prop.FOI.postsym.detected.all.CI <-HDIofMCMC(Prop.FOI.postsym.detected.all[100,], 0.95) 
  Prop.FOI.postsym.detected.low.all.CI <-HDIofMCMC(Prop.FOI.postsym.detected.low.all[100,], 0.95) 
  Prop.FOI.postsym.detected.high.all.CI <-HDIofMCMC(Prop.FOI.postsym.detected.high.all[100,], 0.95) 
  
  Outcomes <- data.frame(Prop.FOI.sym.CI, Prop.FOI.asym.CI, Prop.FOI.inap.CI, Prop.FOI.inap.asym.CI,
                         Prop.FOI.presym.CI, Prop.FOI.postsym.CI, Prop.FOI.inap.presym.CI, Prop.FOI.inap.postsym.CI,
                         Prop.FOI.presym.of.ap.all.CI,Prop.FOI.postsym.detected.all.CI,Prop.FOI.postsym.detected.low.all.CI,Prop.FOI.postsym.detected.high.all.CI)
}

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # CREDITS: John Kruschke
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = floor( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

get.DIC <- function(LL){
  Deviance = -2 * LL
  expectation = mean(Deviance)   # posterior mean of the deviance
  effective.parameters = 2 * var(LL)  # as per Gelman
  print(effective.parameters)
  DIC = expectation + effective.parameters
  return(DIC)
}
