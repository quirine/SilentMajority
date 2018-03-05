
# Samplers ----------------------------------------------------------------

within.host.sampler <- function(dis.cat,n,joinedposterior=NaN){
  Viremia.mod <- numeric()
  if (typeof(joinedposterior) == 'double') {
    Index1 <- sample(6000:dim(joinedposterior)[1], n, replace = TRUE)  
    Index2 <- sample(1:dim(joinedposterior)[2], n, replace = TRUE)  
    Z1.samples = joinedposterior[Index1, Index2, 8] 
    Beta.samples = joinedposterior[Index1, Index2, 1]
    Kappa.samples = joinedposterior[Index1, Index2, 6]
    Nu.samples = joinedposterior[Index1, Index2, 4]
    IIP = joinedposterior[Index1, Index2, 5]
  }
  
  else {
    Z1.samples = positivesampler(n, Z1.mean[dis.cat], sd=Z1.sd[dis.cat], n*100)
    Beta.samples = rnorm(n, mean = Beta[dis.cat], sd = Beta.sd[dis.cat])
    Kappa.samples = rnorm(n, mean = Kappa[dis.cat], sd = Kappa.sd[dis.cat])
    Nu.samples = positivesampler(n, Nu.mean[dis.cat], sd=Nu.sd[dis.cat], n*100)
  }
  for (ii in 1:n){
    z1 <-Z1.samples[ii,ii]   
    Init <- c(X1,Y1,V1,z1)
    Out.1 <- ode(Init,Time,dengue.within.host.model,parms=c(A, gamma, Beta.samples[ii,ii], delta, alpha, omega, Kappa.samples[ii,ii], Nu.samples[ii,ii]))
    Out.1[Out.1<0] <- 0
    Viremia.mod <- cbind(Viremia.mod,Out.1[,4]) 
  }
  Viremia.mod <- log10(Viremia.mod)
  Viremia.mod[is.infinite(Viremia.mod)] <- 0
  
  return(rbind(Viremia.mod, 
               diag(IIP), diag(Z1.samples), diag(Beta.samples),diag(Kappa.samples), diag(Nu.samples)))
}

remove.reignition <- function(Viremia){
  # browser()
  for (ii in 1:n){
    idx <- which(Viremia[,ii] < -5) # find the peak 
    if (length(idx) >0 ){
      Viremia[idx[1]:nrow(Viremia),ii] <- 0
    }
  }
  return(Viremia)
}

state.ratios.sampler <- function(type,n){
  Ratios <- numeric()
  ind <- which(Duong.viremia$Type==type)
  rtype <- rnorm(n = n, mean = Duong.viremia$Mean[ind], sd = Duong.viremia$SE[ind])
  rref <- rnorm(n = n, mean = Duong.viremia$Mean[4], sd = Duong.viremia$SE[4])
  Ratios <- rtype / rref

  return(Ratios)
}

state.viremia.sampler <- function(type,dis.cat,n){
  Viremia <- numeric()
  if (dis.cat == 1) Viremia.mod = Viremia.mod1
  else if (dis.cat == 2) Viremia.mod = Viremia.mod2
  else if (dis.cat == 3) Viremia.mod = Viremia.mod3
  if (type == 'asym') Ratios = Ratios.asym
  else if (type =='presym') Ratios = rep(1,n)
  else if (type =='sym') Ratios = rep(1,n)
  for (ii in 1:n){
  Viremia <- cbind(Viremia, Ratios[ceiling(runif(1,1,n))] * Viremia.mod[,ceiling(runif(1,1,n))])
  }
  return(Viremia)
}

reg.coefs.sampler <- function(type,n,Duong.infprob){
  Coefs <- Ints <- numeric();
  ind <- which(Duong.infprob$Type==type)
  Ints <- rnorm(n = n, mean = Duong.infprob$Int[ind], sd = Duong.infprob$Int.se[ind])
  Coefs <-rnorm(n = n, mean = Duong.infprob$Coefs[ind], sd = Duong.infprob$Coefs.se[ind])
  return(data.frame(Ints,Coefs))
}

infectiousness.sampler <- function(type,dis.cat,n,gender='no'){
  Prob <- numeric()
  Viremia <- state.viremia.sampler(type,dis.cat,n)
  if (type == 'sym' & gender == 'no') Reg = Sym.reg
  else if (type =='asym' & gender == 'no') Reg = Asym.reg
  else if (type =='presym' & gender == 'no') Reg = Presym.reg
  else if (type == 'sym' & gender == 'male') Reg = Sym.reg.male
  else if (type =='asym' & gender == 'male') Reg = Asym.reg.male
  else if (type =='presym' & gender == 'male') Reg = Presym.reg.male
  else if (type == 'sym' & gender == 'female') Reg = Sym.reg.female
  else if (type =='asym' & gender == 'female') Reg = Asym.reg.female
  else if (type =='presym' & gender == 'female') Reg = Presym.reg.female
  
  for (ii in 1:n){
  Prob  <- cbind(Prob, logistic.function.infectiousness(Reg[ii,1],Reg[ii,2],Viremia[,ii]))
  }
  return(Prob)
}

IIP.sampler <- function(n, params = Chan.IIP, dist = 'log.norm'){
  Prob.IIP <- numeric()  
  Betas <- rnorm(n, mean = params$beta.mean, sd = (params$beta.up - params$beta.low) / 3.9)
  Taus <- rnorm(n, mean = params$tau.mean, sd = (params$tau.up - params$tau.low) / 3.9)
  if (dist == 'log.norm'){
    for (ii in 1:n){
      tau <- Taus[ii];     mu <- exp(Betas[ii])
      Prob.IIP <-cbind(Prob.IIP, dlnorm(Time,meanlog=mu, sdlog = 1/sqrt(tau)))
    }
  }
  else if (dist == 'weib'){
    for (ii in 1:n){
      tau <- Taus[ii];     lambda <- exp(Betas[ii])
      Prob.IIP <-cbind(Prob.IIP, dweibull(Time, scale = lambda^(-1/tau), shape = tau))
    }
  }
  else if (dist == 'gamma'){
    for (ii in 1:n){
      tau <- Taus[ii];     lambda <- exp(Betas[ii])
      Prob.IIP <-cbind(Prob.IIP, dgamma(Time, shape = tau, rate = tau/lambda))
    }
  }
  else if (dist == 'exp'){
    for (ii in 1:n){
      mu <- exp(Betas[ii])
      Prob.IIP <-cbind(Prob.IIP, dexp(Time, rate = 1/mu))
    }
  }
  
  return(Prob.IIP)
}

Weighting.symptos.iip <- function(n,dis.cat,Time,gender='no',IIP.vector = NaN) {
  if (dis.cat == 1 & gender == 'no') { Prob.sympto.mod = Prob.sympto.mod1;  Prob.presympto.mod = Prob.presympto.mod1; }
  else if (dis.cat == 2 & gender == 'no') { Prob.sympto.mod = Prob.sympto.mod2;   Prob.presympto.mod = Prob.presympto.mod2; }
  else if (dis.cat == 3 & gender == 'no') { Prob.sympto.mod = Prob.sympto.mod3;  Prob.presympto.mod = Prob.presympto.mod3; }
  else if (dis.cat == 1 & gender == 'male') { Prob.sympto.mod = Prob.sympto.mod1.male;  Prob.presympto.mod = Prob.presympto.mod1.male; }
  else if (dis.cat == 2 & gender == 'male') { Prob.sympto.mod = Prob.sympto.mod2.male;  Prob.presympto.mod = Prob.presympto.mod2.male; }
  else if (dis.cat == 3 & gender == 'male') { Prob.sympto.mod = Prob.sympto.mod3.male;  Prob.presympto.mod = Prob.presympto.mod3.male; }
  else if (dis.cat == 1 & gender == 'female') { Prob.sympto.mod = Prob.sympto.mod1.female;  Prob.presympto.mod = Prob.presympto.mod1.female; }
  else if (dis.cat == 2 & gender == 'female') { Prob.sympto.mod = Prob.sympto.mod2.female;  Prob.presympto.mod = Prob.presympto.mod2.female; }
  else if (dis.cat == 3 & gender == 'female') { Prob.sympto.mod = Prob.sympto.mod3.female;  Prob.presympto.mod = Prob.presympto.mod3.female; }
  if (is.nan(IIP.vector)) {
    Prop.prior.to.symptoms <- numeric()
    Weighted_Sympto <- matrix(data = NA, nrow = length(Time), ncol = n)
    someData <- rep(NaN, length(Time)*(length(Time)-1));  
    Sympto  <- array(someData,c(length(Time),(length(Time)-1)))
    for (jj in 1:n) {
      Temp <- matrix(data = NA, nrow = length(Time), ncol = length(Time)-1);
      for (ii in 1:(length(Time)-1)){
        Temp[,ii] <- c(Prob.presympto.mod[1:ii,jj], Prob.sympto.mod[(ii+1):length(Time),jj])
      }
      Sympto <- abind(Sympto,Temp,along=3)
      
    }
    Sympto <- Sympto[,,-1]
    for (jj in 1:n){
      Temp <- numeric()
      idx <- sample(seq_len(length(Time)-1), 1, prob=Prob.IIP[1:length(Time)-1,jj])
      Temp <- Sympto[,idx,jj]
      Weighted_Sympto[,jj] <- Temp
      Prop.prior.to.symptoms <- c(Prop.prior.to.symptoms, sum(diff(Time[1:idx])*rollmean(Temp[1:idx],2)) / sum(diff(Time)*rollmean(Temp,2)) )
    }
  }
  if (typeof(IIP.vector) == 'double'){
    Prop.prior.to.symptoms <- numeric()
    Weighted_Sympto <- matrix(data = NA, nrow = length(Time), ncol = n)
    for (jj in 1:n){
      idx <- round(IIP.vector[jj]/diff(Time)[1]) + 1
      Temp <- c(Prob.presympto.mod[1:idx,jj], Prob.sympto.mod[(idx+1):length(Time),jj])
      Weighted_Sympto[,jj] <- Temp
      Prop.prior.to.symptoms <- c(Prop.prior.to.symptoms, sum(diff(Time[1:idx])*rollmean(Temp[1:idx],2)) / sum(diff(Time)*rollmean(Temp,2)) )
    }
  }
  if (IIP.vector == 'use.peaks'){
    Prop.prior.to.symptoms <- numeric()
    Weighted_Sympto <- matrix(data = NA, nrow = length(Time), ncol = n)
    Order <- order(colSums(Prob.sympto.mod))
    for (jj in 1:n){
      idx <- which(Prob.sympto.mod[,jj]==max(Prob.sympto.mod[,jj]))
      Temp <- c(Prob.presympto.mod[1:idx,Order[jj]], Prob.sympto.mod[(idx+1):length(Time),jj])
      Weighted_Sympto[,jj] <- Temp
      Prop.prior.to.symptoms <- c(Prop.prior.to.symptoms, sum(diff(Time[1:idx])*rollmean(Temp[1:idx],2)) / sum(diff(Time)*rollmean(Temp,2)) )
    }
  }
  return(rbind(Weighted_Sympto,Prop.prior.to.symptoms))

}   


# Nested functions --------------------------------------------------------

dengue.within.host.model <- function (t, x, params) {
  
  X <- x[1];  Y <- x[2];  V <- x[3];   Z <- x[4];   
  A <- params[1];      gamma <- params[2];  beta <- params[3];     delta <- params[4]
  alpha <- params[5];  omega <- params[6];  kappa <- params[7];    nu <- params[8]
  
  # differential equations
  dX <- A - gamma*X - beta*X*V
  dY <- beta*X*V - delta*Y - alpha*Z*Y
  dV <- omega*Y - kappa*V
  dZ <- nu*Y*Z
  
  list(c(dX,dY,dV,dZ))
}


logistic.function.infectiousness <- function(tweek.parameter1, tweek.parameter2, viremia) {
    Prob <- 1 / (1 + exp(-1*(tweek.parameter1 + tweek.parameter2*viremia) ) ) 
    return(Prob)
}


positivesampler <- function(n, mean, sd, nnorm){
  samp <- rnorm(nnorm, mean, sd)
  samp <- samp[samp >= 0]
  if (length(samp) >= n) {
    return(sample(samp, n))
  }  
  stop(simpleError("Not enough values to sample from. Try increasing nnorm."))
}

# Miscelaneous ------------------------------------------------------------
# source:
# Wan, Xiang, et al. "Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range." 
# BMC medical research methodology 14.1 (2014): 135.

sd.calculator <- function(IQ.low, IQ.up){
  SD <- (IQ.up - IQ.low) / 1.35 
  return(SD)
}

mean.calculator <- function(Median, IQ.low, IQ.up){
  Mean <- (IQ.low + Median + IQ.up) / 3
  return(Mean)
}


