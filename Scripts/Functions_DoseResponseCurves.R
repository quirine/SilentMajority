Ferguson = function(log_V, inf.dose, slope){
  exponent = - (log_V / inf.dose)^slope
  prob_inf = 1 - exp(exponent)
  return(prob_inf)
}

Hill = function(log_V, dis.const, hill.coef){
  prob_inf = (log_V^hill.coef) / (dis.const + log_V^hill.coef)
  return(prob_inf)
}

LogFunction = function(log_V, tweek.parameter1, tweek.parameter2){
  prob_inf = 1 / (1 + exp(-1*(tweek.parameter1 + tweek.parameter2*log_V) ) ) 
  return(prob_inf)
}


FergusonNLL = function(params,k,N,log_V, Serotypes = NULL){ 
  # browser()
  if (is.null(Serotypes)) { # if this vector is not given, asume all to be serotype 1
    Serotypes = rep(1,length(k))
    k.1 = k
    N.1 = N
  } 
  slope <- params[length(params)]
  
  prob_inf <- k.all <- N.all <- numeric()
  DENVs = sort(unique(Serotypes)) 
  for (ii in DENVs){
    ind = which(Serotypes == ii)
    to.eval = paste('inf.dose.',ii,' <- params[',ii,']',sep='')
    eval(parse(text = to.eval))
    to.eval = paste('prob_inf.',ii,' = Ferguson(log_V[ind], inf.dose.',ii,', slope)',sep='')  
    eval(parse(text = to.eval))
    to.eval = paste('prob_inf = c(prob_inf, prob_inf.',ii,')',sep='')
    eval(parse(text = to.eval))
    to.eval = paste('k.',ii,' = k[ind]; k.all = c(k.all, k.',ii,')',sep='')
    eval(parse(text = to.eval))
    to.eval = paste('N.',ii,' = N[ind]; N.all = c(N.all, N.',ii,')',sep='')
    eval(parse(text = to.eval))
  }
  
  
  NLL = -sum(dbinom(k.all, prob = prob_inf, size = N.all, log = TRUE)) 
  print(NLL)
  if(is.na(NLL)) { NLL = 2000}
  if(NLL == -Inf) { NLL = 2000}
  return(NLL)
}

parnames(FergusonNLL) <- c("inf.dose.1","inf.dose.2","inf.dose.3","inf.dose.4", "slope")

FergusonBetaBinomNLL = function(params,k,N,log_V){ 
  
  inf.dose <- params[1]
  slope <- params[2]
  overdispersion <- params[3]
  print(params)
  prob_inf = Ferguson(log_V, inf.dose, slope)
  
  NLL = -sum(dbetabinom(k, prob = prob_inf, size = N, theta = overdispersion, log = TRUE))
  # print(NLL) 
  if (is.nan(NLL) ) { NLL = 2000 } #{browser()} #
  return (NLL)
}

parnames(FergusonBetaBinomNLL) <- c("inf.dose", "slope","overdispersion")

HillNLL = function(params,k,N,log_V){ 
  # browser()
  
  dis.const <- params[1]
  hill.coef <- params[2]
  
  
  prob_inf = Hill(log_V, dis.const, hill.coef)
  
  -sum(dbinom(k, prob = prob_inf, size = N, log = TRUE)) 
  
}

parnames(HillNLL) <- c("dis.const", "hill.coef")

HillBetaBinomNLL = function(params,k,N,log_V){ 
  # browser()
  dis.const <- params[1]
  hill.coef <- params[2]
  overdispersion <- params[3]
  
  prob_inf = Hill(log_V, dis.const, hill.coef)
  
  NLL = -sum(dbetabinom(k, prob = prob_inf, size = N,  theta = overdispersion, log = TRUE))
  if (is.nan(NLL) ) { NLL = 2000 }
  print(NLL)
  return (NLL)
  
}

parnames(HillBetaBinomNLL) <- c("dis.const", "hill.coef", "overdispersion")

LogNLL = function(params,k,N,log_V){ 
  # browser()
  
  tweek.parameter1 <- params[1]
  tweek.parameter2 <- params[2]

  prob_inf = LogFunction(log_V, tweek.parameter1, tweek.parameter2)
  
  -sum(dbinom(k, prob = prob_inf, size = N, log = TRUE)) 
  
}

parnames(LogNLL) <- c("tweek.parameter1", "tweek.parameter2")

LogNLL.ST = function(params,k,N,log_V,Serotypes = NULL){ 
  # browser()
  if (is.null(Serotypes)) { # if this vector is not given, asume all to be serotype 1
    Serotypes = rep(1,length(k))
    k.1 = k
    N.1 = N
  } 
  intercept <- params[length(params)]
  
  prob_inf <- k.all <- N.all <- numeric()
  DENVs = sort(unique(Serotypes)) 
  for (ii in DENVs){
    ind = which(Serotypes == ii)
    to.eval = paste('coef.',ii,' <- params[',ii,']',sep='')
    eval(parse(text = to.eval))
    to.eval = paste('prob_inf.',ii,' = LogFunction(log_V[ind], intercept, coef.',ii,')',sep='')  
    eval(parse(text = to.eval))
    to.eval = paste('prob_inf = c(prob_inf, prob_inf.',ii,')',sep='')
    eval(parse(text = to.eval))
    to.eval = paste('k.',ii,' = k[ind]; k.all = c(k.all, k.',ii,')',sep='')
    eval(parse(text = to.eval))
    to.eval = paste('N.',ii,' = N[ind]; N.all = c(N.all, N.',ii,')',sep='')
    eval(parse(text = to.eval))
  }
  
  
  NLL = -sum(dbinom(k.all, prob = prob_inf, size = N.all, log = TRUE)) 
  print(NLL)
  if(is.na(NLL)) { NLL = 2000}
  if(NLL == -Inf) { NLL = 2000}
  return(NLL)
  
}

parnames(LogNLL.ST) <- c( "coef.1", "coef.2", "coef.3", "coef.4","intercept")

LogBetaBinomNLL = function(params,k,N,log_V){ 
  # browser()
  
  tweek.parameter1 <- params[1]
  tweek.parameter2 <- params[2]
  overdispersion <- params[3]
  
  prob_inf = LogFunction(log_V, tweek.parameter1, tweek.parameter2)
  
  NLL = -sum(dbetabinom(k, prob = prob_inf, size = N,  theta = overdispersion, log = TRUE))
  if (is.nan(NLL) ) { NLL = 2000 }
  return (NLL) 
  
}

parnames(LogBetaBinomNLL) <- c("tweek.parameter1", "tweek.parameter2", "overdispersion")

get.binom.model.characteristics <- function(model, type = 'binom'){
  p1.m = model@coef[1]
  p1.up = model@coef[1] + 1.96 * summary(model)@coef[1,2]
  p1.low = model@coef[1] - 1.96 * summary(model)@coef[1,2]
  p1.pv = summary(model)@coef[1,4]
  p2.m = model@coef[2]
  p2.up = model@coef[2] + 1.96 * summary(model)@coef[2,2]
  p2.low = model@coef[2] - 1.96 * summary(model)@coef[2,2]
  p2.pv = summary(model)@coef[2,4]
  LL = summary(model)@m2logL
  if (type == 'betabinom'){
    od.m  = model@coef[3]
    od.up = model@coef[3] + 1.96 * summary(model)@coef[3,2]
    od.low = model@coef[3] - 1.96 * summary(model)@coef[3,2]
    od.pv = summary(model)@coef[1,4]
    return(data.frame(p1.m, p1.low, p1.up, p1.pv,p2.m, p2.low, p2.up, p2.pv, od.m, od.low, od.up, od.pv, LL ))
  } else {
    return(data.frame(p1.m, p1.low, p1.up, p1.pv, p2.m, p2.low, p2.up, p2.pv, LL ))    
  }
}

AIC.bb <- function(model, data.set, observed.data){
  Residuals = model$fitted.values - observed.data
  AIC = nrow(data.set)*(log(2*pi)+1+log((sum(Residuals^2)/nrow(data.set)))) + ((length(model$coefficients)+1)*2)
  return(AIC)
}

lrt <- function (LL0, LL1, df0, df1) {
    L01 <- as.vector( 2 * (LL0 - LL1)) # assuming Negative log likelihoods
    df <- df1 - df0
    list(L01 = L01, df = df,
                  "p-value" = pchisq(L01, df, lower.tail = FALSE))
}

lrt.dev <- function (dev0, dev1, df0, df1) {
  L01 <- as.vector(dev0 - dev1)
  df <- df1 - df0
  list(L01 = L01, df = df,
       "p-value" = pchisq(L01, df, lower.tail = FALSE))
}
