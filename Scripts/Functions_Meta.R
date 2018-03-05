
# Functions ---------------------------------------------------------------
run.mcmc <- function(initial, nbIteration, k, N, sdProposal, verbose = 0){
  if (verbose == 1) { browser() }
  logLik<-rep(0,nbIteration)
  accept<-matrix(0,ncol=length(initial),nrow=nbIteration)
  parameters<-matrix(0,ncol=length(initial),nrow=nbIteration)
  
  iteration<-1
  parameters[iteration,]<-initial
  
  logLik[iteration]<-calculateLL(parameters[iteration,],k, N, verbose )
  
  for(iteration in 2:nbIteration)
  {
    parameters[iteration,]<-parameters[iteration-1,]
    logLik[iteration]<-logLik[iteration-1]
    for(iParam in 1:length(initial))
    {
      oldParam<-parameters[iteration,iParam]
      newParam<-oldParam*exp(sdProposal[iParam]*rnorm(1))
      if (iParam == 1 & (newParam <=0 | newParam >1) ) {
        parameters[iteration,iParam]<-oldParam; break();
      } else if (iParam == 2 & newParam <= 0) {
        parameters[iteration,iParam]<-oldParam; break();
      }
      parameters[iteration,iParam]<-newParam
      newLogLik <- calculateLL(parameters[iteration,], k, N) # + calculatePriorLL(parameters[iteration,])
      if(log(runif(1))<newLogLik-logLik[iteration]+log(newParam)-log(oldParam))
      {
        logLik[iteration]<-newLogLik
        accept[iteration,iParam]<-1
      }
      else { parameters[iteration,iParam]<-oldParam     }
    }
  }
  DIC <- getDIC(logLik,parameters = apply(parameters,2,mean), k, N)
  return(list(logLik,parameters,accept,DIC[[1]],DIC[[2]]))
}

calculateLL <- function(param, k, N, verbose = 0){
  if(verbose==1) { browser() }
  prob.inapp <- param[1]
  if (length(param)==1){
  LOG.LIK <- sum(dbinom(k, prob = prob.inapp, size = N, log = TRUE)) 
  } else if (length(param)==2) {
    overdispersion <- param[2]
    LOG.LIK <- sum(dbetabinom(k, prob = prob.inapp, size = N,  theta = overdispersion, log = TRUE))
  }
  # if (is.na(LOG.LIK)) { LOG.LIK = -2000}
  return(LOG.LIK)
}

getDIC <- function(LL, parameters = NULL, k = NULL, N = NULL){
  #   browser()
  Deviance = -2 * LL
  expectation = mean(Deviance)   # posterior mean of the deviance
  if (is.null(parameters)){
    effective.parameters = 2 * var(LL)  # as per Gelman
    print(effective.parameters)
    DIC = expectation + effective.parameters
  } else {
    LL.hat = calculateLL(parameters,k, N)
    deviance.hat = -2 * LL.hat
    effective.parameters = expectation - deviance.hat  # as per Spiegelhalter
    print(effective.parameters)
    DIC = -2 * LL.hat + 2 * effective.parameters
  } 
  return(list(round(DIC,2),round(effective.parameters,2)))
}

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # CREDITS: 
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


