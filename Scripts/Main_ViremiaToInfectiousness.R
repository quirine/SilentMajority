# Code to estimate net-infectiousness by infection class
# Based on Duong 2015, Clapham 2013 and Chan 2012

rm(list=ls())

# Load packages -----------------------------------------------------------
library(deSolve)
library(data.table)
library(abind)
library(zoo)
library(foreach)
library(doParallel)
library(mgcv)

# Load source files -------------------------------------------------------
source('Functions_ViremiaToInfectiousness.R')
source('Parameters.R')

# Load workspaces ---------------------------------------------------------
load("../Data/WithinHostPosteriors.RData")
load('../Data/DoseResponseCurves.RData')

# set running parameters ------------------------------------------------
n = 3000   # number of samples to run for each uncertainty exercise
MaxTime <- 20 #15
Time <- seq(0,MaxTime,by=0.05)
set.seed(36)

# Derive viremia of symptomatic cases as a function of time ------------------------------------
Viremia.mod1 <- within.host.sampler(1,n,posteriors.primary) # for primary DF
NU.1 <- Viremia.mod1[nrow(Viremia.mod1),]
KAPPA.1 <- Viremia.mod1[nrow(Viremia.mod1)-1,]
BETA.1 <- Viremia.mod1[nrow(Viremia.mod1)-2,]
Z1.1 <- Viremia.mod1[nrow(Viremia.mod1)-3,]
IIP.1 <- Viremia.mod1[nrow(Viremia.mod1)-4,]
Viremia.mod1 <- Viremia.mod1[-(nrow(Viremia.mod1)-seq(0,4)),]

Viremia.mod2 <- within.host.sampler(2,n,posteriors.secondary.mild)   # secondary DF
NU.2 <- Viremia.mod2[nrow(Viremia.mod2),]
KAPPA.2 <- Viremia.mod2[nrow(Viremia.mod2)-1,]
BETA.2 <- Viremia.mod2[nrow(Viremia.mod2)-2,]
Z1.2 <- Viremia.mod2[nrow(Viremia.mod2)-3,]
IIP.2 <- Viremia.mod2[nrow(Viremia.mod2)-4,]
Viremia.mod2 <- Viremia.mod2[-(nrow(Viremia.mod2)-seq(0,4)),]

Viremia.mod3 <- within.host.sampler(3,n,posteriors.secondary.severe)   # secondary DHF
NU.3 <- Viremia.mod3[nrow(Viremia.mod3),]
KAPPA.3 <- Viremia.mod3[nrow(Viremia.mod3)-1,]
BETA.3 <- Viremia.mod3[nrow(Viremia.mod3)-2,]
Z1.3 <- Viremia.mod3[nrow(Viremia.mod3)-3,]
IIP.3 <- Viremia.mod3[nrow(Viremia.mod3)-4,]
Viremia.mod3 <- Viremia.mod3[-(nrow(Viremia.mod3)-seq(0,4)),]

# Sample viremia scaling factors for asymptomatic and presymptomatic cases ----------------
Ratios.asym <- state.ratios.sampler('asym',n)

# Derive viremia of asymptomatic and presymptomatic cases (not neccessary, because not in nested function, but helpful for plotting)-------
Viremia.mod1.asym <- state.viremia.sampler('asym',1,n)
Viremia.mod2.asym <- state.viremia.sampler('asym',2,n)

# Sample parameters viremia to infectiousness regression model--------------------------
MV.Coefs.sym = rmvn(n=n, mu = model.symp$coefficients, V <- vcov(model.symp))
MV.Coefs.asym = rmvn(n=n, mu = model.asymp$coefficients, V <- vcov(model.asymp))
MV.Coefs.presym = rmvn(n=n, mu = model.presymp$coefficients, V <- vcov(model.presymp))
Sym.reg <- cbind(MV.Coefs.sym[,1], MV.Coefs.sym[,2]); colnames(Sym.reg) <- c('Ints','Coefs')
Asym.reg <- cbind(MV.Coefs.asym[,1], MV.Coefs.asym[,2]); colnames(Asym.reg) <- c('Ints','Coefs')
Presym.reg <- cbind(MV.Coefs.presym[,1], MV.Coefs.presym[,2]); colnames(Presym.reg) <- c('Ints','Coefs')

# Derive infectiousness curves --------------------------------------------
# uses a nested function to derive viremic curves for asymptomatics and presymptomatics based on the ratios acquired earlier
Prob.sympto.mod1 <- infectiousness.sampler('sym',1,n)
Prob.sympto.mod2 <- infectiousness.sampler('sym',2,n)
Prob.sympto.mod3 <- infectiousness.sampler('sym',3,n)
Prob.asympto.mod1 <- infectiousness.sampler('asym',1,n)
Prob.asympto.mod2 <- infectiousness.sampler('asym',2,n)
Prob.presympto.mod1 <- infectiousness.sampler('presym',1,n)
Prob.presympto.mod2 <- infectiousness.sampler('presym',2,n)
Prob.presympto.mod3 <- infectiousness.sampler('presym',3,n)

# Get IIP Probabilistic curve from Chan and Johansson -------------------------------------------------
Prob.IIP <- IIP.sampler(n)  # redundant if in 'Weighting.symptos.iip' an IIP.vector is defined

# Derive symptomatic infectiousness curve from presympto and sympto --------
Prob.sympto.weighted.1 <- Weighting.symptos.iip(n,1,Time,IIP.vector = IIP.1 )
Prob.sympto.weighted.2 <- Weighting.symptos.iip(n,2,Time,IIP.vector = IIP.2 )
Prob.sympto.weighted.3 <- Weighting.symptos.iip(n,3,Time,IIP.vector = IIP.3 )

Prop.infectivity.prior.to.IIP.1 <- Prob.sympto.weighted.1[nrow(Prob.sympto.weighted.1),]
Prop.infectivity.prior.to.IIP.2 <- Prob.sympto.weighted.2[nrow(Prob.sympto.weighted.2),]
Prop.infectivity.prior.to.IIP.3 <- Prob.sympto.weighted.3[nrow(Prob.sympto.weighted.3),]

Prob.sympto.weighted.1 <- Prob.sympto.weighted.1[-nrow(Prob.sympto.weighted.1),]
Prob.sympto.weighted.2 <- Prob.sympto.weighted.2[-nrow(Prob.sympto.weighted.2),]
Prob.sympto.weighted.3 <- Prob.sympto.weighted.3[-nrow(Prob.sympto.weighted.3),]

# AUCs --------------------------------------------------------------------
id <- order(Time)

AUC.asympto.1 <- AUC.asympto.2 <- AUC.sympto.w.1 <- AUC.sympto.w.2 <- AUC.sympto.w.3 <- numeric()
for (ii in 1:n) {
 AUC.asympto.1 <- c(AUC.asympto.1, sum(diff(Time[id])*rollmean(Prob.asympto.mod1[id,ii],2)))
 AUC.asympto.2 <- c(AUC.asympto.2, sum(diff(Time[id])*rollmean(Prob.asympto.mod2[id,ii],2)))
 AUC.sympto.w.1 <- c(AUC.sympto.w.1, sum(diff(Time[id])*rollmean(Prob.sympto.weighted.1[id,ii],2)))
 AUC.sympto.w.2 <- c(AUC.sympto.w.2, sum(diff(Time[id])*rollmean(Prob.sympto.weighted.2[id,ii],2)))
 AUC.sympto.w.3 <- c(AUC.sympto.w.3, sum(diff(Time[id])*rollmean(Prob.sympto.weighted.3[id,ii],2)))
}

Ratio.1 <- AUC.asympto.1/AUC.sympto.w.1
Ratio.2 <- AUC.asympto.2/(AUC.sympto.w.2)  

# save workspace ----------------------------------------------------------
save(list=ls(),file='Workspace.lowerAsym.RData')


