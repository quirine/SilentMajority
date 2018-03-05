# Load packages -----------------------------------------------------------
library(deSolve)
library(data.table)
library(abind)
library(coda)
library(foreach)
library(doParallel)

# Load source files -------------------------------------------------------
source('Functions_FOICalculations.R')
source('Parameters_ViremiatoInfectiousness.R')

# set running parameters --------------------------------------------------
Ca.4 = c(median(AUC.asympto.1), rep(median(AUC.asympto.2,na.rm=T),3) )   # note, we assume infectiousness of tertiary and quartenary infections to equal secondary ones
Cm.4 = c(median(AUC.sympto.w.1), rep(median(AUC.sympto.w.2,na.rm=T),3))
Cs.4 = c(median(AUC.sympto.w.3), rep(median(AUC.sympto.w.3,na.rm=T),3)) 
Ca.all.4 = cbind(AUC.asympto.1, AUC.asympto.2, AUC.asympto.2, AUC.asympto.2) 
Cm.all.4 = cbind(AUC.sympto.w.1, AUC.sympto.w.2, AUC.sympto.w.2, AUC.sympto.w.2) 
Cs.all.4 = cbind(AUC.sympto.w.3, AUC.sympto.w.3, AUC.sympto.w.3, AUC.sympto.w.3) 

# Derive Proportions disease of total prevalence per immune history ---------------------------
X.brazil.4 <- Disease.prop.calculator(infec.hist.brazil,Prim.dis, Sec.dis, Tert.dis)
X.thailand.4 <- Disease.prop.calculator(infec.hist.thailand,Prim.dis, Sec.dis, Tert.dis)

# Derive contribution symptomatics to total foi as a function of f means --------
Prop.FOI.brazil.4 <- Sympto.FOI.calculator(prev = infec.hist.brazil$prev.vec4,X = X.brazil.4,Ca.4,Cm.4,Cs.4,4,foi.vec)
Prop.FOI.thailand.4 <- Sympto.FOI.calculator(prev = infec.hist.thailand$prev.vec4,X = X.thailand.4,Ca.4,Cm.4,Cs.4,4,foi.vec)

# Derive contribution symptomatics to total foi with uncertainty ----------
Prop.FOI.All.brazil.4 <- Sympto.FOI.sampler(prev = infec.hist.brazil$prev.vec4,X = X.brazil.4,4,foi.vec,Ca.all.4,Cm.all.4,Cs.all.4)
Prop.FOI.All.thailand.4 <- Sympto.FOI.sampler(prev = infec.hist.thailand$prev.vec4,X = X.thailand.4,4,foi.vec,Ca.all.4,Cm.all.4,Cs.all.4)




