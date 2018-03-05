# Script to calculate the contribution to the total FOI by infection class given an foi
# includes different scenarios for sensitivity analysis

# Load packages -----------------------------------------------------------
library(deSolve)
library(data.table)
library(abind)
library(coda)
library(foreach)
library(doParallel)

# Load source files -------------------------------------------------------
source('Functions_FOICalculations.R')
source('Parameters.R')

# set running parameters --------------------------------------------------
Ca = c(median(AUC.asympto.1), median(AUC.asympto.2))   # note, we assume infectiousness of tertiary and quartenary infections to equal secondary ones
Cm = c(median(AUC.sympto.w.1), median(AUC.sympto.w.2))
Cs.sup = c(mean(AUC.sympto.w.3,na.rm=T), mean(AUC.sympto.w.3)) # for supplementary that includes DHF cases
Cs = Cm # ignoring differential effects of DHF or DF cases. 

Ca.all = cbind(AUC.asympto.1, AUC.asympto.2) 
Cm.all = cbind(AUC.sympto.w.1, AUC.sympto.w.2) 
Cs.all.sup = cbind(AUC.sympto.w.3, AUC.sympto.w.3) 
Cs.all = Cm.all

Prim.dis.inap <- c(inapp.ratio.rev[1], (1-inapp.ratio.rev[1])*(1-prop.DHF[1]), (1-inapp.ratio.rev[1])*prop.DHF[1]) 
Sec.dis.inap <- c(inapp.ratio.rev[2], (1-inapp.ratio.rev[2])*(1-prop.DHF[2]), (1-inapp.ratio.rev[2])*prop.DHF[2])
Tert.dis.inap <- c(inapp.ratio.rev[3], (1-inapp.ratio.rev[3])*(1-prop.DHF[3]), (1-inapp.ratio.rev[3])*prop.DHF[3])  

# Derive Proportions disease of total prevalence per immune history ---------------------------
X.brazil <- Disease.prop.calculator(infec.hist.brazil,Prim.dis, Sec.dis, Tert.dis)
X.thailand <- Disease.prop.calculator(infec.hist.thailand,Prim.dis, Sec.dis, Tert.dis)
X.brazil.low <- Disease.prop.calculator(infec.hist.brazil,Prim.dis.low, Sec.dis.low, Tert.dis)
X.thailand.low <- Disease.prop.calculator(infec.hist.thailand,Prim.dis.low, Sec.dis.low, Tert.dis)
X.brazil.high <- Disease.prop.calculator(infec.hist.brazil,Prim.dis.high, Sec.dis.high, Tert.dis)
X.thailand.high <- Disease.prop.calculator(infec.hist.thailand,Prim.dis.high, Sec.dis.high, Tert.dis)
X.brazil.emerging <- Disease.prop.calculator(infec.hist.brazil.emerging,Prim.dis, Sec.dis, Tert.dis)
X.thailand.emerging <- Disease.prop.calculator(infec.hist.thailand.emerging,Prim.dis, Sec.dis, Tert.dis)
X.brazil.inap <- Disease.prop.calculator(infec.hist.brazil,Prim.dis.inap, Sec.dis.inap, Tert.dis.inap)
X.thailand.inap <- Disease.prop.calculator(infec.hist.thailand,Prim.dis.inap, Sec.dis.inap, Tert.dis.inap)

# Derive contribution symptomatics to total foi as a function of f means --------
Prop.FOI.brazil <- Sympto.FOI.calculator(prev = infec.hist.brazil$prev.vec2,X = X.brazil,Ca,Cm,Cs,2, foi.vec)
Prop.FOI.thailand <- Sympto.FOI.calculator(prev = infec.hist.thailand$prev.vec2,X = X.thailand,Ca,Cm,Cs,2, foi.vec)
Prop.FOI.brazil.low <- Sympto.FOI.calculator(prev = infec.hist.brazil$prev.vec2,X = X.brazil.low,Ca,Cm,Cs,2, foi.vec)
Prop.FOI.thailand.low <- Sympto.FOI.calculator(prev = infec.hist.thailand$prev.vec2,X = X.thailand.low,Ca,Cm,Cs,2, foi.vec)
Prop.FOI.brazil.high <- Sympto.FOI.calculator(prev = infec.hist.brazil$prev.vec2,X = X.brazil.high,Ca,Cm,Cs,2, foi.vec)
Prop.FOI.thailand.high <- Sympto.FOI.calculator(prev = infec.hist.thailand$prev.vec2,X = X.thailand.high,Ca,Cm,Cs,2, foi.vec)
Prop.FOI.brazil.emerging <- Sympto.FOI.calculator(prev = infec.hist.brazil.emerging$prev.vec1,X = X.brazil.emerging,Ca,Cm,Cs,1, foi.vec)
Prop.FOI.thailand.emerging <- Sympto.FOI.calculator(prev = infec.hist.thailand.emerging$prev.vec1,X = X.thailand.emerging,Ca,Cm,Cs,1, foi.vec)
Prop.FOI.brazil.DHF <- Sympto.FOI.calculator(prev = infec.hist.brazil$prev.vec2,X = X.brazil,Ca,Cm,Cs.sup,2, foi.vec)
Prop.FOI.thailand.DHF <- Sympto.FOI.calculator(prev = infec.hist.thailand$prev.vec2,X = X.thailand,Ca,Cm,Cs.sup,2, foi.vec)
Prop.FOI.brazil.inap <- Sympto.FOI.calculator(prev = infec.hist.brazil$prev.vec2,X = X.brazil.inap,Ca,Cm,Cs.sup,2, foi.vec)
Prop.FOI.thailand.inap <- Sympto.FOI.calculator(prev = infec.hist.thailand$prev.vec2,X = X.thailand.inap,Ca,Cm,Cs.sup,2, foi.vec)

# Derive contribution symptomatics to total foi with uncertainty ----------
Prop.FOI.All.brazil <- Sympto.FOI.sampler(prev = infec.hist.brazil$prev.vec2,X = X.brazil, 2, foi.vec, Ca.all, Cm.all, Cs.all, cores = 3)
Prop.FOI.All.thailand <- Sympto.FOI.sampler(prev = infec.hist.brazil$prev.vec2,X = X.thailand, 2, foi.vec, Ca.all, Cm.all, Cs.all, cores = 3)

# Derive contribution symptomatics to total foi in emerging setting with uncertainty ----------
Prop.FOI.All.brazil.emerging <- Sympto.FOI.sampler(prev = infec.hist.brazil.emerging$prev.vec1,X = X.brazil.emerging, 1, foi.vec, Ca.all, Cm.all, Cs.all, cores = 3)
Prop.FOI.All.thailand.emerging <- Sympto.FOI.sampler(prev = infec.hist.brazil.emerging$prev.vec1,X = X.thailand.emerging, 1, foi.vec, Ca.all, Cm.all, Cs.all, cores = 3)

# Derive contribution sympomatics to total foi including DHF contribution--------------
Prop.FOI.All.brazil.DHF <- Sympto.FOI.sampler(prev = infec.hist.brazil$prev.vec2,X = X.brazil, 2, foi.vec, Ca.all, Cm.all, Cs.all.sup, cores = 3)
Prop.FOI.All.thailand.DHF <- Sympto.FOI.sampler(prev = infec.hist.brazil$prev.vec2,X = X.thailand, 2, foi.vec, Ca.all, Cm.all, Cs.all.sup, cores = 3)

Prop.FOI.All.brazil.DHF.severe <- Sympto.FOI.sampler(prev = infec.hist.brazil$prev.vec2,X = X.brazil, 2, foi.vec, Ca.all, Cm.all, Cs.all.sup, cores = 3, severe = TRUE)
Prop.FOI.All.thailand.DHF.severe <- Sympto.FOI.sampler(prev = infec.hist.thailand$prev.vec2,X = X.thailand, 2, foi.vec, Ca.all, Cm.all, Cs.all.sup, cores = 3, severe = TRUE)

# Derive contribution sympomatics to total foi assuming inapparent infections are similar to asymptomatics--------------
Prop.FOI.All.brazil.inap <- Sympto.FOI.sampler(prev = infec.hist.brazil$prev.vec2,X = X.brazil.inap, 2, foi.vec, Ca.all, Cm.all, Cs.all.sup, cores = 3)
Prop.FOI.All.thailand.inap <- Sympto.FOI.sampler(prev = infec.hist.brazil$prev.vec2,X = X.thailand.inap, 2, foi.vec, Ca.all, Cm.all, Cs.all.sup, cores = 3)







