## Parameters used to estimate infectiousness over time given disease status

# plotting parameters for logistic regression 
Viremia = seq(0,10,by=0.1)
Gender = seq(0,1,by=1)
Asymptomatic = seq(0,1,by=1)

# model parameters for within host model (as per Clapham 2013)
A <-    1.4e6   #1.4e7                 # target cell production per ml per day
gamma <- 0.14 #0.14                         # daily death rate uninfected cells 
Beta <- c(1.72e-10, 2.30e-10, 2.82e-10)   # infection rate (primary DF, secondary DF, secondary DHF)
Beta.low <- c(1.51e-10, 2.17e-10, 2.62e-10)   # infection rate, lower bounds of 95% Cred Int (primary DF, secondary DF, secondary DHF)
Beta.up <- c(2.04e-10, 2.43e-10, 3.01e-10)   # infection rate, upper bounds of 95% Cred Int (primary DF, secondary DF, secondary DHF)
delta <- 0.14                         # daily death rate infected cells    
alpha <- 0.001                        # removal rate of infected cells per immune cell per day 
omega <-  1e4                         # production rate of virions per infected cell per day
Kappa <- c(3.48, 5.29, 6.07)          # virion clearance rate per day
Kappa.low <- c(3.30, 5.08, 5.71)      # virion clearance rate per day, lower bounds of 95% Cred Int (primary DF, secondary DF, secondary DHF)
Kappa.up <- c(3.67, 5.52, 6.43)       # virion clearance rate per day, higher bounds of 95% Cred Int (primary DF, secondary DF, secondary DHF)
Nu <- c(1.29e-5, 2.95e-5, 2.71e-5)    # proliferation rate of immune cells per infected cell per day
Nu.low <- c(8.86e-6, 6.77e-6, 9.48e-6) # proliferation rate of immune cells per infected cell per day, lower bounds of IQ (primary DF, secondary DF, secondary DHF)
Nu.up <- c(6.24e-5, 9.48e-5, 8.07e-5) # proliferation rate of immune cells per infected cell per day, higher bounds of IQ (primary DF, secondary DF, secondary DHF)
Nu.min <- c(7.96e-7, 7.36e-7, 5.01e-7) # proliferation rate of immune cells per infected cell per day, min (primary DF, secondary DF, secondary DHF)
Nu.max <- c(1.09e-3, 3.95e-3, 0.224)  # proliferation rate of immune cells per infected cell per day, max (primary DF, secondary DF, secondary DHF) 
IP <- c(5.80, 5.77, 5.69)             # Incubation Period (not used in initial model)
IP.low <- c(5.48, 5.42, 5.09)         # Incubation Period (not used in initial model), lower bounds of IQ (primary DF, secondary DF, secondary DHF)
IP.up <- c(6.33, 6.02, 6.72)          # Incubation Period (not used in initial model), higher bounds of IQ (primary DF, secondary DF, secondary DHF)
IP.min <- c(4.79, 4.54, 1.72)          # Incubation Period (not used in initial model), min (primary DF, secondary DF, secondary DHF)
IP.max <- c(6.76, 8.29, 9.32)         # Incubation Period (not used in initial model), max (primary DF, secondary DF, secondary DHF)
  
# initial conditions
X1 <- A/gamma                         # uninfected cells
Y1 <- 0                               # infected cells
V1 <- 1                               # virus per ml
Z1 <- c(0.347, 0.411, 0.380)          # immune effector population per ml
Z1.low <- c(0.326, 0.338, 0.326)      # immune effector population per ml, lower bounds of IQ (primary DF, secondary DF, secondary DHF)
Z1.up <- c(0.402, 0.466, 0.439)       # immune effector population per ml, higher bounds of IQ (primary DF, secondary DF, secondary DHF) 
Z1.min <- c(0.285, 0.248, 0.0298)     # immune effector population per ml, min (primary DF, secondary DF, secondary DHF) 
Z1.max <- c(0.848, 0.594, 0.526)      # immune effector population per ml, max (primary DF, secondary DF, secondary DHF) 
  
# Duong viremia descriptives data 
Type <- c('sym', 'asym','presym', 'sym.compiled'); Mean <- c(6.12, 4.75, 6.74, 6.27); SE <- c(0.17, 0.39, 0.25, 0.14); Sample <- c(126, 13, 42, 126 + 42) # sym.compiled is weighted mean of presym and sym
Duong.viremia <- data.frame(Type, Mean, SE, Sample)
Duong.viremia$SD <- Duong.viremia$SE * sqrt(Duong.viremia$Sample)
rm(Mean,SE,Sample)

# Duong viremia vs infectiousness model parameters
Int <- c(-6.5753024, -11.461235, -5.5863807 ); Int.se <- c(0.3850144, 3.4450358, 0.7561838)
Coefs <- c(0.91234842, 2.15818436, 0.98330334); Coefs.se <- c(0.0539375, 0.6420758, 0.1187793)
Duong.infprob <- data.frame(Type[1:3], Int, Int.se, Coefs, Coefs.se)
Duong.infprob$Int.SD <- Duong.infprob$Int.se * sqrt(Duong.viremia$Sample[1:3])
Duong.infprob$Coefs.SD <- Duong.infprob$Coefs.se * sqrt(Duong.viremia$Sample[1:3])
rm(Int,Int.se,Coefs,Coefs.se)

# Duong viremia vs infectiousness model parameters - Indirect feeding
Int <- c(-6.3261364, -6.5198138, -5.1595991 ); Int.se <- c(0.3960377, 1.9240025, 0.7581733)
Coefs <- c(0.8135159, 1.14746273, 0.86484999); Coefs.se <- c(0.054424, 0.3511563, 0.1175872)
Duong.infprob.indirect <- data.frame(Type[1:3], Int, Int.se, Coefs, Coefs.se)
Duong.infprob.indirect$Int.SD <- Duong.infprob$Int.se * sqrt(Duong.viremia$Sample[1:3])
Duong.infprob.indirect$Coefs.SD <- Duong.infprob$Coefs.se * sqrt(Duong.viremia$Sample[1:3])
rm(Int,Int.se,Coefs,Coefs.se)

# Duong viremia vs infectiousness model parameters male
Int <- c(-6.0403259, -11.461235, -11.830976); Int.se <- c(0.5107691, 3.4450358, 2.4400995)  # NOTE: relationship  is assumed to be the same between gender for asymptomatics, because of lack of female infections in Duong
Coefs <- c(0.87463296, 2.15818436, 2.12621479); Coefs.se <- c(0.0710145, 0.6420758, 0.4321683)
Duong.infprob.male <- data.frame(Type[1:3], Int, Int.se, Coefs, Coefs.se)
Duong.infprob.male$Int.SD <- Duong.infprob.male$Int.se * sqrt(Duong.viremia$Sample[1:3])  # WATCH OUT!! THIS IS NOT ACTUALLY THE NUMBER OF MALES, BUT TOTAL. SD IS NOT USED HOWEVER
Duong.infprob.male$Coefs.SD <- Duong.infprob.male$Coefs.se * sqrt(Duong.viremia$Sample[1:3])
rm(Int,Int.se,Coefs,Coefs.se)

# Duong viremia vs infectiousness model parameters female
Int <- c(-6.8040019, -11.461235, -4.1083211); Int.se <- c(0.5914096, 3.4450358, 0.9043418)
Coefs <- c(0.89077643, 2.15818436, 0.72526898); Coefs.se <- c(0.0838419, 0.6420758, 0.1340457)
Duong.infprob.female <- data.frame(Type[1:3], Int, Int.se, Coefs, Coefs.se)
Duong.infprob.female$Int.SD <- Duong.infprob.female$Int.se * sqrt(Duong.viremia$Sample[1:3])
Duong.infprob.female$Coefs.SD <- Duong.infprob.female$Coefs.se * sqrt(Duong.viremia$Sample[1:3])
rm(Int,Int.se,Coefs,Coefs.se)

# IIP parameters from Chan and Johansson
tau.mean <- 13.7; tau.low <- 10.9; tau.up <- 16.9;      # 95% credibility intervals
beta.mean <- 0.56; beta.low <- 0.51; beta.up <- 0.60;   # 95% credibility intervals
Chan.IIP <- data.frame(tau.mean, tau.low, tau.up, beta.mean, beta.low, beta.up) 
rm(tau.mean, tau.low, tau.up, beta.mean, beta.low, beta.up)

# weibull
tau.mean <- 4.0; tau.low <- 3.5; tau.up <- 4.4;      # 95% credibility intervals
beta.mean <- -7.1; beta.low <- -8.0; beta.up <- -6.2;   # 95% credibility intervals
Chan.IIP.weib <- data.frame(tau.mean, tau.low, tau.up, beta.mean, beta.low, beta.up) 
rm(tau.mean, tau.low, tau.up, beta.mean, beta.low, beta.up)

# gamma
tau.mean <- 16; tau.low <- 13; tau.up <- 20;      # 95% credibility intervals
beta.mean <- 1.78; beta.low <- 1.64; beta.up <- 1.92;   # 95% credibility intervals
Chan.IIP.gamma <- data.frame(tau.mean, tau.low, tau.up, beta.mean, beta.low, beta.up) 
rm(tau.mean, tau.low, tau.up, beta.mean, beta.low, beta.up)

# exponential
tau.mean <- NaN; tau.low <- NaN; tau.up <- NaN;      # 95% credibility intervals
beta.mean <- 1.8; beta.low <- 1.6; beta.up <- 2.0;   # 95% credibility intervals
Chan.IIP.exp <- data.frame(tau.mean, tau.low, tau.up, beta.mean, beta.low, beta.up) 
rm(tau.mean, tau.low, tau.up, beta.mean, beta.low, beta.up)


# Nishiura infectiousness over time estimates (Nishiura 2007)
Nish.time <- c(-2, -1, 0, 1, 2, 3)
Nish.time.frominf <- Nish.time + 5.7 - 0.5 # first is the mean IIP, the second is to adjust for censoring
Nish.Inf <- c(0.248691183213025, 0.803568191753655, 1.00642413127619, 0.967604785311397, 0.499712113619158, -0.00321206563809486)

# FOI part of the story
inf.period = 0.01
inv.cross.im = 1/2                           # 2 years cross-immune, as per Reich et al. 
foi.vec = seq(.001, .1, by = .001)
ga = 1 / 4
# rho = 5
Xs = seq(1e-5,5e-2,1e-4)

inapp.ratio.rev = c(1-0.24, 1-0.18, 0.86)                   # @@ HARD CODED meta-analysis results!!!! (based on Clapham et al 2017)
Prim.dis = c(0.092,0.901,0.007) #c(0.74,0.257,0.003)        # Meta-analysis + Montoya for DHF proportions (0.8 % of sympto)     
Sec.dis = c(0.092,0.88,0.028)  #c(0.30,0.682,0.018)                                                      (3 % of sympto)
Tert.dis = c(0.092,0.88,0.028) #c(0.93,0.068,0.002)                                                     (3 % of sympto)
detection.rate = 0.08           # as per Stanaway et al. 2016
detection.rate.CI = c(0.05,0.15)  # as per Stanaway et al. 2016
prop.DHF = c(0.008,0.03,0.03) # as per Montoya, prim, sec, post-sec
# NOTE: this assumes asympto is the same in post-sec, just inapparent is increased. (only part we have data on)    

Prim.dis.low = c(0.05,0.9405,0.0095) #c(0.74,0.257,0.003)        # Meta-analysis + Montoya for DHF proportions    
Sec.dis.low = c(0.05,0.92625,0.02375)  #c(0.30,0.682,0.018)

Prim.dis.high = c(0.20,0.792,0.008) #c(0.74,0.257,0.003)        # Meta-analysis + Montoya for DHF proportions    
Sec.dis.high = c(0.20,0.78,0.02)  #c(0.30,0.682,0.018)


