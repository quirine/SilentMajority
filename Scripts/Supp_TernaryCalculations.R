
# load library ------------------------------------------------------------
library(vcd)
library(ggplot2)
library(ggtern)
library(proto)
library(deSolve)
library(data.table)
library(abind)
library(coda)
library(foreach)
library(doParallel)
library(grid)
library(Rmisc)

# Load data ---------------------------------------------------------------
load('Workspace.lowerAsym.RData')

# Load source files -------------------------------------------------------
source('Functions_FOICalculations.R')
source('Parameters.R')
source('Functions_ProportionImmuneHistory.R')

# set parameters ----------------------------------------------------------
SA.foi.vec = c(0.001, 0.05, 0.1)
stepsize = 10

# Create vectors for proportions in ternary plot --------------------------
Inap <- seq(0,100,by=stepsize)
Ap <- As <- numeric()

for ( ii in Inap ) {
  temp <- seq(0, 100 - ii, by = stepsize)
  for ( jj in temp ) {  
    As <- c(As, jj)
    Ap <- c(Ap, 100 - ii - jj)  
  }
}

count = length(Inap)
Ia <- numeric()
for (ii in 1 : length(Inap) ) {
  Ia <- c(Ia, rep(Inap[ii],count))
  count <- count - 1
}

df <- data.frame(As,Ia,Ap)

# Estimate proportions immune given force of infection and total prevalence--------------------
infec.hist.brazil.small <- infec.hist.calculator(SA.foi.vec,agepop.brazil, inv.cross.im = inv.cross.im, ode = 'yes')

# Idem using full distribution --------------------------------------------
PROPS.brazil.low.all <- PROPS.brazil.mid.all <- PROPS.brazil.high.all <- numeric()

for (ii in 1:length(Ia)){
  X <- Disease.prop.calculator(infec.hist.brazil.small,c(df$As[ii],sum(df$Ia[ii],df$Ap[ii]),0), c(df$As[ii],sum(df$Ia[ii],df$Ap[ii]),0), Tert.dis) 
  PROPS.brazil.low.all<- rbind(PROPS.brazil.low.all, Sympto.FOI.sampler(prev = infec.hist.brazil$prev.vec2,X = X[1,], 2, SA.foi.vec[1], Ca.all, Cm.all, Cs.all, cores = 3))
}

# derive contributions and residuals --------------------------------------
PROP.inap.all <- PROP.ap.all <- PROP.as.all <- PROP.inap.all.up <- PROP.ap.all.up <- PROP.as.all.up <-PROP.inap.all.low <- PROP.ap.all.low <- PROP.as.all.low <- numeric()
Res.as <-  Res.inap <-  Res.ap <-  Res.inap.as <-numeric()
Res.as.low <-  Res.inap.low <-  Res.ap.low <-  Res.inap.as.low <-numeric()
Res.as.up <-  Res.inap.up <-  Res.ap.up <-  Res.inap.as.up <-numeric()
for (ii in 1:length(Ia)) {
  # browser()
  ia.prop <- df$Ia[ii]/100
  as.prop <- df$As[ii]/100
  ap.prop <- df$Ap[ii]/100
  rel.ia.vs.ap <- ia.prop/ (ia.prop + ap.prop)
  if (is.nan(rel.ia.vs.ap)) { rel.ia.vs.ap = 0; } 
  prop.inap.all <- PROPS.brazil.low.all[ii,] * rel.ia.vs.ap 
  PROP.inap.all <- c(PROP.inap.all, mean(prop.inap.all) )
  prop.ap.all <- PROPS.brazil.low.all[ii,] * (1 - rel.ia.vs.ap)
  PROP.ap.all <- c(PROP.ap.all, mean(prop.ap.all ) )
  prop.as.all <- 1 - PROPS.brazil.low.all[ii,]
  PROP.as.all <- c(PROP.as.all, mean(prop.as.all))
  
  Res.as <- c(Res.as, as.prop - PROP.as.all[ii]) 
  Res.inap <- c(Res.inap, ia.prop - PROP.inap.all[ii]) 
  Res.ap <- c(Res.ap, ap.prop - PROP.ap.all[ii])
  Res.inap.as <- c(Res.inap.as, (as.prop+ia.prop) - (PROP.inap.all[ii] + PROP.as.all[ii] ))
}

df <- data.frame(Ia,As,Ap,PROP.as.all, PROP.inap.all, PROP.ap.all, Res.as, Res.inap, Res.ap, Res.inap.as)
colnames(df) <- c('Ia','As','Ap','PROP.as.all', 'PROP.inap.all', 'PROP.ap.all', 'Res.as', 'Res.inap', 'Res.ap', 'Res.inap.as')


# save --------------------------------------------------------------------

save(df,file='TempTernary.RData')

