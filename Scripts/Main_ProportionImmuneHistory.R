# Script to estimate the proportion of the population immune to x-number of serotypes

# load packages ----------------------------------------------------------
library(deSolve)

# Load source files -------------------------------------------------------
source('Functions_ProportionImmuneHistory.R')
source('Parameters.R')

# Set parameters and output matrices --------------------------------------
Filenames <- c('pop_age_brazil.csv', 'pop_age_thailand.csv')

# Load data ---------------------------------------------------------------
path <- getwd()
setwd('../data')
agepop.brazil = read.csv(Filenames[1],col.names=c('age','pop'), header = FALSE)
agepop.thailand = read.csv(Filenames[2],col.names=c('age','pop'), header = FALSE )
setwd(path)

# Estimate proportions immune given force of infection and total prevalence--------------------
infec.hist.brazil <- infec.hist.calculator(foi.vec,agepop.brazil,inv.cross.im = inv.cross.im,  ode = 'yes')
infec.hist.thailand <- infec.hist.calculator(foi.vec,agepop.thailand,inv.cross.im = inv.cross.im, ode = 'yes')
infec.hist.brazil.emerging <- infec.hist.calculator(foi.vec,agepop.brazil,endemic = 'no',inv.cross.im = inv.cross.im, ode = 'yes')
infec.hist.thailand.emerging <- infec.hist.calculator(foi.vec,agepop.thailand,endemic = 'no',inv.cross.im = inv.cross.im, ode = 'yes')




