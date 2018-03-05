# # Code to run meta-analysis on the proportion of apparent 
# in post-secondary infections 

rm(list=ls())
# libraries ---------------------------------------------------------------
library(tictoc)
library(binom)
library(bbmle)
library(emdbook)
library(BayesianTools)

# source files ------------------------------------------------------------
source('Functions_Meta.R')

# set parameters ----------------------------------------------------------
nbParam<-1
names = c('inapp_prop')
param<-c(0.5)
param.bb <- c(param, 1)
sdProposal<-c(.03)
sdProposal.bb <- c(0.4, 1)
nbIteration<-6e3
burnIn<-1e3

# define data -------------------------------------------------------------
asymptomatic_dat_cohorts<- list(s_iq = 2,          
                                apppostsec_iq=c(7,72),
                                inapppostsec_iq = c(468,1035),
                                Npostsec_iq=c(7+468, 72+1035),      
                                
                                y_nic=1,                  
                                app_postsec_nic=c(58),
                                inapp_postsec_nic=c(116),
                                Npostsec_nic=c(58+116) )
k_inapp = unlist(c(asymptomatic_dat_cohorts[c('inapppostsec_iq','inapp_postsec_nic')]))
N_total = unlist(c(asymptomatic_dat_cohorts[c('Npostsec_iq','Npostsec_nic')]))

# run MCMC ----------------------------------------------------------------
tic(); out.full <- run.mcmc(initial = param, nbIteration = nbIteration,
                            sdProposal = sdProposal, k = k_inapp, N = N_total, 
                            verbose = 0 ); 
acceptance = sum(out.full[[3]][-c(1:burnIn)]) / nbIteration; acceptance; toc()
median(out.full[[2]][,1])
out.full[[4]]

tic(); out.full.bb <- run.mcmc(initial = param.bb, nbIteration = nbIteration, 
                            sdProposal = sdProposal.bb, k = k_inapp, N = N_total, 
                            verbose = 0 ); 
acceptance = colSums(out.full.bb[[3]][-c(1:burnIn),]) / nbIteration; acceptance; toc()
median(out.full.bb[[2]][,1])
HDIofMCMC(out.full.bb[[2]][,1])
correlationPlot(out.full.bb[[2]][])
out.full.bb[[4]]
# Compile results ------------------------------------------------------------
