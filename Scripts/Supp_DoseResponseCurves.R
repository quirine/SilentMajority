# Script to compare dose-response curves of mosquito infections between Duong et al, Ferguson et al, and Nguyen et al
rm(list=ls())

# load libraries ----------------------------------------------------------
library(bbmle)
library(emdbook)
library(ggplot2)
library(lme4)
library(mgcv)
library(HRQoL)

# load source files -------------------------------------------------------
source('Functions_DoseResponseCurves.R')

# load data files ---------------------------------------------------------
prep_data_Ferg = read.table('Ferguson.csv', sep = ',', dec = '.', header = TRUE)
prep_data_Nguyen = read.table('Nguyen.csv', sep = ',', dec = '.', header = TRUE)
prep_data_Duong = read.table('Duong_symp.csv', sep = ',', dec = '.', header = TRUE)
prep_data_Duong_alldisease = read.table('Duong.csv', sep = ',', dec = '.', header = TRUE)

# set parameters ----------------------------------------------------------
tolerance = 1e-4
Duong.Int <- -6.5753024
Duong.Coef <- 0.91234842
Ferg.inf.dose.1 <- 5.90
Ferg.inf.dose.2 <- 6.78
Ferg.inf.dose.3 <- 8.41
Ferg.inf.dose.4 <- 9.50
Ferg.slope <- 2.88
Ferg.overdisp <- 0.46

# prepare data (Ferguson)------------------------------------------------------------
log_V_Ferg = prep_data_Ferg$log10_virus
N_Ferg = prep_data_Ferg$N_Ferg
k_Ferg = prep_data_Ferg$k_Ferg
serotypes_Ferg = as.numeric(prep_data_Ferg$Serotype)

# prepare data (Duong) ----------------------------------------------------
# just symptomatic
log_V_Duong = prep_data_Duong$log10_viremia
N_Duong = prep_data_Duong$N_Duong
k_Duong = prep_data_Duong$k_Duong
serotypes_Duong = as.numeric(prep_data_Duong$serotype)
# all disease classes
N_Duong_alldisease = prep_data_Duong_alldisease$N_Duong
k_Duong_alldisease = prep_data_Duong_alldisease$k_Duong
log_V_Duong_alldisease = prep_data_Duong_alldisease$log10_viremia
serotypes_Duong_alldisease = prep_data_Duong_alldisease$serotype
disease_class_Duong_alldisease = prep_data_Duong_alldisease$disease_class

# prepare data (Nguyen) ---------------------------------------------------
log_V_Nguyen = prep_data_Nguyen$Viremia
N_Nguyen = prep_data_Nguyen$N_Nguyen
k_Nguyen = prep_data_Nguyen$k_Nguyen
serotypes_Nguyen = as.numeric(prep_data_Nguyen$Serotype)

# Fit to Ferguson model (binomial)---------------------------------------------------
binom_Ferg <- mle2(minuslogl = FergusonNLL, 
                   start = c(inf.dose.1 = 7, slope = 3),
           data = list(N = N_Ferg, k = k_Ferg, log_V = log_V_Ferg))

dif = 1
while(dif>tolerance) {
  old = binom_Ferg@min
  binom_Ferg <- mle2(minuslogl = FergusonNLL, 
                     start = c(inf.dose.1 = as.numeric(binom_Ferg@coef[1]), 
                               slope = as.numeric(binom_Ferg@coef[2])),
                     data = list(N = N_Ferg, k = k_Ferg, log_V = log_V_Ferg))
  new = old = binom_Ferg@min
  dif = abs(old - new)
}

binom_Ferg_ST <- mle2(minuslogl = FergusonNLL, 
                   start = c(inf.dose.1 = 7,inf.dose.2 = 7,inf.dose.3 = 7,inf.dose.4 = 7,
                             slope = 3),
                   data = list(N = N_Ferg, k = k_Ferg, log_V = log_V_Ferg, Serotypes = serotypes_Ferg))
anova(binom_Ferg, binom_Ferg_ST)

binom_Ferg_Duong <- mle2(minuslogl = FergusonNLL, 
                   start = c(inf.dose.1 = 7, slope = 3),
                   data = list(N = N_Duong, k = k_Duong, log_V = log_V_Duong))

dif = 1
while(dif>tolerance) {
  old = binom_Ferg_Duong@min
  binom_Ferg_Duong <- mle2(minuslogl = FergusonNLL, 
                     start = c(inf.dose.1 = as.numeric(binom_Ferg_Duong@coef[1]), 
                               slope = as.numeric(binom_Ferg_Duong@coef[2])),
                     data = list(N = N_Duong, k = k_Duong, log_V = log_V_Duong))
  new = old = binom_Ferg_Duong@min
  dif = abs(old - new)
}

binom_Ferg_Duong_ST <- mle2(minuslogl = FergusonNLL, 
                         start = c(inf.dose.1 = 7,inf.dose.2 = 7,inf.dose.3 = 7,inf.dose.4 = 7,
                                   slope = 3),
                         data = list(N = N_Duong, k = k_Duong, log_V = log_V_Duong, Serotypes = serotypes_Duong))

anova(binom_Ferg_Duong, binom_Ferg_Duong_ST)

binom_Ferg_Nguyen <- mle2(minuslogl = FergusonNLL, 
                   start = c(inf.dose.1 = 7, slope = 3),
                   data = list(N = N_Nguyen, k = k_Nguyen, log_V = log_V_Nguyen))

binom_Ferg_Nguyen_ST <- mle2(minuslogl = FergusonNLL, 
                      start = c(inf.dose.1 = 7,inf.dose.2 = 7,inf.dose.3 = 7,inf.dose.4 = 7,
                                slope = 3),
                      data = list(N = N_Nguyen, k = k_Nguyen, log_V = log_V_Nguyen, Serotypes = serotypes_Nguyen))
anova(binom_Ferg_Nguyen, binom_Ferg_Nguyen_ST)

binom_Ferg_all <- mle2(minuslogl = FergusonNLL, 
                          start = c(inf.dose.1 = 7, slope = 3),
                          data = list(N = c(N_Ferg,N_Duong,N_Nguyen), k = c(k_Ferg,k_Duong,k_Nguyen), 
                                      log_V = c(log_V_Ferg,log_V_Duong,log_V_Nguyen)))

binom_Ferg_all_ST <- mle2(minuslogl = FergusonNLL, 
                             start = c(inf.dose.1 = 7,inf.dose.2 = 7,inf.dose.3 = 7,inf.dose.4 = 7,
                                       slope = 3),
                          data = list(N = c(N_Ferg,N_Duong,N_Nguyen), k = c(k_Ferg,k_Duong,k_Nguyen), 
                                      log_V = c(log_V_Ferg,log_V_Duong,log_V_Nguyen), Serotypes = c(serotypes_Ferg, serotypes_Duong, serotypes_Nguyen)))
anova(binom_Ferg_all, binom_Ferg_all_ST)

# Fit to Ferguson model (beta-binomial)---------------------------------------------------
betabinom_Ferg <- mle2(minuslogl = FergusonBetaBinomNLL, 
                       start = c(inf.dose = 9, slope = 2.88, overdispersion = 0.46),
                   data = list(N = N_Ferg, k = k_Ferg, log_V = log_V_Ferg))

anova(binom_Ferg, betabinom_Ferg)

betabinom_Ferg_Duong <- mle2(minuslogl = FergusonBetaBinomNLL, 
                       start = c(inf.dose = 9, slope = 2.88, overdispersion = 0.46),
                       data = list(N = N_Duong, k = k_Duong, log_V = log_V_Duong))

betabinom_Ferg_Nguyen <- mle2(minuslogl = FergusonBetaBinomNLL, 
                       start = c(inf.dose = 9, slope = 2.88, overdispersion = 0.46),
                       data = list(N = N_Nguyen, k = k_Nguyen, log_V = log_V_Nguyen))

# Fit to Hill function (binomial) ----------------------------------------------------
binom_Hill <- mle2(minuslogl = HillNLL, 
                   start = c(dis.const = 0.01, hill.coef = 0.01),
                   data = list(N = N_Ferg, k = k_Ferg, log_V = log_V_Ferg))

dif = 1
while(dif>tolerance) {
  old = binom_Hill@min
  binom_Hill <- mle2(minuslogl = HillNLL, 
                     start = c(dis.const = as.numeric(binom_Hill@coef[1]), 
                               hill.coef = as.numeric(binom_Hill@coef[2])),
                     data = list(N = N_Ferg, k = k_Ferg, log_V = log_V_Ferg))
  new = old = binom_Hill@min
  dif = abs(old - new)
}

binom_Hill_Duong <- mle2(minuslogl = HillNLL, 
                   start = c(dis.const = 0.01, hill.coef = 0.01),
                   data = list(N = N_Duong, k = k_Duong, log_V = log_V_Duong))

dif = 1
while(dif>tolerance) {
  old = binom_Hill_Duong@min
  binom_Hil_Duong <- mle2(minuslogl = HillNLL, 
                     start = c(dis.const = as.numeric(binom_Hill_Duong@coef[1]), 
                               hill.coef = as.numeric(binom_Hill_Duong@coef[2])),
                     data = list(N = N_Duong, k = k_Duong, log_V = log_V_Duong))
  new = old = binom_Hill_Duong@min
  dif = abs(old - new)
}

binom_Hill_Nguyen <- mle2(minuslogl = HillNLL, 
                         start = c(dis.const = 0.01, hill.coef = 0.01),
                         data = list(N = N_Nguyen, k = k_Nguyen, log_V = log_V_Nguyen))

# Fit to Hill function (beta - binomial) ----------------------------------------------------
betabinom_Hill <- mle2(minuslogl = HillBetaBinomNLL, 
                  start = c(dis.const = 50, hill.coef = 0.1, overdispersion = 0.1),
                  data = list(N = N_Ferg, k = k_Ferg, log_V = log_V_Ferg))

betabinom_Hill_Duong <- mle2(minuslogl = HillBetaBinomNLL, 
                       start = c(dis.const = 50, hill.coef = 0.1, overdispersion = 0.1),
                       data = list(N = N_Duong, k = k_Duong, log_V = log_V_Duong))

betabinom_Hill_Nguyen <- mle2(minuslogl = HillBetaBinomNLL, 
                       start = c(dis.const = 50, hill.coef = 0.1, overdispersion = 0.1),
                       data = list(N = N_Nguyen, k = k_Nguyen, log_V = log_V_Nguyen))


# Fit to logistic function ------------------------------------------------
binom_Log <- mle2(minuslogl = LogNLL, 
                   start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01),
                   data = list(N = N_Ferg, k = k_Ferg, log_V = log_V_Ferg))

dif = 1
while(dif>tolerance) {
  old = binom_Log@min
  binom_Log <- mle2(minuslogl = LogNLL, 
                     start = c(tweek.parameter1 = as.numeric(binom_Log@coef[1]), 
                               tweek.parameter1 = as.numeric(binom_Log@coef[2])),
                     data = list(N = N_Ferg, k = k_Ferg, log_V = log_V_Ferg))
  new = old = binom_Log@min
  dif = abs(old - new)
}

binom_Log_Duong <- mle2(minuslogl = LogNLL, 
                  start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01),
                  data = list(N = N_Duong, k = k_Duong, log_V = log_V_Duong))

dif = 1
while(dif>tolerance) {
  old = binom_Log_Duong@min
  binom_Log_Duong <- mle2(minuslogl = LogNLL, 
                    start = c(tweek.parameter1 = as.numeric(binom_Log_Duong@coef[1]), 
                              tweek.parameter1 = as.numeric(binom_Log_Duong@coef[2])),
                    data = list(N = N_Duong, k = k_Duong, log_V = log_V_Duong))
  new = old = binom_Log_Duong@min
  dif = abs(old - new)
}

binom_Log_Nguyen <- mle2(minuslogl = LogNLL, 
                        start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01),
                        data = list(N = N_Nguyen, k = k_Nguyen, log_V = log_V_Nguyen))

binom_Log_all <- mle2(minuslogl = LogNLL, 
                         start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01),
                         data = list(N = c(N_Ferg,N_Duong,N_Nguyen), k = c(k_Ferg,k_Duong,k_Nguyen), log_V = c(log_V_Ferg,log_V_Duong,log_V_Nguyen)))

binom_Log_DF <- mle2(minuslogl = LogNLL, 
                      start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01),
                      data = list(N = c(N_Ferg,N_Duong), k = c(k_Ferg,k_Duong), log_V = c(log_V_Ferg,log_V_Duong)))

binom_Log_DN <- mle2(minuslogl = LogNLL, 
                      start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01),
                      data = list(N = c(N_Duong,N_Nguyen), k = c(k_Duong,k_Nguyen), log_V = c(log_V_Duong,log_V_Nguyen)))

binom_Log_FN <- mle2(minuslogl = LogNLL, 
                      start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01),
                      data = list(N = c(N_Ferg,N_Nguyen), k = c(k_Ferg,k_Nguyen), log_V = c(log_V_Ferg,log_V_Nguyen)))

# and by serotype
binom_Log_ST <- mle2(minuslogl = LogNLL.ST, 
                  start = c(coef.1 = 1, coef.2 = 1, coef.3 = 1, coef.4 = 1, intercept = -7),
                  data = list(N = N_Ferg, k = k_Ferg, log_V = log_V_Ferg, Serotypes = serotypes_Ferg))

binom_Log_Duong_ST <- mle2(minuslogl = LogNLL.ST, 
                     start = c(coef.1 = 1, coef.2 = 1, coef.3 = 1, coef.4 = 1,intercept = -7),
                     data = list(N = N_Duong, k = k_Duong, log_V = log_V_Duong, Serotypes = serotypes_Duong))

binom_Log_Nguyen_ST <- mle2(minuslogl = LogNLL.ST, 
                     start = c( coef.1 = 1, coef.2 = 1, coef.3 = 1, coef.4 = 1, intercept = -7),
                     data = list(N = N_Nguyen, k = k_Nguyen, log_V = log_V_Nguyen, Serotypes = serotypes_Nguyen))

binom_Log_all_ST <- mle2(minuslogl = LogNLL.ST, 
                            start = c( coef.1 = 1, coef.2 = 1, coef.3 = 1, coef.4 = 1, intercept = -7),
                            data = list(N = c(N_Ferg,N_Duong,N_Nguyen), k = c(k_Ferg,k_Duong,k_Nguyen), 
                                        log_V = c(log_V_Ferg,log_V_Duong,log_V_Nguyen), 
                                        Serotypes = c(serotypes_Ferg, serotypes_Duong,serotypes_Nguyen) ) )


# idem for data stratified by disease class. Not ideal because of the small sample size for asymp
ind = which(prep_data_Duong_alldisease$disease_class=='symp')
binom_Log_Duong.symp <- mle2(minuslogl = LogNLL, 
                                 start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01),
                                 data = list(N = prep_data_Duong_alldisease$N_Duong[ind], k = prep_data_Duong_alldisease$k_Duong[ind], log_V = prep_data_Duong_alldisease$log10_viremia[ind]))

ind = which(prep_data_Duong_alldisease$disease_class=='presymp')
binom_Log_Duong.presymp <- mle2(minuslogl = LogNLL, 
                                    start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01),
                                    data = list(N = prep_data_Duong_alldisease$N_Duong[ind], k = prep_data_Duong_alldisease$k_Duong[ind], log_V = prep_data_Duong_alldisease$log10_viremia[ind]))

ind = which(prep_data_Duong_alldisease$disease_class=='asymp')
binom_Log_Duong.asymp <- mle2(minuslogl = LogNLL, 
                                  start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01),
                                  data = list(N = prep_data_Duong_alldisease$N_Duong[ind], k = prep_data_Duong_alldisease$k_Duong[ind], log_V = prep_data_Duong_alldisease$log10_viremia[ind]))


# Fit to logistic function (beta - binomial) ------------------------------------------------
betabinom_Log <- mle2(minuslogl = LogBetaBinomNLL, 
                  start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01, overdispersion = 0.1),
                  data = list(N = N_Ferg, k = k_Ferg, log_V = log_V_Ferg))


betabinom_Log_Duong <- mle2(minuslogl = LogBetaBinomNLL, 
                        start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01, overdispersion = 0.1),
                        data = list(N = N_Duong, k = k_Duong, log_V = log_V_Duong))


betabinom_Log_Nguyen <- mle2(minuslogl = LogBetaBinomNLL, 
                         start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01, overdispersion = 0.1),
                         data = list(N = N_Nguyen, k = k_Nguyen, log_V = log_V_Nguyen))


ind = which(prep_data_Duong_alldisease$disease_class=='symp')
betabinom_Log_Duong.symp <- mle2(minuslogl = LogBetaBinomNLL, 
                            start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01, overdispersion = 0.1),
                            data = list(N = prep_data_Duong_alldisease$N_Duong[ind], k = prep_data_Duong_alldisease$k_Duong[ind], log_V = prep_data_Duong_alldisease$log10_viremia[ind]))

ind = which(prep_data_Duong_alldisease$disease_class=='presymp')
betabinom_Log_Duong.presymp <- mle2(minuslogl = LogBetaBinomNLL, 
                                    start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01, overdispersion = 0.1),
                                    data = list(N = prep_data_Duong_alldisease$N_Duong[ind], k = prep_data_Duong_alldisease$k_Duong[ind], log_V = prep_data_Duong_alldisease$log10_viremia[ind]))

ind = which(prep_data_Duong_alldisease$disease_class=='asymp')
betabinom_Log_Duong.asymp <- mle2(minuslogl = LogBetaBinomNLL, 
                                  start = c(tweek.parameter1 = 0.01, tweek.parameter2 = 0.01, overdispersion = 0.1),
                                  data = list(N = prep_data_Duong_alldisease$N_Duong[ind], k = prep_data_Duong_alldisease$k_Duong[ind], log_V = prep_data_Duong_alldisease$log10_viremia[ind]))

# Likelihood ratio tests --------------------------------------------------
# the log-models are significantly different between datasets: 
lrt(LL0 = binom_Log_all@details$value, 
    LL1 = sum(binom_Log_Duong@details$value,binom_Log@details$value,binom_Log_Nguyen@details$value),
    df0 = 2, df1 = 6) 
# also between Duong and Ferguson:
lrt(LL0 = binom_Log_DF@details$value, 
    LL1 = sum(binom_Log_Duong@details$value,binom_Log@details$value),
    df0 = 2, df1 = 4) 
# and Duong and Nguyen
lrt(LL0 = binom_Log_DN@details$value, 
    LL1 = sum(binom_Log_Duong@details$value,binom_Log_Nguyen@details$value),
    df0 = 2, df1 = 4) 
# and Nguyen and Ferguson
lrt(LL0 = binom_Log_FN@details$value, 
    LL1 = sum(binom_Log@details$value,binom_Log_Nguyen@details$value),
    df0 = 2, df1 = 4) 

lrt(LL0 = binom_Log_all@details$value, 
    LL1 = sum(binom_Log_DN@details$value,binom_Log@details$value),
    df0 = 2, df1 = 4) 

# betabinom vs binom for log model (Note that we're comparing the logistic-regression model however)
lrt(LL0 = binom_Log_Duong@details$value, LL1 = betabinom_Log_Duong@details$value, df0 = 2, df1 = 3)
anova(binom_Log_Duong, betabinom_Log_Duong)

# get AICs for each model -------------------------------------------------
AICs_Ferg = AIC(binom_Ferg,binom_Ferg_ST, binom_Hill, binom_Log,
                betabinom_Ferg, betabinom_Hill, betabinom_Log)
row.names(AICs_Ferg) = c('binom_Ferg','binom_Ferg_ST', 'binom_Hill', 'binom_Log',
                'betabinom_Ferg', 'betabinom_Hill', 'betabinom_Log')
AICs_Duong = AIC(binom_Ferg_Duong,binom_Ferg_Duong_ST, binom_Hill_Duong, binom_Log_Duong,
                 betabinom_Ferg_Duong, betabinom_Hill_Duong, betabinom_Log_Duong)
row.names(AICs_Duong) = c('binom_Ferg_Duong','binom_Ferg_Duong_ST', 'binom_Hill_Duong', 'binom_Log_Duong',
                        'betabinom_Ferg_Duong', 'betabinom_Hill_Duong', 'betabinom_Log_Duong')
AICs_Nguyen = AIC(binom_Ferg_Nguyen,binom_Ferg_Nguyen_ST, binom_Hill_Nguyen, binom_Log_Nguyen,
                  betabinom_Ferg_Nguyen, betabinom_Hill_Nguyen, betabinom_Log_Nguyen)
row.names(AICs_Nguyen) = c('binom_Ferg_Nguyen','binom_Ferg_Nguyen_ST', 'binom_Hill_Nguyen', 'binom_Log_Nguyen',
                  'betabinom_Ferg_Nguyen', 'betabinom_Hill_Nguyen', 'betabinom_Log_Nguyen')


# Get coefficients and CIs ------------------------------------------------
get.binom.model.characteristics(binom_Log_Duong)
get.binom.model.characteristics(binom_Hill_Duong)
get.binom.model.characteristics(binom_Ferg_Duong)

get.binom.model.characteristics(binom_Log)
get.binom.model.characteristics(binom_Hill)
get.binom.model.characteristics(binom_Ferg)

get.binom.model.characteristics(binom_Log_Nguyen)
get.binom.model.characteristics(binom_Hill_Nguyen)
get.binom.model.characteristics(binom_Ferg_Nguyen)

get.binom.model.characteristics(betabinom_Log_Duong, type = 'betabinom')
get.binom.model.characteristics(betabinom_Hill_Duong, type = 'betabinom')
get.binom.model.characteristics(betabinom_Ferg_Duong, type = 'betabinom')

get.binom.model.characteristics(betabinom_Log, type = 'betabinom')
get.binom.model.characteristics(betabinom_Hill, type = 'betabinom')
get.binom.model.characteristics(betabinom_Ferg, type = 'betabinom')

get.binom.model.characteristics(betabinom_Log_Nguyen, type = 'betabinom')
get.binom.model.characteristics(betabinom_Hill_Nguyen, type = 'betabinom')
get.binom.model.characteristics(betabinom_Ferg_Nguyen, type = 'betabinom')

# multivariate logistic regression using probabililties --------------------------------
prep_data_Duong_alldisease <- within(prep_data_Duong_alldisease, disease_class <- relevel(disease_class, ref = 'symp'))
prep_data_Duong_alldisease <- within(prep_data_Duong_alldisease, serotype <- relevel(serotype, ref = 'D2'))

log.reg.mod.Duong.2.ng <- glm(k_Duong / N_Duong ~ log10_viremia + disease_class, weights = N_Duong,
                           family = binomial(link = "logit"), data = prep_data_Duong_alldisease)
summary(log.reg.mod.Duong.2.ng)

# serotype model 
log.reg.mod.Duong.uniST <- glm(k_Duong / N_Duong ~ log10_viremia + serotype, weights = N_Duong,
              family = binomial(link = "logit"), data = prep_data_Duong_alldisease)
summary(log.reg.mod.Duong.uniST) 

log.reg.mod.Duong.ST.DC <- glm(k_Duong / N_Duong ~ log10_viremia + serotype + disease_class, weights = N_Duong,
                               family = binomial(link = "logit"), data = prep_data_Duong_alldisease)
summary(log.reg.mod.Duong.ST.DC) 
step(log.reg.mod.Duong.ST.DC) 

log.reg.mod.Duong.base <- glm(k_Duong / N_Duong ~ log10_viremia, weights = N_Duong,
                               family = binomial(link = "logit"), data = prep_data_Duong_alldisease)
summary(log.reg.mod.Duong.base)
anova(log.reg.mod.Duong.base, log.reg.mod.Duong.uniST)
pchisq(27.434, 2, lower.tail = FALSE)

anova(log.reg.mod.Duong.2.ng, log.reg.mod.Duong.ST.DC)
pchisq(15.681, 2, lower.tail = FALSE)

# by disease class fits
ind = which(prep_data_Duong_alldisease$disease_class == 'symp')
model.symp <- glm(k_Duong / N_Duong ~ log10_viremia, weights = N_Duong,
              family = binomial(link = "logit"), data = prep_data_Duong_alldisease[ind,])
summary(model.symp)
# including the serotype model 
log.reg.mod.Duong.uniST.symp <- glm(k_Duong / N_Duong ~ log10_viremia + serotype, weights = N_Duong,
                  family = binomial(link = "logit"), data = prep_data_Duong_alldisease[ind,])
summary(log.reg.mod.Duong.uniST.symp)
anova(model.symp, log.reg.mod.Duong.uniST.symp)
pchisq(10.456, 2, lower.tail = FALSE)

ind = which(prep_data_Duong_alldisease$disease_class == 'asymp')
model.asymp <- glm(k_Duong / N_Duong ~ log10_viremia , weights = N_Duong,
             family = binomial(link = "logit"), data = prep_data_Duong_alldisease[ind,])
summary(model.asymp)

ind = which(prep_data_Duong_alldisease$disease_class == 'presymp')
model.presymp <- glm(k_Duong / N_Duong ~ log10_viremia, weights = N_Duong,
             family = binomial(link = "logit"), data = prep_data_Duong_alldisease[ind,])
summary(model.presymp)

# Likelihood ratio test 
lrt.dev(log.reg.mod.Duong.2.ng$deviance, sum(model.symp$deviance, model.presymp$deviance, model.asymp$deviance), df0 = 4, df1 = 6)

# by serotype fits
ind = which(prep_data_Duong_alldisease$serotype == 'D1')
model.d1 <- glm(k_Duong / N_Duong ~ log10_viremia, weights = N_Duong,
                  family = binomial(link = "logit"), data = prep_data_Duong_alldisease[ind,])
summary(model.d1)

ind = which(prep_data_Duong_alldisease$serotype == 'D2')
model.d2 <- glm(k_Duong / N_Duong ~ log10_viremia, weights = N_Duong,
                family = binomial(link = "logit"), data = prep_data_Duong_alldisease[ind,])
summary(model.d2)

ind = which(prep_data_Duong_alldisease$serotype == 'D4')
model.d4 <- glm(k_Duong / N_Duong ~ log10_viremia, weights = N_Duong,
                family = binomial(link = "logit"), data = prep_data_Duong_alldisease[ind,])
summary(model.d4)

# Likelihood ratio test 
lrt.dev(log.reg.mod.Duong.uniST$deviance, sum(model.d1$deviance, model.d2$deviance, model.d4$deviance), df0 = 4, df1 = 6)

# logistic regression with beta-binomial distribution ---------------------
model <- BBreg(k_Duong~log_V_Duong + factor(serotypes_Duong),N_Duong)
summary(model)

model <- BBreg(k_Duong~log10_viremia + factor(disease_class) + factor(serotype), prep_data_Duong_alldisease$N_Duong, data = prep_data_Duong_alldisease)
summary(model)

log.reg.mod.Duong.bb <- BBreg(k_Duong~log10_viremia + factor(disease_class), prep_data_Duong_alldisease$N_Duong, data = prep_data_Duong_alldisease)
summary(log.reg.mod.Duong.bb)

AIC.Duong.bb = AIC.bb(model = log.reg.mod.Duong.bb, data.set = prep_data_Duong_alldisease, observed.data = k_Duong_alldisease/N_Duong_alldisease )
AIC.Duong = AIC.bb(model = log.reg.mod.Duong.2.ng, data.set = prep_data_Duong_alldisease, observed.data = k_Duong_alldisease/N_Duong_alldisease )

lrt.dev(log.reg.mod.Duong.2.ng$deviance, log.reg.mod.Duong.bb$deviance, df0 = 4, df1 = 5)

# add serotype again, for serotype specific analysis
log.reg.mod.Duong.bb.st <- BBreg(k_Duong~log10_viremia + factor(disease_class) + factor(serotype), prep_data_Duong_alldisease$N_Duong, data = prep_data_Duong_alldisease)
summary(log.reg.mod.Duong.bb.st)

ind = which(prep_data_Duong_alldisease$disease_class == 'symp')
model.bb.symp <- BBreg(k_Duong~log10_viremia, prep_data_Duong_alldisease$N_Duong[ind], data = prep_data_Duong_alldisease[ind,])
summary(model.bb.symp)

ind = which(prep_data_Duong_alldisease$disease_class == 'presymp')
model.bb.presymp <- BBreg(k_Duong~log10_viremia, prep_data_Duong_alldisease$N_Duong[ind], data = prep_data_Duong_alldisease[ind,])
summary(model.bb.presymp)

ind = which(prep_data_Duong_alldisease$disease_class == 'asymp')
model.bb.asymp <- BBreg(k_Duong~log10_viremia, prep_data_Duong_alldisease$N_Duong[ind], data = prep_data_Duong_alldisease[ind,])
summary(model.bb.asymp)

lrt.dev(log.reg.mod.Duong.bb$deviance, sum(model.bb.symp$deviance, model.bb.presymp$deviance, model.bb.asymp$deviance), df0 = 5, df1 = 9)
lrt.dev(model.symp$deviance, model.bb.symp$deviance, df0 = 2, df1 = 3)
lrt.dev(model.asymp$deviance, model.bb.asymp$deviance, df0 = 2, df1 = 3)
lrt.dev(model.presymp$deviance, model.bb.presymp$deviance, df0 = 2, df1 = 3)

# save final model for further analysis -----------------------------------

# save('log.reg.mod.Duong.2.ng', file = 'doseresponsecurves.RData')
save('model.symp','model.presymp','model.asymp', file = 'doseresponsecurves.RData')
# save('log.reg.mod.Duong.bb', file = 'doseresponsecurves.bb.RData')
save('model.bb.symp','model.bb.presymp','model.bb.asymp', file = 'doseresponsecurves.bb.RData')

# plot the Log fits and data (Ferguson, Duong, and Nguyen data) ------------------------------------

pdf('FigS2.1.pdf',width=4.75,height=4.75)
par(mfrow = c(1,1),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,2)) 

LogV_vector = seq(0,12,by = 0.1)
Prob_inf_Log = LogFunction(LogV_vector, binom_Log@coef[1], binom_Log@coef[2] )
Prob_inf_Log_Duong = LogFunction(LogV_vector, binom_Log_Duong@coef[1], binom_Log_Duong@coef[2] )
Prob_inf_Log_Nguyen = LogFunction(LogV_vector, binom_Log_Nguyen@coef[1], binom_Log_Nguyen@coef[2] )
Prob_inf_Log_all = LogFunction(LogV_vector, binom_Log_all@coef[1], binom_Log_all@coef[2] )

plot(LogV_vector, Prob_inf_Log, type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1))
lines(LogV_vector, Prob_inf_Log_Duong, lty = 2, lwd = 2)
lines(LogV_vector, Prob_inf_Log_Nguyen, lty = 3, lwd = 2)
lines(LogV_vector, Prob_inf_Log_all, lty = 4, lwd = 2)

points(log_V_Ferg, k_Ferg/N_Ferg,cex=0.5)
points(log_V_Duong, k_Duong/N_Duong,pch=3,cex=0.5)
points(log_V_Nguyen, k_Nguyen/N_Nguyen,pch=5,cex=0.5)

mtext('Viremia (log10)', side = 1, line = 1, outer=TRUE)
mtext('Probability of human-mosquito transmission', side = 2, line = 1, outer=TRUE)

h=legend(x=-0.4, y=0.95, c('Ferguson fit','Duong fit','Nguyen fit','Combined fit','Ferguson data','Duong data','Nguyen data'),
         col = c('black','black','black','black','black','black'),
         lwd = c(1,1,1), lty = c(1,2,3,4,NA,NA,NA),pch=c(NA,NA,NA,NA,1,3,5), bty = 'n',cex=.8)

dev.off()

# Compare the fits of different functional forms for each data set --------
Prob_inf_Hill_Nguyen = Hill(LogV_vector, binom_Hill_Nguyen@coef[1], binom_Hill_Nguyen@coef[2] )
Prob_inf_Ferg_Nguyen = Ferguson(LogV_vector, binom_Ferg_Nguyen@coef[1], binom_Ferg_Nguyen@coef[2] )
Prob_inf_Ferg_Duong = Ferguson(LogV_vector, binom_Ferg_Duong@coef[1], binom_Ferg_Duong@coef[2] )
Prob_inf_Hill_Duong = Hill(LogV_vector, binom_Hill_Duong@coef[1], binom_Hill_Duong@coef[2] )
Prob_inf_Ferg = Ferguson(LogV_vector, binom_Ferg@coef[1], binom_Ferg@coef[2] )
Prob_inf_Hill = Hill(LogV_vector, binom_Hill@coef[1], binom_Hill@coef[2] )

pdf('FigS2.2.pdf',width=7.5,height=4.00)
par(mfrow = c(1,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,2))

plot(LogV_vector, Prob_inf_Log_Duong, type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1))
lines(LogV_vector, Prob_inf_Hill_Duong, lty = 2, lwd = 2)
lines(LogV_vector, Prob_inf_Ferg_Duong, lty = 3, lwd = 2)
points(log_V_Duong, k_Duong/N_Duong,pch=3,cex=0.5)
mtext(text='Duong et al', side = 3, line = 0.1)

LogV_vector = seq(0,12,by = 0.1)
plot(LogV_vector, Prob_inf_Log, type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1))
lines(LogV_vector, Prob_inf_Hill, lty = 2, lwd = 2)
lines(LogV_vector, Prob_inf_Ferg, lty = 3, lwd = 2)
points(log_V_Ferg, k_Ferg/N_Ferg,cex=0.5)
mtext(text='Ferguson et al', side = 3, line = 0.1)

plot(LogV_vector, Prob_inf_Log_Nguyen, type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1))
lines(LogV_vector, Prob_inf_Hill_Nguyen, lty = 2, lwd = 2)
lines(LogV_vector, Prob_inf_Ferg_Nguyen, lty = 3, lwd = 2)
points(log_V_Nguyen, k_Nguyen/N_Nguyen,pch=5,cex=0.5)
mtext(text='Nguyen et al', side = 3, line = 0.1)

h=legend(x=-0.4, y=0.9, c('Logistic fit','Hill fit','Ferguson fit'),
         col = c('black','black','black'),
         lwd = c(2,2,2), lty = c(1,2,3), bty = 'n',cex=.8)

mtext('Viremia (log10)', side = 1, line = 1, outer=TRUE)
mtext('Probability of human-mosquito transmission', side = 2, line = 1, outer=TRUE)


dev.off()

# plot the logistic fits by Study by Serotype ------------------------------------
cols = c('red', 'blue', 'green', 'cyan')

pdf('FigS2.3.pdf',width=5.75,height=4.75)
par(mfrow = c(2,2),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,2))  

LogV_vector = seq(0,12,by = 0.1)
samples = 1e3

plot(LogV_vector, LogFunction(LogV_vector, binom_Log_Duong_ST@coef[5], binom_Log_Duong_ST@coef[1])
     ,type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1), col = cols[1])
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_Duong_ST@coef[5], binom_Log_Duong_ST@coef[2] ), lty = 1, lwd = 2, col = cols[2])
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_Duong_ST@coef[5], binom_Log_Duong_ST@coef[4] ), lty = 1, lwd = 2, col = cols[4])
mtext('Duong et al', side = 3, line = 0.5)

plot(LogV_vector, LogFunction(LogV_vector, binom_Log_ST@coef[5], binom_Log_ST@coef[1])
     ,type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1), col = cols[1])
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_ST@coef[5], binom_Log_ST@coef[2] ), lty = 1, lwd = 2, col = cols[2])
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_ST@coef[5], binom_Log_ST@coef[3] ), lty = 1, lwd = 2, col = cols[3])
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_ST@coef[5], binom_Log_ST@coef[4] ), lty = 1, lwd = 2, col = cols[4])
mtext('Ferguson et al', side = 3, line = 0.5)

plot(LogV_vector, LogFunction(LogV_vector, binom_Log_Nguyen_ST@coef[5], binom_Log_Nguyen_ST@coef[1])
     ,type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1), col = cols[1])
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_Nguyen_ST@coef[5], binom_Log_Nguyen_ST@coef[2] ), lty = 1, lwd = 2, col = cols[2])
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_Nguyen_ST@coef[5], binom_Log_Nguyen_ST@coef[3] ), lty = 1, lwd = 2, col = cols[3])
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_Nguyen_ST@coef[5], binom_Log_Nguyen_ST@coef[4] ), lty = 1, lwd = 2, col = cols[4])
mtext('Nguyen et al', side = 3, line = 0.5)

plot(LogV_vector, LogFunction(LogV_vector, binom_Log_all_ST@coef[5], binom_Log_all_ST@coef[1])
     ,type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1), col = cols[1])
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_all_ST@coef[5], binom_Log_all_ST@coef[2] ), lty = 1, lwd = 2, col = cols[2])
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_all_ST@coef[5], binom_Log_all_ST@coef[3] ), lty = 1, lwd = 2, col = cols[3])
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_all_ST@coef[5], binom_Log_all_ST@coef[4] ), lty = 1, lwd = 2, col = cols[4])
mtext('Combined', side = 3, line = 0.5)

h=legend(x=-0.5, y=1, c('DENV1','DENV2','DENV3','DENV4'),
         col = cols,
         lwd = c(1,1,1,1), lty = 1, bty = 'n',cex=.8)

mtext('Viremia (log10)', side = 1, line = 1, outer=TRUE)
mtext('Probability of human-mosquito transmission', side = 2, line = 1, outer=TRUE)

dev.off()

# Compare binom vs beta-binom: uncertainty  -------------------------------
LogV_vector = seq(0,12,by = 0.1)
samples = 1e3
Prob_inf_Log_Duong_symp <- Prob_inf_Log_Duong_presymp <- Prob_inf_Log_Duong_asymp <- matrix(NA,ncol = length(LogV_vector), nrow = samples)
Prob_inf_Log_bb_Duong_symp <- Prob_inf_Log_bb_Duong_presymp <- Prob_inf_Log_bb_Duong_asymp <- matrix(NA,ncol = length(LogV_vector), nrow = samples)

MV.Coefs = rmvn(n=samples, mu = log.reg.mod.Duong.2.ng$coefficients, V <- vcov(log.reg.mod.Duong.2.ng))
for (ii in 1:dim(MV.Coefs)[1]){
  Prob_inf_Log_Duong_symp[ii,] = LogFunction(LogV_vector, MV.Coefs[ii,1], MV.Coefs[ii,2] )
}
Prob_inf_Log_Duong_symp = Prob_inf_Log_Duong_symp[order(Prob_inf_Log_Duong_symp[,60]),]
for (ii in 1:dim(MV.Coefs)[1]){
  Prob_inf_Log_Duong_presymp[ii,] = LogFunction(LogV_vector, MV.Coefs[ii,1] + MV.Coefs[ii,4], MV.Coefs[ii,2] )
}
Prob_inf_Log_Duong_presymp = Prob_inf_Log_Duong_presymp[order(Prob_inf_Log_Duong_presymp[,60]),]
for (ii in 1:dim(MV.Coefs)[1]){
  Prob_inf_Log_Duong_asymp[ii,] = LogFunction(LogV_vector, MV.Coefs[ii,1] + MV.Coefs[ii,3], MV.Coefs[ii,2] )
}
Prob_inf_Log_Duong_asymp = Prob_inf_Log_Duong_asymp[order(Prob_inf_Log_Duong_asymp[,60]),]

MV.Coefs = rmvn(n=samples, mu = t(log.reg.mod.Duong.bb$coefficients)[1,], log.reg.mod.Duong.bb$vcov)
for (ii in 1:dim(MV.Coefs)[1]){
  Prob_inf_Log_bb_Duong_symp[ii,] = LogFunction(LogV_vector, MV.Coefs[ii,1], MV.Coefs[ii,2] )
}
Prob_inf_Log_bb_Duong_symp = Prob_inf_Log_bb_Duong_symp[order(Prob_inf_Log_bb_Duong_symp[,60]),]
for (ii in 1:dim(MV.Coefs)[1]){
  Prob_inf_Log_bb_Duong_presymp[ii,] = LogFunction(LogV_vector, MV.Coefs[ii,1] +  MV.Coefs[ii,4], MV.Coefs[ii,2] )
}
Prob_inf_Log_bb_Duong_presymp = Prob_inf_Log_bb_Duong_presymp[order(Prob_inf_Log_bb_Duong_presymp[,60]),]
for (ii in 1:dim(MV.Coefs)[1]){
  Prob_inf_Log_bb_Duong_asymp[ii,] = LogFunction(LogV_vector, MV.Coefs[ii,1] +  MV.Coefs[ii,3], MV.Coefs[ii,2] )
}
Prob_inf_Log_bb_Duong_asymp = Prob_inf_Log_bb_Duong_asymp[order(Prob_inf_Log_bb_Duong_asymp[,60]),]


pdf('FigS2.4.pdf',width=7.5,height=3.00)
par(mfrow = c(1,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,2))  

plot(LogV_vector, LogFunction(LogV_vector, log.reg.mod.Duong.2.ng$coefficients[1], log.reg.mod.Duong.2.ng$coefficients[2] )
     ,type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1), col = 'black')
matplot(replicate(samples - 50, LogV_vector),t(data.frame(Prob_inf_Log_bb_Duong_symp[26:975,])), type = 'l', col = 'lightblue', add = T) # Hardcoded indices: ADJUST
matplot(replicate(samples - 50, LogV_vector),t(data.frame(Prob_inf_Log_Duong_symp[26:975,])), type = 'l', col = 'blue', add = T)
lines(LogV_vector,LogFunction(LogV_vector, log.reg.mod.Duong.2.ng$coefficients[1], log.reg.mod.Duong.2.ng$coefficients[2]), lty = 1, lwd = 2)
lines(LogV_vector,LogFunction(LogV_vector, log.reg.mod.Duong.bb$coefficients[1], log.reg.mod.Duong.bb$coefficients[2]), lty = 2, lwd = 2)
mtext('symptomatic', side = 3, line = 1)

h=legend(x=-.5, y= 1.05, c('binomial fit (mle)','beta-binomial fit (mle)','binomial fit (95% CI)','beta-binomial fit (95% CI)'),
         col = c('black','black','blue','lightblue'),
         lwd = c(1,1,8,8), lty = c(1,2,1,1), bty = 'n',cex=.8)

plot(LogV_vector, LogFunction(LogV_vector, log.reg.mod.Duong.2.ng$coefficients[1] + log.reg.mod.Duong.2.ng$coefficients[4], log.reg.mod.Duong.2.ng$coefficients[2] )
     ,type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1), col = 'black')
matplot(replicate(samples- 50, LogV_vector),t(data.frame(Prob_inf_Log_bb_Duong_presymp[26:975,])), type = 'l', col = 'lightblue', add = T)
matplot(replicate(samples- 50, LogV_vector),t(data.frame(Prob_inf_Log_Duong_presymp[26:975,])), type = 'l', col = 'blue', add = T)
lines(LogV_vector,LogFunction(LogV_vector, log.reg.mod.Duong.2.ng$coefficients[1] + log.reg.mod.Duong.2.ng$coefficients[4], log.reg.mod.Duong.2.ng$coefficients[2]), lty = 1, lwd = 2)
lines(LogV_vector,LogFunction(LogV_vector, log.reg.mod.Duong.bb$coefficients[1] + log.reg.mod.Duong.bb$coefficients[4], log.reg.mod.Duong.bb$coefficients[2]), lty = 2, lwd = 2)
mtext('presymptomatic', side = 3, line = 1)

plot(LogV_vector, LogFunction(LogV_vector, log.reg.mod.Duong.2.ng$coefficients[1] + log.reg.mod.Duong.2.ng$coefficients[3], log.reg.mod.Duong.2.ng$coefficients[2] )
     ,type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1), col = 'black')
matplot(replicate(samples-50, LogV_vector),t(data.frame(Prob_inf_Log_bb_Duong_asymp[26:975,])), type = 'l', col = 'lightblue', add = T)
matplot(replicate(samples-50, LogV_vector),t(data.frame(Prob_inf_Log_Duong_asymp[26:975,])), type = 'l', col = 'blue', add = T)
lines(LogV_vector,LogFunction(LogV_vector, log.reg.mod.Duong.2.ng$coefficients[1] + log.reg.mod.Duong.2.ng$coefficients[3], log.reg.mod.Duong.2.ng$coefficients[2]), lty = 1, lwd = 2)
lines(LogV_vector,LogFunction(LogV_vector, log.reg.mod.Duong.bb$coefficients[1] + log.reg.mod.Duong.bb$coefficients[3], log.reg.mod.Duong.bb$coefficients[2]), lty = 2, lwd = 2)
mtext('asymptomatic', side = 3, line = 1)

mtext('Viremia (log10)', side = 1, line = 1, outer=TRUE, cex = .8)
mtext('Probability of human-mosquito transmission', side = 2, line = 1, outer=TRUE, cex = .8)

dev.off()

# plot curves by infection class + uncertainty ----------------------------

pdf('FigS2.5.pdf',width=7.5,height=4.75)
par(mfrow = c(1,1),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,2))  

plot(LogV_vector, LogFunction(LogV_vector, log.reg.mod.Duong.2.ng$coefficients[1], log.reg.mod.Duong.2.ng$coefficients[2] )
     ,type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1), col = 'black')
matplot(replicate(samples - 50, LogV_vector),t(data.frame(Prob_inf_Log_Duong_symp[26:975,])), type = 'l', col = alpha(colour = 'lightblue', alpha = 0.1), add = T)
matplot(replicate(samples- 50, LogV_vector),t(data.frame(Prob_inf_Log_Duong_presymp[26:975,])), type = 'l', col = alpha(colour = 'lightgreen', alpha = 0.1), add = T)
matplot(replicate(samples-50, LogV_vector),t(data.frame(Prob_inf_Log_Duong_asymp[26:975,])), type = 'l', col = alpha(colour = 'pink', alpha = 0.1), add = T)
lines(LogV_vector,LogFunction(LogV_vector, log.reg.mod.Duong.2.ng$coefficients[1], log.reg.mod.Duong.2.ng$coefficients[2]), lty = 1, lwd = 2)
lines(LogV_vector,LogFunction(LogV_vector, log.reg.mod.Duong.2.ng$coefficients[1] + log.reg.mod.Duong.2.ng$coefficients[3], log.reg.mod.Duong.2.ng$coefficients[2]), lty = 3, lwd = 2)
lines(LogV_vector,LogFunction(LogV_vector, log.reg.mod.Duong.2.ng$coefficients[1] + log.reg.mod.Duong.2.ng$coefficients[4], log.reg.mod.Duong.2.ng$coefficients[2]), lty = 2, lwd = 2)

mtext('Viremia (log10)', side = 1, line = 1, outer=TRUE)
mtext('Probability of human-mosquito transmission', side = 2, line = 1, outer=TRUE)

h=legend(x=6.5, y=0.3, c('symp','presym','asym'),
         col = c('black','black','black'),
         lwd = c(2,2,2), lty = c(1,2,3), bty = 'n',cex=.8)

dev.off()

# plot curves by infection class + uncertainty (random intercepts and binomial dist) ----------------------------
LogV_vector = seq(0,12,by = 0.1)
samples = 1e3

Prob_inf_Log_Duong_symp <- Prob_inf_Log_Duong_presymp <- Prob_inf_Log_Duong_asymp <- matrix(NA,ncol = length(LogV_vector), nrow = samples)

MV.Coefs = rmvn(n=samples, mu = binom_Log_Duong.symp@coef, V = vcov(binom_Log_Duong.symp))
for (ii in 1:dim(MV.Coefs)[1]){
  Prob_inf_Log_Duong_symp[ii,] = LogFunction(LogV_vector, MV.Coefs[ii,1], MV.Coefs[ii,2] )
}
Prob_inf_Log_Duong_symp = Prob_inf_Log_Duong_symp[order(Prob_inf_Log_Duong_symp[,60]),]

MV.Coefs = rmvn(n=samples, mu = binom_Log_Duong.presymp@coef, V = vcov(binom_Log_Duong.presymp))
for (ii in 1:dim(MV.Coefs)[1]){
  Prob_inf_Log_Duong_presymp[ii,] = LogFunction(LogV_vector, MV.Coefs[ii,1], MV.Coefs[ii,2] )
}
Prob_inf_Log_Duong_presymp = Prob_inf_Log_Duong_presymp[order(Prob_inf_Log_Duong_presymp[,60]),]

MV.Coefs = rmvn(n=samples, mu = binom_Log_Duong.asymp@coef, V = vcov(binom_Log_Duong.asymp))
for (ii in 1:dim(MV.Coefs)[1]){
  Prob_inf_Log_Duong_asymp[ii,] = LogFunction(LogV_vector, MV.Coefs[ii,1], MV.Coefs[ii,2] )
}
Prob_inf_Log_Duong_asymp = Prob_inf_Log_Duong_asymp[order(Prob_inf_Log_Duong_asymp[,60]),]

pdf('FigS2.6.pdf',width=7.5,height=4.75)
par(mfrow = c(1,1),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,2))  

plot(LogV_vector, LogFunction(LogV_vector, log.reg.mod.Duong.2.ng$coefficients[1], log.reg.mod.Duong.2.ng$coefficients[2] )
     ,type = 'l', xlab = "", ylab = "", lwd = 2, ylim=c(0,1), col = 'white')
matplot(replicate(samples - 50, LogV_vector),t(data.frame(Prob_inf_Log_Duong_symp[26:975,])), type = 'l', col = alpha(colour = 'lightblue', alpha = 0.1), add = T)
matplot(replicate(samples- 50, LogV_vector),t(data.frame(Prob_inf_Log_Duong_presymp[26:975,])), type = 'l', col = alpha(colour = 'lightgreen', alpha = 0.1), add = T)
matplot(replicate(samples-50, LogV_vector),t(data.frame(Prob_inf_Log_Duong_asymp[26:975,])), type = 'l', col = alpha(colour = 'pink', alpha = 0.1), add = T)
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_Duong.symp@coef[1], binom_Log_Duong.symp@coef[2]), lty = 1, lwd = 2)
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_Duong.presymp@coef[1] , binom_Log_Duong.presymp@coef[2]), lty = 3, lwd = 2)
lines(LogV_vector,LogFunction(LogV_vector, binom_Log_Duong.asymp@coef[1] , binom_Log_Duong.asymp@coef[2]), lty = 2, lwd = 2)

mtext('Viremia (log10)', side = 1, line = 1, outer=TRUE)
mtext('Probability of human-mosquito transmission', side = 2, line = 1, outer=TRUE)

h=legend(x=6.5, y=0.3, c('symp','presym','asym'),
         col = c('black','black','black'),
         lwd = c(2,2,2), lty = c(1,2,3), bty = 'n',cex=.8)

dev.off()







