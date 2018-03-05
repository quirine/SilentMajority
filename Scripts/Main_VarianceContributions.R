# decomposition of variance

# Get variances of the full model -----------------------------------------
load('Workspace.lowerAsym.RData')
Y.Variances = array(data = NA, dim = 5)
Y.Variances[1] = var(AUC.asympto.1); Y.Variances[2] = var(AUC.asympto.2); 
Y.Variances[3] = var(AUC.sympto.w.1); Y.Variances[4] = var(AUC.sympto.w.2); Y.Variances[5] = var(AUC.sympto.w.3)

Y.Variances.PPS = array(data = NA, dim = 3)
Y.Variances.PPS[1] = var(Prop.infectivity.prior.to.IIP.1); Y.Variances.PPS[2] = var(Prop.infectivity.prior.to.IIP.2); Y.Variances.PPS[3] = var(Prop.infectivity.prior.to.IIP.3)


# Get variances for all but viremia ---------------------------------------

load('Workspace.Uncertainty.viremia.RData')
Vir.Variances = array(data = NA, dim = 5)
Vir.Variances[1] = var(AUC.asympto.1); Vir.Variances[2] = var(AUC.asympto.2); 
Vir.Variances[3] = var(AUC.sympto.w.1); Vir.Variances[4] = var(AUC.sympto.w.2); Vir.Variances[5] = var(AUC.sympto.w.3)

Vir.Variances.PPS = array(data = NA, dim = 3)
Vir.Variances.PPS[1] = var(Prop.infectivity.prior.to.IIP.1); Vir.Variances.PPS[2] = var(Prop.infectivity.prior.to.IIP.2); Vir.Variances.PPS[3] = var(Prop.infectivity.prior.to.IIP.3)


# Get variances for all but infectiousness ---------------------------------------

load('Workspace.Uncertainty.infectiousness.RData')
Inf.Variances = array(data = NA, dim = 5)
Inf.Variances[1] = var(AUC.asympto.1); Inf.Variances[2] = var(AUC.asympto.2)
Inf.Variances[3] = var(AUC.sympto.w.1); Inf.Variances[4] = var(AUC.sympto.w.2); Inf.Variances[5] = var(AUC.sympto.w.3)

Inf.Variances.PPS = array(data = NA, dim = 3)
Inf.Variances.PPS[1] = var(Prop.infectivity.prior.to.IIP.1); Inf.Variances.PPS[2] = var(Prop.infectivity.prior.to.IIP.2); Inf.Variances.PPS[3] = var(Prop.infectivity.prior.to.IIP.3)


# Get variances for all but ratios ---------------------------------------

load('Workspace.Uncertainty.ratios.RData')
Ratio.Variances = array(data = NA, dim = 5)
Ratio.Variances[1] = var(AUC.asympto.1); Ratio.Variances[2] = var(AUC.asympto.2)
Ratio.Variances[3] = var(AUC.sympto.w.1); Ratio.Variances[4] = var(AUC.sympto.w.2); Ratio.Variances[5] = var(AUC.sympto.w.3)

Ratio.Variances.PPS = array(data = NA, dim = 3)
Ratio.Variances.PPS[1] = var(Prop.infectivity.prior.to.IIP.1); Ratio.Variances.PPS[2] = var(Prop.infectivity.prior.to.IIP.2); Ratio.Variances.PPS[3] = var(Prop.infectivity.prior.to.IIP.3)

# Get variances for all but iip ---------------------------------------

load('Workspace.Uncertainty.iip.lognorm.RData')
IIP.Variances = array(data = NA, dim = 5)
IIP.Variances[1] = var(AUC.asympto.1); IIP.Variances[2] = var(AUC.asympto.2); 
IIP.Variances[3] = var(AUC.sympto.w.1); IIP.Variances[4] = var(AUC.sympto.w.2); IIP.Variances[5] = var(AUC.sympto.w.3)

IIP.Variances.PPS = array(data = NA, dim = 3)
IIP.Variances.PPS[1] = var(Prop.infectivity.prior.to.IIP.1); IIP.Variances.PPS[2] = var(Prop.infectivity.prior.to.IIP.2); IIP.Variances.PPS[3] = var(Prop.infectivity.prior.to.IIP.3)

# Get variances for all but iip - mean ---------------------------------------

load('Workspace.Uncertainty.iip.RData')
IIP.mean.Variances = array(data = NA, dim = 5)
IIP.mean.Variances[1] = var(AUC.asympto.1); IIP.mean.Variances[2] = var(AUC.asympto.2); 
IIP.mean.Variances[3] = var(AUC.sympto.w.1); IIP.mean.Variances[4] = var(AUC.sympto.w.2); IIP.mean.Variances[5] = var(AUC.sympto.w.3)

IIP.mean.Variances.PPS = array(data = NA, dim = 3)
IIP.mean.Variances.PPS[1] = var(Prop.infectivity.prior.to.IIP.1); IIP.mean.Variances.PPS[2] = var(Prop.infectivity.prior.to.IIP.2); IIP.mean.Variances.PPS[3] = var(Prop.infectivity.prior.to.IIP.3)


# Derive first order variance ---------------------------------------------
S.Vir = 1- (Vir.Variances / Y.Variances)
S.Inf = 1 - (Inf.Variances / Y.Variances)
S.Ratio = 1 - (Ratio.Variances / Y.Variances)
S.IIP = 1 - (IIP.Variances / Y.Variances)
S.mean.IIP = 1 - (IIP.mean.Variances / Y.Variances)

S.Vir.PPS = 1- (Vir.Variances.PPS / Y.Variances.PPS)
S.Inf.PPS = 1 - (Inf.Variances.PPS / Y.Variances.PPS)
S.Ratio.PPS = 1 - (Ratio.Variances.PPS / Y.Variances.PPS)
S.IIP.PPS = 1 - (IIP.Variances.PPS / Y.Variances.PPS)
S.mean.IIP.PPS = 1 - (IIP.mean.Variances.PPS / Y.Variances.PPS)


# Set directory to save figures -------------------------------------------
path <- '../figures'
cd <- getwd()

# Plot total effect index decomposition of variance (pie chart) ------------------------------
setwd(path)
pdf('SFig4.pdf',width=4.5,height=7.5)
par(mfrow = c(4,2), mar=c(0,0,1.5,0),oma=c(0,2,3,1))

pie(c(S.Vir[1], S.Inf[1], S.mean.IIP[1]),labels = c('','',''), col = c('white','cyan','deeppink'),radius = 1) 
Lines <- list("Asymptomatic","primary")
mtext(do.call(expression, Lines),2,line=c(0,-1.5),cex = 1)
mtext('Contribution to uncertainty about',3,line=.5, outer = TRUE,cex = 1)
mtext('Net infectiousness',3,line=.5, outer = FALSE,cex = 1)
plot(1, type="n", axes=F, xlab="", ylab="")
Lines <- list("Proportion of infectiousness"," prior to symptoms")
mtext(do.call(expression, Lines),3,line=c(0.2,-1), outer = FALSE, cex = 1)
legend('bottom',inset = c(-0,.15),c('Viremia','Infectiousness','IIP'),
       horiz = FALSE, fill = c('white','cyan','deeppink'), 
       bty = 'n',cex=1.8, pt.cex = 1, density = c(NA,NA,NA))

pie(c(S.Vir[2], S.Inf[2], S.mean.IIP[2]), labels = c('','',''), col = c('white','cyan','deeppink'),radius = 1) 
Lines <- list("Asymptomatic","secondary")
mtext(do.call(expression, Lines),2,line=c(0,-1.5),cex = 1)
plot(1, type="n", axes=F, xlab="", ylab="")

pie(c(S.Vir[3], S.Inf[3], S.mean.IIP[3]), labels = c('','',''), col = c('white','cyan','deeppink'),radius = 1) 
Lines <- list("Symptomatic","primary")
mtext(do.call(expression, Lines),2,line=c(0,-1.5),cex = 1)
pie(c(S.Vir.PPS[1], S.Inf.PPS[1], S.mean.IIP.PPS[1]), labels = c('','',''), col = c('white','cyan','deeppink'),radius = 1) 

pie(c(S.Vir[4], S.Inf[4], S.mean.IIP[4]), 
        labels = c('', '', ''), col = c('white','cyan','deeppink'),radius = 1) 
Lines <- list("Symptomatic","secondary")
mtext(do.call(expression, Lines),2,line=c(0,-1.5),cex = 1)
pie(c(S.Vir.PPS[2], S.Inf.PPS[2], S.mean.IIP.PPS[2]),labels = c('','',''), col = c('white','cyan','deeppink'),radius = 1)

dev.off()
setwd(cd)
