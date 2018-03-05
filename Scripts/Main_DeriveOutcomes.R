
# load packages -----------------------------------------------------------
library(coda)

# Load source files -------------------------------------------------------
source('Functions_Outcomes.R')

# Derive the relevant proportions of the prevalence -----------------------
XProps.brazil <- Derive.outcomes.X(X = X.brazil, Prop.infectivity.prior.to.IIP.1, Prop.infectivity.prior.to.IIP.2, 
                              Prop.infectivity.prior.to.IIP.3, foi.vec)
XProps.brazil.4 <- Derive.outcomes.X(X = X.brazil.4, Prop.infectivity.prior.to.IIP.1, Prop.infectivity.prior.to.IIP.2, 
                              Prop.infectivity.prior.to.IIP.3, foi.vec, inf = 4)
XProps.brazil.emerging <- Derive.outcomes.X(X = X.brazil.emerging, Prop.infectivity.prior.to.IIP.1, Prop.infectivity.prior.to.IIP.2, 
                              Prop.infectivity.prior.to.IIP.3, foi.vec)
XProps.thailand <- Derive.outcomes.X(X = X.thailand, Prop.infectivity.prior.to.IIP.1, Prop.infectivity.prior.to.IIP.2, 
                                       Prop.infectivity.prior.to.IIP.3, foi.vec)
XProps.thailand.4 <- Derive.outcomes.X(X = X.thailand.4, Prop.infectivity.prior.to.IIP.1, Prop.infectivity.prior.to.IIP.2, 
                                Prop.infectivity.prior.to.IIP.3, foi.vec, inf = 4)
XProps.thailand.emerging <- Derive.outcomes.X(X = X.thailand.emerging, Prop.infectivity.prior.to.IIP.1, Prop.infectivity.prior.to.IIP.2, 
                                       Prop.infectivity.prior.to.IIP.3, foi.vec)

# Calculate these in contribution space -----------------------------------

FOIProps.brazil <- Derive.outcomes.FOI(Prop.FOI = Prop.FOI.All.brazil, XProps = XProps.brazil)
FOIProps.brazil.4 <- Derive.outcomes.FOI(Prop.FOI = Prop.FOI.All.brazil.4, XProps = XProps.brazil.4)
FOIProps.brazil.emerging <- Derive.outcomes.FOI(Prop.FOI = Prop.FOI.All.brazil.emerging, XProps = XProps.brazil.emerging)
FOIProps.brazil.DHF <- Derive.outcomes.FOI(Prop.FOI = Prop.FOI.All.brazil.DHF, XProps = XProps.brazil, Prop.FOI.DHF = Prop.FOI.All.brazil.DHF.severe )

FOIProps.thailand <- Derive.outcomes.FOI(Prop.FOI = Prop.FOI.All.thailand, XProps = XProps.thailand)
FOIProps.thailand.4 <- Derive.outcomes.FOI(Prop.FOI = Prop.FOI.All.thailand.4, XProps = XProps.thailand.4)
FOIProps.thailand.emerging <- Derive.outcomes.FOI(Prop.FOI = Prop.FOI.All.thailand.emerging, XProps = XProps.thailand.emerging)

# Estimate CI's -----------------------------------------------------------
FOIProps.brazil.CIs <- Derive.outcomes.FOI.CIs(Prop.FOI = Prop.FOI.All.brazil, XProps = XProps.brazil)
FOIProps.brazil.4.CIs <- Derive.outcomes.FOI.CIs(Prop.FOI = Prop.FOI.All.brazil.4, XProps = XProps.brazil.4)
FOIProps.brazil.emerging.CIs <- Derive.outcomes.FOI.CIs(Prop.FOI = Prop.FOI.All.brazil.emerging, XProps = XProps.brazil.emerging)
FOIProps.brazil.DHF.CIs <- Derive.outcomes.FOI.CIs(Prop.FOI = Prop.FOI.All.brazil.DHF, XProps = XProps.brazil,Prop.FOI.DHF = Prop.FOI.All.brazil.DHF.severe)
FOIProps.thailand.CIs <- Derive.outcomes.FOI.CIs(Prop.FOI = Prop.FOI.All.thailand, XProps = XProps.thailand)
FOIProps.thailand.4.CIs <- Derive.outcomes.FOI.CIs(Prop.FOI = Prop.FOI.All.thailand.4, XProps = XProps.thailand.4)
FOIProps.thailand.emerging.CIs <- Derive.outcomes.FOI.CIs(Prop.FOI = Prop.FOI.All.thailand.emerging, XProps = XProps.thailand.emerging)

# and pre-sym CI for default ---------------------------------------------------------
Prop.sym.prim.presym.all <- XProps.brazil$Prop.sym.prim[length(foi.vec)] * Prop.infectivity.prior.to.IIP.1
Prop.sym.sec.presym.all <- XProps.brazil$Prop.sym.sec[length(foi.vec)] * Prop.infectivity.prior.to.IIP.2 
Prop.sym.presym.CI <- HPDinterval(mcmc(Prop.sym.prim.presym.all + Prop.sym.sec.presym.all), 0.95)
# tests -------------------------------------------------------------------

A.1vsS.2 <- wilcox.test(AUC.asympto.1, AUC.sympto.w.2)  
A.2vsS.2 <- wilcox.test(AUC.asympto.2, AUC.sympto.w.2)
A.2vsS.1 <- wilcox.test(AUC.asympto.2, AUC.sympto.w.1)
A.1vsS.1 <- wilcox.test(AUC.asympto.1, AUC.sympto.w.1)
A.1vsA.2 <- wilcox.test(AUC.asympto.1, AUC.asympto.2)
S.1vsS.2 <- wilcox.test(AUC.sympto.w.1, AUC.sympto.w.2)

WSR.W <- data.frame(A.1vsS.2$statistic, A.2vsS.2$statistic, A.2vsS.1$statistic,
                    A.1vsS.1$statistic, A.1vsA.2$statistic, S.1vsS.2$statistic)

WSR.P <- data.frame(A.1vsS.2$p.value, A.2vsS.2$p.value ,A.2vsS.1$p.value,
                   A.1vsS.1$p.value, A.1vsA.2$p.value, S.1vsS.2$p.value)


ks.A.1vsS.2 <- ks.test(AUC.asympto.1, AUC.sympto.w.2,alternative='two.sided')  
ks.A.2vsS.2 <- ks.test(AUC.asympto.2, AUC.sympto.w.2,alternative='two.sided')
ks.A.2vsS.1 <- ks.test(AUC.asympto.2, AUC.sympto.w.1,alternative='two.sided')
ks.A.1vsS.1 <- ks.test(AUC.asympto.1, AUC.sympto.w.1,alternative='two.sided')
ks.A.1vsA.2 <- ks.test(AUC.asympto.1, AUC.asympto.2,alternative='two.sided')
ks.S.1vsS.2 <- ks.test(AUC.sympto.w.1, AUC.sympto.w.2,alternative='two.sided')

KS.D <- data.frame(ks.A.1vsS.2$statistic , ks.A.2vsS.2$statistic ,ks.A.2vsS.1$statistic,
                   ks.A.1vsS.1$statistic, ks.A.1vsA.2$statistic, ks.S.1vsS.2$statistic)

KS.P <- data.frame(ks.A.1vsS.2$p.value, ks.A.2vsS.2$p.value ,ks.A.2vsS.1$p.value,
                   ks.A.1vsS.1$p.value, ks.A.1vsA.2$p.value, ks.S.1vsS.2$p.value)


