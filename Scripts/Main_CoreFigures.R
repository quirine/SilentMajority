# script to create figures for core part of the paper

# Load packages -----------------------------------------------------------
library(scales)
library(plotrix)
library(vioplot)
library(deSolve)

# Set directory to save figures -------------------------------------------
figure.path <- '../figures'
cd <- getwd()

# Load source files -------------------------------------------------------
source('Functions_ProportionImmuneHistory.R')
source('Functions_ProportionImmuneFigures.R')
source('Functions_InfectiousnessFigures.R')
source('Functions_FOIFigures.R')

# Nishiura data for grounding of infectiousness results -------------------

Nish.Inf <- c(0.248691183213025, 0.803568191753655, 1.00, 0.967604785311397, 0.499712113619158, -0.00321206563809486)
Nish.Inf.low <- c(0, 0.4489, 1.00, 0.8904, 0.0951, 0)
Nish.Inf.high <- c(0.6774, 1.00, 1.00, 1.00, 0.90, 0 )
Nish.Ints <- c(sum(Nish.Inf), sum(Nish.Inf.low), sum(Nish.Inf.high))

# Figure 2 ----------------------------------------------------------------

setwd(figure.path)
tiff('Fig2.tiff',width=4.75,height=4.75,units='in',res=300)
par(mfrow = c(4,3),mai=c(0.1,0.3,0.1,0.1) ,oma=c(3,4,1,1))  #mar=c(1,3,1,1.9)

ref <- median(AUC.sympto.w.1)

par(xaxt='n')
par(yaxt='s')
Figure.Viralload('asym',single='prim')
mtext(expression(bold("A")), side = 3, adj = 0.05, 
      line = -1.3,cex=.8)
mtext(expression("Asymptomatic - 1" ^"o"),2,line=4, cex=.7)

Figure.Infectiousness.Unweighted('asym',1)
mtext(expression(bold("B")), side = 3, adj = 0.05, 
      line = -1.3,cex=.8)

d <- hist((AUC.asympto.1), breaks = 30, plot = FALSE) 
par(yaxt='s')
d$density <- d$density*diff(d$breaks)[1]
plot(d,main='', xlim=c(0,15), col = 'lightskyblue', border = FALSE, freq = FALSE)
mtext(expression(bold("C")), side = 3, adj = 0.05, 
      line = -1.3,cex=.8)
box()
abline(v = ref, col = 'navyblue', lty = 2, lwd = 2)
abline(v = median(AUC.asympto.1), col = 'navyblue', lty = 1, lwd = 2)

par(yaxt='s')
Figure.Viralload('asym',single='sec')
mtext(expression(bold("D")), side = 3, adj = 0.05, 
      line = -1.3,cex=.8)
mtext(expression("Asymptomatic - 2" ^"o"),2,line=4, cex=.7)
Figure.Infectiousness.Unweighted('asym',2)
mtext(expression(bold("E")), side = 3, adj = 0.05, 
      line = -1.3,cex=.8)

d <- hist(AUC.asympto.2, breaks = 30, plot = FALSE) 
par(yaxt='s')
d$density <- d$density*diff(d$breaks)[1]
plot(d,main='', xlim=c(0,15), col = 'lightskyblue', border = FALSE, freq = FALSE)
mtext(expression(bold("F")), side = 3, adj = 0.05, 
      line = -1.3,cex=.8)
box()
abline(v = ref, col = 'navyblue', lty = 2, lwd = 2)
abline(v = median(AUC.asympto.2), col = 'navyblue', lty = 1, lwd = 2)


par(yaxt='s')
Figure.Viralload('sym',single='prim')
mtext(expression(bold("G")), side = 3, adj = 0.05, 
      line = -1.3,cex=.8)
mtext(expression("Symptomatic - 1" ^"o"),2,line=4, cex=.7)
Figure.Infectiousness.Weighted(1)
mtext(expression(bold("H")), side = 3, adj = 0.05, 
      line = -1.3,cex=.8)
d <- hist(AUC.sympto.w.1, breaks = 20, plot = FALSE) 
par(yaxt='s')
d$density <- d$density*diff(d$breaks)[1]
plot(d,main='', xlim=c(0,15), type = 'n', xlab = "", ylab = "", border = FALSE, freq = FALSE)
plot(d,main='', xlim=c(0,15), col = 'lightskyblue', border = FALSE, freq = FALSE, add = TRUE)
mtext(expression(bold("I")), side = 3, adj = 0.05, 
      line = -1.3,cex=.8)
box()
abline(v = Nish.Ints[2], col = 'darkred', lty = 3, lwd = 1.5)
abline(v = Nish.Ints[3], col = 'darkred', lty = 3, lwd = 1.5)
abline(v = ref, col = 'navyblue', lty = 2, lwd = 2)
abline(v = median(AUC.sympto.w.1), col = 'navyblue', lty = 1, lwd = 2)
abline(v = Nish.Ints[1], col = 'darkred', lty = 1, lwd = 2)

par(yaxt='s')
par(xaxt='s')
Figure.Viralload('sym',single='sec')
mtext(expression(bold("J")), side = 3, adj = 0.05, 
      line = -1.3,cex=.8)
mtext(expression("Symptomatic - 2" ^"o"),2,line=4, cex=.7)
mtext('Time (days)',1,line=2.25,cex=.7)
Figure.Infectiousness.Weighted(2)
mtext(expression(bold("K")), side = 3, adj = 0.05, 
      line = -1.3,cex=.8)
mtext('Time (days)',1,line=2.25,cex=.7)
d <- hist(AUC.sympto.w.2, breaks = 20, plot = FALSE) 
par(yaxt='s')
d$density <- d$density*diff(d$breaks)[1]
plot(d,main='', xlim=c(0,15), col = 'lightskyblue', border = FALSE, freq = FALSE)
mtext(expression(bold("L")), side = 3, adj = 0.05, 
      line = -1.3,cex=.8)
box()
abline(v = ref, col = 'navyblue', lty = 2, lwd = 2)
abline(v = mean(AUC.sympto.w.2), col = 'navyblue', lty = 1, lwd = 2)
mtext(expression("Net Infectiousness"),1,line=2.25,cex=.7)

mtext(expression("Viremia (log"[10]*"cDNA copies/mL)" ), side = 2, outer = TRUE, line = 0.1, cex=.7)
mtext('Probability of mosquito infection', side = 2, outer = TRUE, line = -10.7, cex=.7)  
mtext('Probability density', side = 2, outer = TRUE, line = -21, cex=.7)  

dev.off()
setwd(cd)

# Figure 3 ----------------------------------------------------------------
setwd(figure.path)
tiff('Fig3.tiff',width=10,height=5,units='in',res=300)
par(mfrow = c(2,2), mar=c(0.5,4,1,1),oma=c(3,2,.5,1))

foi = .01
seroage.brazil <- sero.age.calculator(agepop.brazil, foi,inv.cross.im, ode = 'yes')
seroprev.brazil = matrix(data=0,100,5)
seroprev.brazil[,1] = (seroage.brazil[,1] * agepop.brazil$pop) / sum(agepop.brazil$pop)
seroprev.brazil[,2] = (seroage.brazil[,2] * agepop.brazil$pop) / sum(agepop.brazil$pop)
seroprev.brazil[,3] = (seroage.brazil[,3] * agepop.brazil$pop) / sum(agepop.brazil$pop)
seroprev.brazil[,4] = (seroage.brazil[,4] * agepop.brazil$pop) / sum(agepop.brazil$pop)
seroprev.brazil[,5] = (seroage.brazil[,5] * agepop.brazil$pop) / sum(agepop.brazil$pop)
infec.hist.by.age.figure(agepop.brazil, seroprev.brazil)
Axis(side=1, at = NULL ,labels=FALSE, cex=.6)
Axis(side=2, at = seq(0,.03,.01), labels=seq(0,.03,.01))
mtext('FoI = 0.01', side = 3, line =  0, cex = 1, at = 15)
mtext(expression(bold("A")), side = 3, adj = 0.05, at = -1,
      line = 0.1,cex=1)
floating.pie(
  85,1.5e6/ sum(agepop.brazil$pop),
  colSums(seroprev.brazil),
  radius = 10,  col = c(rgb(0,1,0,seq(0,0.75,by = 0.25)),'black'))
h=legend(x=52, y=0.0205,title = 'Pre-exposed to',c('0 Serotypes','1 Serotype','2 Serotypes','3 Serotypes','4 Fully Immune'), 
         fill = c(rgb(0,1,0,seq(0,0.75,by = 0.25)),'black'),bty = 'n',cex=.8)

A.age.brazil = dis.per.age.calculator(foi, agepop.brazil, 'asym', inv.cross.im, ode = 'yes' ) / sum(agepop.brazil$pop)
I.age.brazil = dis.per.age.calculator(foi, agepop.brazil, 'inapparent', inv.cross.im, ode = 'yes' ) / sum(agepop.brazil$pop)
M.age.brazil = dis.per.age.calculator(foi, agepop.brazil, 'mild', inv.cross.im, ode = 'yes' ) / sum(agepop.brazil$pop)
S.age.brazil = dis.per.age.calculator(foi, agepop.brazil, 'severe', inv.cross.im, ode = 'yes' ) / sum(agepop.brazil$pop)
Immune.age.brazil = immune.per.age.calculator(agepop.brazil, foi, inv.cross.im) / sum(agepop.brazil$pop) 
dis.by.age.figure.detected(A.age.brazil,I.age.brazil, M.age.brazil, S.age.brazil, agepop.brazil,seroage.brazil,Immune.age.brazil,detection.rate = detection.rate)
Axis(side=1, at = NULL ,labels=FALSE)
Axis(side=2, at = seq(0,.03,.01), labels=seq(0,.03,.01))
mtext('FoI = 0.01', side = 3, line = 0, cex = 1, at = 15)
floating.pie(90,1.5e6/ sum(agepop.brazil$pop),c(sum(A.age.brazil),sum(I.age.brazil) - sum(A.age.brazil),
                                                ((sum(M.age.brazil)+sum(S.age.brazil)) - (sum(I.age.brazil) - sum(A.age.brazil)))*0.92,  
                                                ((sum(M.age.brazil)+sum(S.age.brazil)) - (sum(I.age.brazil) - sum(A.age.brazil)))*0.08,  
                                                sum(rowSums(Immune.age.brazil[,c(1,2,3)])*agepop.brazil),sum(Immune.age.brazil[,4]*agepop.brazil)), 
             radius = 10,  col = c(rgb(1,0,0,c(.1,.2,.5,1)),'gray','black'))
h=legend(x=50, y=0.0205,c('Asymptomatic','Inapparent symptomatic','Undetected apparent symptomatic','Detected apparent symptomatic','Cross-immune','Fully immune'), 
         fill = c(rgb(1,0,0,c(.1,.2,.5,1)),'gray','black'),bty = 'n',cex=.8)
mtext(expression(bold("B")), side = 3, adj = 0.05, at = -1,
      line = 0.1,cex=1)

foi = .1
seroage.brazil <- sero.age.calculator(agepop.brazil, foi)
seroprev.brazil = matrix(data=0,100,5)
seroprev.brazil[,1] = (seroage.brazil[,1] * agepop.brazil$pop) / sum(agepop.brazil$pop)
seroprev.brazil[,2] = (seroage.brazil[,2] * agepop.brazil$pop) / sum(agepop.brazil$pop)
seroprev.brazil[,3] = (seroage.brazil[,3] * agepop.brazil$pop) / sum(agepop.brazil$pop)
seroprev.brazil[,4] = (seroage.brazil[,4] * agepop.brazil$pop) / sum(agepop.brazil$pop)
seroprev.brazil[,5] = (seroage.brazil[,5] * agepop.brazil$pop) / sum(agepop.brazil$pop)
infec.hist.by.age.figure(agepop.brazil, seroprev.brazil)
Axis(side=1, at = NULL ,labels=TRUE, cex=.6)
Axis(side=2, at = seq(0,.03,.01), labels=seq(0,.03,.01))
mtext('FoI = 0.1', side = 3, line = 0, cex = 1, at = 14)
mtext('Age (years)', side = 1, line = 2.25, cex = 1.2)
floating.pie(
  85,2.5e6/ sum(agepop.brazil$pop),
  colSums(seroprev.brazil),
  radius = 10,  col = c(rgb(0,1,0,seq(0,0.75,by = 0.25)),'black'))
mtext("% of total population", side = 3, line = c(-1),at=85, cex=.8)
mtext(expression(bold("C")), side = 3, adj = 0.05, at = -1,
      line = 0.1,cex=1)

A.age.brazil = dis.per.age.calculator(foi, agepop.brazil, 'asym', inv.cross.im, ode = 'yes' )  / sum(agepop.brazil$pop)
I.age.brazil = dis.per.age.calculator(foi, agepop.brazil, 'inapparent', inv.cross.im, ode = 'yes' )  / sum(agepop.brazil$pop)
M.age.brazil = dis.per.age.calculator(foi, agepop.brazil, 'mild', inv.cross.im, ode = 'yes' )  / sum(agepop.brazil$pop)
S.age.brazil = dis.per.age.calculator(foi, agepop.brazil, 'severe', inv.cross.im, ode = 'yes' )  / sum(agepop.brazil$pop)
Immune.age.brazil = immune.per.age.calculator(agepop.brazil, foi, inv.cross.im)  / sum(agepop.brazil$pop)
dis.by.age.figure.detected(A.age.brazil,I.age.brazil, M.age.brazil, S.age.brazil, agepop.brazil,seroage.brazil,Immune.age.brazil,detection.rate = detection.rate)
Axis(side=1, at = NULL ,labels=TRUE)
Axis(side=2, at = seq(0,.03,.01), labels=seq(0,.03,.01))
mtext('FoI = 0.1', side = 3, line = 0, cex = 1, at = 14)
mtext('Age (years)', side = 1, line = 2.25, cex = 1.2)
Lines <- list("Proportion of total infection events","by age and clinical outcome")
mtext(do.call(expression, Lines),2,line = c(-29.5,-30.5),outer = TRUE)
Lines <- list("Proportion of total infection events","by age and pre-exposure history")
mtext(do.call(expression, Lines),2,line = c(-1,-2),outer = TRUE)
mtext(expression(bold("D")), side = 3, adj = 0.05, at = -1,
      line = 0.1,cex=1)

floating.pie(90,2.5e6/ sum(agepop.brazil$pop),c(sum(A.age.brazil),sum(I.age.brazil) - sum(A.age.brazil),
                                                ((sum(M.age.brazil)+sum(S.age.brazil)) - (sum(I.age.brazil) - sum(A.age.brazil) ))*0.92, 
                                                ((sum(M.age.brazil)+sum(S.age.brazil)) - (sum(I.age.brazil) - sum(A.age.brazil) ))*0.08, 
                                                sum(rowSums(Immune.age.brazil[,c(1,2,3)])*agepop.brazil),
                                                sum(Immune.age.brazil[,4]*agepop.brazil)), 
             radius = 10,  col = c(rgb(1,0,0,c(.1,.2,.5,1)),'gray','black'))
mtext("% of total population", side = 3, line = c(-1),at=90, cex=.8)

dev.off()
setwd(cd)


# Figure 4 ---------------------------------------------------------------

setwd(path)
tiff('Fig4.tiff',width=7.5,height=4.5,units='in',res=300)
par( mar=rep(0,4),oma=c(4,4,.5,.5))
par(cex.axis = 1.5)
layout(matrix(c(1,1,1,2,1,1,1,2), 2, 4, byrow = TRUE))

Figure.FOIcontribution.detected(FOIProps.brazil$Prop.FOI.asym, FOIProps.brazil$Prop.FOI.inap, FOIProps.brazil$Prop.FOI.presym, FOIProps.brazil$Prop.FOI.postsym, 
                                FOIProps.brazil$Prop.FOI.inap.presym, FOIProps.brazil$Prop.FOI.inap.postsym)
lines(foi.vec, FOIProps.brazil$Prop.FOI.asym, lwd = 3)
lines(foi.vec, FOIProps.brazil$Prop.FOI.inap + FOIProps.brazil$Prop.FOI.asym, lwd = 3)
mtext('FoI',1,line=2.25)
mtext('Contribution to FoI',2,line=2.25)
yh <- hist(1-Prop.FOI.All.brazil, breaks = 40,plot = FALSE)
barplot((yh$density),space=0,horiz=T, axes = FALSE,col=rgb(1,0,0,.1),border = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

boxed.labels(-0.1827692, -0.6939367, 'Asymptomatic',bg="white",border=NA,cex = .9,xpad = 1.1, ypad = 2)
boxed.labels(-0.1827692, -0.4691403, 'Inapparent pre-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.1, ypad = 2)
boxed.labels(-0.1827692, 0.09285068, 'Inapparent post-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.1, ypad = 2)
boxed.labels(-0.1827692, 0.5561034, 'Undetected apparent pre-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.1, ypad = 1.5)
boxed.labels(-0.1827692, 0.7794036, 'Undetected apparent post-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.1, ypad = 2)
Lines <- list("Detected apparent","pre-symptomatic")
mtext(do.call(expression, Lines),3,line = c(-5.7,-6.7),at = 0.64, adj = 0, outer = FALSE, cex = 0.9)
Lines <- list("Detected apparent","post-symptomatic")
mtext(do.call(expression, Lines),3,line = c(-2.7,-3.7),at = 0.64, adj = 0, outer = FALSE, cex = 0.9)
arrows(x1 = 0.504, y1 = 0.950, x0 = 0.62, y0 = 0.8, length = 0.03, lwd = 2)
arrows(x1 = 0.504, y1 = 0.97, x0 = 0.62, y0 = 0.97, length = 0.03, lwd = 2 )

dev.off()
setwd(cd)


