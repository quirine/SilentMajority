# script to create figures for the supplementary

# Load packages -----------------------------------------------------------
library(scales)
library(vioplot)
library(fields)
library(plotrix)

# Load source files -------------------------------------------------------
source('ProportionImmune.Functions.R')
source('ProportionImmune.Figures.R')
source('InfectiousnessFigures.Functions.R')
source('FOI.Figures.R')

# Set directory to save figures -------------------------------------------
path <- '../figures'
cd <- getwd()

# SI Figure 2a: pre-exposure history Brazil---------------------------------------

setwd(path)
pdf('SFig2a.pdf',width=4.5,height=4.5)
par(mfrow = c(1,1), mar=rep(0,4),oma=c(4,4,1,1))

infect.hist.figure(foi.vec, infec.hist.brazil)
mtext( bquote(bold("a")), side = 0, line = 0 ,at =0.015,cex=1) 
mtext('Force of infection',1,line=2.25)
mtext('Proportion of susceptible population ',2,line=2.5)

dev.off()
setwd(cd)


# SI Figure 2b idem Thailand------------------------------------------------------------

setwd(path)
pdf('SFig2b.pdf',width=4.5,height=4.5)
par(mfrow = c(1,1), mar=rep(0,4),oma=c(4,4,1,1))

infect.hist.figure(foi.vec, infec.hist.thailand)
mtext( bquote(bold("b")), side = 3, line = 0 ,at =0,cex=1) 
mtext('Force of infection',1,line=2.25)
legend('topright',title = 'susceptible to', c('1 serotype','2 serotypes','3 serotypes','4 serotypes'), fill = c(rgb(0,1,0,c(.75,.5,.25,0))),bty='o', bg = 'white')

dev.off()
setwd(cd)

# SI Figure 2 a and b: pre-exposure history Brazil---------------------------------------

setwd(path)
pdf('SFig2.pdf',width=7.5,height=4.5)
par(mfrow = c(1,2), mar=rep(0,4),oma=c(4,4,1,1))

par(yaxt='s')
infect.hist.figure(foi.vec, infec.hist.brazil)
mtext( bquote(bold("a")), side = 3, line = 0 ,at =0,cex=1) 
mtext('Force of infection',1,line=2.25)
mtext('Proportion of susceptible population ',2,line=2.5)
par(yaxt='n')
infect.hist.figure(foi.vec, infec.hist.thailand)
mtext( bquote(bold("b")), side = 3, line = 0 ,at =0,cex=1) 
mtext('Force of infection',1,line=2.25)
legend('topright',title = 'susceptible to', c('1 serotype','2 serotypes','3 serotypes','4 serotypes'), fill = c(rgb(0,1,0,c(.75,.5,.25,0))),bty='o', bg = 'white')

dev.off()
setwd(cd)

# SI Figure 3: FOI contributions Thailand ----------------------------------------

setwd(path)
pdf('SFig3.pdf',width=7.5,height=4.5)
par( mar=rep(0,4),oma=c(4,4,.5,.5))
layout(matrix(c(1,1,1,2,1,1,1,2), 2, 4, byrow = TRUE))

Figure.FOIcontribution.detected(FOIProps.thailand$Prop.FOI.asym, FOIProps.thailand$Prop.FOI.inap, FOIProps.thailand$Prop.FOI.presym, FOIProps.thailand$Prop.FOI.postsym,
                       FOIProps.thailand$Prop.FOI.inap.presym, FOIProps.thailand$Prop.FOI.inap.postsym)
lines(foi.vec, FOIProps.thailand$Prop.FOI.asym, lwd = 3)
lines(foi.vec, FOIProps.thailand$Prop.FOI.inap + FOIProps.thailand$Prop.FOI.asym, lwd = 3)
mtext('FoI',1,line=2.25)
mtext('Contribution to FoI',2,line=2.25)
yh <- hist(1-Prop.FOI.All.thailand, breaks = 40,plot = FALSE)
barplot((yh$density),space=0,horiz=T, axes = FALSE,col=rgb(1,0,0,.1),border = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
boxed.labels(-0.1827692, -0.6939367, 'Asymptomatic',bg="white",border=NA,cex = .9)
boxed.labels(-0.1827692, -0.4691403, 'Inapparent pre-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.2, ypad = 2)
boxed.labels(-0.1827692, 0.09285068, 'Inapparent post-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.2, ypad = 2)
boxed.labels(-0.1827692, 0.5766516, 'Undetected apparent pre-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.2, ypad = 1.5)
boxed.labels(-0.1827692, 0.8258824, 'Undetected apparent post-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.2, ypad = 2)
Lines <- list("Detected apparent","pre-symptomatic")
mtext(do.call(expression, Lines),3,line = c(-5.7,-6.7),at = 0.64, adj = 0, outer = FALSE, cex = 0.9)
Lines <- list("Detected apparent","post-symptomatic")
mtext(do.call(expression, Lines),3,line = c(-2.7,-3.7),at = 0.64, adj = 0, outer = FALSE, cex = 0.9)
arrows(x1 = 0.504, y1 = 0.950, x0 = 0.62, y0 = 0.8, length = 0.03, lwd = 2)
arrows(x1 = 0.504, y1 = 0.97, x0 = 0.62, y0 = 0.97, length = 0.03, lwd = 2 )

dev.off()
setwd(cd)

# SI Figure 5: FOI contribution in emerging setting ------------------------------

setwd(path)
pdf('SFig5.pdf',width=7.5,height=4.5)
par( mar=rep(0,4),oma=c(4,4,.5,.5))
layout(matrix(c(1,1,1,2,1,1,1,2), 2, 4, byrow = TRUE))

Figure.FOIcontribution.detected(FOIProps.brazil.emerging$Prop.FOI.asym, FOIProps.brazil.emerging$Prop.FOI.inap, FOIProps.brazil.emerging$Prop.FOI.presym, FOIProps.brazil.emerging$Prop.FOI.postsym,
                       FOIProps.brazil.emerging$Prop.FOI.inap.presym, FOIProps.brazil.emerging$Prop.FOI.inap.postsym)
lines(foi.vec, FOIProps.brazil.emerging$Prop.FOI.asym, lwd = 3)
lines(foi.vec, FOIProps.brazil.emerging$Prop.FOI.inap + FOIProps.brazil.emerging$Prop.FOI.asym, lwd = 3)
mtext('FoI',1,line=2.25)
mtext('Contribution to FoI',2,line=2.25)
yh <- hist(1-Prop.FOI.All.brazil.emerging, breaks = 40,plot = FALSE)
barplot((yh$density),space=0,horiz=T, axes = FALSE,col=rgb(1,0,0,.1),border = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

boxed.labels(-0.1827692, -0.6908508, 'Asymptomatic',bg="white",border=NA,cex = .9)
boxed.labels(-0.1827692, -0.4691403, 'Inapparent pre-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.2, ypad = 2)
boxed.labels(-0.1827692, 0.09285068, 'Inapparent post-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.2, ypad = 2)
boxed.labels(-0.1827692,  0.5582238, 'Undetected apparent pre-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.2, ypad = 1.5)
boxed.labels(-0.1827692, 0.8258824, 'Detected apparent post-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.2, ypad = 2)
Lines <- list("Detected apparent","pre-symptomatic")
mtext(do.call(expression, Lines),3,line = c(-5.7,-6.7),at = 0.64, adj = 0, outer = FALSE, cex = 0.9)
Lines <- list("Detected apparent","post-symptomatic")
mtext(do.call(expression, Lines),3,line = c(-2.7,-3.7),at = 0.64, adj = 0, outer = FALSE, cex = 0.9)
arrows(x1 = 0.504, y1 = 0.950, x0 = 0.62, y0 = 0.8, length = 0.03, lwd = 2)
arrows(x1 = 0.504, y1 = 0.97, x0 = 0.62, y0 = 0.97, length = 0.03, lwd = 2 )


dev.off()
setwd(cd)

# SI Figure 1: 4-infections FOI contribution ------------------------------

setwd(path)
pdf('SFig1.pdf',width=7.5,height=4.5)
par( mar=rep(0,4),oma=c(4,4,.5,.5))
layout(matrix(c(1,1,1,2,1,1,1,2), 2, 4, byrow = TRUE))

Figure.FOIcontribution.detected(FOIProps.brazil.4$Prop.FOI.asym, FOIProps.brazil.4$Prop.FOI.inap, FOIProps.brazil.4$Prop.FOI.presym, FOIProps.brazil.4$Prop.FOI.postsym,
                       FOIProps.brazil.4$Prop.FOI.inap.presym, FOIProps.brazil.4$Prop.FOI.inap.postsym)
lines(foi.vec, FOIProps.brazil.4$Prop.FOI.asym, lwd = 3)
lines(foi.vec, FOIProps.brazil.4$Prop.FOI.inap + FOIProps.brazil.4$Prop.FOI.asym, lwd = 3)
mtext('FoI',1,line=2.25)
mtext('Contribution to FoI',2,line=2.25)
yh <- hist(1-Prop.FOI.All.brazil.4, breaks = 40,plot = FALSE)
barplot((yh$density),space=0,horiz=T, axes = FALSE,col=rgb(1,0,0,.1),border = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
boxed.labels(-0.1827692, -0.6808703, 'Asymptomatic',bg="white",border=NA,cex = .9)
boxed.labels(-0.1827692, -0.4691403, 'Inapparent pre-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.2, ypad = 2)
boxed.labels(-0.1827692, 0.3126821, 'Inapparent post-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.2, ypad = 2)
boxed.labels(-0.1827692, 0.7692362, 'Undetected apparent pre-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.1, ypad = 1.2)
boxed.labels(-0.1827692, 0.8926998, 'Undetected apparent post-symptomatic',bg="white",border=NA,cex = .9,xpad = 1.2, ypad = 1.5)
Lines <- list("Detected apparent","pre-symptomatic")
mtext(do.call(expression, Lines),3,line = c(-5.7,-6.7),at = 0.64, adj = 0, outer = FALSE, cex = 0.9)
Lines <- list("Detected apparent","post-symptomatic")
mtext(do.call(expression, Lines),3,line = c(-2.7,-3.7),at = 0.64, adj = 0, outer = FALSE, cex = 0.9)
arrows(x1 = 0.504, y1 = 0.950, x0 = 0.62, y0 = 0.8, length = 0.03, lwd = 2)
arrows(x1 = 0.504, y1 = 0.97, x0 = 0.62, y0 = 0.97, length = 0.03, lwd = 2 )

dev.off()
setwd(cd)

