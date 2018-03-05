library(fields)
library(akima)

# Get data ----------------------------------------------------------------
load('../Output/Ternary.RData')

# Set directory to save figures -------------------------------------------
path <- '../figures'
cd <- getwd()

# prepare data for plotting -----------------------------------------------
a<- df$Ia
b<- df$As
c<- df$Ap
d<- 1 - df$PROP.ap.all
e<- 100-c
f<- df$Res.inap.as
data.f<- data.frame (b, a, c)

tern2cart <- function(coord){
  coord[1]->x
  coord[2]->y
  coord[3]->z
  x+y+z -> tot
  x/tot -> x  
  y/tot -> y
  z/tot -> z
  (2*y + z)/(2*(x+y+z)) -> x1 
  sqrt(3)*z/(2*(x+y+z)) -> y1
  return(c(x1,y1))
}


t(apply(data.f,1,tern2cart)) -> tern


resolution <- 0.0005
interp(tern[,1],tern[,2],z=d*100, xo=seq(0,1,by=resolution), yo=seq(0,1,by=resolution)) -> tern.grid
interp(tern[,1],tern[,2],z=e, xo=seq(0,1,by=resolution), yo=seq(0,1,by=resolution)) -> tern.grid.exp
interp(tern[,1],tern[,2],z=f*100, xo=seq(0,1,by=resolution), yo=seq(0,1,by=resolution)) -> tern.grid.res


# Plot A:IS:CS ratios --------------------------------------------------------------------
setwd(path)
pdf('Fig5.rev.pdf',width=20,height=20)
par(mfrow = c(1,1), mar=c(0,0.5,1,0.5),oma=c(4,4,.5,.5))

plot(NA,NA,xlim=c(0,1),ylim=c(0,sqrt(3)/2),asp=1,bty="n",axes=F,xlab="",ylab="")

Lines <- list(bquote(underline("(As,IS,AS)")),bquote("As = asymptomatic"), bquote("IS = inapparent symptomatic"), bquote("AS = apparent symptomatic"))
mtext(do.call(expression, Lines),side=3,line= seq(-11,-18.5,by=-2.5),at = .18,.9, cex=3)

segments(0,0,0.5,sqrt(3)/2)
segments(0.5,sqrt(3)/2,1,0)
segments(1,0,0,0)

image(tern.grid,breaks=seq(0,100,by=1),col=rev(tim.colors(100)),add=T) 
contour(tern.grid,add=T,levels=seq(0,100,by=20),labcex=3,col = 'black')
contour(tern.grid.exp,add=T,levels=seq(0,100,by=20),labcex=3,col = 'white',lty=2,lwd = 3, drawlabels = F)

text(0.5,(sqrt(3)/2),"(0,0,1)", pos=3,cex=3)
text(0,0,"(1,0,0)", pos=1,cex=3)
text(1,0,"(0,1,0)", pos=1,cex=3)

coor <- tern2cart(c(9,64,27))
points(coor[1],coor[2],pch=19,  cex = 2)
text(coor[1],coor[2],'(0.09,0.64,0.27)',pos = 1, cex = 3)
coor <- tern2cart(c(100/3,100/3,100/3))
points(coor[1],coor[2],pch=19,  cex = 2)
text(coor[1],coor[2],'(0.33,0.33,0.33)',pos = 1, cex = 3)
coor <- tern2cart(c(50,0,50))
text(coor[1],coor[2],'(0.5,0,0.5)',pos = 4,  cex = 3)
coor <- tern2cart(c(0,50,50))
text(coor[1],coor[2],'(0,0.5,0.5)', pos = 2,  cex = 3)
coor <- tern2cart(c(50,50,0))
text(coor[1],coor[2],'(0.5,0.5,0)', pos = 1,  cex = 3)

dev.off()
setwd(cd)


# color bar ---------------------------------------------------------------
setwd(path)
pdf('Fig5.colorbar.pdf',width=7.5,height=7.5)
par(mfrow = c(1,1), mar=c(0,0.5,1,0.5),oma=c(4,4,.5,4))

clrs <- rev(tim.colors(100))
clrs <- clrs[1:100]
ticks<- seq(0,100,by=20)
plot.new()
image.plot(tern.grid,breaks=seq(0,100,by=1),col=clrs,legend.only = TRUE,nlevel=24, legend.width = 2, 
           axis.args=list( at=(ticks), labels=ticks),
           legend.args=list(text='Contribution to overall FoI (%) ', side=4, font=2, line=3, cex=1.3) )

dev.off()
setwd(cd)


