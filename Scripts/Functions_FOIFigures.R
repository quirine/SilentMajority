

Figure.FOIcontribution.single = function(Prop.FOI.all) {
  plot(foi.vec,Prop.FOI.all[,1],ylim=c(0,1), col = alpha('lightskyblue',0.05), type = 'l', lwd = 1, bty='n')
  for (ii in 2:dim(Prop.FOI.all)[2]){
    lines(foi.vec,Prop.FOI.all[,ii],col = alpha('lightskyblue',0.05), type = 'l', lwd = 1)
  }
  lines(foi.vec,rowMeans(Prop.FOI.all,na.rm=TRUE),col = 'navyblue', type = 'l', lwd = 3)
}

Figure.FOIcontribution = function(Prop.FOI.asym, Prop.FOI.inap, Prop.FOI.presym, Prop.FOI.postsym, Prop.FOI.inap.presym, Prop.FOI.inap.postsym) {
  plot(NULL,xlim=range(foi.vec),ylim=c(0,1))  
  polygon(
    c(foi.vec,tail(foi.vec,1),foi.vec[1]),
    c(Prop.FOI.asym,0,0),
    col=rgb(1,0,0,.1))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym,rev(Prop.FOI.asym + Prop.FOI.inap.presym)),
    col=rgb(1,0,0,.2))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap.presym,rev(Prop.FOI.asym + Prop.FOI.inap.presym + Prop.FOI.inap.postsym)),
    col=rgb(1,0,0,.2))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap.presym,rev(Prop.FOI.asym + Prop.FOI.inap.presym + Prop.FOI.inap.postsym)),
    col='black',density = 25)
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap,rev(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym)),
    col=rgb(1,0,0,.5))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym,rev(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym + Prop.FOI.postsym)),
    col=rgb(1,0,0,.5))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym,rev(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym + Prop.FOI.postsym)),
    col='black',density = 25)
}

Figure.FOIcontribution.detected = function(Prop.FOI.asym, Prop.FOI.inap, Prop.FOI.presym, Prop.FOI.postsym, Prop.FOI.inap.presym, Prop.FOI.inap.postsym) {
  plot(NULL,xlim=range(foi.vec),ylim=c(0,1))  
  polygon(
    c(foi.vec,tail(foi.vec,1),foi.vec[1]),
    c(Prop.FOI.asym,0,0),
    col=rgb(1,0,0,.1))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym,rev(Prop.FOI.asym + Prop.FOI.inap.presym)),
    col=rgb(1,0,0,.2))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap.presym,rev(Prop.FOI.asym + Prop.FOI.inap.presym + Prop.FOI.inap.postsym)),
    col=rgb(1,0,0,.2))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap.presym,rev(Prop.FOI.asym + Prop.FOI.inap.presym + Prop.FOI.inap.postsym)),
    col='black',density = 25)
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap,rev(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym*0.92)),
    col=rgb(1,0,0,.5))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym*0.92,rev(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym*0.92 + Prop.FOI.postsym*0.92)),
    col=rgb(1,0,0,.5))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym*0.92,rev(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym*0.92 + Prop.FOI.postsym*0.92)),
    col='black',density = 25)
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym*0.92 + Prop.FOI.postsym*0.92,rev(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym*0.92 + Prop.FOI.postsym*0.92 + Prop.FOI.presym*0.08)),
    col=rgb(1,0,0,1))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym*0.92 + Prop.FOI.postsym*0.92 + Prop.FOI.presym*0.08,rev(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym*0.92 + Prop.FOI.postsym*0.92 + Prop.FOI.presym*0.08 + Prop.FOI.postsym*0.08)),
    col=rgb(1,0,0,1))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym*0.92 + Prop.FOI.postsym*0.92 + Prop.FOI.presym*0.08,rev(Prop.FOI.asym + Prop.FOI.inap + Prop.FOI.presym*0.92 + Prop.FOI.postsym*0.92 + Prop.FOI.presym*0.08 + Prop.FOI.postsym*0.08)),
    col='black',density = 25)
}

