# Figure functions infectiousness over time 


# Compare primary vs secondary infection including presymptomatics ----------------------------------

Figure.Viralload <- function(type,single='no',IIP='no') {
  if (is.double(IIP)){
    XLIM = c(-5,10)
  }
  if (IIP == 'no') { 
    IIP = rep(0,n) 
    XLIM = c(0,20)
  }
  if (type == 'sym') { Viremia.1 <- Viremia.mod1; Viremia.2 <- Viremia.mod2 }
  else if (type == 'asym') { Viremia.1 <- Viremia.mod1.asym; Viremia.2 <- Viremia.mod2.asym }
  else if (type == 'presym') { Viremia.1 <- Viremia.mod1.presym; Viremia.2 <- Viremia.mod2.presym }
  else if (type == 'severe') { Viremia.1 <- matrix(data=NA,nrow=dim(Viremia.mod1)[1], ncol=dim(Viremia.mod1)[2]); Viremia.2 <- Viremia.mod3 }
  else if (type == 'severe presym') { Viremia.1 <- matrix(data=NA,nrow=dim(Viremia.mod1)[1], ncol=dim(Viremia.mod1)[2]); Viremia.2 <- Viremia.mod3.presym }
  if (single == 'no' || single == 'prim'){
  plot(Time-IIP[1],Viremia.1[,1], ylim=c(0,12), xlim=XLIM, col = alpha('lightskyblue',0.01), type = 'l', lwd = 1, bty='o', ylab = '')
  for (ii in 2:dim(Viremia.1)[2]) {
    lines(Time-IIP[ii],Viremia.1[,ii], col = alpha('lightskyblue',0.01), type = 'l', lwd = 1, bty='n')
  }
  lines(Time-mean(IIP),rowMeans(Viremia.1),col = 'navyblue', type = 'l', lwd = 3, bty='n')
  }
  if (single == 'no'){
  lines(Time-IIP[1],Viremia.2[,1], ylim=c(0,12), xlim=XLIM, col = alpha('springgreen',0.01), type = 'l', lwd = 1, bty='n')
  for (ii in 2:dim(Viremia.1)[2]) {
    lines(Time-IIP[ii],Viremia.2[,ii], col = alpha('springgreen',0.01), type = 'l', lwd = 1, bty='n')
  }
  lines(Time-mean(IIP),rowMeans(Viremia.2),col = 'darkgreen', type = 'l', lwd = 3, bty='n')
  }
  if (single == 'sec'){
    plot(Time--IIP[1],Viremia.2[,1], ylim=c(0,12), xlim=XLIM, col = alpha('lightskyblue',0.01), type = 'l', lwd = 1, bty='o', ylab = '')
    for (ii in 2:dim(Viremia.1)[2]) {
      lines(Time-IIP[ii],Viremia.2[,ii], col = alpha('lightskyblue',0.01), type = 'l', lwd = 1, bty='n')
    }
    lines(Time-mean(IIP),rowMeans(Viremia.2),col = 'navyblue', type = 'l', lwd = 3, bty='n') 
    
  }
}

# Plot infectiousness as a function of viremia for symptomatics an --------

Figure.InfvsViremia <- function(type){
  if (type == 'sym') Reg = Sym.reg
  else if (type =='asym') Reg = Asym.reg
  else if (type =='presym') Reg = Presym.reg
  plot(Viremia,logistic.function.infectiousness(Reg[1,1],Reg[1,2],Viremia), ylim=c(0,1), type='l', lwd=2, col='lightskyblue',bty='n', ylab = '',xaxs='i',yaxs='i')
  for (ii in 2:length(Viremia)) {
    lines(Viremia,logistic.function.infectiousness(Reg[ii,1],Reg[ii,2],Viremia), col = alpha('lightskyblue',0.5), type = 'l', lwd = 1, bty='n')
  }
  lines(Viremia,logistic.function.infectiousness(colMeans(Reg)[1],colMeans(Reg)[2],Viremia),col = 'navyblue', type = 'l', lwd = 3, bty='n')
}

# Plot infectiousness as a function of time --------

Figure.Infectiousness.Unweighted <- function(type, dis.cat, gender='no',IIP='no'){
  if (is.double(IIP)){
    XLIM = c(-5,10)
  }
  if (IIP == 'no') { 
    IIP = rep(0,n) 
    XLIM = c(0,20)
  }
  if (type == 'sym' & dis.cat == 1 & gender == 'no') Prob <- Prob.sympto.mod1
  else if (type == 'sym' & dis.cat == 2 & gender == 'no') Prob <- Prob.sympto.mod2
  else if (type == 'asym' & dis.cat == 1 & gender == 'no') Prob <- Prob.asympto.mod1
  else if (type == 'asym' & dis.cat == 2 & gender == 'no') Prob <- Prob.asympto.mod2
  else if (type == 'presym' & dis.cat == 1 & gender == 'no') Prob <- Prob.presympto.mod1
  else if (type == 'presym' & dis.cat == 2 & gender == 'no') Prob <- Prob.presympto.mod2
  else if (type == 'sym' & dis.cat == 1 & gender == 'male') Prob <- Prob.sympto.mod1.male
  else if (type == 'sym' & dis.cat == 2 & gender == 'male') Prob <- Prob.sympto.mod2.male
  else if (type == 'asym' & dis.cat == 1 & gender == 'male') Prob <- Prob.asympto.mod1.male
  else if (type == 'asym' & dis.cat == 2 & gender == 'male') Prob <- Prob.asympto.mod2.male
  else if (type == 'presym' & dis.cat == 1 & gender == 'male') Prob <- Prob.presympto.mod1.male
  else if (type == 'presym' & dis.cat == 2 & gender == 'male') Prob <- Prob.presympto.mod2.male
  else if (type == 'sym' & dis.cat == 1 & gender == 'female') Prob <- Prob.sympto.mod1.female
  else if (type == 'sym' & dis.cat == 2 & gender == 'female') Prob <- Prob.sympto.mod2.female
  else if (type == 'asym' & dis.cat == 1 & gender == 'female') Prob <- Prob.asympto.mod1.female
  else if (type == 'asym' & dis.cat == 2 & gender == 'female') Prob <- Prob.asympto.mod2.female
  else if (type == 'presym' & dis.cat == 1 & gender == 'female') Prob <- Prob.presympto.mod1.female
  else if (type == 'presym' & dis.cat == 2 & gender == 'female') Prob <- Prob.presympto.mod2.female
  
  plot(Time - IIP[1],Prob[,1], ylim=c(0,1), xlim=XLIM, col = alpha('lightskyblue',0.01), type = 'l', lwd = 1, bty='o', ylab = '')
  for (ii in 2:n) {
    lines(Time- IIP[ii],Prob[,ii], col = alpha('lightskyblue',0.01), type = 'l', lwd = 1, bty='n')
  }
  lines(Time- mean(IIP),rowMeans(Prob),col = 'navyblue', type = 'l', lwd = 3, bty='n')
  
}

# Plot Probability Density Distribution of IIP ----------------------------
# As per Chan et al

Figure.IIP <- function(Prob.IIP = Prob.IIP) {
  plot(Time,Prob.IIP[,1], xlim=c(0,12),ylim=c(0,0.4), col = alpha('lightskyblue',0.5), type='l', lwd=1, bty='n', ylab = '')
  for (ii in 2:n) {
    lines(Time,Prob.IIP[,ii], col = alpha('lightskyblue',0.5), type = 'l', lwd = 1, bty='n')
  }
  lines(Time,rowMeans(Prob.IIP),col = 'navyblue', type = 'l', lwd = 3, bty='n')
}


# Plot Infectiousness vs Time: weighted average over different dur --------

Figure.Infectiousness.Weighted <- function(dis.cat, gender = 'no',IIP='no'){
  if (is.double(IIP)){
    XLIM = c(-5,10)
  }
  if (IIP == 'no') { 
    IIP = rep(0,n) 
    XLIM = c(0,20)
  }
  if (dis.cat == 1 & gender == 'no') Prob <- Prob.sympto.weighted.1
  else if ( dis.cat == 2 & gender == 'no') Prob <- Prob.sympto.weighted.2
  else if ( dis.cat == 3 & gender == 'no') Prob <- Prob.sympto.weighted.3
  else if (dis.cat == 1 & gender == 'male') Prob <- Prob.sympto.weighted.1.male
  else if ( dis.cat == 2 & gender == 'male') Prob <- Prob.sympto.weighted.2.male
  else if ( dis.cat == 3 & gender == 'male') Prob <- Prob.sympto.weighted.3.male
  else if (dis.cat == 1 & gender == 'female') Prob <- Prob.sympto.weighted.1.female
  else if ( dis.cat == 2 & gender == 'female') Prob <- Prob.sympto.weighted.2.female
  else if ( dis.cat == 3 & gender == 'female') Prob <- Prob.sympto.weighted.3.female
  
  plot(Time - IIP[1],Prob[,1], ylim=c(0,1), xlim=XLIM, col = alpha('lightskyblue',0.01), type = 'l', lwd = 1, bty='o', ylab = '')
  for (ii in 2:n) {
    lines(Time - IIP[ii],Prob[,ii], col = alpha('lightskyblue',0.01), type = 'l', lwd = 1, bty='n')
  }
  lines(Time - mean(IIP),rowMeans(Prob),col = 'navyblue', type = 'l', lwd = 3, bty='n')
}

# small code for error bars -----------------------------------------------

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

