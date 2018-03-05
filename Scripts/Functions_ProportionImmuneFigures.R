

infect.hist.figure <- function(foi.vec, infec.hist){
  plot(NULL,xlim=range(foi.vec),ylim=c(0,1),xaxs='i',yaxs='i',las=1,
       xlab='',
       ylab='')
  polygon(
    c(foi.vec,tail(foi.vec,1),foi.vec[1]),
    c(infec.hist[,1],0,0),
    col=rgb(0,1,0,0))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(infec.hist[,1],rev(rowSums(infec.hist[,1:2]))),
    col=rgb(0,1,0,.25))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(rowSums(infec.hist[,1:2]),rev(rowSums(infec.hist[,1:3]))),
    col=rgb(0,1,0,.5))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(rowSums(infec.hist[,1:3]),rev(rowSums(infec.hist[,1:4]))),
    col=rgb(0,1,0,.75))
}

sero.by.age.figure <- function(seroage, agepop){
  barplot(
    t(seroage * agepop$pop),
    col = c(rgb(0,1,0,c(.1,.4,.7,1)),'black'),border=NA)
  lines(agepop$age*1.2, (seroage * agepop$pop)[,1])
  lines(agepop$age*1.2, rowSums((seroage * agepop$pop)[,1:2]))
  lines(agepop$age*1.2, rowSums((seroage * agepop$pop)[,1:3]))
  lines(agepop$age*1.2, rowSums((seroage * agepop$pop)[,1:4]))
}

dis.prop.figure <- function(foi.vec, A, M, S){
  plot(NULL,xlim=range(foi.vec),ylim=c(0,1),xaxs='i',yaxs='i',las=1,
       xlab='',
       ylab='')
  polygon(
    c(foi.vec,tail(foi.vec,1),foi.vec[1]),
    c(A,0,0),
    col=rgb(1,0,0,.1))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(A,rev(A+M)),
    col=rgb(1,0,0,.55))
  polygon(
    c(foi.vec,rev(foi.vec)),
    c(A+M,rev(A+M+S)),
    col=rgb(1,0,0,1))
}

dis.by.age.figure <- function(A.age, I.age, M.age, S.age, agepop, seroage, Immune.age){
  plot(NULL,xlim=c(0,100),ylim=c(0,.02),xaxs='i',yaxs='i',las=1,
       xlab='',
       ylab='', axes=FALSE)
  polygon(
    c(agepop$age,tail(agepop$age,1),agepop$age[1]),
    c(A.age,0,0),
    col=rgb(1,0,0,.1))
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(A.age,rev(A.age+(I.age-A.age))),
    col=rgb(1,0,0,.2))
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(A.age+(I.age-A.age),rev(A.age+(I.age-A.age) + ((M.age+S.age)-(I.age-A.age)))),
    col=rgb(1,0,0,.5))
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(A.age+(I.age-A.age) + ((M.age+S.age)-(I.age-A.age)),rev(A.age+(I.age-A.age) + ((M.age+S.age)-(I.age-A.age)) + rowSums(Immune.age[,c(1,2,3)])*agepop$pop)),
    col = 'gray')
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(A.age+(I.age-A.age) + ((M.age+S.age)-(I.age-A.age)+ rowSums(Immune.age[,c(1,2,3)])*agepop$pop),rev(A.age+(I.age-A.age) + ((M.age+S.age)-(I.age-A.age)) + rowSums(Immune.age[,c(1,2,3,4)])*agepop$pop)),
    col = 'black')
  
}


dis.by.age.figure.detected <- function(A.age, I.age, M.age, S.age, agepop, seroage, Immune.age, detection.rate=0.08){
  # browser()
  plot(NULL,xlim=c(0,100),ylim=c(0,.02),xaxs='i',yaxs='i',las=1,
       xlab='',
       ylab='', axes=FALSE)
  polygon(
    c(agepop$age,tail(agepop$age,1),agepop$age[1]),
    c(A.age,0,0),
    col=rgb(1,0,0,.1))
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(A.age,rev(A.age+(I.age-A.age))),
    col=rgb(1,0,0,.2))
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(A.age+(I.age-A.age),rev(A.age+(I.age-A.age) + ((M.age*(1-detection.rate)+S.age*(1-detection.rate))-(I.age-A.age)))),
    col=rgb(1,0,0,.5))
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(A.age+(I.age-A.age) + (M.age*(1-detection.rate)+S.age*(1-detection.rate))-(I.age-A.age),rev(A.age+(I.age-A.age) + (M.age*(1-detection.rate)+S.age*(1-detection.rate)) -(I.age-A.age) + (M.age*detection.rate+S.age*detection.rate))),
    col=rgb(1,0,0,1))
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(A.age+(I.age-A.age) + ((M.age+S.age)-(I.age-A.age)),rev(A.age+(I.age-A.age) + ((M.age+S.age)-(I.age-A.age)) + rowSums(Immune.age[,c(1,2,3)])*agepop$pop)),
    col = 'gray')
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(A.age+(I.age-A.age) + ((M.age+S.age)-(I.age-A.age)+ rowSums(Immune.age[,c(1,2,3)])*agepop$pop),rev(A.age+(I.age-A.age) + ((M.age+S.age)-(I.age-A.age)) + rowSums(Immune.age[,c(1,2,3,4)])*agepop$pop)),
    col = 'black')
  
}

infec.hist.by.age.figure <- function(agepop, seroprev){
  plot(NULL,xlim=c(0,100),ylim=c(0,0.02),xaxs='i',yaxs='i',las=1,
       xlab='',
       ylab='', axes=FALSE)
  polygon(
    c(agepop$age,tail(agepop$age,1),agepop$age[1]),
    c(seroprev[,1],0,0),
    col=rgb(0,1,0,0))
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(seroprev[,1],rev(rowSums(seroprev[,1:2]))),
    col=rgb(0,1,0,.25))
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(rowSums(seroprev[,1:2]),rev(rowSums(seroprev[,1:3]))),
    col=rgb(0,1,0,.5))
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(rowSums(seroprev[,1:3]),rev(rowSums(seroprev[,1:4]))),
    col=rgb(0,1,0,.75))
  polygon(
    c(agepop$age,rev(agepop$age)),
    c(rowSums(seroprev[,1:4]),rev(rowSums(seroprev[,1:5]))),
    col='black')
}
