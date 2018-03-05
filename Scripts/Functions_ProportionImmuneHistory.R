

prop.dis.calculator <- function(infec.hist,type) {
  if (type == 'asym') ind = 1
  else if (type == 'mild') ind = 2
  else if (type == 'severe') ind = 3
  X =
    infec.hist[,1] * Prim.dis[ind] +
    infec.hist[,2] * Sec.dis[ind] +
    infec.hist[,3] * Tert.dis[ind] +
    infec.hist[,4] * Tert.dis[ind]
  
  return(X)
}

infec.hist.calculator <- function (foi.vec, agepop, gender = 'no', totalpop = 'NaN', 
                                   endemic = 'yes', inv.cross.im = 'NA', ode = 'yes'){
# browser()
  infec.hist = matrix(0,length(foi.vec),4)
prev.vec4 <- prev.vec2 <- prev.vec1 <- numeric()
for(ff in 1:length(foi.vec)){
  foi = foi.vec[ff]
  seroage = sero.age.calculator(agepop, foi, inv.cross.im, ode)
  if (endemic == 'no') {
    seroage = matrix(data = 0, nrow(agepop),5)
    seroage[,1] = rep(1,nrow(agepop))
  }
  
  if (gender == 'no') {
  seroprev = rep(0,5)
  seroprev[1] = sum(seroage[,1] * agepop$pop)
  seroprev[2] = sum(seroage[,2] * agepop$pop)
  seroprev[3] = sum(seroage[,3] * agepop$pop)
  seroprev[4] = sum(seroage[,4] * agepop$pop)
  seroprev[5] = sum(seroage[,5] * agepop$pop)
  seroprev = seroprev / sum(seroprev)
  
  infec.hist[ff,] = seroprev[1:4] / sum(seroprev[1:4])
  if (ode == 'no') {
  prev.vec4 = c(prev.vec4, sum((1 - dpois(0, seq(4,0,-1) * foi * inf.period)) * colSums(seroage * agepop$pop)) / sum(agepop$pop))
  prev.vec2 = c(prev.vec2, sum((1 - dpois(0, seq(4,3,-1) * foi * inf.period)) * colSums(seroage[,1:2] * agepop$pop)) / sum(agepop$pop))
  prev.vec1 = c(prev.vec1, sum((1 - dpois(0, 4 * foi * inf.period)) * (seroage[,1] * agepop$pop)) / sum(agepop$pop))
  }
  else {
    seroprob <- ode(c(1,rep(0,8)),agepop$age,proportion.exposed.model.forprev,parms=c(foi,inv.cross.im))
    if (inv.cross.im == 'NA'){  # DOESN'T WORK YET!!! NON-NUMERIC ARGUMENT IN NO-CI MODEL
      prev.vec4 = c(prev.vec4, sum((seroprob[,c(3, 4, 5, 6)]) * (agepop$pop)) / sum(agepop$pop) * ((inf.period)/(1/seq(4,1,-1) * foi)) )
      prev.vec2 = c(prev.vec2, sum((seroprob[,c(3, 4)]) * (agepop$pop)) / sum(agepop$pop) * ((inf.period)/(1/seq(4,3,-1) * foi)) )
      prev.vec1 = c(prev.vec1, sum((seroprob[,c(3)]) * (agepop$pop)) / sum(agepop$pop) * ((inf.period)/(1/4 * foi)) )
      
    }
    else {
      prev.vec4 = c(prev.vec4, sum((seroprob[,c(3, 5, 7, 9)]) * (agepop$pop)) / sum(agepop$pop) * ((inf.period)/(1/inv.cross.im)) )
      prev.vec2 = c(prev.vec2, sum((seroprob[,c(3, 5)]) * (agepop$pop)) / sum(agepop$pop) * ((inf.period)/(1/inv.cross.im)) )
      prev.vec1 = c(prev.vec1, sum((seroprob[,c(3)]) * (agepop$pop)) / sum(agepop$pop) * ((inf.period)/(1/inv.cross.im)) )
      
      
    } 
      
  }
 
  }
  else if (gender == 'male' | gender == 'female') {
    seroprev <- serotot <- rep(0,5)
    seroprev[1] = sum(seroage[,1] * agepop$pop  * (agepop$pop/totalpop))
    seroprev[2] = sum(seroage[,2] * agepop$pop  * (agepop$pop/totalpop))
    seroprev[3] = sum(seroage[,3] * agepop$pop  * (agepop$pop/totalpop))
    seroprev[4] = sum(seroage[,4] * agepop$pop  * (agepop$pop/totalpop))
    seroprev[5] = sum(seroage[,5] * agepop$pop  * (agepop$pop/totalpop))
    serotot[1] = sum(seroage[,1] * totalpop)
    serotot[2] = sum(seroage[,2] * totalpop)
    serotot[3] = sum(seroage[,3] * totalpop)
    serotot[4] = sum(seroage[,4] * totalpop)
    serotot[5] = sum(seroage[,5] * totalpop)
    seroprev = seroprev / sum(serotot) 
    
    infec.hist[ff,] = seroprev[1:4] #/ sum(serotot[1:4] ) 
    prev.vec4 = c(prev.vec4, sum((1 - dpois(0, seq(4,0,-1) * foi * inf.period)) * colSums(seroage * totalpop)) / sum(totalpop))
    prev.vec2 = c(prev.vec2, sum((1 - dpois(0, seq(4,3,-1) * foi * inf.period)) * colSums(seroage[,1:2] * totalpop)) / sum(totalpop))
    prev.vec1 = c(prev.vec1, sum((1 - dpois(0, 4 * foi * inf.period)) * (seroage[,1] * totalpop)) / sum(totalpop))
  }
  
}
return(data.frame(infec.hist,prev.vec4,prev.vec2,prev.vec1))
}

dis.per.age.calculator <- function(foi, agepop, type, inv.cross.im = 'NA', ode = 'yes' ){
  seroage = sero.age.calculator(agepop,foi,inv.cross.im,ode)
  if (type == 'inapparent') {
    A =
      seroage[,1] * agepop$pop * inapp.ratio.rev[1] +
      seroage[,2] * agepop$pop * inapp.ratio.rev[2] +
      seroage[,3] * agepop$pop * inapp.ratio.rev[3] +
      seroage[,4] * agepop$pop * inapp.ratio.rev[3]  
  }
  else {
    if (type == 'asym') ind = 1
    else if (type == 'mild') ind = 2
    else if (type == 'severe') ind = 3
    A =
      seroage[,1] * agepop$pop * Prim.dis[ind] +
      seroage[,2] * agepop$pop * Sec.dis[ind] +
      seroage[,3] * agepop$pop * Tert.dis[ind] +
      seroage[,4] * agepop$pop * Tert.dis[ind]
  }
  
  return(A)
}

sero.age.calculator <- function(agepop, foi, inv.cross.im = 'NA', ode = 'yes'){
  if (ode == 'yes'){
    seroage = matrix(0,nrow(agepop),5)
    if(inv.cross.im == 'NA'){
      seroprob <- ode(c(1,rep(0,4)),agepop$age,proportion.exposed.model.noci,parms=c(foi))
      seroage[,1] = (seroprob[,2])  # prim
      seroage[,2] = (seroprob[,3])  # sec
      seroage[,3] = (seroprob[,4])  # tert
      seroage[,4] = (seroprob[,5])  # quart
      seroage[,5] = (seroprob[,6])  # immmune
    }
    else {
      seroprob <- ode(c(1,rep(0,7)),agepop$age,proportion.exposed.model,parms=c(foi,inv.cross.im))
      seroage[,1] = (seroprob[,2] ) #
      seroage[,2] = (seroprob[,4] ) #+ seroprob[,3]) 
      seroage[,3] = (seroprob[,6] ) # + seroprob[,5]) 
      seroage[,4] = (seroprob[,8] ) #+ seroprob[,7]) 
      seroage[,5] = (seroprob[,9] + seroprob[,3] + seroprob[,5] + seroprob[,7]) 
    }
  }
    else if (ode == 'no' ){  # original function with poisson draws
  seroage = matrix(0,nrow(agepop),5)
  seroage[,1] = dpois(0,foi*agepop$age)^4
  seroage[,2] = 4 * (1 - dpois(0,foi*agepop$age)) * dpois(0,foi*agepop$age)^3
  seroage[,3] = choose(4,2) * (1 - dpois(0,foi*agepop$age))^2 * dpois(0,foi*agepop$age)^2 
  seroage[,4] = choose(4,3) * (1 - dpois(0,foi*agepop$age))^3 * dpois(0,foi*agepop$age) 
  seroage[,5] = (1 - dpois(0,foi*agepop$age))^4
  }
  return(seroage)
}

immune.per.age.calculator = function(agepop, foi, inv.cross.im){
  immuneage = matrix(0,nrow(agepop),4)
  seroprob <- ode(c(1,rep(0,7)),agepop$age,proportion.exposed.model,parms=c(foi,inv.cross.im))
  immuneage[,1] = seroprob[,3]
  immuneage[,2] = seroprob[,5]
  immuneage[,3] = seroprob[,7]
  immuneage[,4] = seroprob[,9]
  return(immuneage)
}

# nested models -----------------------------------------------------------
proportion.exposed.model <- function (t, x, params) {
  
  Zero <- x[1];  OneCI <- x[2];  One <- x[3];   TwoCI <- x[4];     Two <- x[5];    ThreeCI <- x[6];  Three <- x[7];   Four <- x[8];   
  FOI <- params[1];      inv.cross.im <- params[2];  
  
  # differential equations
  dZero    <-              - 4*FOI*Zero
  dOneCI   <- 4*FOI*Zero   - inv.cross.im*OneCI
  dOne     <- inv.cross.im*OneCI    - 3*FOI*One
  dTwoCI   <- 3*FOI*One    - inv.cross.im*TwoCI
  dTwo     <- inv.cross.im*TwoCI    - 2*FOI*Two
  dThreeCI <- 2*FOI*Two    - inv.cross.im*ThreeCI
  dThree   <- inv.cross.im*ThreeCI  - 1*FOI*Three
  dFour    <- 1*FOI*Three
  
  list(c(dZero,dOneCI,dOne,dTwoCI,dTwo,dThreeCI,dThree,dFour))
}

proportion.exposed.model.forprev <- function (t, x, params) {
  
  Zero <- x[1];  OneCI <- x[2];  One <- x[3];   TwoCI <- x[4];     Two <- x[5];    ThreeCI <- x[6];  Three <- x[7];   FourCI <- x[8];     Four <-x[9];   
  FOI <- params[1];      inv.cross.im <- params[2];  
  
  # differential equations
  dZero    <-              - 4*FOI*Zero
  dOneCI   <- 4*FOI*Zero   - inv.cross.im*OneCI
  dOne     <- inv.cross.im*OneCI    - 3*FOI*One
  dTwoCI   <- 3*FOI*One    - inv.cross.im*TwoCI
  dTwo     <- inv.cross.im*TwoCI    - 2*FOI*Two
  dThreeCI <- 2*FOI*Two    - inv.cross.im*ThreeCI
  dThree   <- inv.cross.im*ThreeCI  - 1*FOI*Three
  dFourCI  <- 1*FOI*Three - inv.cross.im*FourCI
  dFour    <- inv.cross.im*FourCI
  
  list(c(dZero,dOneCI,dOne,dTwoCI,dTwo,dThreeCI,dThree,dFourCI,dFour))
}

proportion.exposed.model.noci <- function (t, x, params) {
  
  Zero <- x[1];  One <- x[2];   Two <- x[3];    Three <- x[4];   Four <- x[5];   
  FOI <- params[1];      
  
  # differential equations
  dZero    <-              - 4*FOI*Zero
  dOne     <- 4*FOI*Zero   - 3*FOI*One
  dTwo     <- 3*FOI*One    - 2*FOI*Two
  dThree   <- 2*FOI*Two    - 1*FOI*Three
  dFour    <- 1*FOI*Three
  
  list(c(dZero,dOne,dTwo,dThree,dFour))
}


