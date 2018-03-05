# Code to create Table S2: proportion of run higher or lower than a specific reference

# Load packages -----------------------------------------------------------
library(zoo)

# Load Workspace ----------------------------------------------------------
load('Workspace.lowerAsym.RData')

# Using samples ---------------------------------------------

AUCs <- cbind(AUC.asympto.1, AUC.asympto.2, AUC.sympto.w.1, AUC.sympto.w.2)

IntProbs <- IntProbs.20 <- IntProbs.200 <- matrix(nrow=4,ncol=4)
for (kk in 1:4) {
  Reference <- AUCs[,kk]
  for (jj in 1:4) {
    Comparer <- AUCs[,jj]
    print(kk,jj)
    test.1 <- test.2 <- test.2.20 <- test.2.200 <- numeric()
    for (ii in 1:n){
      Reference[ii]
      Ind <- which(Comparer>Reference[ii])
      Ind.20 <- which(Comparer>0.2*Reference[ii])
      Ind.200 <- which(Comparer<2*Reference[ii])
      test.1 <- c(test.1,1/n)
      test.2 <- c(test.2,length(Ind)/n)
      test.2.20 <- c(test.2.20,length(Ind.20)/n)
      test.2.200 <- c(test.2.200,length(Ind.200)/n)
      
    }
    IntProbs[jj,kk]<-sum(test.1 *  (1-test.2))
    IntProbs.20[jj,kk]<-sum(test.1 *  (1-test.2.20))
    IntProbs.200[jj,kk]<-sum(test.1 *  (1-test.2.200))
  }
}

IntProbs <- round(IntProbs, 2)
IntProbs.20 <- round(IntProbs.20, 2)
IntProbs.200 <- round(IntProbs.200, 2)

