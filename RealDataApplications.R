rm(list = ls()) # clear the global environment
source("Implementations.R")


##### 4. Real data applications ------------------------------------------------

### 4.1. Physician data --------------------------------------------------------

rm(list = ls()) # clear the global environment
source("Implementations.R")
library(AER)
data(NMES1988)
data <- na.omit(data.matrix(NMES1988[,c(1,6,8,11,15,16)]))

par(mgp=c(2.2,1,0),mar=c(1,3.5,2,1),mfrow=c(1,3))
boxplot(data[,1],pch=16,cex.lab=1,cex.main = 1.3,main="Visits", col = "dodgerblue4",ylab = "Number of physician visits")
boxplot(data[,2],pch=16,cex.lab=1,cex.main = 1.3,main="Hospital", col = "dodgerblue4",ylab = "Number of hospital stays")
boxplot(data[,6],pch=16,cex.lab=1,cex.main = 1.3,main="Income", col = "dodgerblue4",ylab = "Family income in 10,000 USD")
par(mfrow=c(1,1))

DL <- DirectLingam(data) 
TSL <- TSLiNGAM(data,indep_meas = 1,slope = 1)
colnames(data)[DL$K] # causal order found by DirectLiNGAM
colnames(data)[TSL$K] # causal order found by TSLiNGAM
B_DL <- Bprune_adlasso(data,k=DL$K,B=DL$B); colnames(B_DL) <- colnames(data); rownames(B_DL) <- colnames(data)
B_TSL <- Bprune_adlasso(data,k=TSL$K,B=TSL$B); colnames(B_TSL) <- colnames(data); rownames(B_TSL) <- colnames(data)
B_DL; B_TSL # Pruned adjacency matrices B


### 4.2. GAGurine data ---------------------------------------------------------
# age is a cause of the concentration of the chemical GAG => causal order: (1 2)

rm(list = ls()) # clear the global environment
source("Implementations.R")
library(MASS)
X <- data.matrix(GAGurine)
DirectLingam(X)$K
TSLiNGAM(X,indep_meas = 1, slope=1)$K
# both methods find the correct causal order on the full data set

sizes <- seq(5,50,5)
correctDL <- rep(0,length(sizes))
correctTSL <- rep(0,length(sizes))
for(s in 1:length(sizes)){
  print(sizes[s])
  for(i in 1:1000){
    set.seed(i)
    set <- sample(1:314,sizes[s])
    if ( (DirectLingam(X[set,])$K)[2] == 2){
      correctDL[s] = correctDL[s] + 1
    }
    if ((TSLiNGAM(X[set,],indep_meas = 1, slope=1)$K)[2] == 2){
      correctTSL[s] = correctTSL[s] + 1
    }
  }
}

par(mgp=c(2.2,1,0),mar=c(3.5,3.5,1,1))
plot(sizes, correctTSL, type = "b", pch = 18, 
     col = "dodgerblue4", xlab = "sample size", ylab = "correct",ylim = c(450,800),lwd=2,cex.axis=0.9)
lines(sizes, correctDL, pch = 18, col = "darkolivegreen4", type = "b", lty = 2,lwd=2)
legend("bottomright", legend=c("TSLiNGAM","DirectLiNGAM"),
       col=c("dodgerblue4","darkolivegreen4"), lty = 1:2, cex=1,pch=18)


### 4.3. FMRI data -------------------------------------------------------------
# 5 variables, causal order is (1 2 3 4 5)

rm(list = ls()) # clear the global environment
source("Implementations.R")
datatotal <- data.matrix(read.table('fMRI_sim1.txt',header=FALSE,sep=",",dec="."))
# https://www.fmrib.ox.ac.uk/datasets/netsim/index.html

# NO OUTLIERS
data <- datatotal
DirectLingam(data)$K 
TSLiNGAM(data,indep_meas = 1, slope = 1)$K
# both methods find the correct causal order on the original data set

# 20 OUTLIERS ADDED: 0.2% contamination
DLcorrect = 0
TSLcorrect = 0
n_cont <- 20
for (i in 1:100){
  set.seed(i)
  data <- datatotal
  data[(10000-n_cont+1):10000,1] <- 25 + matrix(rt(n_cont,3),nrow = n_cont,ncol = 1) # contaminate the data
  if (sum(DirectLingam(data)$K == c(1,2,3,4,5)) == 5){
    DLcorrect = DLcorrect + 1
  }
  if (sum(TSLiNGAM(data,indep_meas = 1, slope=1)$K == c(1,2,3,4,5)) == 5){
    TSLcorrect = TSLcorrect + 1
  }
  print(i)
}
DLcorrect # number of times DirectLiNGAM finds the right causal order on the contaminated data
TSLcorrect # number of times TSLiNGAM finds the right causal order on the contaminated data


### 4.4. Sociological data -----------------------------------------------------

rm(list = ls()) 
source("Implementations.R")
# Data needs to be downloaded on their website
# https://gssdataexplorer.norc.org/gss_data 

DL <- DirectLingam(X); DL$K # DirectLiNGAM
TSL <- TSLiNGAM(X,indep_meas = 1,slope=1); TSL$K # TSLiNGAM
# Adaptive Lasso
B_DL <- Bprune_adlasso(X, k = DL$K, B = DL$B); B_DL
B_TSL <- Bprune_adlasso(X, k = TSL$K, B = TSL$B); B_TSL
# different independence measure: distance correlation
TSL_dcor <- TSLiNGAM(X,indep_meas = 2,slope=1); TSL_dcor$K
B_TSL_dcor <- Bprune_adlasso(X, k = TSL_dcor$K, B = TSL_dcor$B); B_TSL_dcor
