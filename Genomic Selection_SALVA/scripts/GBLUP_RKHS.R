############################################
##                                        ##
##          Marker x ENVironment          ##
##                                        ##
############################################

#Within-Environment (i.e., stratified) GBLUP (model fitting
library('BGLR')
Pheno_rust <- as.matrix(read.table(file = "BLUP_GS_rust.txt", header = T)) #load phenotypes
Pheno_rust<- mapply(Pheno_rust, FUN=as.numeric) #convert matrix to numeric
Pheno_rust <- matrix(data=Pheno_rust, ncol=15, nrow=320) # convert matrix to numeric 2

#traits:
R19 <- 3 # choose any number in 1:ncol(Pheno_rust)
AUDPC <- 5
IF <- 7
IT <- 8
DS <- 9
Index <- 14

G <- as.matrix(read.xlsx(xlsxFile = "GenPea_SilDArT_Kinship_rust.xlsx", sheet = "Sheet 3", colNames = F)) #G matrix
prefix <- paste(colnames(Pheno_rust)[R19],"_",sep="")
# Fitting the model
ETA <- list(G=list(K=G,model='RKHS'))
fm_R19 <- BGLR(y=Pheno_rust[,R19],ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix) 
fm_AUDPC <- BGLR(y=Pheno_rust[,AUDPC],ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
fm_IF <- BGLR(y=Pheno_rust[,IF],ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
fm_IT <- BGLR(y=Pheno_rust[,IT],ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
fm_DS <- BGLR(y=Pheno_rust[,DS],ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
fm_Index <- BGLR(y=Pheno_rust[,Index],ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix) 

#Within-R19environment GBLUP (post-hoc)
# Extracting some estimates & predictions in R19
fm_R19$varE # residual variance
fm_R19$ETA[[1]]$varU # genomic variance
cor(Pheno_rust[,R19], fm_R19$ETA[[1]]$u) #Predictive ability
cor(Pheno_rust[,R19], fm_R19$yHat)
# Extracting some estimates & predictions in AUDPC
fm_AUDPC$varE # residual variance
fm_AUDPC$ETA[[1]]$varU # genomic variance
cor(Pheno_rust[,AUDPC], fm_AUDPC$ETA[[1]]$u) #Predictive ability
# Extracting some estimates & predictions in IF
fm_IF$varE # residual variance
fm_IF$ETA[[1]]$varU # genomic variance
cor(Pheno_rust[,IF], fm_IF$ETA[[1]]$u) #Predictive ability
# Extracting some estimates & predictions in IT
fm_IT$varE # residual variance
fm_IT$ETA[[1]]$varU # genomic variance
cor(Pheno_rust[,IT], fm_IT$ETA[[1]]$u) #Predictive ability
# Extracting some estimates & predictions in DS
fm_DS$varE # residual variance
fm_DS$ETA[[1]]$varU # genomic variance
cor(Pheno_rust[,DS], fm_DS$ETA[[1]]$u) #Predictive ability
# Extracting some estimates & predictions in Index
fm_Index$varE # residual variance
fm_Index$ETA[[1]]$varU # genomic variance
cor(Pheno_rust[,Index], fm_Index$ETA[[1]]$u) #Predictive ability

#Box 2a. Across-Environment Model (model fitting)
env_2 <- c(AUDPC,R19) # choose any set of environments from 1:ncol(Y)
nEnv <- length(env_2)
prefix <- paste(c('Across',colnames(Pheno_rust)[env_2],''),collapse='_')
y <- as.vector(Pheno_rust[,env_2])
# Fixed effect (env-intercepts)
envID <- rep(env_2,each=nrow(Pheno_rust))
ETA <- list(list(~factor(envID)-1,model="FIXED"))
# Effects of markers
G0 <- kronecker(matrix(nrow=nEnv,ncol=nEnv,1),G)
ETA[[2]] <- list(K=G0,model='RKHS')
# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
#Box 2b. Across-Environment Model (post-hoc)
# Extracting estimates of variance parameters
fm$varE # residual variance
fm$ETA[[2]]$varU # genomic variance
# Predictions (this is all within training)
tmpEnv <- 2
plot(y[envID==env_2[tmpEnv]]~fm$yHat[envID==env_2[tmpEnv]])
cor(y[envID==env_2[tmpEnv]],fm$yHat[envID==env_2[tmpEnv]])
# Samples
varE <- scan(paste(prefix,'varE.dat',sep=''))
plot(varE,type='o',cex=.5,col=4)
varU0 <- scan(paste(prefix,'ETA_2_varU.dat',sep=''))
plot(varU0,type='o',cex=.5,col=4)

#Box 3a. Marker-by-Environment Interaction Model (model fitting)
env_2 <- c(AUDPC,R19) # choose any set of environments from 1:ncol(Y)
nEnv <- length(env_2)
prefix <- paste(c('MxE',colnames(Pheno_rust)[env_2],''),collapse='_')
y <- as.vector(Pheno_rust[,env_2])
# Fixed effect (env-intercepts)
envID <- rep(env_2,each=nrow(Pheno_rust))
ETA <- list(list(~factor(envID)-1,model="FIXED"))
# Main effects of markers
G0 <- kronecker(matrix(nrow=nEnv,ncol=nEnv,1),G)
ETA[[2]] <- list(K=G0,model='RKHS')
# Adding interaction terms
for(i in 1:nEnv){
  tmp <- rep(0,nEnv) ; tmp[i] <- 1
  G1 <- kronecker(diag(tmp),G)
  ETA[[(i+2)]] <- list(K=G1, model='RKHS')
}
# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)


   # same with loop:
  folds=sample(1:10,size=n,replace=T)
  yHatCV=rep(NA,n)

    timeIn=proc.time()[3]
      for(i in 1:max(folds)){
        tst=which(folds==i)
        yNA=y
        yNA[tst]=NA
        fm=BGLR(y=yNA,ETA=ETA,nIter=6000,burnIn=1000)
        yHatCV[tst]=fm$yHat[tst]
      }
    proc.time()[3]-timeIn

#Box 3b. Marker-by-Environment Interaction Model (post-hoc)
# Extracting estimates of variance parameters
fm$varE # residual variance
fm$ETA[[2]]$varU # genomic variance (main effect)
vGInt <- rep(NA,nEnv)
for(i in 1:nEnv){ # interaction variances
  vGInt[i] <- fm$ETA[[(i+2)]]$varU
}
vGInt
# Predictions (this is all within training)
tmpEnv <- 2
plot(y[envID==env_2[tmpEnv]]~fm$yHat[envID==env_2[tmpEnv]])
cor(y[envID==env_2[tmpEnv]],fm$yHat[envID==env_2[tmpEnv]])

# Samples
varE <- scan(paste(prefix,'varE.dat',sep=''))
plot(varE,type='o',cex=.5,col=4)
varU0 <- scan(paste(prefix,'ETA_2_varU.dat',sep=''))
plot(varU0,type='o',cex=.5,col=4)

varU1 <- matrix(nrow=length(varU0),ncol=nEnv,NA)
for(i in 1:nEnv){
  varU1[,i] <- scan(paste(prefix,'ETA_',i+2,'_varU.dat',sep=''))
}

tmpEnv <- 2
plot(varU1[,tmpEnv],type='o',col=4,cex=.5)

#Box 4a. Creating a Testing Sets for CV1
env <- R19 # choose any set of environments from 1:ncol(Y)
nEnv <- length(env)
Y <- Pheno_rust[,env]
n <- nrow(Y)
percTST<-0.1
nTST <- round(percTST*320)
tst<-sample(1:320,size=nTST,replace=FALSE)
YNA <- Y
YNA[tst]<-NA

fm_R19 <- BGLR(y=YNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix) 
fm_R19$varE # residual variance
fm_R19$ETA[[1]]$varU # genomic variance
cor(Pheno_rust[,R19], fm_R19$ETA[[1]]$u) #Predictive ability

  #Vamos a usarlo, primero con single ENV
  YHatSE <- matrix(nrow=nrow(Y),ncol=ncol(Y),NA)
  ETA <- list(G=list(K=G,model='RKHS'))
    for(i in 1:nEnv){
      prefix <- paste(colnames(Y)[i],"_",sep="")
      fm <-BGLR(y=YNA[,i],ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
      YHatSE[,i] <- fm$yHat
    } 
  
  ## Across environment model (ignoring GxE) #######################
  yNA <- as.vector(YNA)
  # Fixed effect (env-intercepts)
  envID <- rep(env_2,each=nrow(Y))
  ETA <- list(list(~factor(envID)-1,model="FIXED"))
  # Main effects of markers
  G0 <- kronecker(matrix(nrow=nEnv,ncol=nEnv,1),G)
  ETA[[2]] <- list(K=G0,model='RKHS')
  # Model Fitting
  prefix <- paste(c('Across',colnames(Y),''),collapse='_')
  fm <- BGLR(y=yNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
  YHatAcross <- matrix(fm$yHat,ncol=nEnv)
  cor(YHatAcross)






#Box 4b. Creating a Testing Sets for CV2
env_2 <- c(AUDPC,R19) # choose any set of environments from 1:ncol(Y)
nEnv <- length(env_2)
Y <- Pheno_rust[,env_2]
n <- nrow(Y)
percTST<-0.1
nTST <- round(percTST*n)
nNA <- nEnv*nTST
if(nNA<n){ indexNA <- sample(1:n,nNA,replace=FALSE) }
if(nNA>=n){
  nRep <- floor(nNA/n)
  remain <- sample(1:n,nNA%%n,replace=FALSE)
  a0 <- sample(1:n,n,replace=FALSE)
  indexNA <- rep(a0,nRep)
  if(length(remain)>0){
    a1 <- floor(length(indexNA)/nTST)*nTST
    a2 <- nNA - a1 - length(remain)
    bb <- sample(a0[!a0%in%remain],a2,replace=FALSE)
    noInIndexNA <- c(rep(a0,nRep-1),a0[!a0%in%bb])
    indexNA <- c(noInIndexNA,bb,remain)
  }
}
indexEnv <- rep(1:nEnv,each=nTST)
YNA <- Y
for(j in 1:nEnv) YNA[indexNA[indexEnv==j],j] <- NA

fm$varE # residual variance
fm$ETA[[1]]$varU # genomic variance
cor(Pheno_rust[,R19], fm$ETA[[1]]$u) #Predictive ability

#Box 5. Fitting Models to TRN-TST Partitions (continues from Box 4b)
## Single environments models #####################################
YHatSE <- matrix(nrow=nrow(Y),ncol=ncol(Y),NA)
ETA <- list(G=list(K=G,model='RKHS'))
for(i in 1:nEnv){
  prefix <- paste(colnames(Y)[i],"_",sep="")
  fm <-BGLR(y=YNA[,i],ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
  YHatSE[,i] <- fm$yHat
}
## Across environment model (ignoring GxE) #######################
yNA <- as.vector(YNA)
# Fixed effect (env-intercepts)
envID <- rep(env,each=nrow(Y))
ETA <- list(list(~factor(envID)-1,model="FIXED"))
# Main effects of markers
G0 <- kronecker(matrix(nrow=nEnv,ncol=nEnv,1),G)
ETA[[2]] <- list(K=G0,model='RKHS')
# Model Fitting
prefix <- paste(c('Across',colnames(Y),''),collapse='_')
fm <- BGLR(y=yNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
YHatAcross <- matrix(fm$yHat,ncol=nEnv)
cor(YHatSE)
## MxE Interaction Model #########################################
# Adding interaction terms
for(i in 1:nEnv){
  tmp <- rep(0,nEnv) ; tmp[i] <- 1; G1 <- kronecker(diag(tmp),G)
  ETA[[(i+2)]] <- list(K=G1,model='RKHS')
}
# Model Fitting
prefix <- paste(c('MxE',colnames(Y),''),collapse='_')
fm <- BGLR(y=yNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
YHatInt <- matrix(fm$yHat,ncol=nEnv)


#Single environment with NA [CV0]####
#Single ENV R19 with NA [CV0]:
{
env <- R19 # choose any set of environments from 1:ncol(Y)
nEnv <- length(env)
Y <- Pheno_rust[,env]
n <- nrow(as.matrix(Y))
percTST<-0.1
nTST <- round(percTST*n)
rep = 50
mymat_R19 <- matrix(nrow = rep, ncol = 1)
prefix_R19 <- paste(colnames(Pheno_rust)[R19],"_",sep="")
for (i in 1:rep) {
  tst<-sample(1:n,size=nTST,replace=FALSE)
  YNA <- Y
  YNA[tst]<-NA
  fm_R19_CV1 <- BGLR(y=YNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix_R19, verbose = F) 
  mymat_R19[i,] <- cor(Pheno_rust[tst,R19], fm_R19_CV1$yHat[tst])
}
}
mean(mymat_R19)
#Single ENV AUDPC with NA [CV0]:
{
  env <- AUDPC # choose any set of environments from 1:ncol(Y)
  nEnv <- length(env)
  Y <- Pheno_rust[,env]
  n <- nrow(as.matrix(Y))
  percTST<-0.1
  nTST <- round(percTST*n)
  rep = 50
  mymat_AUDPC <- matrix(nrow = rep, ncol = 1)
  prefix_AUDPC <- paste(colnames(Pheno_rust)[AUDPC],"_",sep="")
  for (i in 1:rep) {
    tst<-sample(1:n,size=nTST,replace=FALSE)
    YNA <- Y
    YNA[tst]<-NA
    fm_AUDPC_CV1 <- BGLR(y=YNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix_AUDPC, verbose = F) 
    mymat_AUDPC[i,] <- cor(Pheno_rust[tst,AUDPC], fm_AUDPC_CV1$yHat[tst])
  }
}
mean(mymat_AUDPC)
#Single ENV IF with NA [CV0]:
{
  env <- IF # choose any set of environments from 1:ncol(Y)
  nEnv <- length(env)
  Y <- Pheno_rust[,env]
  n <- nrow(as.matrix(Y))
  percTST<-0.1
  nTST <- round(percTST*n)
  rep = 50
  mymat_IF <- matrix(nrow = rep, ncol = 1)
  prefix_IF <- paste(colnames(Pheno_rust)[IF],"_",sep="")
  for (i in 1:rep) {
    tst<-sample(1:n,size=nTST,replace=FALSE)
    YNA <- Y
    YNA[tst]<-NA
    fm_IF_CV1 <- BGLR(y=YNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix_IF, verbose = F) 
    mymat_IF[i,] <- cor(Pheno_rust[tst,IF], fm_IF_CV1$yHat[tst])
  }
}
mean(mymat_IF)
#Single ENV IT with NA [CV0]:
{
  env <- IT # choose any set of environments from 1:ncol(Y)
  nEnv <- length(env)
  Y <- Pheno_rust[,env]
  n <- nrow(as.matrix(Y))
  percTST<-0.1
  nTST <- round(percTST*n)
  rep = 50
  mymat_IT <- matrix(nrow = rep, ncol = 1)
  prefix_IT <- paste(colnames(Pheno_rust)[IT],"_",sep="")
  for (i in 1:rep) {
    tst<-sample(1:n,size=nTST,replace=FALSE)
    YNA <- Y
    YNA[tst]<-NA
    fm_IT_CV1 <- BGLR(y=YNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix_IT, verbose = F) 
    mymat_IT[i,] <- cor(Pheno_rust[tst,IT], fm_IT_CV1$yHat[tst])
  }
}
mean(mymat_IT)
#Single ENV DS with NA [CV0]:
{
  env <- DS # choose any set of environments from 1:ncol(Y)
  nEnv <- length(env)
  Y <- Pheno_rust[,env]
  n <- nrow(as.matrix(Y))
  percTST<-0.1
  nTST <- round(percTST*n)
  rep = 50
  mymat_DS <- matrix(nrow = rep, ncol = 1)
  prefix_DS <- paste(colnames(Pheno_rust)[DS],"_",sep="")
  for (i in 1:rep) {
    tst<-sample(1:n,size=nTST,replace=FALSE)
    YNA <- Y
    YNA[tst]<-NA
    fm_DS_CV1 <- BGLR(y=YNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix_DS, verbose = F) 
    mymat_DS[i,] <- cor(Pheno_rust[tst,DS], fm_DS_CV1$yHat[tst])
  }
}
mean(mymat_DS)
#Single ENV Index with NA [CV0]:
{
  env <- Index # choose any set of environments from 1:ncol(Y)
  nEnv <- length(env)
  Y <- Pheno_rust[,env]
  n <- nrow(as.matrix(Y))
  percTST<-0.1
  nTST <- round(percTST*n)
  rep = 50
  mymat_Index <- matrix(nrow = rep, ncol = 1)
  prefix_Index <- paste(colnames(Pheno_rust)[Index],"_",sep="")
  for (i in 1:rep) {
    tst<-sample(1:n,size=nTST,replace=FALSE)
    YNA <- Y
    YNA[tst]<-NA
    fm_Index_CV1 <- BGLR(y=YNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix_Index, verbose = F) 
    mymat_Index[i,] <- cor(Pheno_rust[tst,Index], fm_Index_CV1$yHat[tst])
  }
}
mean(mymat_Index)

GBLUP_singlENV <- cbind(mymat_R19, mymat_AUDPC, mymat_IF, mymat_IT, mymat_DS, mymat_Index)
write.xlsx(GBLUP_singlENV, "GBLUP_singlENV.xlsx")

#Across ENV (CC -AUDPC- and R19 -DS-)####
#Across-Environment Model (model fitting)
env_2 <- c(AUDPC,R19) # choose any set of environments from 1:ncol(Y)
nEnv <- length(env_2)
prefix <- paste(c('Across',colnames(Pheno_rust)[env_2],''),collapse='_')
y <- as.vector(Pheno_rust[,env_2])
# Fixed effect (env-intercepts)
envID <- rep(env_2,each=nrow(Pheno_rust))
ETA <- list(list(~factor(envID)-1,model="FIXED"))
# Effects of markers
G0 <- kronecker(matrix(nrow=nEnv,ncol=nEnv,1),G)
ETA[[2]] <- list(K=G0,model='RKHS')
# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
#Across-Environment Model (post-hoc)
# Extracting estimates of variance parameters
fm$varE # residual variance
fm$ETA[[2]]$varU # genomic variance
# Predictions (this is all within training)
tmpEnv <- 1
plot(y[envID==env_2[tmpEnv]]~fm$yHat[envID==env_2[tmpEnv]])
cor(y[envID==env_2[tmpEnv]],fm$yHat[envID==env_2[tmpEnv]])

#Across ENV with NA ignoring MxE. CV2:
env_2 <- c(AUDPC,R19) # choose any set of environments from 1:ncol(Y)
nEnv <- length(env_2)
Y <- Pheno_rust[,env_2]
n <- nrow(Y)
percTST<-0.1
nTST <- round(percTST*n)
nNA <- nEnv*nTST
if(nNA<n){ indexNA <- sample(1:n,nNA,replace=FALSE) }
if(nNA>=n){
  nRep <- floor(nNA/n)
  remain <- sample(1:n,nNA%%n,replace=FALSE)
  a0 <- sample(1:n,n,replace=FALSE)
  indexNA <- rep(a0,nRep)
  if(length(remain)>0){
    a1 <- floor(length(indexNA)/nTST)*nTST
    a2 <- nNA - a1 - length(remain)
    bb <- sample(a0[!a0%in%remain],a2,replace=FALSE)
    noInIndexNA <- c(rep(a0,nRep-1),a0[!a0%in%bb])
    indexNA <- c(noInIndexNA,bb,remain)
  }
}
indexEnv <- rep(1:nEnv,each=nTST)
YNA <- Y
for(j in 1:nEnv) YNA[indexNA[indexEnv==j],j] <- NA
# Fixed effect (env-intercepts)
envID <- rep(env_2,each=nrow(Y))
ETA <- list(list(~factor(envID)-1,model="FIXED"))
# Main effects of markers
G0 <- kronecker(matrix(nrow=nEnv,ncol=nEnv,1),G)
ETA[[2]] <- list(K=G0,model='RKHS')
# Model Fitting
prefix <- paste(c('Across',colnames(Y),''),collapse='_')
fm_acrossCV2 <- BGLR(y=YNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
YHatAcross <- matrix(fm_acrossCV2$yHat,ncol=nEnv)
cor(YHatAcross[,1], Pheno_rust[,3]) #R19 predicted vs R19 valid
cor(YHatAcross[,1], Pheno_rust[,5]) #R19 predicted vs AUDPC valid
cor(YHatAcross[,2], Pheno_rust[,3]) #AUDPC predicted vs R19 valid
cor(YHatAcross[,2], Pheno_rust[,5]) #AUDPC predicted vs AUDPC valid

#Across ENV with NA with MxE interaction. CV2:
env_2 <- c(AUDPC,R19) # choose any set of environments from 1:ncol(Y)
nEnv <- length(env_2)
Y <- Pheno_rust[,env_2]
n <- nrow(Y)
percTST<-0.1
nTST <- round(percTST*n)
nNA <- nEnv*nTST
if(nNA<n){ indexNA <- sample(1:n,nNA,replace=FALSE) }
if(nNA>=n){
  nRep <- floor(nNA/n)
  remain <- sample(1:n,nNA%%n,replace=FALSE)
  a0 <- sample(1:n,n,replace=FALSE)
  indexNA <- rep(a0,nRep)
  if(length(remain)>0){
    a1 <- floor(length(indexNA)/nTST)*nTST
    a2 <- nNA - a1 - length(remain)
    bb <- sample(a0[!a0%in%remain],a2,replace=FALSE)
    noInIndexNA <- c(rep(a0,nRep-1),a0[!a0%in%bb])
    indexNA <- c(noInIndexNA,bb,remain)
  }
}
indexEnv <- rep(1:nEnv,each=nTST)
YNA <- Y
for(j in 1:nEnv) YNA[indexNA[indexEnv==j],j] <- NA

yNA <- as.vector(YNA)
# Adding interaction terms
for(i in 1:nEnv){
  tmp <- rep(0,nEnv)
  tmp[i] <- 1
  G1 <- kronecker(diag(tmp),G)
  ETA[[(i+2)]] <- list(K=G1,model='RKHS')
}
# Model Fitting
prefix <- paste(c('MxE',colnames(Y),''),collapse='_')
fm_MxE_CV2 <- BGLR(y=yNA,ETA=ETA,nIter=12000,burnIn=2000,saveAt=prefix)
YHatInt <- matrix(fm_MxE_CV2$yHat,ncol=nEnv)
cor(YHatInt[,1], Pheno_rust[,3]) #R19 predicted vs R19 valid
cor(YHatInt[,1], Pheno_rust[,5]) #R19 predicted vs AUDPC valid
cor(YHatInt[,2], Pheno_rust[,3]) #AUDPC predicted vs R19 valid
cor(YHatInt[,2], Pheno_rust[,5]) #AUDPC predicted vs AUDPC valid

#computing correlations:
COR <- matrix(nrow=length(env),ncol=3,NA)
colnames(COR) <- c('SingleEnv', 'AcrossEnv', 'MxE')
rownames(COR) <- colnames(Y)
for(i in 1:nEnv){
  tst <- which(is.na(YNA[,i]))
  COR[i,1] <- cor(Y[tst,i],YHatSE[tst,i])
  COR[i,2] <- cor(Y[tst,i],YHatAcross[tst,i])
  COR[i,3] <- cor(Y[tst,i],YHatInt[tst,i])
}
COR
