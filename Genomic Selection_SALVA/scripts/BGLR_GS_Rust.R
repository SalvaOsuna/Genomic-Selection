##########################################################
##                                                      ##
##          GS models in RUST pea G-BLUP models         ##
##                                                      ##
##########################################################
#From book GS machine learning
#Appendix 3: R Code Example 1
rm(list=ls())
library(BGLR)
load('dat_ls_E1.RData',verbose=TRUE)

#Phenotypic data
dat_F = as.matrix(read.xlsx(xlsxFile = "BLUP_field.xlsx", sep= "\t", rowNames = F, colNames = T, sheet = "BLUP_GS_rust"))
head(dat_F)

#Marker data
dat_M = as.matrix(read.table("data/DArT_rrBLUP_withNA.txt", header = T))
dat_M <- dat_M[-c(288, 294, 300, 320, 325), ] #Estas entradas no están evaluadas en CC así que las quito

dat_SVD <-
  impute.svd(dat_M, # data matrix with missing values
             k = 4, #the rank of the SVD approximation, I use k = 4 following Nazzicari, N. 2016
             #tol = max(24279, 325) * 1e-10, #the convergence tolerance for the EM algorithm
             maxiter = 100 #the maximum number of EM steps to take
  )$x

dim(dat_M)
dat_F = transform(dat_F, GEN = as.character(GEN))
head(dat_F,5)

#Matrix design of markers
Pos = match(dat_F$GEN,row.names(dat_M))
XM = dat_M[Pos,]
XM = scale(XM)
dim(XM)
n = dim(dat_F)[1]
y = dat_F$y

#10 random partitions
K = 10
set.seed(1)
PT = replicate(K,sample(n,0.20*n))

#BRR
ETA_BRR = list(list(model='BRR', X = XM))
Tab = data.frame(PT = 1:K, MSEP = NA)
set.seed(1)
for(k in 1:K)
{
  Pos_tst = PT[,k]
  y_NA = y
  y_NA[Pos_tst] = NA
  A = BGLR(y= y_NA, ETA= ETA_BRR, nIter = 1000, burnIn = 500, verbose = FALSE)
  yp_ts = A$yHat[Pos_tst]
  Tab$MSEP[k] = mean((y[Pos_tst]-yp_ts)^2)
}
#GBLUP
dat_M = scale(dat_M)
G = tcrossprod(XM)/dim(XM)[2]
dim(G)

#Matrix design of GIDs
Z = model.matrix(~0+GID,data=dat_F,xlev = list(GID=unique
                                               (dat_F$GID)))
K_L = Z%*%G%*%t(Z)
ETA_GB = list(list(model='RKHS',K = K_L))

#Tab = data.frame(PT = 1:K,MSEP = NA)
set.seed(1)
for(k in 1:K)
{
  Pos_tst = PT[,k]
  y_NA = y
  y_NA[Pos_tst] = NA
  A = BGLR(y=y_NA,ETA=ETA_GB,nIter = 1e4,burnIn = 1e3,verbose = FALSE)
  yp_ts = A$yHat[Pos_tst]
  Tab$MSEP_GB[k] = mean((y[Pos_tst]-yp_ts)^2)
}

#BA
ETA_BA = list(list(model='BayesA',X=XM))

#Tab = data.frame(PT = 1:K,MSEP = NA)
set.seed(1)
for(k in 1:K)
{
  Pos_tst = PT[,k]
  y_NA = y
  y_NA[Pos_tst] = NA
  A = BGLR(y=y_NA,ETA=ETA_BA,nIter = 1e4,burnIn = 1e3,verbose = FALSE)
  yp_ts = A$yHat[Pos_tst]
  Tab$MSEP_BA[k] = mean((y[Pos_tst]-yp_ts)^2)
}

#BB
ETA_BB = list(list(model='BayesB',X=XM))

#Tab = data.frame(PT = 1:K,MSEP = NA)
set.seed(1)
for(k in 1:K)
{
  Pos_tst = PT[,k]
  y_NA = y
  y_NA[Pos_tst] = NA
  A = BGLR(y=y_NA,ETA=ETA_BB,nIter = 1e4,burnIn = 1e3,verbose = FALSE)
  Appendix 3: R Code Example 1 201yp_ts = A$yHat[Pos_tst]
  Tab$MSEP_BB[k] = mean((y[Pos_tst]-yp_ts)^2)
}

#BC
ETA_BC = list(list(model='BayesC',X=XM))

#Tab = data.frame(PT = 1:K,MSEP = NA)
set.seed(1)
for(k in 1:K)
{
  Pos_tst = PT[,k]
  y_NA = y
  y_NA[Pos_tst] = NA
  A = BGLR(y=y_NA,ETA=ETA_BC,nIter = 1e4,burnIn = 1e3,verbose = FALSE)
  yp_ts = A$yHat[Pos_tst]
  Tab$MSEP_BC[k] = mean((y[Pos_tst]-yp_ts)^2)
}

#BL
ETA_BL = list(list(model='BL',X=XM))

#Tab = data.frame(PT = 1:K,MSEP = NA)
set.seed(1)
for(k in 1:K)
{
  Pos_tst = PT[,k]
  y_NA = y
  y_NA[Pos_tst] = NA
  A = BGLR(y=y_NA,ETA=ETA_BL,nIter = 1e4,burnIn = 1e3,verbose = FALSE)
  yp_ts = A$yHat[Pos_tst]
  Tab$MSEP_BL[k] = mean((y[Pos_tst]-yp_ts)^2)
}
Tab

#Mean and SD across the five partitions
apply(Tab[,-1],2,function(x)c(mean(x),sd(x)))

