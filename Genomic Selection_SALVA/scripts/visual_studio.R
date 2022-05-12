#read marker DArT data set from file "Dart.xlsx"
DArT <- as.matrix(read.table("data/Dart.txt", sep = "\t", head = TRUE))
head(DArT[1:10,1:10])
library(openxlsx)
#read phenotype file
pheno <- read.xlsx("BLUPs_scaled.xlsx", colNames = T)
install.packages("GROAN")
install.packages("rrBLUP")
install.packages("bcv")
library(GROAN); library(rrBLUP); library(bcv)
#set wb
wb <- createWorkbench(folds = 10, reps = 20, regressor.name = "rrBLUP")
wb
addRegressor(wb, regressor = phenoRegressor.BGLR, regressor.name = "BL")
