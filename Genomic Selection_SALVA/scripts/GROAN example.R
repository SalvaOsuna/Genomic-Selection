#####################################################################################
###                                                                               ###
###                      GENOMIC SELECTION FOR PEA RUST                           ###
###                                                                               ###
#####################################################################################

#ISTALL AND ACTIVATE GROAN PACKAGE

  install.packages("GROAN")
  install.packages("BGLR")
  library(GROAN)
  library(BGLR)

  ## Example data
#arrays of phenotypes
yield.KI <- GROAN.KI$yield
yield.AI <- GROAN.AI$yield
#dataframes of SNP genotypes
df.SNP.KI <- GROAN.KI$SNPs
df.SNP.AI <- GROAN.AI$SNPs
#dataframes of realized genotypic kinship
df.kinship.KI <- GROAN.KI$kinship
df.kinship.AI <- GROAN.AI$kinship


#creating a GROAN.NoisyDataset without any extra noise injected
nds.no_noise = createNoisyDataset(
  name = 'PEA KI, no noise',
  genotypes = df.SNP.KI, 
  phenotypes = yield.KI
)

#creating a GROAN.NoisyDataset adding noise sampled from a normal distribution
nds.normal_noise = createNoisyDataset(
  name = 'PEA KI, normal noise',
  genotypes = df.SNP.KI, 
  phenotypes = yield.KI,
  noiseInjector = noiseInjector.norm,
  mean = 0,
  sd = sd(yield.KI) * 0.5
)

#creating a third dataset, this time with data from the AI lines
nds.no_noise.AI = createNoisyDataset(
  name = 'PEA AI, no noise',
  genotypes = df.SNP.AI, 
  phenotypes = yield.AI
)

#plotting the original phenotypes
plot(yield.KI, pch=20, main = 'True (black) vs. Noisy (red)', xlab = 'Samples', ylab = 'Phenotypes')
#plotting an instance of the phenotypes with noise injected 
points(getNoisyPhenotype(nds.normal_noise), col='red')

#average correlation oscillates around 0.89
cor(GROAN.KI$yield, getNoisyPhenotype(nds.normal_noise))
cor(GROAN.KI$yield, getNoisyPhenotype(nds.normal_noise))
cor(GROAN.KI$yield, getNoisyPhenotype(nds.normal_noise))

#Obviously in absence of noise no variability is present.
#no noise injector ==> the original phenotypes are returned
all(GROAN.KI$yield == getNoisyPhenotype(nds.no_noise))

#It is possible to invoke both print and summary methods to quickly inspect the created objects
print(nds.no_noise)
print(nds.no_noise.AI)

#creating a Workbench with default values 
#creating a GROAN.Workbench with default values explicitly assigned 
wb = createWorkbench(
  #parameters defining crossvalidation
  folds = 5, reps = 10, stratified = FALSE, 
  
  #parameters defining save-on-hard-disk policy
  outfolder = NULL, saveHyperParms = FALSE, saveExtraData = FALSE,
  
  #a regressor
  regressor = phenoRegressor.rrBLUP, regressor.name = 'rrBLUP'
)

#It is possible to update the GROAN.Workbench object by adding other regressors using the addRegressor function:
#adding a regressor to an existing Workbench
wb = addRegressor(
  #the Workbench to be updater
  wb,
  #the new regressor
  regressor = phenoRegressor.BGLR, regressor.name = 'Bayesian Lasso',
  
  #regressor-specific parameters
  type = 'BL'
)
print(wb)

#For reporting purposes it is possible to collate together (and filter) several result data frames, as shown here:
#executing two GROAN test, same workbench, different datasets
res.no_noise     = GROAN.run(nds.no_noise, wb)
res.normal_noise = GROAN.run(nds.normal_noise, wb)

summary(res.no_noise)
summary(res.normal_noise)

plotResult(
  res.no_noise,
  variable = c("pearson"),
  x.label = c("both"),
  plot.type = c("box"),
  strata = c("no_strata")
)


#putting the results together for further analysis
res.total = rbind(res.no_noise, res.normal_noise)

#defaults is a boxplot of Pearson's correlations
p = plotResult(res.total)
print(p)

#a barplot with 95% confidence interval of Pearson's correlations
p = plotResult(res.total, plot.type = 'bar_conf95')
print(p)

#a barplot of execution times per fold, in seconds
p = plotResult(res.total, plot.type = 'bar', variable = 'time')
print(p)

##Working with strata
#In the following code datasets GROAN.KI and GROAN.AI are put together, and strata are used to separate the two crosses.
#collating the two example datasets
nds.double = createNoisyDataset(
  name = 'KI and AI', 
  genotypes = rbind(GROAN.KI$SNPs, GROAN.AI$SNPs), 
  phenotypes = c(GROAN.KI$yield, GROAN.AI$yield),
  strata = c(rep('KI', 103), rep('AI', ,105)) #we have 103 KI and 105 AI
)

#the workbench will take into account strata
wb = createWorkbench(stratified = TRUE)

#ready to go
res = GROAN.run(nds.double, wb)
plotResult(res, strata = 'avg_strata', plot.type = 'bar')
summary(res)
