library(rrBLUP)
library(optparse)


option_list = list(
  make_option(c("-G", "--genoFile"), type="character", default=NULL,
              help="genotype file", metavar="character"),
  make_option(c("-P", "--phenoFile"), type="character", default=NULL,
              help="phenotype file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="OUT",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-I", "--iterations"), type="integer", default="100",
              help="the number of iterations for cross validation [default= %default]", metavar="integer"),
  make_option(c("ts", "--trainSize"), type="integer", default="70",
              help="proportion of population to use for training (value between 0 and 100), [default= %default]", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$genoFile)){
  print_help(opt_parser)
  stop("Let me know what the genotype file is called (-G).\n", call.=FALSE)
}

if (is.null(opt$phenoFile)){
  print_help(opt_parser)
  stop("Let me know what the phenotype file is called (-P).\n", call.=FALSE)
}


#read in the genotype file and impute
geno = as.matrix(read.table(opt$genoFile, header=FALSE))
impute = A.mat(geno, max.missing=0.5, impute.method="mean", return.imputed=T)
markers_impute = impute$imputed

#establish training and testing sets
N = nrow(markers_impute)
trainset = N*(opt$trainSize/100)
testset = N-trainset

#read in pheno data
pheno = as.matrix(read.table(opt$phenoFile, header=TRUE))

#establish a matrix to hold the results of each iteration
res<-matrix(NA,opt$iterations,2)

#################################

#perform cross validation using rrBLUP for predictive modelling


for (iter in 1:opt$iterations)
{

  train = as.matrix(sample(1:N,size=trainset,replace=FALSE))
  test = setdiff(1:N, train)


  #use above to pull out the records correpsonding to these (pheno and geno)
  geno_train = markers_impute[train,]
  geno_test = markers_impute[test,]
  pheno_train = pheno[train,]
  pheno_test = pheno[test,]

###############################################################################
## Using rrBLUP for prediction ##
###############################################################################
  #generate marker effects
  mxslv = mixed.solve(pheno_train, Z=geno_train, K=NULL, SE=FALSE, return.Hinv=FALSE)
  TPr = mxslv$u
  eff = as.matrix(TPr)
  #predict in the test set
  pred_test = geno_test %*% eff
  predictions = (pred_test[,1])+mxslv$beta
################################################################################


  #predictive ability
  acc = cor(predictions, pheno_test, use="complete")

  #caluclate the bias by regressing the phenotypes (y) on the GEBVs (x)
  reg = lm(pheno_test ~ predictions)
  #extract the predictionsHD coefficeint from the reg object
  coef = summary(reg)$coefficients[2,1]


  # Store the results
  res[iter,1] = acc
  res[iter,2] = coef

  write.table(res, file=paste(opt$out, ".its.txt", sep=""))
}

colnames(res) = c("PA", "Beta")
write.table(res, file=paste(opt$out, ".its.txt", sep=""))
write.table(summary(res), file=paste(opt$out, "summary.txt", sep="_"))