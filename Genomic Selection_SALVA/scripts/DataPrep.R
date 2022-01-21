############################################################################
###                                                                      ###
###                         DATA PREPARATION                             ###
###                                                                      ###
############################################################################
  library(openxlsx)
# 1) Read it in R as text data frame (or tibble)
  #I've open the DArT file (G,A and N) with TASSEL, I convert them into numerical and save the file as "Markers.txt":
  #where 1 = presence marker, 0 = absence markers, and NA = missing value. 
  DArT <- as.matrix(read.table("data/Markers.txt", header = T))
  DArT2 <- DArT[,-1]
  rownames(DArT2) <- DArT[,1]
  write.table(x = DArT2, file = "data/DArT.txt" ,sep= "\t")
  DArT <- as.matrix(read.table("data/DArT.txt", header = T))
  
# 2) Make sure NAs are properly read it.
  # Calculate the % of missing value within the matrix n x m:
  dim(DArT)
  (sum(is.na(DArT))/(length(DArT)))*100  # there is a 5.38 % of missing data into the matrix n x m

# 3) From letters to numbers:
  # By TASSEL. I've convert the hapmap format (C,A, N) to numerical format (AA = 1, aa = 0, aA = 0,5 = NA),
    # TASSEL didn't recognice het. = N as 0.5, instead is NA value. 
    # So here I import the TASSEL markers and convert it to rrBLUP format (coded as {-1,0,1} = {aa,Aa,AA})
    # Substitute 0 to -1 (homozygous minor aa)
    DArT[DArT == 0] <- -1
    write.table(x = DArT, file = "data/DArT_rrBLUP_withNA.txt" ,sep= "\t")
    
# 4) Filter data:
  # Missing per marker and missing per sample
    # It seems that there is not missing values in the marker/sample
    
  # MAF (5%)
    # Already done 

  # Heterozygosity (remove > 10%)
    # Already done
      
# 5) Imput your data "filling" all the NA's (Random forest, KNN imputation, BEAGLE...)
    #Random Forest:
    install.packages("missForest")
    library(missForest)
    DArT <- t(DArT)
    DArT <- as.data.frame(DArT)
    missForest(xmis =  DArT, # data matrix with missing values
               maxiter = 10, #m aximum number of iterations 
               ntree = 100, # number of trees to grow in each forest.
               variablewise = FALSE, # If 'TRUE' the OOB error is returned for each variable separately
               decreasing = FALSE, #If F then the variables are sorted w.r.t. increasing amount of missing entries during computation
               verbose = TRUE, #If 'TRUE' the user is supplied with additional output between iterations
               mtry = floor(sqrt(ncol(xmis))), #number of variables randomly sampled at each split
               replace = TRUE,
               classwt = NULL, 
               cutoff = NULL, 
               strata = NULL,
               sampsize = NULL, 
               nodesize = NULL, 
               maxnodes = NULL,
               xtrue = NA 
               #parallelize = 'variables'
               )
    #SVDI
    install.packages("bcv")
    library(bcv)  
    DArT_noNA <-
      impute.svd(DArT, # data matrix with missing values
               k = 4, #the rank of the SVD approximation, I use k = 4 following Nazzicari, N. 2016
               #tol = max(24279, 325) * 1e-10, #the convergence tolerance for the EM algorithm
               maxiter = 100 #the maximum number of EM steps to take
               )$x
    #Sustitute the values to get -1, 0, 1:
    DArT_noNA[DArT_noNA >= 0.5] <- 1
    DArT_noNA[DArT_noNA <= -0.5] <- -1
    DArT_noNA[DArT_noNA < 0.5 & DArT_noNA > -0.5]<- 0
    
    write.table(DArT_noNA, "DArT_noNA_SVDmethod.txt",sep= "\t" )
    
    

    

    