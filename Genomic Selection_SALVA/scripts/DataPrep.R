################################################################################
###                                                                          ###
###                      DATA PREPARATION                                    ###
###                                                                          ###
################################################################################

# 1) Read it in R as text data frame (or tibble)

  myG <- read.table(file="GenPea_SilDArT_sort_def.hmp.txt",sep= "\t",head=FALSE)

# 2) Make sure NAs are properly read it.
  
      # It seems that there is not missing values in the marker/sample

# 3) From letters to numbers:
  
  #I've donde this with a GAPPIT function and save it as a .txt, now the head is the "rs#" and the first column "taxa
  # "rs#" contains the SNP identifier

  myGAPIT <- GAPIT(G=myG, output.numerical=TRUE)
  myGD= myGAPIT$GD
  myGM= myGAPIT$GM

  write.table(myGD, "myGD.txt", sep = "\t")
  
# 4) Filter data:
    #Missing per marker and missing per sample
      
      # It seems that there is not missing values in the marker/sample
    
    #MAF (5%)
  
      
  
    #Heterozygosity (remove > 10%)
  