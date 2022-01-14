############################################################################
###                                                                      ###
###                         DATA PREPARATION                             ###
###                                                                      ###
############################################################################

# 1) Read it in R as text data frame (or tibble)

  myG <- read.table(file="GenPea_SilDArT_sort_def.hmp.txt",sep= "\t",head=FALSE)

# 2) Make sure NAs are properly read it.
  
      # SilicoDarT markers don't have NA values. 

# 3) From letters to numbers:
  
  #I've donde this with a GAPPIT function and save it as a .txt, now the head is the "rs#" and the first column "taxa
  # "rs#" contains the marker identifier

    myGAPIT <- GAPIT(G=myG, output.numerical=TRUE)
    myGD= myGAPIT$GD
    myGM= myGAPIT$GM
  
    write.table(myGD, "myGD.txt", sep = "\t")
    #Here C = 2 = pres.; A = 0 = aus.; N = 1 = het. 
  
  #By TASSEL. I've convert the hapmap format (C,A, N) to numerical format (AA = 1, aa = 0, aA = 0,5 = NA),
    # TASSEL didn't recognice het. = N as 0.5, instead is Na value. But SilicoDarT markers doesn't contain Na's
    # So here I import the TASSEL markers and convert it to rrBLUP format (coded as {-1,0,1} = {aa,Aa,AA})
    
    DarT <- read.table("Markers.txt",sep= "\t",head=F)
    DarT <- data.frame(DarT)
      #Substitute 0 to -1 (homozygous minor aa) and  NA to 0 (heterozigosis aA or Aa)
    DarT[DarT == 0] <- -1
    DarT[is.na(DarT)] <- 0
    DarT <-t(DarT)
    write.table(x = DarT, file = "Dart.txt" ,sep= "\t")
    
    DarT1 <- read.xlsx("DarT.xlsx", sep= "\t" ,colNames = T)
    DarT1 %>%
      count()
    
    #set first row as header:
    install.packages("janitor")
    library(janitor)
    DarT <- row_to_names(dat = DarT, row_number = 1)
     #Count presence/absence (1/-1) of a marker ()
    DarT <- DarT %>%
      arrange(X1)
    

# 4) Filter data:
    #Missing per marker and missing per sample
      
      # It seems that there is not missing values in the marker/sample
    
    #MAF (5%)

    names(DarT) <- as.matrix(DarT[1, ])
    DarT <- DarT[-1, ]
    DarT[] <- lapply(DarT, function(x) type.convert(as.character(x)))
    dat

    DarT <- t(DarT)
      DarT %>%
        count(X.Marker.)
  
    #Heterozygosity (remove > 10%)
  