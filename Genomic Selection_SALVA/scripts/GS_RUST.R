############################################
##                                        ##
##          GS models in RUST pea         ##
##                                        ##
############################################

# 1) Load Genotypic data, DArT markers:
    DArT <- as.matrix(read.table("Dart.txt", header = T))
    DArT[DArT == 0] <- -1 #change 0 to -1
    DArT_rust <- DArT[-c(288, 294, 300, 320, 325), ] #Estas entradas no están evaluadas en CC así que las quito
    ##rrBLUP necesita la matriz en formato {-1, 0, 1} así que cambio los 0 por -1:

    #Also load SNP markers:
    myG <- as.matrix(read_excel("data/SNP_numeric2.xlsx", col_types = num, col_names = T))
    num <- sample("numeric", 11512, replace = T)
    
    myG_rust <- myG[-c(288, 294, 300, 320, 325), -c(1)] #Estas entradas no están evaluadas en CC así que las quito
    
    myG_rust[myG_rust == 1] <- 2
    myG_rust[myG_rust == 0.5] <- 1
    myG_rust <- as.data.frame(myG_rust)
    

    # 2) Make sure NAs are properly read it.
    # Calculate the % of missing value within the matrix n x m:
    dim(DArT_rust)
    (sum(is.na(DArT_rust))/(length(DArT_rust)))*100  # there is a 5.31 % of missing data into the matrix n x m
    
# 3) Impute your data "filling" all the NA's (SVDI)
    library(bcv)
    DArT_rust_SVDI <-
      impute.svd(DArT_rust, # data matrix with missing values
                 k = 4, #the rank of the SVD approximation, I use k = 4 following Nazzicari, N. 2016
                 #tol = max(24279, 325) * 1e-10, #the convergence tolerance for the EM algorithm
                 maxiter = 100 #the maximum number of EM steps to take
      )$x

# 4) Input phenotype. Traits in a matrix format. rows = GEN; column = trait
    Pheno_rust <- as.matrix(read.xlsx(xlsxFile = "BLUP_field.xlsx", sep= "\t", rowNames = F, colNames = T, sheet = "BLUP_GS_rust"))
    head(Pheno_rust)
    dim(Pheno_rust)
    dim(DArT_rust_SVDI)
    dim(myG_rust)
  
# 5) Define the training (70 % = 224 genotypes) and validation (30 % = 96 genotypes) populations
    train = as.matrix(sample(1:320, 224))
    test <- setdiff(1:320, train)      
    
    # Pheno_train and m_train are the phenotype and marker matrices for the values in the training population
    # Pheno_valid and m_valid will be the validation populations
    Pheno_train = Pheno_rust[train, ]
    m_train = DArT_rust_SVDI[train, ]
    
    Pheno_valid = Pheno_rust[test, ]
    m_valid = DArT_rust_SVDI[test, ]

# Run the function to calculate Correlation accuracy with 500 iterations:
    {
    #R18:
    traits = 1
    cycles = 500
    R18_accuracy = matrix(nrow = cycles, ncol = traits)
    for (r in 1:cycles) {
      train = as.matrix(sample(1:320, 224))
      test = setdiff(1:320, train)
      Pheno_train = Pheno_rust[train, ]
      m_train = DArT_rust_SVDI[train, ]
      Pheno_valid = Pheno_rust[test, ]
      m_valid = DArT_rust_SVDI[test, ]
      
      R18 = (Pheno_train[, 1])
      R18_answer <- mixed.solve(R18, Z = m_train, K = NULL, SE = F, return.Hinv = F)
      u = R18_answer$u
      e = as.matrix(u)
      pred_R18_valid = m_valid %*% e
      pred_R18 = pred_R18_valid[, 1] + R18_answer$beta  
      pred_R18
      R18_valid = Pheno_valid[, 1]
      R18_accuracy[r, 1] <-  cor(pred_R18_valid, R18_valid, use = "complete")
    }
    mean(R18_accuracy)
    
    #R19:
    traits = 1
    cycles = 500
    R19_accuracy = matrix(nrow = cycles, ncol = traits)
    for (r in 1:cycles) {
      train = as.matrix(sample(1:320, 224))
      test = setdiff(1:320, train)
      Pheno_train = Pheno_rust[train, ]
      m_train = DArT_rust_SVDI[train, ]
      Pheno_valid = Pheno_rust[test, ]
      m_valid = DArT_rust_SVDI[test, ]
      
      R19 = (Pheno_train[, 2])
      R19_answer <- mixed.solve(R19, Z = m_train, K = NULL, SE = F, return.Hinv = F)
      u = R19_answer$u
      e = as.matrix(u)
      pred_R19_valid = m_valid %*% e
      pred_R19 = pred_R19_valid[, 1] + R19_answer$beta  
      pred_R19
      R19_valid = Pheno_valid[, 2]
      R19_accuracy[r, 1] <-  cor(pred_R19_valid, R19_valid, use = "complete")
    }
    mean(R19_accuracy)
    
    #R20: 
    traits = 1
    cycles = 500
    R20_accuracy = matrix(nrow = cycles, ncol = traits)
    for (r in 1:cycles) {
      train = as.matrix(sample(1:320, 224))
      test = setdiff(1:320, train)
      Pheno_train = Pheno_rust[train, ]
      m_train = DArT_rust_SVDI[train, ]
      Pheno_valid = Pheno_rust[test, ]
      m_valid = DArT_rust_SVDI[test, ]
      
      R20 = (Pheno_train[, 3])
      R20_answer <- mixed.solve(R20, Z = m_train, K = NULL, SE = F, return.Hinv = F)
      u = R20_answer$u
      e = as.matrix(u)
      pred_R20_valid = m_valid %*% e
      pred_R20 = pred_R20_valid[, 1] + R20_answer$beta  
      pred_R20
      R20_valid = Pheno_valid[, 3]
      R20_accuracy[r, 1] <-  cor(pred_R20_valid, R20_valid, use = "complete")
    }
    mean(R20_accuracy)
    
    #AUDPC:
    traits = 1
    cycles = 500
    AUDPC_accuracy = matrix(nrow = cycles, ncol = traits)
    for (r in 1:cycles) {
      train = as.matrix(sample(1:320, 224))
      test = setdiff(1:320, train)
      Pheno_train = Pheno_rust[train, ]
      m_train = DArT_rust_SVDI[train, ]
      Pheno_valid = Pheno_rust[test, ]
      m_valid = DArT_rust_SVDI[test, ]
      
      AUDPC = (Pheno_train[, 4])
      AUDPC_answer <- mixed.solve(AUDPC, Z = m_train, K = NULL, SE = F, return.Hinv = F)
      u = AUDPC_answer$u
      e = as.matrix(u)
      pred_AUDPC_valid = m_valid %*% e
      pred_AUDPC = pred_AUDPC_valid[, 1] + AUDPC_answer$beta  
      pred_AUDPC
      AUDPC_valid = Pheno_valid[, 4]
      AUDPC_accuracy[r, 1] <-  cor(pred_AUDPC_valid, AUDPC_valid, use = "complete")
    }
    mean(AUDPC_accuracy)
    
    #LP50:
    traits = 1
    cycles = 500
    LP50_accuracy = matrix(nrow = cycles, ncol = traits)
    for (r in 1:cycles) {
      train = as.matrix(sample(1:320, 224))
      test = setdiff(1:320, train)
      Pheno_train = Pheno_rust[train, ]
      m_train = DArT_rust_SVDI[train, ]
      Pheno_valid = Pheno_rust[test, ]
      m_valid = DArT_rust_SVDI[test, ]
      
      LP50 = (Pheno_train[, 5])
      LP50_answer <- mixed.solve(LP50, Z = m_train, K = NULL, SE = F, return.Hinv = F)
      u = LP50_answer$u
      e = as.matrix(u)
      pred_LP50_valid = m_valid %*% e
      pred_LP50 = pred_LP50_valid[, 1] + LP50_answer$beta  
      pred_LP50
      LP50_valid = Pheno_valid[, 5]
      LP50_accuracy[r, 1] <-  cor(pred_LP50_valid, LP50_valid, use = "complete")
    }
    mean(LP50_accuracy)
    
    #IF:
    traits = 1
    cycles = 500
    IF_accuracy = matrix(nrow = cycles, ncol = traits)
    for (r in 1:cycles) {
      train = as.matrix(sample(1:320, 224))
      test = setdiff(1:320, train)
      Pheno_train = Pheno_rust[train, ]
      m_train = DArT_rust_SVDI[train, ]
      Pheno_valid = Pheno_rust[test, ]
      m_valid = DArT_rust_SVDI[test, ]
      
      IF = (Pheno_train[, 6])
      IF_answer <- mixed.solve(IF, Z = m_train, K = NULL, SE = F, return.Hinv = F)
      u = IF_answer$u
      e = as.matrix(u)
      pred_IF_valid = m_valid %*% e
      pred_IF = pred_IF_valid[, 1] + IF_answer$beta  
      pred_IF
      IF_valid = Pheno_valid[, 6]
      IF_accuracy[r, 1] <-  cor(pred_IF_valid, IF_valid, use = "complete")
    }
    mean(IF_accuracy)
    
    #IT:
    traits = 1
    cycles = 500
    IT_accuracy = matrix(nrow = cycles, ncol = traits)
    for (r in 1:cycles) {
      train = as.matrix(sample(1:320, 224))
      test = setdiff(1:320, train)
      Pheno_train = Pheno_rust[train, ]
      m_train = DArT_rust_SVDI[train, ]
      Pheno_valid = Pheno_rust[test, ]
      m_valid = DArT_rust_SVDI[test, ]
      
      IT = (Pheno_train[, 7])
      IT_answer <- mixed.solve(IT, Z = m_train, K = NULL, SE = F, return.Hinv = F)
      u = IT_answer$u
      e = as.matrix(u)
      pred_IT_valid = m_valid %*% e
      pred_IT = pred_IT_valid[, 1] + IT_answer$beta  
      pred_IT
      IT_valid = Pheno_valid[, 7]
      IT_accuracy[r, 1] <-  cor(pred_IT_valid, IT_valid, use = "complete")
    }
    mean(IT_accuracy)
    
    #DS:
    traits = 1
    cycles = 500
    DS_accuracy = matrix(nrow = cycles, ncol = traits)
    for (r in 1:cycles) {
      train = as.matrix(sample(1:320, 224))
      test = setdiff(1:320, train)
      Pheno_train = Pheno_rust[train, ]
      m_train = DArT_rust_SVDI[train, ]
      Pheno_valid = Pheno_rust[test, ]
      m_valid = DArT_rust_SVDI[test, ]
      
      DS = (Pheno_train[, 8])
      DS_answer <- mixed.solve(DS, Z = m_train, K = NULL, SE = F, return.Hinv = F)
      u = DS_answer$u
      e = as.matrix(u)
      pred_DS_valid = m_valid %*% e
      pred_DS = pred_DS_valid[, 1] + DS_answer$beta  
      pred_DS
      DS_valid = Pheno_valid[, 8]
      DS_accuracy[r, 1] <-  cor(pred_DS_valid, DS_valid, use = "complete")
    }
    mean(DS_accuracy)
    }
    accuracies.rust <- data.frame(R18_accuracy, R19_accuracy, R20_accuracy, AUDPC_accuracy, LP50_accuracy, IF_accuracy, IT_accuracy, DS_accuracy)
    write.xlsx(accuracies.rust, "df_accuracies_rust.xlsx", sep = "/t")
    boxplot(accuracies.rust,horizontal = T, las = 1,
            names = c("R18", "R19", "R20", "AUDPC", "LP50", "IF", "IT", "DS"),
            main = "Accuracies by rust trait with 500 itinerations"
            )
    #los datos de campo no salen demasiado bien, así que voy a probar haciendo arc.sin transform. sobre los blups a ver:
    #IC1:
    traits = 1
    cycles = 500
    IC1_accuracy = matrix(nrow = cycles, ncol = traits)
    for (r in 1:cycles) {
      train = as.matrix(sample(1:320, 224))
      test = setdiff(1:320, train)
      Pheno_train = Pheno_rust[train, ]
      m_train = DArT_rust_SVDI[train, ]
      Pheno_valid = Pheno_rust[test, ]
      m_valid = DArT_rust_SVDI[test, ]
      
      IC1 = (Pheno_train[, 9])
      IC1_answer <- mixed.solve(IC1, Z = m_train, K = NULL, SE = F, return.Hinv = F)
      u = IC1_answer$u
      e = as.matrix(u)
      pred_IC1_valid = m_valid %*% e
      pred_IC1 = pred_IC1_valid[, 1] + IC1_answer$beta  
      pred_IC1
      IC1_valid = Pheno_valid[, 9]
      IC1_accuracy[r, 1] <-  cor(pred_IC1_valid, IC1_valid, use = "complete")
    }
    mean(IC1_accuracy)
    
    #R19_t:
    traits = 1
    cycles = 500
    R19_t_accuracy = matrix(nrow = cycles, ncol = traits)
    for (r in 1:cycles) {
      train = as.matrix(sample(1:320, 224))
      test = setdiff(1:320, train)
      Pheno_train = Pheno_rust[train, ]
      m_train = DArT_rust_SVDI[train, ]
      Pheno_valid = Pheno_rust[test, ]
      m_valid = DArT_rust_SVDI[test, ]
      
      R19_t = (Pheno_train[, 10])
      R19_t_answer <- mixed.solve(R19_t, Z = m_train, K = NULL, SE = F, return.Hinv = F)
      u = R19_t_answer$u
      e = as.matrix(u)
      pred_R19_t_valid = m_valid %*% e
      pred_R19_t = pred_R19_t_valid[, 1] + R19_t_answer$beta  
      pred_R19_t
      R19_t_valid = Pheno_valid[, 10]
      R19_t_accuracy[r, 1] <-  cor(pred_R19_t_valid, R19_t_valid, use = "complete")
    }
    mean(R19_t_accuracy)
    
    #R20_t: 
    traits = 1
    cycles = 500
    R20_t_accuracy = matrix(nrow = cycles, ncol = traits)
    for (r in 1:cycles) {
      train = as.matrix(sample(1:320, 224))
      test = setdiff(1:320, train)
      Pheno_train = Pheno_rust[train, ]
      m_train = DArT_rust_SVDI[train, ]
      Pheno_valid = Pheno_rust[test, ]
      m_valid = DArT_rust_SVDI[test, ]
      
      R20_t = (Pheno_train[, 11])
      R20_t_answer <- mixed.solve(R20_t, Z = m_train, K = NULL, SE = F, return.Hinv = F)
      u = R20_t_answer$u
      e = as.matrix(u)
      pred_R20_t_valid = m_valid %*% e
      pred_R20_t = pred_R20_t_valid[, 1] + R20_t_answer$beta  
      pred_R20_t
      R20_t_valid = Pheno_valid[, 11]
      R20_t_accuracy[r, 1] <-  cor(pred_R20_t_valid, R20_t_valid, use = "complete")
    }
    mean(R20_t_accuracy)

    
#Hoy voy a probar prediction accuracy con multi trait index en lugar de solo DS en cámara y ver si se puede aplicar a campo:
    {
    #Markers:
     dim(DArT_rust_SVDI)
    #Pheno:
     dim(Pheno_rust)
    #rrBLUP 500 itinerancies:
     #I_cc_FAI_LP as IC1:
     traits = 1
     cycles = 500
     IC1_accuracy = matrix(nrow = cycles, ncol = traits)
     for (r in 1:cycles) {
       train = as.matrix(sample(1:320, 224))
       test = setdiff(1:320, train)
       Pheno_train = Pheno_rust[train, ]
       m_train = DArT_rust_SVDI[train, ]
       Pheno_valid = Pheno_rust[test, ]
       m_valid = DArT_rust_SVDI[test, ]
       
       IC1 = (Pheno_train[, 12])
       IC1_answer <- mixed.solve(IC1, Z = m_train, K = NULL, SE = F, return.Hinv = F)
       u = IC1_answer$u
       e = as.matrix(u)
       pred_IC1_valid = m_valid %*% e
       pred_IC1 = pred_IC1_valid[, 1] + IC1_answer$beta  
       pred_IC1
       IC1_valid = Pheno_valid[, 12]
       IC1_accuracy[r, 1] <-  cor(pred_IC1_valid, IC1_valid, use = "complete")
     }
     mean(IC1_accuracy)
     #I_cc_FAI as IC2:
     traits = 1
     cycles = 500
     IC2_accuracy = matrix(nrow = cycles, ncol = traits)
     for (r in 1:cycles) {
       train = as.matrix(sample(1:320, 224))
       test = setdiff(1:320, train)
       Pheno_train = Pheno_rust[train, ]
       m_train = DArT_rust_SVDI[train, ]
       Pheno_valid = Pheno_rust[test, ]
       m_valid = DArT_rust_SVDI[test, ]
       
       IC2 = (Pheno_train[, 13])
       IC2_answer <- mixed.solve(IC2, Z = m_train, K = NULL, SE = F, return.Hinv = F)
       u = IC2_answer$u
       e = as.matrix(u)
       pred_IC2_valid = m_valid %*% e
       pred_IC2 = pred_IC2_valid[, 1] + IC2_answer$beta  
       pred_IC2
       IC2_valid = Pheno_valid[, 13]
       IC2_accuracy[r, 1] <-  cor(pred_IC2_valid, IC2_valid, use = "complete")
     }
     mean(IC2_accuracy)
     #I_cc_SH as IC3:
     traits = 1
     cycles = 500
     IC3_accuracy = matrix(nrow = cycles, ncol = traits)
     for (r in 1:cycles) {
       train = as.matrix(sample(1:320, 224))
       test = setdiff(1:320, train)
       Pheno_train = Pheno_rust[train, ]
       m_train = DArT_rust_SVDI[train, ]
       Pheno_valid = Pheno_rust[test, ]
       m_valid = DArT_rust_SVDI[test, ]
       
       IC3 = (Pheno_train[, 14])
       IC3_answer <- mixed.solve(IC3, Z = m_train, K = NULL, SE = F, return.Hinv = F)
       u = IC3_answer$u
       e = as.matrix(u)
       pred_IC3_valid = m_valid %*% e
       pred_IC3 = pred_IC3_valid[, 1] + IC3_answer$beta  
       pred_IC3
       IC3_valid = Pheno_valid[, 14]
       IC3_accuracy[r, 1] <-  cor(pred_IC3_valid, IC3_valid, use = "complete")
     }
     mean(IC3_accuracy)
    }
    
    df_index <- data.frame(IC1_accuracy, IC2_accuracy, IC3_accuracy)
    write.xlsx(df_index, "results/accuracies_indixes_CC.xlsx", sep = "/t")
    df_index <- as.data.frame(read.xlsx(xlsxFile = "results/accuracies_indixes_CC.xlsx", colNames = T))
    ggplot(df_index, aes(x= Trait, y = Accuracy)) +
      geom_boxplot(aes(fill = Trait)) +
      coord_flip() +
      labs(title = "Accuracy by single trait or indices", 
           subtitle = "500 itinerances, 70% training 30% testing",
           caption = "FAI 1 = With LP50, FAI 2 = Without LP50") +
      theme_bw() +
      theme(legend.position = "none")
    
    #Multi trait index FAI-BLUP mejora los resultados de accuracy usando los traits AUDPC, IF, IT e DS

#Ahora voy a hacer lo mismo con GROAN intentando 50 rep. de 10 fold cross-validation:
    library(GROAN)
  Pheno_rust_df <- as.data.frame(Pheno_rust)
  head(Pheno_rust_df)
  #parece que hay que dejar los datos en formato {2, 1, 0}
  DArT <- as.matrix(read.table("DArT.txt", header = T))
  DArT[DArT == 1] <- 2 #change 1 to 2
  DArT_GROAN <- DArT[-c(288, 294, 300, 320, 325), ] #Estas entradas no están evaluadas en CC así que las quito
  DArT_GROAN_SVDI <-
    impute.svd(DArT_GROAN, # data matrix with missing values
               k = 4, #the rank of the SVD approximation, I use k = 4 following Nazzicari, N. 2016
               #tol = max(24279, 325) * 1e-10, #the convergence tolerance for the EM algorithm
               maxiter = 100 #the maximum number of EM steps to take
    )$x
  
  DArT_GROAN_SVDI[DArT_GROAN_SVDI >= 1.5] <- 2
  DArT_GROAN_SVDI[DArT_GROAN_SVDI <= 0.5] <- 0
  DArT_GROAN_SVDI[DArT_GROAN_SVDI > 0.5 & DArT_GROAN_SVDI < 1.5]<- 1
  
  DArT_GROAN_SVDI <- as.data.frame(DArT_GROAN_SVDI)
  DArT_rust_SVDI <- as.matrix(read.table("DArT_noNA_SVDmethod.txt",header = T))
    #creating a dataset for DS
    nds.no_noise.DS <- createNoisyDataset(
      name = 'PEA DS, no noise',
      genotypes = DArT_rust_SVDI, 
      phenotypes = Pheno_rust_df$DS)
    #create wb
    wb <- createWorkbench(regressor = phenoRegressor.rrBLUP, 
      regressor.name = 'rrBLUP',
      stratified = T,
      reps = 50,
      folds = 10)
    #Run GROAN 
    res.no_noise <- GROAN.run(nds.no_noise.DS, wb)
    plotResult(res.no_noise)
    
    
    #Create others dataset:
    #I_cc_FAI_LP
    nds.no_noise.FAI1<- createNoisyDataset(
      name = 'PEA FAI, no noise',
      genotypes = DArT_GROAN_SVDI, 
      phenotypes = Pheno_rust_df$I_cc_FAI_LP)
    
    #I_cc_FAI
    nds.no_noise.FAI2<- createNoisyDataset(
      name = 'PEA FAI2, no noise',
      genotypes = DArT_rust_SVDI, 
      phenotypes = Pheno_rust_df$I_cc_FAI)
    
    #I_cc_SH
    nds.no_noise.SH<- createNoisyDataset(
      name = 'PEA SH, no noise',
      genotypes = DArT_GROAN_SVDI, 
      phenotypes = Pheno_rust_df$I_cc_SH)
    
    #RUN GROAN for those three with the same wb:
    mean(res.no_noise$pearson)
    res.no_noise.FAI1 <- GROAN.run(nds.no_noise.FAI1, wb)
    plotResult(res.no_noise.FAI1)
    mean(res.no_noise.FAI1$pearson)
    res.no_noise.FAI2 <- GROAN.run(nds.no_noise.FAI2, wb)
    plotResult(res.no_noise.FAI2)
    mean(res.no_noise.FAI2$pearson)
    res.no_noise.SH <- GROAN.run(nds.no_noise.SH, wb)
    plotResult(res.no_noise.SH)
    mean(res.no_noise.SH$pearson)
    res.total1 <- rbind(res.no_noise, res.no_noise.FAI1, res.no_noise.FAI2, res.no_noise.SH)
    plotResult(res.total1) +
      theme_bw() +
      labs(title = "Accuracies between diferents Rust Indices",
           subtitle = "10 folds - 50 repetitions",
           caption = "FAI-BLUP without LP50 is better than the others")
    write.xlsx(res.total1, "Accuracies_RustIndices_rrBLUP.xlsx")
    
    #creating the BL regressor:
    library(BGLR)
    wb2 <- createWorkbench(phenoRegressor.BGLR, 
                           regressor.name = 'Bayesian Lasso',
                           outfolder = "/results",
                           saveHyperParms=F,
                           saveExtraData = F,
                          type = 'BL',
                          stratified = T,
                          reps = 50,
                          folds = 10)
    print(wb2)
    #re run the code with the new wb:
    res.no_noise.FAI2.bl <- GROAN.run(nds.no_noise.FAI2, wb2)
    plotResult(res.no_noise.FAI2.bl)
    write.xlsx(res.total, "BL_vs_rrBLUP_IndexDS.xlsx")
    res.total <- read.xlsx(xlsxFile = "/results/BL_vs_rrBLUP_IndexDS.xlsx")
    res.total %>%
      group_by(regressor) %>%
      summarise("meanP" = median(pearson))
    
    plotResult(res.total, variable = "pearson") #No differences between regressor
    ggbetweenstats(data = res.total, x = regressor, y = pearson, type = "p", pairwise.comparisons = T, 
                   p.adjust.method = "holm", k = 4)
    
    
    #putting the results best index together for further analysis
    res.total = rbind(res.no_noise.FAI2.bl, res.no_noise.FAI2)
    
    #defaults is a boxplot of Pearson's correlations
    p = plotResult(res.total)
    p +
      theme_bw() +
      labs(title = "Accuracies between models in the same Rust Index",
           subtitle = "10 folds - 50 repetitions",
           caption = "BL seems slightly better")

    #a barplot with 95% confidence interval of Pearson's correlations
    p = plotResult(res.total, plot.type = 'bar_conf95')
    print(p)
    
    plotResult(res.total, plot.type = 'bar', variable = 'time')

#train model with the index and test in field:
    #R18
    nds.no_noise.R18<- createNoisyDataset(
      name = 'R18, no noise',
      genotypes = DArT_GROAN_SVDI, 
      phenotypes = Pheno_rust_df$R18)
    #R19
    nds.no_noise.R19<- createNoisyDataset(
      name = 'R19, no noise',
      genotypes = DArT_rust_SVDI, 
      phenotypes = Pheno_rust_df$R19)
    #R20
    nds.no_noise.R20<- createNoisyDataset(
      name = 'R20, no noise',
      genotypes = DArT_GROAN_SVDI, 
      phenotypes = Pheno_rust_df$R20)
    #a new GROAN.Workbench with NO crossvalidation and only one repetition
    wb3 = createWorkbench(
      folds = 10, reps = 15, 
      regressor.name = 'rrBLUP', regressor = phenoRegressor.rrBLUP)

    #training on PEA.KI, testing on PEA.AI
    res = GROAN.run(nds = nds.no_noise.FAI2, wb = wb3, nds.test = 
                      list(nds.no_noise.R18, nds.no_noise.R19, nds.no_noise.R20))
    
    print(res[,c('dataset.train', 'dataset.test', 'pearson')])
    res.1 = GROAN.run(nds = nds.no_noise.FAI2, wb = wb3, nds.test = 
      nds.no_noise.DS)
    print(res.1[,c('dataset.train', 'dataset.test', 'pearson')])
    res %>%
      group_by(dataset.train, dataset.test) %>%
      summarise("meanPA" = abs(mean(pearson)))

  #Cual es el mejor predictor sobre DS2019 de los datos de cámara y campo?
    res.2 = GROAN.run(nds = nds.no_noise.DS, wb = wb3, nds.test = nds.no_noise.R19)
    print(res.2[,c('dataset.train', 'dataset.test', 'pearson')])
    
    res.3 = GROAN.run(nds = nds.no_noise.FAI2, wb = wb3, nds.test = nds.no_noise.R19)
    print(res.3[,c('dataset.train', 'dataset.test', 'pearson')])
    predictors19 <- rbind(res.2,res.3)
    predictors19$pearson <- abs(predictors19$pearson)
    head(predictors19)
    plotResult(predictors19, variable = "pearson")
    predictors19 %>%
      group_by(dataset.train, dataset.test) %>%
      summarise("meanPA" = mean(pearson))
    
    #AUDPC
    nds.no_noise.AUDPC<- createNoisyDataset(
      name = 'AUDPC, no noise',
      genotypes = DArT_rust_SVDI, 
      phenotypes = Pheno_rust_df$AUDPC)
    res.AUDPC = GROAN.run(nds = nds.no_noise.AUDPC, wb = wb3, nds.test = 
                      list(nds.no_noise.R18, nds.no_noise.R19, nds.no_noise.R20))
    #LP50
    nds.no_noise.LP50<- createNoisyDataset(
      name = 'LP50, no noise',
      genotypes = DArT_GROAN_SVDI, 
      phenotypes = Pheno_rust_df$LP50)
    res.LP50 = GROAN.run(nds = nds.no_noise.LP50, wb = wb3, nds.test = 
                            list(nds.no_noise.R18, nds.no_noise.R19, nds.no_noise.R20))
    #IF
    nds.no_noise.IF<- createNoisyDataset(
      name = 'IF, no noise',
      genotypes = DArT_rust_SVDI, 
      phenotypes = Pheno_rust_df$IF)
    res.IF = GROAN.run(nds = nds.no_noise.IF, wb = wb3, nds.test = 
                            list(nds.no_noise.R18, nds.no_noise.R19, nds.no_noise.R20))
    #IT
    nds.no_noise.IT<- createNoisyDataset(
      name = 'IT, no noise',
      genotypes = DArT_rust_SVDI, 
      phenotypes = Pheno_rust_df$IT)
    res.IT = GROAN.run(nds = nds.no_noise.IT, wb = wb3, nds.test = 
                            list(nds.no_noise.R18, nds.no_noise.R19, nds.no_noise.R20))
    #DS
    nds.no_noise.DS<- createNoisyDataset(
      name = 'DS, no noise',
      genotypes = DArT_GROAN_SVDI, 
      phenotypes = Pheno_rust_df$DS)
    res.DS = GROAN.run(nds = nds.no_noise.DS, wb = wb3, nds.test = 
                            list(nds.no_noise.R18, nds.no_noise.R19, nds.no_noise.R20))
    #FAI2
    nds.no_noise.FAI2<- createNoisyDataset(
      name = 'FAI2, no noise',
      genotypes = DArT_GROAN_SVDI, 
      phenotypes = Pheno_rust_df$I_cc_FAI)
    res.FAI2 = GROAN.run(nds = nds.no_noise.FAI2, wb = wb3, nds.test = 
                            list(nds.no_noise.R18, nds.no_noise.R19, nds.no_noise.R20))
    #bind and plot
    test1 <- rbind(res.AUDPC, res.LP50, res.IF, res.IT, res.DS, res.FAI2)
    write.xlsx(test1,"results/CCtrain_fieldtest.xlsx")
    test.rrblup <- test1 %>%
      group_by(dataset.train, dataset.test) %>%
      summarise("meanPA" = (mean(pearson)))
    plotResult(test1, variable = "pearson")
    write.xlsx(test.rrblup, "results/CCtrain_fieldtest.xlsx")
    ggplot(test1, aes(x= dataset.train, y = abs(pearson)))+
      geom_boxplot(aes(fill = dataset.test)) +
      theme_bw()
    
#ahora voy a hacer un wb que incluya los cuatro regresores más usados: rrBLUP, GBLUP, BL y RKHS
    #sabiendo que FAI2 es el mejor predictor para DS2019, ahora voy a ver qué modelo es mejor:
    wb4 = createWorkbench(
      #parameters defining crossvalidation
      folds = 10, reps = 50, stratified = F, 
      
      #parameters defining save-on-hard-disk policy
      outfolder = NULL, saveHyperParms = F, saveExtraData = F,
      
      #a regressor
      regressor = phenoRegressor.rrBLUP, regressor.name = 'rrBLUP'
    )
    wb5 = addRegressor(
      #the Workbench to be updater
      wb4,
      #the new regressor
      regressor = phenoRegressor.BGLR, regressor.name = 'Bayesian Lasso',
      
      #regressor-specific parameters
      type = 'BL'
    )
    print(wb5)
    dwb6 = addRegressor(
      wb5,
      regressor = phenoRegressor.BGLR, regressor.name = 'G-BLUP'
    )
    wb7 = addRegressor(
      wb6,
      regressor = phenoRegressor.BGLR, regressor.name = 'RKHS'
    )
print(wb7)    
res.multimodel = GROAN.run(nds = nds.no_noise.FAI2, wb = wb7, nds.test = nds.no_noise.R19)
write.xlsx(res.multimodel, "res.multimodel_multitrait_DS2019.xlsx")
    plotResult(res.multimodel)
    res.multimodel %>%
      group_by(dataset.train, dataset.test, regressor) %>%
      summarise("meanPA" = mean(pearson))

res.BL.rrBLUP = GROAN.run(nds = nds.no_noise.FAI2, wb = wb5, nds.test = nds.no_noise.R19)
plotResult(res.BL.rrBLUP)
res.BL.rrBLUP %>%
  group_by(dataset.train, dataset.test, regressor) %>%
  summarise("meanPA" = mean(pearson))

 my.pheno = Pheno_rust_df$I_cc_FAI
my.pheno[160:320] = NA
res1 = phenoRegressor.rrBLUP(phenotypes = my.pheno, genotypes = DArT_GROAN_SVDI)
print(names(res1))
plot(
   x = Pheno_rust_df$R19, xlab = "Real values",
   y = res1$predictions, ylab = "Predicted values")
abline(a=0, b=1) #adding first quadrant bisector, for reference
cor(Pheno_rust_df$R19,res1$predictions)
cor(Pheno_rust_df$DS,Pheno_rust_df$R19)

#Voy a repetir un GROAN con los dos modelos más usados rrBLUP y GBLUP para cada trait campo y cámara:
res.no_noise.AUDPC <- GROAN.run(nds.no_noise.AUDPC, wb5)
res.no_noise.DS <- GROAN.run(nds.no_noise.DS, wb5)
res.no_noise.FAI2 <- GROAN.run(nds.no_noise.FAI2, wb5)
res.no_noise.IF <- GROAN.run(nds.no_noise.IF, wb5)
res.no_noise.IT <- GROAN.run(nds.no_noise.IT, wb5)
res.no_noise.LP50 <- GROAN.run(nds.no_noise.LP50, wb5)
res.no_noise.R18 <- GROAN.run(nds.no_noise.R18, wb5)
res.no_noise.R19 <- GROAN.run(nds.no_noise.R19, wb5)
res.no_noise.R20 <- GROAN.run(nds.no_noise.R20, wb5)

res.total2 <- rbind(res.no_noise.AUDPC, res.no_noise.DS, res.no_noise.FAI2, res.no_noise.IF, res.no_noise.IT,
                         res.no_noise.LP50, res.no_noise.R18, res.no_noise.R19, res.no_noise.R20)
write.xlsx(res.total2, "results/field&CC_traits_twomodels.xlsx")

plotResult(res.total2)
res.total2 %>%
  group_by(dataset.train, dataset.test, regressor) %>%
  summarise("meanPA" = mean(pearson))
print(wb2)

res.DStoR19 <- GROAN.run(nds = nds.no_noise.DS, wb = wb2, nds.test = nds.no_noise.R19)
write.xlsx(res.DStoR19, "results/res.DStoR19_BL.xlsx")
plotResult(res.DStoR19)
res.DStoR19 %>%
  group_by(dataset.train, dataset.test, regressor) %>%
  summarise("meanPA" = mean(pearson))

#Bayesian Genomic single-trait, multi-environment Linear Regression Model
#Function description:
BMORS_Env(
  data = NULL, #(data.frame) Phenotypic response where each column is a different trait and the first column are the name of the environment where it was evaluated.
  testingEnv = "", #(string) Name of the Environment to test.
  ETA = NULL, #(matrix) This is a two-level list used to specify the regression function (or linear predictor).
  covModel = "BRR", # (string) Name of the covariates model to implement (BRR, BayesA, BayesB, BayesC)
  predictor_Sec_complete = FALSE, #(Logical) FALSE by default
  nIter = 2500, #(integer) Number of iterations to fit the model.
  burnIn = 500, #(integer) Number of items to burn at the beginning of the model
  thin = 5, #(integer) Number of items to thin the model.
  progressBar = TRUE, #(Logical) Show the progress bar.
  digits = 4 #(integer) Number of digits of accuracy in the results.
)

pheno <- data.frame(read.xlsx("multiENV_rust.xlsx", sheet = "Sheet 1", colNames = T))
#Matrix design
geno <- data.frame(read.xlsx("data/GenPea_SilDArT_Kinship_rust.xlsx", sheet = "Sheet 1", colNames = T, rowNames = T))
LG <- cholesky(geno)
write.csv(LG, "cholesky_geno.txt", sep = "_")
ZG <- model.matrix(~0 + as.factor(pheno$GEN))

Z.G <- ZG %*% LG

#Pheno data
Y <- as.matrix(pheno[, -c(1:3)])
# Check fitting
fm <- BME(Y = Y, Z1 = Z.G, nIter = 10000, burnIn = 5000, thin = 2, bs = 50)
#Linear Predictor
ETA <- list(Gen = list(X = Z.G, model = 'BRR'))
dataset <- pheno[, 2:4] #Must Include in the first column the environments
dataset$Rust <- as.numeric(dataset$Rust)
#Check predictive capacities of the model
pm <- BMORS_Env(dataset, testingEnv = 'R19', ETA = ETA, covModel = "BRR", nIter = 10000,
                burnIn = 5000, thin = 2, progressBar = FALSE, digits = 3)
dataset <-length(phenoMaizeToy[, 2:5]) == 1

length(dataset) == 1

data('MaizeToy')
phenoMaizeToy <- phenoMaizeToy[order(phenoMaizeToy$Env, phenoMaizeToy$Line),]
#Matrix design
LG <- cholesky(genoMaizeToy)
ZG <- model.matrix(~0 + as.factor(phenoMaizeToy$Line))
Z.G <- ZG %*% LG
length(ETA)
#Linear Predictor
ETA <- list(Gen = list(X = Z.G, model = 'BL'))
dataset <- phenoMaizeToy[, 2:5] #Must Include in the first column the environments
#Check predictive capacities of the model
pm <- BMORS_Env(dataset, testingEnv = c("KTI",'EBU'), ETA = ETA, covModel = 'BL', nIter = 10000,
                burnIn = 5000, thin = 2, progressBar = FALSE, digits = 3)
names(pm)
length(dataset) == 1
class(dataset)
summary(pm)


#Voy a repetir un GROAN con los dos modelos más usados rrBLUP y GBLUP para cada trait campo y cámara, esta vez con SNPs:
dim(myG_rust)
dim(Pheno_rust_df)
Pheno_rust_df$AUDPC <- as.numeric(Pheno_rust_df$AUDPC)
nds.no_noise.AUDPC<- createNoisyDataset(
  name = 'AUDPC, no noise',
  genotypes = myG_rust, 
  phenotypes = Pheno_rust_df$AUDPC)
res.no_noise.AUDPC <- GROAN.run(nds.no_noise.AUDPC, wb5)

Pheno_rust_df$IF <- as.numeric(Pheno_rust_df$IF)
nds.no_noise.IF<- createNoisyDataset(
  name = 'IF, no noise',
  genotypes = myG_rust, 
  phenotypes = Pheno_rust_df$IF)
res.no_noise.IF <- GROAN.run(nds.no_noise.IF, wb5)

Pheno_rust_df$IT <- as.numeric(Pheno_rust_df$IT)
nds.no_noise.IT<- createNoisyDataset(
  name = 'IT, no noise',
  genotypes = myG_rust, 
  phenotypes = Pheno_rust_df$IT)
res.no_noise.IT <- GROAN.run(nds.no_noise.IT, wb5)

Pheno_rust_df$DS <- as.numeric(Pheno_rust_df$DS)
nds.no_noise.DS<- createNoisyDataset(
  name = 'DS, no noise',
  genotypes = myG_rust, 
  phenotypes = Pheno_rust_df$DS)
res.no_noise.DS <- GROAN.run(nds.no_noise.DS, wb5)

Pheno_rust_df$R19 <- as.numeric(Pheno_rust_df$R19)
nds.no_noise.R19<- createNoisyDataset(
  name = 'R19, no noise',
  genotypes = myG_rust, 
  phenotypes = Pheno_rust_df$R19)
res.no_noise.R19 <- GROAN.run(nds.no_noise.R19, wb5)

Pheno_rust_df$I_cc_FAI <- as.numeric(Pheno_rust_df$I_cc_FAI)
nds.no_noise.I_cc_FAI<- createNoisyDataset(
  name = 'I_cc_FAI, no noise',
  genotypes = myG_rust, 
  phenotypes = Pheno_rust_df$I_cc_FAI)
res.no_noise.I_cc_FAI <- GROAN.run(nds.no_noise.I_cc_FAI, wb5)

res.total3 <- rbind(res.no_noise.AUDPC, res.no_noise.DS, res.no_noise.I_cc_FAI, res.no_noise.IF, res.no_noise.IT,
                    res.no_noise.R19)
write.xlsx(res.total3, "results/field&CC_traits_twomarkers.xlsx")

plotResult(res.total3)
res.total3 %>%
  group_by(dataset.train, dataset.test, regressor) %>%
  summarise("meanPA" = mean(pearson))

#Cual es el mejor predictor sobre DS2019 de los datos de cámara y campo con los marcadores SNP?
res.AUDPC = GROAN.run(nds = nds.no_noise.AUDPC, wb = wb5, nds.test = nds.no_noise.R19)
print(res.AUDPC[,c('dataset.train', 'dataset.test', 'pearson')])

res.IF = GROAN.run(nds = nds.no_noise.IF, wb = wb5, nds.test = nds.no_noise.R19)
print(res.IF[,c('dataset.train', 'dataset.test', 'pearson')])

res.IT = GROAN.run(nds = nds.no_noise.IT, wb = wb5, nds.test = nds.no_noise.R19)
print(res.IT[,c('dataset.train', 'dataset.test', 'pearson')])

res.DS = GROAN.run(nds = nds.no_noise.DS, wb = wb5, nds.test = nds.no_noise.R19)
print(res.DS[,c('dataset.train', 'dataset.test', 'pearson')])

res.I_cc_FAI = GROAN.run(nds = nds.no_noise.I_cc_FAI, wb = wb5, nds.test = nds.no_noise.R19)
print(res.I_cc_FAI[,c('dataset.train', 'dataset.test', 'pearson')])

predictorsSNP <- rbind(res.AUDPC,res.IF,res.IT,res.DS,res.I_cc_FAI)
predictorsSNP$pearson <- abs(predictorsSNP$pearson)
head(predictorsSNP)
plotResult(predictorsSNP, variable = "pearson")
predictorsSNP %>%
  group_by(dataset.train, dataset.test) %>%
  summarise("meanPA" = mean(pearson))

write.xlsx(predictorsSNP, "results/field&CC_traits_predictorsSNP.xlsx")

#Cual es el mejor predictor sobre DS2019 de los datos de cámara y campo con los marcadores DArT?
res.AUDPC = GROAN.run(nds = nds.no_noise.AUDPC, wb = wb5, nds.test = nds.no_noise.R19)
print(res.AUDPC[,c('dataset.train', 'dataset.test', 'pearson')])

res.IF = GROAN.run(nds = nds.no_noise.IF, wb = wb5, nds.test = nds.no_noise.R19)
print(res.IF[,c('dataset.train', 'dataset.test', 'pearson')])

res.IT = GROAN.run(nds = nds.no_noise.IT, wb = wb5, nds.test = nds.no_noise.R19)
print(res.IT[,c('dataset.train', 'dataset.test', 'pearson')])

res.DS = GROAN.run(nds = nds.no_noise.DS, wb = wb5, nds.test = nds.no_noise.R19)
print(res.DS[,c('dataset.train', 'dataset.test', 'pearson')])

res.I_cc_FAI = GROAN.run(nds = nds.no_noise.I_cc_FAI, wb = wb5, nds.test = nds.no_noise.R19)
print(res.I_cc_FAI[,c('dataset.train', 'dataset.test', 'pearson')])

predictorsDArT <- rbind(res.AUDPC,res.IF,res.IT,res.DS,res.I_cc_FAI)
predictorsDArT$pearson <- abs(predictorsDArT$pearson)
head(predictorsDArT)
plotResult(predictorsDArT, variable = "pearson")
predictorsDArT %>%
  group_by(dataset.train, dataset.test) %>%
  summarise("meanPA" = mean(pearson))

write.xlsx(predictorsDArT, "results/field&CC_traits_predictorsDArT.xlsx")

#Hasta ahora tengo los resultados del escenario [CV1], rehacer en el escenario [CV2]: new env new lines:
  # Define the training (90 % = 224 genotypes) and validation (10 % = 96 genotypes) populations
train <- as.matrix(sample(1:320, 288))
test <- setdiff(1:320, train)      

# Pheno_train and m_train are the phenotype and marker matrices for the values in the training population
# Pheno_valid and m_valid will be the validation populations
Pheno_train = Pheno_rust_df[train, ]
m_train = DArT_rust_SVDI[train, ]

Pheno_valid = Pheno_rust_df[test, ]
m_valid = DArT_rust_SVDI[test, ]

traits = 1
cycles = 500

AUDPC_R19_accuracy = matrix(nrow = cycles, ncol = traits)
for (r in 1:cycles) {
  train = as.matrix(sample(1:320, 288))
  test = setdiff(1:320, train)
  Pheno_train = Pheno_rust_df[train, ]
  m_train = DArT_rust_SVDI[train, ]
  Pheno_valid = Pheno_rust_df[test, ]
  m_valid = DArT_rust_SVDI[test, ]

  AUDPC = (Pheno_train[, 5])
  AUDPC_answer <- mixed.solve(AUDPC, Z = m_train, K = NULL, SE = F, return.Hinv = F)
  u = AUDPC_answer$u
  e = as.matrix(u)
  pred_AUDPC_valid = m_valid %*% e
  pred_AUDPC = pred_AUDPC_valid[, 1] + AUDPC_answer$beta  
  pred_AUDPC
  AUDPC_R19_valid = as.numeric(Pheno_valid[, 3])
  AUDPC_R19_accuracy[r, 1] <-  cor(pred_AUDPC_valid, AUDPC_R19_valid, use = "complete")
}
mean(AUDPC_R19_accuracy)

DS_R19_accuracy = matrix(nrow = cycles, ncol = traits)
for (r in 1:cycles) {
  train = as.matrix(sample(1:320, 288))
  test = setdiff(1:320, train)
  Pheno_train = Pheno_rust_df[train, ]
  m_train = DArT_rust_SVDI[train, ]
  Pheno_valid = Pheno_rust_df[test, ]
  m_valid = DArT_rust_SVDI[test, ]
  
  DS = (Pheno_train[, 9])
  DS_answer <- mixed.solve(DS, Z = m_train, K = NULL, SE = F, return.Hinv = F)
  u = DS_answer$u
  e = as.matrix(u)
  pred_DS_valid = m_valid %*% e
  pred_DS = pred_DS_valid[, 1] + DS_answer$beta  
  pred_DS
  DS_R19_valid = as.numeric(Pheno_valid[, 3])
  DS_R19_accuracy[r, 1] <-  cor(pred_DS_valid, DS_R19_valid, use = "complete")
}
mean(DS_R19_accuracy)

IF_R19_accuracy = matrix(nrow = cycles, ncol = traits)
for (r in 1:cycles) {
  train = as.matrix(sample(1:320, 288))
  test = setdiff(1:320, train)
  Pheno_train = Pheno_rust_df[train, ]
  m_train = DArT_rust_SVDI[train, ]
  Pheno_valid = Pheno_rust_df[test, ]
  m_valid = DArT_rust_SVDI[test, ]
  
  IF = (Pheno_train[, 7])
  IF_answer <- mixed.solve(IF, Z = m_train, K = NULL, SE = F, return.Hinv = F)
  u = IF_answer$u
  e = as.matrix(u)
  pred_IF_valid = m_valid %*% e
  pred_IF = pred_IF_valid[, 1] + IF_answer$beta  
  pred_IF
  IF_R19_valid = as.numeric(Pheno_valid[, 3])
  IF_R19_accuracy[r, 1] <-  cor(pred_IF_valid, IF_R19_valid, use = "complete")
}
mean(IF_R19_accuracy)

IT_R19_accuracy = matrix(nrow = cycles, ncol = traits)
for (r in 1:cycles) {
  train = as.matrix(sample(1:320, 288))
  test = setdiff(1:320, train)
  Pheno_train = Pheno_rust_df[train, ]
  m_train = DArT_rust_SVDI[train, ]
  Pheno_valid = Pheno_rust_df[test, ]
  m_valid = DArT_rust_SVDI[test, ]
  
  IT = (Pheno_train[, 8])
  IT_answer <- mixed.solve(IT, Z = m_train, K = NULL, SE = F, return.Hinv = F)
  u = IT_answer$u
  e = as.matrix(u)
  pred_IT_valid = m_valid %*% e
  pred_IT = pred_IT_valid[, 1] + IT_answer$beta  
  pred_IT
  IT_R19_valid = as.numeric(Pheno_valid[, 3])
  IT_R19_accuracy[r, 1] <-  cor(pred_IT_valid, IT_R19_valid, use = "complete")
}
mean(IT_R19_accuracy)

Index_R19_accuracy = matrix(nrow = cycles, ncol = traits)
for (r in 1:cycles) {
  train = as.matrix(sample(1:320, 288))
  test = setdiff(1:320, train)
  Pheno_train = Pheno_rust_df[train, ]
  m_train = DArT_rust_SVDI[train, ]
  Pheno_valid = Pheno_rust_df[test, ]
  m_valid = DArT_rust_SVDI[test, ]
  
  Index = (Pheno_train[, 14])
  Index_answer <- mixed.solve(Index, Z = m_train, K = NULL, SE = F, return.Hinv = F)
  u = Index_answer$u
  e = as.matrix(u)
  pred_Index_valid = m_valid %*% e
  pred_Index = pred_Index_valid[, 1] + Index_answer$beta  
  pred_Index
  Index_R19_valid = as.numeric(Pheno_valid[, 3])
  Index_R19_accuracy[r, 1] <-  cor(pred_Index_valid, Index_R19_valid, use = "complete")
}
mean(Index_R19_accuracy)

CV2_DArT_rrBLUP <- cbind(AUDPC_R19_accuracy,DS_R19_accuracy,IF_R19_accuracy,IT_R19_accuracy,Index_R19_accuracy, )
CV2_DArT_rrBLUP <- data.frame(CV2_DArT_rrBLUP)
write.xlsx(x = CV2_DArT_rrBLUP, file = "results/CV2_DArT_rrBLUP.xlsx", overwrite = T)
