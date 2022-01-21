############################################
##                                        ##
##          GS models in RUST pea         ##
##                                        ##
############################################

# 1) Load Genotypic data, DArT markers:
    DArT <- as.matrix(read.table("data/DArT.txt", header = T))
    DArT_rust <- DArT[-c(288, 294, 300, 320, 325), ] #Estas entradas no están evaluadas en CC así que las quito
    ##rrBLUP necesita la matriz en formato {-1, 0, 1} así que cambio los 0 por -1:

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
    Pheno_rust <- as.matrix(read.xlsx(xlsxFile = "BLUP_field.xlsx", sep= "\t", rowNames = T, colNames = T, sheet = "BLUP_GS_rust"))
    head(Pheno_rust)
    dim(Pheno_rust)
    dim(DArT_rust_SVDI)
  
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
    #R18_t:
    traits = 1
    cycles = 500
    R18_t_accuracy = matrix(nrow = cycles, ncol = traits)
    for (r in 1:cycles) {
      train = as.matrix(sample(1:320, 224))
      test = setdiff(1:320, train)
      Pheno_train = Pheno_rust[train, ]
      m_train = DArT_rust_SVDI[train, ]
      Pheno_valid = Pheno_rust[test, ]
      m_valid = DArT_rust_SVDI[test, ]
      
      R18_t = (Pheno_train[, 9])
      R18_t_answer <- mixed.solve(R18_t, Z = m_train, K = NULL, SE = F, return.Hinv = F)
      u = R18_t_answer$u
      e = as.matrix(u)
      pred_R18_t_valid = m_valid %*% e
      pred_R18_t = pred_R18_t_valid[, 1] + R18_t_answer$beta  
      pred_R18_t
      R18_t_valid = Pheno_valid[, 9]
      R18_t_accuracy[r, 1] <-  cor(pred_R18_t_valid, R18_t_valid, use = "complete")
    }
    mean(R18_t_accuracy)
    
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
#Try cross validation:
    