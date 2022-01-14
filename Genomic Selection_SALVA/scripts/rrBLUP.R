# 1) Load required packages and load files:
  library(openxlsx)
  library(rrBLUP)

  #Input Genotype. Markers in a matrix format. rows = GEN; column = marker
   Markers <- as.matrix(read.table("DarTrrBLUP.txt", sep= "\t" ,head = T))
    head(Markers)

  #Input phenotype. Traits in a matrix format. rows = GEN; column = trait
    Pheno <- as.matrix(read.xlsx(xlsxFile = "BLUP_field.xlsx", sep= "\t", rowNames = T, colNames = T, sheet = "BLUP_field"))
    head(Pheno)
  
      #Check Matrix dimensions
      dim(Markers)
      dim(Pheno)  
  
# 2) Impute missing markers using A.mat()
      impute = A.mat(Markers,max.missing = 0.5,impute.method = "mean", return.imputed = T)
      Markers_impute = impute$imputed
      
      dim(Markers_impute) #There's no NA's in the input markers matrix
      dim(Markers)
       
# 3) Define the training (60 % = 195 genotypes) and validation (40 % = 130 genotypes) populations
      train = as.matrix(sample(1:325, 195))
      head(train)
      
      test <- setdiff(1:325, train)      
      test
      
      # Pheno_train and m_train are the phenotype and marker matrices for the values in the training population
      # Pheno_valid and m_valid will be the validation populations
      Pheno_train = Pheno[train, ]
      m_train = Markers[train, ]
      
      Pheno_valid = Pheno[test, ]
      m_valid = Markers[test, ]

# 4) Run mixed.solve() and determine accuracy of predictions
      GDD_F = (Pheno_train[, 1])
      GDD_F_answer <- mixed.solve(GDD_F, Z = m_train, K = NULL, SE = FALSE, return.Hinv = FALSE)
      
      Flowering = as.matrix(GGD_F_answer$u) #output of the marker effects
      head(Flowering) # shows the marker effects for the first six markers
      
      pred_flowering_valid = m_valid %*% Flowering #m_valid * flowering = marker validation matrix times the marker effects
      
      pred_flowering = pred_flowering_valid[, 1] + GGD_F_answer$beta # 
      pred_flowering # pred. flower. based on the markere ffects of the training population with the grand mean added in
      
      #Determine Correlation Accuracy
       #Correlation between the predicted flowering values and the observed flowering values
       flowering_valid = Pheno_valid[, 1]
       FLW_accuracy = cor(pred_flowering_valid, flowering_valid, use = "complete")
       FLW_accuracy
       
#Test the other traits
       #Poding: 
       GDD_P = Pheno_train[, 2]
       GDD_P_answer <- mixed.solve(GDD_P, Z = m_train, K = NULL, SE = F, return.Hinv = F)
       Poding = GDD_P_answer$u
       pred_GDD_P_valid = m_valid %*% Poding
       GDD_P_valid = Pheno_valid[, 2]
       GDD_P_accuracy <- cor(pred_GDD_P_valid, GDD_P_valid, use = "complete")
       GDD_P_accuracy
       
       #Rust (in field): 
       Rust_field = Pheno_train[, 3]
       Rust_field_answer <- mixed.solve(Rust_field, Z = m_train, K = NULL, SE = F, return.Hinv = F)
       Rust = Rust_field_answer$u
       pred_Rust_field_valid = m_valid %*% Rust
       Rust_field_valid = Pheno_valid[, 3]
       Rust_field_accuracy <- cor(pred_Rust_field_valid, Rust_field_valid, use = "complete")
       Rust_field_accuracy
       
       #Biomass: 
       Biomass_ = Pheno_train[, 4]
       Biomass_answer <- mixed.solve(Biomass_, Z = m_train, K = NULL, SE = F, return.Hinv = F)
       Biomass = Biomass_answer$u
       pred_Biomass_valid = m_valid %*% Biomass
       Biomass_valid = Pheno_valid[, 4]
       Biomass_accuracy <- cor(pred_Biomass_valid, Biomass_valid, use = "complete")
       Biomass_accuracy
       
       #Yield: 
       Yield_ = Pheno_train[, 5]
       Yield_answer <- mixed.solve(Yield_, Z = m_train, K = NULL, SE = F, return.Hinv = F)
       Yield = Yield_answer$u
       pred_Yield_valid = m_valid %*% Yield
       Yield_valid = Pheno_valid[, 5]
       Yield_accuracy <- cor(pred_Yield_valid, Yield_valid, use = "complete")
       Yield_accuracy
       
##Correlation accuracy with 500 iterations
       #Rust:
       traits = 1
       cycles = 500
       RUST_accuracy = matrix(nrow = cycles, ncol = traits)
       for (r in 1:cycles) {
         train = as.matrix(sample(1:325,195))
         test = setdiff(1:325, train)
         Pheno_train = Pheno[train, ]
         m_train = Markers[train, ]
         Pheno_valid = Pheno[test, ]
         m_valid = Markers[test, ]
         
         Rust_field = (Pheno_train[, 3])
         Rust_field_answer <- mixed.solve(Rust_field, Z = m_train, K = NULL, SE = F, return.Hinv = F)
         Rust = Rust_field_answer$u
         e = as.matrix(Rust)
         pred_Rust_valid = m_valid %*% e
         pred_Rust = pred_Rust_valid[, 1] + Rust_field_answer$beta  
         pred_Rust
         Rust_valid = Pheno_valid[, 3]
         RUST_accuracy[r, 1] <-  cor(pred_Rust_valid, Rust_valid, use = "complete")
       }
       mean(RUST_accuracy)
       
       #Flowering:
       traits = 1
       cycles = 500
       FLR_accuracy = matrix(nrow = cycles, ncol = traits)
       for (r in 1:cycles) {
        train = as.matrix(sample(1:325,195))
        test = setdiff(1:325, train)
        Pheno_train = Pheno[train, ]
        m_train = Markers[train, ]
        Pheno_valid = Pheno[test, ]
        m_valid = Markers[test, ]
        
        GDD_F = (Pheno_train[, 1])
        GDD_F_answer <- mixed.solve(GDD_F, Z = m_train, K = NULL, SE = F, return.Hinv = F)
        Flowering = GGD_F_answer$u
        e = as.matrix(Flowering)
        pred_flowering_valid = m_valid %*% e
        pred_flowering = pred_flowering_valid[, 1] + GGD_F_answer$beta  
        pred_flowering
        flowering_valid = Pheno_valid[, 1]
        FLR_accuracy[r, 1] <-  cor(pred_flowering_valid, flowering_valid, use = "complete")
       }
       mean(FLR_accuracy)
       
       #Poding:
       traits = 1
       cycles = 500
       POD_accuracy = matrix(nrow = cycles, ncol = traits)
       for (r in 1:cycles) {
         train = as.matrix(sample(1:325,195))
         test = setdiff(1:325, train)
         Pheno_train = Pheno[train, ]
         m_train = Markers[train, ]
         Pheno_valid = Pheno[test, ]
         m_valid = Markers[test, ]
         
         GDD_P = (Pheno_train[, 2])
         GDD_P_answer <- mixed.solve(GDD_P, Z = m_train, K = NULL, SE = F, return.Hinv = F)
         Poding = GDD_P_answer$u
         e = as.matrix(Poding)
         pred_poding_valid = m_valid %*% e
         pred_poding = pred_poding_valid[, 1] + GDD_P_answer$beta  
         pred_poding
         poding_valid = Pheno_valid[, 2]
         POD_accuracy[r, 1] <-  cor(pred_poding_valid, poding_valid, use = "complete")
       }
       mean(POD_accuracy)
       
       #Biomass:
       traits = 1
       cycles = 500
       BMS_accuracy = matrix(nrow = cycles, ncol = traits)
       for (r in 1:cycles) {
         train = as.matrix(sample(1:325,195))
         test = setdiff(1:325, train)
         Pheno_train = Pheno[train, ]
         m_train = Markers[train, ]
         Pheno_valid = Pheno[test, ]
         m_valid = Markers[test, ]
         
         biomass = (Pheno_train[, 4])
         biomass_answer <- mixed.solve(biomass, Z = m_train, K = NULL, SE = F, return.Hinv = F)
         BMS = biomass_answer$u
         e = as.matrix(BMS)
         pred_biomass_valid = m_valid %*% e
         pred_biomass = pred_biomass_valid[, 1] + biomass_answer$beta  
         pred_biomass
         biomass_valid = Pheno_valid[, 4]
         BMS_accuracy[r, 1] <-  cor(pred_biomass_valid, biomass_valid, use = "complete")
       }
       mean(BMS_accuracy)
       
       #Yield:
       traits = 1
       cycles = 500
       YLD_accuracy = matrix(nrow = cycles, ncol = traits)
       for (r in 1:cycles) {
         train = as.matrix(sample(1:325,195))
         test = setdiff(1:325, train)
         Pheno_train = Pheno[train, ]
         m_train = Markers[train, ]
         Pheno_valid = Pheno[test, ]
         m_valid = Markers[test, ]
         
         yield = (Pheno_train[, 4])
         yield_answer <- mixed.solve(yield, Z = m_train, K = NULL, SE = F, return.Hinv = F)
         YLD = yield_answer$u
         e = as.matrix(YLD)
         pred_yield_valid = m_valid %*% e
         pred_yield = pred_yield_valid[, 1] + yield_answer$beta  
         pred_yield
         yield_valid = Pheno_valid[, 4]
         YLD_accuracy[r, 1] <-  cor(pred_yield_valid, yield_valid, use = "complete")
       }
       mean(YLD_accuracy)
      
      
      
      