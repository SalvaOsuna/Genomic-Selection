# 1) Load required packages and load files:
  library(openxlsx)
  library(rrBLUP)

  #Input Genotype. Markers in a matrix format. rows = GEN; column = marker
  Markers <- as.matrix(read.table("data/DarTrrBLUP.txt", sep= "\t" ,head = T)) #here the original data N = 0 as het. (No NA's corrected)
  DArT_SVD <- as.matrix(read.table("data/DArT_noNA_SVDmethod.txt", sep= "\t")) #Load data with NA values corrected trough SVDImput method
  DArT_SVD <- t(DArT_SVD)

  #Input phenotype. Traits in a matrix format. rows = GEN; column = trait
    Pheno <- as.matrix(read.xlsx(xlsxFile = "BLUP_field.xlsx", sep= "\t", rowNames = T, colNames = T, sheet = "BLUP_field"))
    head(Pheno)
  
      #Check Matrix dimensions
      dim(Markers)
      dim(DArT_SVD)
      dim(Pheno)  
  
# 2) Impute missing markers using A.mat() en DArT:
      impute = A.mat(dat_M,max.missing = 0.5,impute.method = "mean", return.imputed = T)
      Markers_impute = impute$imputed
      
      Markers_impute <- t(Markers_impute)
      dim(Markers_impute)
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
        m_train_svd = DArT_SVD[train, ] #Same with NA's corrected
      
      Pheno_valid = Pheno[test, ]
      m_valid = Markers[test, ]
        m_valid_svd = DArT_SVD[test, ] #Same with NA's corrected

# 4) Run mixed.solve() and determine accuracy of predictions
      GDD_F = (Pheno_train[, 1])
      GDD_F_answer <- mixed.solve(GDD_F, Z = m_train, K = NULL, SE = FALSE, return.Hinv = FALSE)
      
      Flowering = as.matrix(GDD_F_answer$u) #output of the marker effects
      head(Flowering) # shows the marker effects for the first six markers
      
      pred_flowering_valid = m_valid %*% Flowering #m_valid * flowering = marker validation matrix times the marker effects
      
      pred_flowering = pred_flowering_valid[, 1] + GDD_F_answer$beta # 
      pred_flowering # pred. flower. based on the markere ffects of the training population with the grand mean added in
      
      #Determine Correlation Accuracy
       #Correlation between the predicted flowering values and the observed flowering values
       flowering_valid = Pheno_valid[, 1]
       FLW_accuracy = cor(pred_flowering_valid, flowering_valid, use = "complete")
       FLW_accuracy
       boxplot(flowering_valid, pred_flowering)

#Test the other traits
       #Poding: 
       GDD_P = Pheno_train[, 2]
       GDD_P_answer <- mixed.solve(GDD_P, Z = m_train, K = NULL, SE = F, return.Hinv = F)
       Poding = GDD_P_answer$u
       pred_GDD_P_valid = m_valid %*% Poding
       GDD_P_valid = Pheno_valid[, 2]
       GDD_P_accuracy <- cor(pred_GDD_P_valid, GDD_P_valid, use = "complete")
       GDD_P_accuracy
       plot(GDD_P_valid, pred_GDD_P_valid)
       
       #Rust (in field):
        #NO NA's corrected:
        Rust_field = Pheno_train[, 3]
        Rust_field_answer <- mixed.solve(Rust_field, Z = m_train, K = NULL, SE = F, return.Hinv = F)
        Rust = Rust_field_answer$u
        pred_Rust_field_valid = m_valid %*% Rust
        Rust_field_valid = Pheno_valid[, 3]
        Rust_field_accuracy <- cor(pred_Rust_field_valid, Rust_field_valid, use = "complete")
        Rust_field_accuracy
        plot(Rust_field_valid, pred_Rust_field_valid)
        
        #NA markers corrected:
        Rust_field = Pheno_train[, 3]
        Rust_field_answer <- mixed.solve(Rust_field, Z = m_train_svd, K = NULL, SE = F, return.Hinv = F)
        Rust = Rust_field_answer$u
        pred_Rust_field_valid = m_valid_svd %*% Rust
        Rust_field_valid = Pheno_valid[, 3]
        Rust_field_accuracy <- cor(pred_Rust_field_valid, Rust_field_valid, use = "complete")
        Rust_field_accuracy
        plot(Rust_field_valid, pred_Rust_field_valid)
       
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
        #with NA's corrected:
        traits = 1
        cycles = 500
        RUST_accuracy_svd = matrix(nrow = cycles, ncol = traits)
        for (r in 1:cycles) {
         train = as.matrix(sample(1:325,195))
         test = setdiff(1:325, train)
         Pheno_train = Pheno[train, ]
         m_train_svd = DArT_SVD[train, ]
         Pheno_valid = Pheno[test, ]
         m_valid_svd = DArT_SVD[test, ]
         
         Rust_field = (Pheno_train[, 3])
         Rust_field_answer_svd <- mixed.solve(Rust_field, Z = m_train_svd, K = NULL, SE = F, return.Hinv = F)
         Rust_svd = Rust_field_answer_svd$u
         e = as.matrix(Rust_svd)
         pred_Rust_valid_svd = m_valid_svd %*% e
         pred_Rust_svd = pred_Rust_valid_svd[, 1] + Rust_field_answer_svd$beta  
         pred_Rust_svd
         Rust_valid_svd = Pheno_valid[, 3]
         RUST_accuracy_svd[r, 1] <-  cor(pred_Rust_valid_svd, Rust_valid_svd, use = "complete")
       }
        mean(RUST_accuracy_svd$V1)
        boxplot(data.frame(x = RUST_accuracy, y = RUST_accuracy_svd))
          
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
        Flowering = GDD_F_answer$u
        e = as.matrix(Flowering)
        pred_flowering_valid = m_valid %*% e
        pred_flowering = pred_flowering_valid[, 1] + GDD_F_answer$beta  
        pred_flowering
        flowering_valid = Pheno_valid[, 1]
        FLR_accuracy[r, 1] <-  cor(pred_flowering_valid, flowering_valid, use = "complete")
       }
       mean(FLR_accuracy)
        #with NA's corrected:
        traits = 1
        cycles = 500
        FLR_accuracy_svd = matrix(nrow = cycles, ncol = traits)
        for (r in 1:cycles) {
         train = as.matrix(sample(1:325,195))
         test = setdiff(1:325, train)
         Pheno_train = Pheno[train, ]
         m_train_svd = DArT_SVD[train, ]
         Pheno_valid = Pheno[test, ]
         m_valid_svd = DArT_SVD[test, ]
         
         GDD_F = (Pheno_train[, 1])
         GDD_F_answer_svd <- mixed.solve(GDD_F, Z = m_train_svd, K = NULL, SE = F, return.Hinv = F)
         GDD_F_svd = GDD_F_answer_svd$u
         e = as.matrix(GDD_F_svd)
         pred_GDD_F_valid_svd = m_valid_svd %*% e
         pred_GDD_F_svd = pred_GDD_F_valid_svd[, 1] + GDD_F_answer_svd$beta  
         pred_GDD_F_svd
         GDD_F_valid_svd = Pheno_valid[, 1]
         FLR_accuracy_svd[r, 1] <-  cor(pred_GDD_F_valid_svd, GDD_F_valid_svd, use = "complete")
       }
        mean(FLR_accuracy_svd$V1)
        boxplot(FLR_accuracy, FLR_accuracy_svd)
       
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
        #with NA's corrected:
        traits = 1
        cycles = 500
        POD_accuracy_svd = matrix(nrow = cycles, ncol = traits)
        for (r in 1:cycles) {
         train = as.matrix(sample(1:325,195))
         test = setdiff(1:325, train)
         Pheno_train = Pheno[train, ]
         m_train_svd = DArT_SVD[train, ]
         Pheno_valid = Pheno[test, ]
         m_valid_svd = DArT_SVD[test, ]
         
         GDD_P = (Pheno_train[, 2])
         GDD_P_answer_svd <- mixed.solve(GDD_P, Z = m_train_svd, K = NULL, SE = F, return.Hinv = F)
         GDD_P_svd = GDD_P_answer_svd$u
         e = as.matrix(GDD_P_svd)
         pred_GDD_P_valid_svd = m_valid_svd %*% e
         pred_GDD_P_svd = pred_GDD_P_valid_svd[, 1] + GDD_P_answer_svd$beta  
         pred_GDD_P_svd
         GDD_P_valid_svd = Pheno_valid[, 2]
         POD_accuracy_svd[r, 1] <-  cor(pred_GDD_P_valid_svd, GDD_P_valid_svd, use = "complete")
       }
        mean(POD_accuracy_svd)
        POD_accuracy_svd <- as.data.frame(POD_accuracy_svd)
       
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
        #with NA's corrected:
        traits = 1
        cycles = 500
        BMS_accuracy_svd = matrix(nrow = cycles, ncol = traits)
        for (r in 1:cycles) {
         train = as.matrix(sample(1:325,195))
         test = setdiff(1:325, train)
         Pheno_train = Pheno[train, ]
         m_train_svd = DArT_SVD[train, ]
         Pheno_valid = Pheno[test, ]
         m_valid_svd = DArT_SVD[test, ]
         
         biomass = (Pheno_train[, 4])
         biomass_answer_svd <- mixed.solve(biomass, Z = m_train_svd, K = NULL, SE = F, return.Hinv = F)
         biomass_svd = biomass_answer_svd$u
         e = as.matrix(biomass_svd)
         pred_biomass_valid_svd = m_valid_svd %*% e
         pred_biomass_svd = pred_biomass_valid_svd[, 1] + biomass_answer_svd$beta  
         pred_biomass_svd
         biomass_valid_svd = Pheno_valid[, 4]
         BMS_accuracy_svd[r, 1] <-  cor(pred_biomass_valid_svd, biomass_valid_svd, use = "complete")
       }
        mean(BMS_accuracy_svd)
       
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
        #with NA's corrected:
        traits = 1
        cycles = 500
        YLD_accuracy_svd = matrix(nrow = cycles, ncol = traits)
        for (r in 1:cycles) {
         train = as.matrix(sample(1:325,195))
         test = setdiff(1:325, train)
         Pheno_train = Pheno[train, ]
         m_train_svd = DArT_SVD[train, ]
         Pheno_valid = Pheno[test, ]
         m_valid_svd = DArT_SVD[test, ]
         
         yield = (Pheno_train[, 5])
         yield_answer_svd <- mixed.solve(yield, Z = m_train_svd, K = NULL, SE = F, return.Hinv = F)
         yield_svd = yield_answer_svd$u
         e = as.matrix(yield_svd)
         pred_yield_valid_svd = m_valid_svd %*% e
         pred_yield_svd = pred_yield_valid_svd[, 1] + yield_answer_svd$beta  
         pred_yield_svd
         yield_valid_svd = Pheno_valid[, 5]
         YLD_accuracy_svd[r, 1] <-  cor(pred_yield_valid_svd, yield_valid_svd, use = "complete")
       }
        mean(YLD_accuracy_svd)
        
#save and write accuracies with 500 itinerations in two imput conditions:
df.accuracies.svd <- data.frame(FLR_accuracy_svd,POD_accuracy_svd,RUST_accuracy_svd,BMS_accuracy_svd, YLD_accuracy_svd)
names(df.accuracies.svd) <- c("FLR.svd", "POD.svd", "RUST.svd", "BMS.svd", "YLD.svd")
boxplot(df.accuracies.svd)
df.accuracies <- data.frame(FLR_accuracy,POD_accuracy,RUST_accuracy,BMS_accuracy, YLD_accuracy)
names(df.accuracies) <- c("FLR", "POD", "RUST", "BMS", "YLD")
boxplot(df.accuracies)
write.xlsx(df.accuracies, "df.acuracies.xlsx", sep = "/t")
write.xlsx(df.accuracies.svd, "df.acuraciessvd.xlsx", sep = "/t")
df.accuracies <- read.xlsx("df.acuraciessvd.xlsx", colNames = T)
df.accuracies <- df.accuracies %>%
   rename("Impute_Condition" = "Impute Condition")

ggplot(df.accuracies, aes(x = Trait, y = Accuracy))+
  geom_boxplot(aes(fill = Impute_Condition)) +
  coord_flip() +
  labs(title = "Accuracy by Impute Condition", 
       subtitle = "With 500 itinerances in five traits",
       caption = "SVD Impute reveals better results in Yield") +
  scale_fill_discrete(labels = c("MNI", "N = het.", "SVDI")) +
  theme_bw()


#Voy a probar sustituyendo los NA ahora con la media (el imput que me deja rrBLUP) para cada trait con itinerancia de 500
  # y añadir los resultados al boxplot:
#Rust:
traits = 1
cycles = 500
RUST_accuracy_m = matrix(nrow = cycles, ncol = traits)
for (r in 1:cycles) {
  train = as.matrix(sample(1:325,195))
  test = setdiff(1:325, train)
  Pheno_train = Pheno[train, ]
  m_train_m = Markers_impute[train, ]
  Pheno_valid = Pheno[test, ]
  m_valid_m = Markers_impute[test, ]
  
  Rust_field = (Pheno_train[, 3])
  Rust_field_answer <- mixed.solve(Rust_field, Z = m_train_m, K = NULL, SE = F, return.Hinv = F)
  Rust = Rust_field_answer$u
  e = as.matrix(Rust)
  pred_Rust_valid = m_valid_m %*% e
  pred_Rust = pred_Rust_valid[, 1] + Rust_field_answer$beta  
  pred_Rust
  Rust_valid = Pheno_valid[, 3]
  RUST_accuracy_m[r, 1] <-  cor(pred_Rust_valid, Rust_valid, use = "complete")
}
mean(RUST_accuracy_m)
#Flowering:
traits = 1
cycles = 500
FLR_accuracy_m = matrix(nrow = cycles, ncol = traits)
for (r in 1:cycles) {
  train = as.matrix(sample(1:325,195))
  test = setdiff(1:325, train)
  Pheno_train = Pheno[train, ]
  m_train_m = Markers_impute[train, ]
  Pheno_valid = Pheno[test, ]
  m_valid_m = Markers_impute[test, ]
  
  GDD_F = (Pheno_train[, 1])
  GDD_F_answer <- mixed.solve(GDD_F, Z = m_train_m, K = NULL, SE = F, return.Hinv = F)
  Flowering = GDD_F_answer$u
  e = as.matrix(Flowering)
  pred_flowering_valid = m_valid_m %*% e
  pred_flowering = pred_flowering_valid[, 1] + GDD_F_answer$beta  
  pred_flowering
  flowering_valid = Pheno_valid[, 1]
  FLR_accuracy_m[r, 1] <-  cor(pred_flowering_valid, flowering_valid, use = "complete")
}
mean(FLR_accuracy_m)
#Poding:
traits = 1
cycles = 500
POD_accuracy_m = matrix(nrow = cycles, ncol = traits)
for (r in 1:cycles) {
  train = as.matrix(sample(1:325,195))
  test = setdiff(1:325, train)
  Pheno_train = Pheno[train, ]
  m_train_m = Markers_impute[train, ]
  Pheno_valid = Pheno[test, ]
  m_valid_m = Markers_impute[test, ]
  
  GDD_P = (Pheno_train[, 2])
  GDD_P_answer <- mixed.solve(GDD_P, Z = m_train_m, K = NULL, SE = F, return.Hinv = F)
  Poding = GDD_P_answer$u
  e = as.matrix(Poding)
  pred_poding_valid = m_valid_m %*% e
  pred_poding = pred_poding_valid[, 1] + GDD_P_answer$beta  
  pred_poding
  poding_valid = Pheno_valid[, 2]
  POD_accuracy_m[r, 1] <-  cor(pred_poding_valid, poding_valid, use = "complete")
}
mean(POD_accuracy_m)
#Biomass:
traits = 1
cycles = 500
BMS_accuracy_m = matrix(nrow = cycles, ncol = traits)
for (r in 1:cycles) {
  train = as.matrix(sample(1:325,195))
  test = setdiff(1:325, train)
  Pheno_train = Pheno[train, ]
  m_train_m = Markers_impute[train, ]
  Pheno_valid = Pheno[test, ]
  m_valid_m = Markers_impute[test, ]
  
  biomass = (Pheno_train[, 4])
  biomass_answer <- mixed.solve(biomass, Z = m_train_m, K = NULL, SE = F, return.Hinv = F)
  BMS = biomass_answer$u
  e = as.matrix(BMS)
  pred_biomass_valid = m_valid_m %*% e
  pred_biomass = pred_biomass_valid[, 1] + biomass_answer$beta  
  pred_biomass
  biomass_valid = Pheno_valid[, 4]
  BMS_accuracy_m[r, 1] <-  cor(pred_biomass_valid, biomass_valid, use = "complete")
}
mean_m(BMS_accuracy)
#Yield:
traits = 1
cycles = 500
YLD_accuracy_m = matrix(nrow = cycles, ncol = traits)
for (r in 1:cycles) {
  train = as.matrix(sample(1:325,195))
  test = setdiff(1:325, train)
  Pheno_train = Pheno[train, ]
  m_train_m = Markers_impute[train, ]
  Pheno_valid = Pheno[test, ]
  m_valid_m = Markers_impute[test, ]
  
  yield = (Pheno_train[, 4])
  yield_answer <- mixed.solve(yield, Z = m_train_m, K = NULL, SE = F, return.Hinv = F)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid = m_valid_m %*% e
  pred_yield = pred_yield_valid[, 1] + yield_answer$beta  
  pred_yield
  yield_valid = Pheno_valid[, 4]
  YLD_accuracy_m[r, 1] <-  cor(pred_yield_valid, yield_valid, use = "complete")
}
mean(YLD_accuracy_m)
df.accuracies.m <- data.frame(FLR_accuracy_m,POD_accuracy_m,RUST_accuracy_m,BMS_accuracy_m, YLD_accuracy_m)
names(df.accuracies.m) <- c("FLR.m", "POD.m", "RUST.m", "BMS.m", "YLD.m")
write.xlsx(df.accuracies.m, "df.acuraciesm.xlsx", sep = "/t")


#R code for implementing three default kernels in rrBLUP (linear, Gaussian, and Exponential)
library(rrBLUP) #load rrBLUP
library(BGLR) #load BGLR
    #Voy a hacerlo con todos los traits de roya que tengo, así puedo comparar G-Blup con SNP-blup
X= dat_SVD
dim(X)
#X <- 2*X-1 #recode genotypes
X=scale(X)
Pheno_rust <- as.matrix(read.xlsx(xlsxFile = "data/BLUP_field.xlsx", sep= "\t", rowNames = T, colNames = T, sheet = "BLUP_GS_rust"))
head(Pheno_rust)

#R18 
{
t=1
y <-Pheno_rust[,t] #Rust from R18
yy=y
MSE=function(yobserved,ypredicted){
  MSE=mean((yobserved-ypredicted)^2)
}
n_records=nrow(X)
n_folds=10
set.seed(10)
sets <- findInterval(cut(sample(1:n_records, n_records),
                         breaks=n_folds), 1:n_records)
results=data.frame()
for (i in 1:n_folds){
  # i=1
  trn <- which(sets!=i)
  tst <- which(sets==i)
  ans.RR<-kinship.BLUP(y=y[trn],
                       G.train=X[trn,],G.pred=X[tst,])
  #accuracy with RR
  Cor_RR=cor(ans.RR$g.pred,yy[tst])
  MSE_RR=MSE(ans.RR$g.pred,yy[tst])
  ans.GAUSS<-kinship.BLUP(y=y[trn],
                          G.train=X[trn,],G.pred=X[tst,],
                          K.method="GAUSS")
  #accuracy with GAUSS
  Cor_GAUSS=cor(ans.GAUSS$g.pred,yy[tst])
  MSE_GAUSS=MSE(ans.GAUSS$g.pred,yy[tst])
  ans.EXP<-kinship.BLUP(y=y[trn],
                        G.train=X[trn,],G.pred=X[tst,],
                        K.method="EXP")
  #accuracy with EXponential
  Cor_EXP=cor(ans.EXP$g.pred,yy[tst])
  MSE_EXP=MSE(ans.EXP$g.pred,yy[tst])
  results=rbind(results,data.frame(Fold=i, MSE_RR=MSE_RR,
                                   MSE_GAUSS=MSE_GAUSS, MSE_EXP=MSE_EXP,Cor_RR=Cor_RR,
                                   Cor_GAUSS=Cor_GAUSS, Cor_EXP=Cor_EXP))
}
results
write.xlsx(results,file = "Kernel_Mixed_Example_R18.xlsx", sep = "/t")
}
#R19
{
  t=2
  y <-Pheno_rust[,t] #Rust from R19
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds){
    # i=1
    trn <- which(sets!=i)
    tst <- which(sets==i)
    ans.RR<-kinship.BLUP(y=y[trn],
                         G.train=X[trn,],G.pred=X[tst,])
    #accuracy with RR
    Cor_RR=cor(ans.RR$g.pred,yy[tst])
    MSE_RR=MSE(ans.RR$g.pred,yy[tst])
    ans.GAUSS<-kinship.BLUP(y=y[trn],
                            G.train=X[trn,],G.pred=X[tst,],
                            K.method="GAUSS")
    #accuracy with GAUSS
    Cor_GAUSS=cor(ans.GAUSS$g.pred,yy[tst])
    MSE_GAUSS=MSE(ans.GAUSS$g.pred,yy[tst])
    ans.EXP<-kinship.BLUP(y=y[trn],
                          G.train=X[trn,],G.pred=X[tst,],
                          K.method="EXP")
    #accuracy with EXponential
    Cor_EXP=cor(ans.EXP$g.pred,yy[tst])
    MSE_EXP=MSE(ans.EXP$g.pred,yy[tst])
    results=rbind(results,data.frame(Fold=i, MSE_RR=MSE_RR,
                                     MSE_GAUSS=MSE_GAUSS, MSE_EXP=MSE_EXP,Cor_RR=Cor_RR,
                                     Cor_GAUSS=Cor_GAUSS, Cor_EXP=Cor_EXP))
  }
  results
  write.xlsx(results,file = "data/Kernel_Mixed_Example_R19.xlsx", sep = "/t")
}
#R20
{
  t=3
  y <-Pheno_rust[,t] #Rust from R20
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds){
    # i=1
    trn <- which(sets!=i)
    tst <- which(sets==i)
    ans.RR<-kinship.BLUP(y=y[trn],
                         G.train=X[trn,],G.pred=X[tst,])
    #accuracy with RR
    Cor_RR=cor(ans.RR$g.pred,yy[tst])
    MSE_RR=MSE(ans.RR$g.pred,yy[tst])
    ans.GAUSS<-kinship.BLUP(y=y[trn],
                            G.train=X[trn,],G.pred=X[tst,],
                            K.method="GAUSS")
    #accuracy with GAUSS
    Cor_GAUSS=cor(ans.GAUSS$g.pred,yy[tst])
    MSE_GAUSS=MSE(ans.GAUSS$g.pred,yy[tst])
    ans.EXP<-kinship.BLUP(y=y[trn],
                          G.train=X[trn,],G.pred=X[tst,],
                          K.method="EXP")
    #accuracy with EXponential
    Cor_EXP=cor(ans.EXP$g.pred,yy[tst])
    MSE_EXP=MSE(ans.EXP$g.pred,yy[tst])
    results=rbind(results,data.frame(Fold=i, MSE_RR=MSE_RR,
                                     MSE_GAUSS=MSE_GAUSS, MSE_EXP=MSE_EXP,Cor_RR=Cor_RR,
                                     Cor_GAUSS=Cor_GAUSS, Cor_EXP=Cor_EXP))
  }
  results
  write.xlsx(results,file = "Kernel_Mixed_Example_R20.xlsx", sep = "/t")
}
#AUDPC
{
  t=4
  y <-Pheno_rust[,t] #Rust from AUDPC
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds){
    # i=1
    trn <- which(sets!=i)
    tst <- which(sets==i)
    ans.RR<-kinship.BLUP(y=y[trn],
                         G.train=X[trn,],G.pred=X[tst,])
    #accuracy with RR
    Cor_RR=cor(ans.RR$g.pred,yy[tst])
    MSE_RR=MSE(ans.RR$g.pred,yy[tst])
    ans.GAUSS<-kinship.BLUP(y=y[trn],
                            G.train=X[trn,],G.pred=X[tst,],
                            K.method="GAUSS")
    #accuracy with GAUSS
    Cor_GAUSS=cor(ans.GAUSS$g.pred,yy[tst])
    MSE_GAUSS=MSE(ans.GAUSS$g.pred,yy[tst])
    ans.EXP<-kinship.BLUP(y=y[trn],
                          G.train=X[trn,],G.pred=X[tst,],
                          K.method="EXP")
    #accuracy with EXponential
    Cor_EXP=cor(ans.EXP$g.pred,yy[tst])
    MSE_EXP=MSE(ans.EXP$g.pred,yy[tst])
    results=rbind(results,data.frame(Fold=i, MSE_RR=MSE_RR,
                                     MSE_GAUSS=MSE_GAUSS, MSE_EXP=MSE_EXP,Cor_RR=Cor_RR,
                                     Cor_GAUSS=Cor_GAUSS, Cor_EXP=Cor_EXP))
  }
  results
  write.xlsx(results,file = "Kernel_Mixed_Example_AUDPC.xlsx", sep = "/t")
}
#LP50
{
  t=5
  y <-Pheno_rust[,t] #Rust from LP50
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds){
    # i=1
    trn <- which(sets!=i)
    tst <- which(sets==i)
    ans.RR<-kinship.BLUP(y=y[trn],
                         G.train=X[trn,],G.pred=X[tst,])
    #accuracy with RR
    Cor_RR=cor(ans.RR$g.pred,yy[tst])
    MSE_RR=MSE(ans.RR$g.pred,yy[tst])
    ans.GAUSS<-kinship.BLUP(y=y[trn],
                            G.train=X[trn,],G.pred=X[tst,],
                            K.method="GAUSS")
    #accuracy with GAUSS
    Cor_GAUSS=cor(ans.GAUSS$g.pred,yy[tst])
    MSE_GAUSS=MSE(ans.GAUSS$g.pred,yy[tst])
    ans.EXP<-kinship.BLUP(y=y[trn],
                          G.train=X[trn,],G.pred=X[tst,],
                          K.method="EXP")
    #accuracy with EXponential
    Cor_EXP=cor(ans.EXP$g.pred,yy[tst])
    MSE_EXP=MSE(ans.EXP$g.pred,yy[tst])
    results=rbind(results,data.frame(Fold=i, MSE_RR=MSE_RR,
                                     MSE_GAUSS=MSE_GAUSS, MSE_EXP=MSE_EXP,Cor_RR=Cor_RR,
                                     Cor_GAUSS=Cor_GAUSS, Cor_EXP=Cor_EXP))
  }
  results
  write.xlsx(results,file = "Kernel_Mixed_Example_LP50.xlsx", sep = "/t")
}
#IF
{
  t=6
  y <-Pheno_rust[,t] #Rust from IF
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds){
    # i=1
    trn <- which(sets!=i)
    tst <- which(sets==i)
    ans.RR<-kinship.BLUP(y=y[trn],
                         G.train=X[trn,],G.pred=X[tst,])
    #accuracy with RR
    Cor_RR=cor(ans.RR$g.pred,yy[tst])
    MSE_RR=MSE(ans.RR$g.pred,yy[tst])
    ans.GAUSS<-kinship.BLUP(y=y[trn],
                            G.train=X[trn,],G.pred=X[tst,],
                            K.method="GAUSS")
    #accuracy with GAUSS
    Cor_GAUSS=cor(ans.GAUSS$g.pred,yy[tst])
    MSE_GAUSS=MSE(ans.GAUSS$g.pred,yy[tst])
    ans.EXP<-kinship.BLUP(y=y[trn],
                          G.train=X[trn,],G.pred=X[tst,],
                          K.method="EXP")
    #accuracy with EXponential
    Cor_EXP=cor(ans.EXP$g.pred,yy[tst])
    MSE_EXP=MSE(ans.EXP$g.pred,yy[tst])
    results=rbind(results,data.frame(Fold=i, MSE_RR=MSE_RR,
                                     MSE_GAUSS=MSE_GAUSS, MSE_EXP=MSE_EXP,Cor_RR=Cor_RR,
                                     Cor_GAUSS=Cor_GAUSS, Cor_EXP=Cor_EXP))
  }
  results
  write.xlsx(results,file = "Kernel_Mixed_Example_IF.xlsx", sep = "/t")
}
#IT
{
  t=7
  y <-Pheno_rust[,t] #Rust from IT
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds){
    # i=1
    trn <- which(sets!=i)
    tst <- which(sets==i)
    ans.RR<-kinship.BLUP(y=y[trn],
                         G.train=X[trn,],G.pred=X[tst,])
    #accuracy with RR
    Cor_RR=cor(ans.RR$g.pred,yy[tst])
    MSE_RR=MSE(ans.RR$g.pred,yy[tst])
    ans.GAUSS<-kinship.BLUP(y=y[trn],
                            G.train=X[trn,],G.pred=X[tst,],
                            K.method="GAUSS")
    #accuracy with GAUSS
    Cor_GAUSS=cor(ans.GAUSS$g.pred,yy[tst])
    MSE_GAUSS=MSE(ans.GAUSS$g.pred,yy[tst])
    ans.EXP<-kinship.BLUP(y=y[trn],
                          G.train=X[trn,],G.pred=X[tst,],
                          K.method="EXP")
    #accuracy with EXponential
    Cor_EXP=cor(ans.EXP$g.pred,yy[tst])
    MSE_EXP=MSE(ans.EXP$g.pred,yy[tst])
    results=rbind(results,data.frame(Fold=i, MSE_RR=MSE_RR,
                                     MSE_GAUSS=MSE_GAUSS, MSE_EXP=MSE_EXP,Cor_RR=Cor_RR,
                                     Cor_GAUSS=Cor_GAUSS, Cor_EXP=Cor_EXP))
  }
  results
  write.xlsx(results,file = "Kernel_Mixed_Example_IT.xlsx", sep = "/t")
}
#DS
{
  t=8
  y <-Pheno_rust[,t] #Rust from DS
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds){
    # i=1
    trn <- which(sets!=i)
    tst <- which(sets==i)
    ans.RR<-kinship.BLUP(y=y[trn],
                         G.train=X[trn,],G.pred=X[tst,])
    #accuracy with RR
    Cor_RR=cor(ans.RR$g.pred,yy[tst])
    MSE_RR=MSE(ans.RR$g.pred,yy[tst])
    ans.GAUSS<-kinship.BLUP(y=y[trn],
                            G.train=X[trn,],G.pred=X[tst,],
                            K.method="GAUSS")
    #accuracy with GAUSS
    Cor_GAUSS=cor(ans.GAUSS$g.pred,yy[tst])
    MSE_GAUSS=MSE(ans.GAUSS$g.pred,yy[tst])
    ans.EXP<-kinship.BLUP(y=y[trn],
                          G.train=X[trn,],G.pred=X[tst,],
                          K.method="EXP")
    #accuracy with EXponential
    Cor_EXP=cor(ans.EXP$g.pred,yy[tst])
    MSE_EXP=MSE(ans.EXP$g.pred,yy[tst])
    results=rbind(results,data.frame(Fold=i, MSE_RR=MSE_RR,
                                     MSE_GAUSS=MSE_GAUSS, MSE_EXP=MSE_EXP,Cor_RR=Cor_RR,
                                     Cor_GAUSS=Cor_GAUSS, Cor_EXP=Cor_EXP))
  }
  results
  write.xlsx(results,file = "Kernel_Mixed_Example_DS.xlsx", sep = "/t")
}

  #visualize results:
  GBLUP_bp <- read.xlsx("results/Kernel_Mixed_GBLUP.xlsx", sheet = "bp", colNames = T)
  library(tidyverse)  
  ggplot(GBLUP_bp, aes(x = Trait, y = PC)) +
    geom_boxplot(aes(fill = Model)) +
    coord_flip() +
    theme_classic() +
    labs(fill = "Kernel Model",
         title = "Genomic BLUP accuracies by 10 k-fold cross-validation",
         subtitle = "Accuracy by Rust traits",
         caption = "R19, DS and AUDPC are goods")
  

  GBLUP_bp$Fold = as.character(GBLUP_bp$Fold)  
  ggplot(GBLUP_bp, aes( x = Fold, y = PC))+
    geom_boxplot(aes(fill = Model))+
    coord_flip() +
    theme_classic() +
    labs(fill = "Kernel Model",
         title = "Accuracy in G-BLUP model by Folds",
         subtitle = "10 k-folds cross-validation",
         caption = "Fold 4 always shows worst results, why?")
  
# SNP BLUP vs G BLUP in rust traits:
  BLUP_bp <- read.xlsx("results/RR_SNP&G_BLUP_rust_traits.xlsx", sheet = "Sheet 1", colNames = T)
  ggplot(BLUP_bp, aes( x = Trait, y = Accuracy))+
    geom_boxplot(aes(fill = Model))+
    coord_flip() +
    theme_classic() +
    labs(fill = "BLUP Model",
         title = "SNP-BLUP vs G-BLUP accuracies in Rust traits",
         subtitle = "SNP with 500 itinerancies, G-BLUP 10 k-folds",
         caption = "G-BLUP kinship calculated with Linear Kernel")
  
#He hecho SNP-blup con 500 itinerancias y G BLUP con cross-validation. Voy a probar SNP-BLUP con cross validation:
  {
  #R18:
  X= DArT_rust_SVDI
  t = 1
  n_folds=10
  y <-Pheno_rust[,t]
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds) {
    
    trn <- which(sets!=i)
    tst <- which(sets==i)
    
    ans.R18 <- mixed.solve(y[trn], Z = X[trn,], K = NULL, SE = F, return.Hinv = F)
    u = ans.R18$u
    e = as.matrix(u)
    pred_R18_valid = X[tst,] %*% e
    pred_R18 = pred_R18_valid[, 1] + ans.R18$beta  
    pred_R18
    
    Cor_RR=cor(pred_R18,yy[tst])
    MSE_RR=MSE(pred_R18,yy[tst])
    
    results=rbind(results,data.frame(Fold=i,
                                     Cor_RR=Cor_RR,
                                     MSE_RR=MSE_RR
    ))
  }
  results
  write.xlsx(results,file = "results/Cross_R18.xlsx", sep = "/t")
  
  #R19:
  X= DArT_rust_SVDI
  t = 2
  n_folds=10
  y <-Pheno_rust[,t]
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds) {
    
    trn <- which(sets!=i)
    tst <- which(sets==i)
    
    ans.R18 <- mixed.solve(y[trn], Z = X[trn,], K = NULL, SE = F, return.Hinv = F)
    u = ans.R18$u
    e = as.matrix(u)
    pred_R18_valid = X[tst,] %*% e
    pred_R18 = pred_R18_valid[, 1] + ans.R18$beta  
    pred_R18
    
    Cor_RR=cor(pred_R18,yy[tst])
    MSE_RR=MSE(pred_R18,yy[tst])
    
    results=rbind(results,data.frame(Fold=i,
                                     Cor_RR=Cor_RR,
                                     MSE_RR=MSE_RR
    ))
  }
  results
  write.xlsx(results,file = "results/Cross_R19.xlsx", sep = "/t")
  
  #R20:
  X= DArT_rust_SVDI
  t = 3
  n_folds=10
  y <-Pheno_rust[,t]
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds) {
    
    trn <- which(sets!=i)
    tst <- which(sets==i)
    
    ans.R18 <- mixed.solve(y[trn], Z = X[trn,], K = NULL, SE = F, return.Hinv = F)
    u = ans.R18$u
    e = as.matrix(u)
    pred_R18_valid = X[tst,] %*% e
    pred_R18 = pred_R18_valid[, 1] + ans.R18$beta  
    pred_R18
    
    Cor_RR=cor(pred_R18,yy[tst])
    MSE_RR=MSE(pred_R18,yy[tst])
    
    results=rbind(results,data.frame(Fold=i,
                                     Cor_RR=Cor_RR,
                                     MSE_RR=MSE_RR
    ))
  }
  results
  write.xlsx(results,file = "results/Cross_R20.xlsx", sep = "/t")
  
  #AUDPC:
  X= DArT_rust_SVDI
  t = 4
  n_folds=10
  y <-Pheno_rust[,t]
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds) {
    
    trn <- which(sets!=i)
    tst <- which(sets==i)
    
    ans.R18 <- mixed.solve(y[trn], Z = X[trn,], K = NULL, SE = F, return.Hinv = F)
    u = ans.R18$u
    e = as.matrix(u)
    pred_R18_valid = X[tst,] %*% e
    pred_R18 = pred_R18_valid[, 1] + ans.R18$beta  
    pred_R18
    
    Cor_RR=cor(pred_R18,yy[tst])
    MSE_RR=MSE(pred_R18,yy[tst])
    
    results=rbind(results,data.frame(Fold=i,
                                     Cor_RR=Cor_RR,
                                     MSE_RR=MSE_RR
    ))
  }
  results
  write.xlsx(results,file = "results/Cross_AUDPC.xlsx", sep = "/t")
  
  #LP50:
  X= DArT_rust_SVDI
  t = 5
  n_folds=10
  y <-Pheno_rust[,t]
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds) {
    
    trn <- which(sets!=i)
    tst <- which(sets==i)
    
    ans.R18 <- mixed.solve(y[trn], Z = X[trn,], K = NULL, SE = F, return.Hinv = F)
    u = ans.R18$u
    e = as.matrix(u)
    pred_R18_valid = X[tst,] %*% e
    pred_R18 = pred_R18_valid[, 1] + ans.R18$beta  
    pred_R18
    
    Cor_RR=cor(pred_R18,yy[tst])
    MSE_RR=MSE(pred_R18,yy[tst])
    
    results=rbind(results,data.frame(Fold=i,
                                     Cor_RR=Cor_RR,
                                     MSE_RR=MSE_RR
    ))
  }
  results
  write.xlsx(results,file = "results/Cross_LP50.xlsx", sep = "/t")
  
  #IF:
  X= DArT_rust_SVDI
  t = 6
  n_folds=10
  y <-Pheno_rust[,t]
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds) {
    
    trn <- which(sets!=i)
    tst <- which(sets==i)
    
    ans.R18 <- mixed.solve(y[trn], Z = X[trn,], K = NULL, SE = F, return.Hinv = F)
    u = ans.R18$u
    e = as.matrix(u)
    pred_R18_valid = X[tst,] %*% e
    pred_R18 = pred_R18_valid[, 1] + ans.R18$beta  
    pred_R18
    
    Cor_RR=cor(pred_R18,yy[tst])
    MSE_RR=MSE(pred_R18,yy[tst])
    
    results=rbind(results,data.frame(Fold=i,
                                     Cor_RR=Cor_RR,
                                     MSE_RR=MSE_RR
    ))
  }
  results
  write.xlsx(results,file = "results/Cross_IF.xlsx", sep = "/t")
  
  #IT:
  X= DArT_rust_SVDI
  t = 7
  n_folds=10
  y <-Pheno_rust[,t]
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds) {
    
    trn <- which(sets!=i)
    tst <- which(sets==i)
    
    ans.R18 <- mixed.solve(y[trn], Z = X[trn,], K = NULL, SE = F, return.Hinv = F)
    u = ans.R18$u
    e = as.matrix(u)
    pred_R18_valid = X[tst,] %*% e
    pred_R18 = pred_R18_valid[, 1] + ans.R18$beta  
    pred_R18
    
    Cor_RR=cor(pred_R18,yy[tst])
    MSE_RR=MSE(pred_R18,yy[tst])
    
    results=rbind(results,data.frame(Fold=i,
                                     Cor_RR=Cor_RR,
                                     MSE_RR=MSE_RR
    ))
  }
  results
  write.xlsx(results,file = "results/Cross_IT.xlsx", sep = "/t")
  
  #DS:
  X= DArT_rust_SVDI
  t = 8
  n_folds=10
  y <-Pheno_rust[,t]
  yy=y
  MSE=function(yobserved,ypredicted){
    MSE=mean((yobserved-ypredicted)^2)
  }
  n_records=nrow(X)
  n_folds=10
  set.seed(10)
  sets <- findInterval(cut(sample(1:n_records, n_records),
                           breaks=n_folds), 1:n_records)
  results=data.frame()
  for (i in 1:n_folds) {
    
    trn <- which(sets!=i)
    tst <- which(sets==i)
    
    ans.R18 <- mixed.solve(y[trn], Z = X[trn,], K = NULL, SE = F, return.Hinv = F)
    u = ans.R18$u
    e = as.matrix(u)
    pred_R18_valid = X[tst,] %*% e
    pred_R18 = pred_R18_valid[, 1] + ans.R18$beta  
    pred_R18
    
    Cor_RR=cor(pred_R18,yy[tst])
    MSE_RR=MSE(pred_R18,yy[tst])
    
    results=rbind(results,data.frame(Fold=i,
                                     Cor_RR=Cor_RR,
                                     MSE_RR=MSE_RR
    ))
  }
  results
  write.xlsx(results,file = "results/Cross_DS.xlsx", sep = "/t")
    }
  
#Vuelvo a representar las accuracy de las tres maneras:
  BLUP_bp <- read.xlsx("results/RR_SNP&G_BLUP_rust_traits.xlsx", sheet = "Sheet 1", colNames = T)
  ggplot(BLUP_bp, aes( x = Trait, y = Accuracy))+
    geom_boxplot(aes(fill = Model))+
    coord_flip() +
    theme_classic() +
    labs(fill = "BLUP Model",
         title = "SNP-BLUP vs G-BLUP accuracies in Rust traits",
         subtitle = "SNP-BLUP calculated by 500 it. and 10 k-folds",
         caption = "G-BLUP kinship calculated with Linear Kernel and 10k-folds")
  