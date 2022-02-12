##LOAD DATA##
library(readxl)
library(openxlsx)
library(metan)
library(ggstatsplot)
library(tidyverse)
library(GGEBiplots)
library(multcompView)
  
  #GWAS Pea Panel INFO
  {
  Panel_info <- read_excel("C:/Users/Salva/Documents/GitHub/Genomic-Selection/Genomic Selection_SALVA/data/ListaPeaGWAS.xlsx",
                           sheet = "R",
                           col_types = 
                             c("text", #GEN
                               "text", #Reference
                               "text", #Species
                               "text", #Name
                               "text", #Origen
                               "text", #Material
                               "text"))#Structure
}

  #season 2018
  {
R18_nc <- read_excel("C:/Users/Salva/Documents/GitHub/Genomic-Selection/Genomic Selection_SALVA/data/R18_19_20.xlsx", 
                     sheet = "R18_nc", 
                     col_types = 
                       c("text", #ENV
                         "text", #REP
                         "text", #BLOCK
                         "numeric", #ROW
                         "text", #Specie
                         "text", #Material
                         "text", #Structure
                         "numeric", #GEN
                         "numeric", #GERM
                         "numeric", #DtF
                         "numeric", #DtP
                         "numeric", #GDD_F	
                         "numeric", #GDD_P
                         "numeric", #biomass_pl
                         "numeric", #yield_pl
                         "numeric", #Asco
                         "numeric", #Rust
                         "numeric" #Oidio
                       ))
R18_nc$ENV <- as.factor(R18_nc$ENV)
R18_nc$REP <- as.factor(R18_nc$REP)
R18_nc$BLOCK <- as.factor(R18_nc$BLOCK)
R18_nc$ROW <- as.factor(R18_nc$ROW)
R18_nc$Material <- as.factor(R18_nc$Material)
R18_nc$Species <- as.factor(R18_nc$Species)
R18_nc$Structure <- as.factor(R18_nc$Structure)
R18_nc$GEN <- as.factor(R18_nc$GEN)
  }

  #season 2019
  {
R19_nc <- read_excel("C:/Users/Salva/Documents/GitHub/Genomic-Selection/Genomic Selection_SALVA/data/R18_19_20.xlsx", 
                     sheet = "R19_nc", 
                     col_types = 
                       c("text", #ENV
                         "text", #REP
                         "text", #BLOCK
                         "numeric", #ROW
                         "text", #Specie
                         "text", #Material
                         "text", #Structure
                         "text", #GEN
                         "numeric", #GERM
                         "numeric", #GDD_FF
                         "numeric", #GDD_F
                         "numeric", #GDD_FP	
                         "numeric", #GDD_P
                         "numeric", #GDD_M
                         "numeric", #Rust
                         "numeric", #Biom_pl_g
                         "numeric" #yield_pl_g
                       ))
R19_nc$ENV <- as.factor(R19_nc$ENV)
R19_nc$REP <- as.factor(R19_nc$REP)
R19_nc$BLOCK <- as.factor(R19_nc$BLOCK)
R19_nc$ROW <- as.factor(R19_nc$ROW)
R19_nc$Material <- as.factor(R19_nc$Material)
R19_nc$Species <- as.factor(R19_nc$Species)
R19_nc$Structure <- as.factor(R19_nc$Structure)
R19_nc$GEN <- as.factor(R19_nc$GEN)
}

  #season 2020
  {
R20_nc <- #NoControls
  read_excel("C:/Users/Salva/Documents/GitHub/Genomic-Selection/Genomic Selection_SALVA/data/R18_19_20.xlsx", 
                     sheet = "R20_nc", 
                     col_types = 
                       c("text", #ENV
                         "text", #REP
                         "text", #BLOCK
                         "numeric", #ROW
                         "text", #Specie
                         "text", #Material
                         "text", #Structure
                         "text", #GEN
                         "numeric", #GERM
                         "numeric", #GDD_FF
                         "numeric", #GDD_F
                         "numeric", #GDD_FP	
                         "numeric", #GDD_P
                         "numeric", #GDD_M
                         "numeric", #Rust
                         "numeric", #Asco
                         "numeric", #Oidio
                         "numeric", #Biom_pl_g
                         "numeric" #yield_pl_g
                       )) 
R20_nc$ENV <- as.factor(R20_nc$ENV)
R20_nc$REP <- as.factor(R20_nc$REP)
R20_nc$BLOCK <- as.factor(R20_nc$BLOCK)
R20_nc$ROW <- as.factor(R20_nc$ROW)
R20_nc$Material <- as.factor(R20_nc$Material)
R20_nc$Species <- as.factor(R20_nc$Species)
R20_nc$Structure <- as.factor(R20_nc$Structure)
R20_nc$GEN <- as.factor(R20_nc$GEN)



}
  
  #Growth Chamber. Controlled conditions:
  {
  RCC_nc <- #NoControls
    read_excel("C:/Users/Salva/Documents/GitHub/Rust-collection/data/Rust_controlled_condition.xlsx", 
               sheet = "CC_n", 
               col_types = 
                 c("text", #REP
                   "text", #BLOCK
                   "text", #GEN
                   "numeric", #L_area						
                   "numeric", #L_aphila (0 = NO aphila; 1 = SI aphila)
                   "numeric", #L_border (0 = Liso; 0,5 = Semi; 1 = dentado)
                   "numeric", #AUDPC_n	
                   "numeric", #PL50_n
                   "numeric", #IF_n
                   "numeric", #DS_n
                   "numeric", #IT_n
                   "numeric", #AUDPC_T square root
                   "numeric", #PL50_T log
                   "numeric", #IF_T square root
                   "numeric" #DS_T arc.sin
                 ))
  
  RCC_nc$REP <- as.factor(RCC_nc$REP)
  RCC_nc$BLOCK <- as.factor(RCC_nc$BLOCK)
  RCC_nc$GEN <- as.factor(RCC_nc$GEN)
  }

  #Three field Seasons together. Data from multi-environment trials
  {
Rall_field <- #NoControls
  read_excel("C:/Users/Salva/Documents/GitHub/Genomic-Selection/Genomic Selection_SALVA/data/R18_19_20.xlsx", 
             sheet = "R_all", 
             col_types = 
               c("text", #ENV
                 "text", #REP
                 "text", #BLOCK
                 "numeric", #ROW
                 "text", #Specie
                 "text", #Material
                 "text", #Structure
                 "text", #GEN
                 "numeric", #GERM
                 "numeric", #GDD_F
                 "numeric", #GDD_P
                 "numeric", #Biom_pl_g
                 "numeric", #yield_pl_g
                 "numeric", #Rust
                 "numeric", #Asco
                 "numeric" #Oidio
               )) 
Rall_field$ENV <- as.factor(Rall_field$ENV)
Rall_field$REP <- as.factor(Rall_field$REP)
Rall_field$BLOCK <- as.factor(Rall_field$BLOCK)
Rall_field$ROW <- as.factor(Rall_field$ROW)
Rall_field$Material <- as.factor(Rall_field$Material)
Rall_field$Species <- as.factor(Rall_field$Species)
Rall_field$Structure <- as.factor(Rall_field$Structure)
Rall_field$GEN <- as.factor(Rall_field$GEN)
}
res_ind <- waasb(Rall_field,
                 env = ENV,
                 gen = GEN,
                 rep = REP,
                 resp = Rust,
                 block = BLOCK,
                 mresp = "l",
                 verbose = FALSE)
res_ind$Rust$ESTIMATES
model_indexes <- blup_indexes(res_ind)
model_indexes$Rust
BLUP_field <- gmd(model_indexes)
write.xlsx(BLUP_field, "BLUP_field2.xlsx")

  #Four ENVs with Rust (%) arc.sin transformed
  {
    Rall_T <- #NoControls
    read_excel("C:/Users/Salva/Documents/GitHub/Genomic-Selection/Genomic Selection_SALVA/data/R18_19_20.xlsx", 
               sheet = "RustT", 
               col_types = 
                 c("text", #ENV
                   "text", #REP
                   "text", #BLOCK
                   "text", #GEN
                   "numeric" #RustT
                 ))
  }
      #Distribution
      D18T <-   Rall_T %>%      
         filter(ENV == "R18")
      ggplot(D18T, aes(x = RustT))+
        geom_histogram(bins = 9,color="black", fill="white") +
        theme_classic()
      ggplot(BLUP18T, aes(x = BLUPg))+
        geom_histogram(bins = 10,color="black", fill="white") +
        theme_classic()
      ggplot(BLUP19T, aes(x = BLUPg))+
        geom_histogram(bins = 10,color="black", fill="white") +
        theme_classic()
      ggplot(BLUP20T, aes(x = BLUPg))+
        geom_histogram(bins = 10,color="black", fill="white") +
        theme_classic()
      ggplot(BLUPCCT, aes(x = BLUPg))+
        geom_histogram(bins = 10,color="black", fill="white") +
        theme_classic()
      ggplot(Rall_T, aes(x = RustT))+
        geom_histogram(bins = 10,color="black", fill="white") +
        theme_classic() +
        facet_grid(~ ENV)
      


#Analysis of single environment/experiment using mixed-models

  #Season 2018:
  {
    stat18 <- desc_stat(R18_nc, stats = "all")
    
gen_alphaR18 <- 
  gamem(R18_nc, 
        GEN, 
        REP, 
        resp = c("GERM", "DtF", "DtP", "GDD_F", "GDD_P", "Biom_plant_g", "yield_plant_g", "Asco", "Rust", "Rust_T", "Oidio"), 
        BLOCK)
    
    BLUP18 <- get_model_data(gen_alphaR18, what = "blupg")
    #write.xlsx(BLUP18, "BLUP18.xlsx", sep = "/t")
    get_model_data(gen_alphaR18)
    get_model_data(gen_alphaR18, "lrt")
    get_model_data(gen_alphaR18, "details")
}

  #Season 2019:
  {
    stat19 <- desc_stat(R19_nc, stats = "all")
gen_alphaR19 <- 
  gamem(R19_nc, 
        GEN, 
        REP, 
        resp = c("GERM", "GDD_FF", "GDD_F", "GDD_FP", "GDD_P","GDD_M", "Rust", "Biom_plant_g", "yield_plant_g"), 
        BLOCK)
    
    BLUP19 <- get_model_data(gen_alphaR19, what = "blupg")
    #write.xlsx(BLUP19, "BLUP19.xlsx", sep = "/t")
    get_model_data(gen_alphaR19)
    get_model_data(gen_alphaR19, "lrt")
    get_model_data(gen_alphaR19, "details")
}

  #Season 2020:
  {
    stat20 <- desc_stat(R20_nc, stats = "all")
gen_alphaR20 <- 
  gamem(R20_nc, 
        GEN, 
        REP, 
        resp = c("GERM", "GDD_DTFF", "GDD_DTF", "GDD_DTFP", "GDD_DTP", "GDD_M", "Rust", "Asco", "Oidio", "Biom_plant_g", "yield_plant_g"), 
        BLOCK)
    
    BLUP20 <- get_model_data(gen_alphaR20, what = "blupg")
    get_model_data(gen_alphaR20)
    get_model_data(gen_alphaR20, "lrt")
    get_model_data(gen_alphaR20, "details")
}

  #Growth Chamber. Controlled conditions:
  {
    statCC <- desc_stat(RCC_nc, stats = "all")
    
    gen_alphaRCC <- 
      gamem(RCC_nc, 
            GEN, 
            REP, 
            resp = c("AUDPC_T", "PL50_T", "IF_T", "IT_n", "DS_T"), 
            BLOCK)
    
  gen_alphaRCC_noLP <- 
  gamem(RCC_nc, 
        GEN, 
        REP, 
        resp = c("AUDPC_n", "IF_n", "IT_n", "DS_n"))
    
    get_model_data(gen_alphaRCC_noLP)
    gen_alphaRCC$AUDPC_T$BLUPgen
    BLUPcc <- get_model_data(gen_alphaRCC, what = "blupg")
    #write.xlsx(BLUPcc, file = "BLUPcc.xlsx", sep = "/t")
    BLUPcc1 <- read.xlsx("BLUPcc.xlsx")
    BLUPcc_noLP <- get_model_data(gen_alphaRCC_noLP, what = "blupg")
    corr_plot(RCC_nc, c("AUDPC_T", "PL50_T", "IF_T", "DS_T", "IT_n"))
    corr_plot(BLUPcc, c("AUDPC_T", "PL50_T", "IF_T", "DS_T", "IT_n"))
    
    get_model_data(gen_alphaRCC)
    get_model_data(gen_alphaRCC, "lrt")
    
    b18.19.20.CCplus <- inner_join(b18.19.20, BLUPcc, by = "GEN")
    corr_plot(b18.19.20.CCplus)
  
    #Genotypes Selection for rust resistance:
    #FAI-BLUP:
    FAI <- fai_blup(
      gen_alphaRCC,
       use_data = "blup",
        DI = c("min", "max", "min", "min", "min"),
        UI = c("max", "min", "max", "max", "max"),
        SI = 15,
        mineval = 1,
       verbose = TRUE)

    FAI_noLP <- fai_blup(
      gen_alphaRCC_noLP,
      use_data = "blup",
      DI = c("max", "max", "max", "max"),
      UI = c("min","min", "min", "min"),
      SI = 15,
      mineval = 1,
      verbose = TRUE)

      FAI_selection <- FAI$FAI
      FAI_selection1 <- FAI_noLP$FAI
      FAI_noLP$FAI
      #write.xlsx(FAI_noLP$FAI, "FAI_noLP$FAI.xlsx", sep = "/t")

    #MGIDI index selection:  

      mgidi_index <- mgidi(gen_alphaRCC,
                          ideotype = c("l", "h", "l", "l","l"),
                     SI = 10) # Selection intensity
      mgidi_index_noLP <- mgidi(gen_alphaRCC_noLP,
                           ideotype = c("l", "l", "l","l"),
                           SI = 15) # Selection intensity
      MGIDI_selection <- mgidi_index$MGIDI
      #write.xlsx(MGIDI_selection, "MGIDI_selection.xlsx", sep = "/t")
      MGIDI_selection1 <- mgidi_index_noLP$MGIDI
      #write.xlsx(MGIDI_selection1, "MGIDI_selection1.xlsx", sep = "/t")
      
    #Smith-Hazel index selection:
      
     covcor1 <-  covcor_design(RCC_nc, 
                    gen = GEN, 
                    rep = REP , 
                    resp = c("AUDPC_n", "DS_n", "IF_n", "IT_n"), 
                    design = "RCBD")
      
      SH <- Smith_Hazel(
        gen_alphaRCC,
        use_data = "blup",
        pcov = covcor1$phen_cov,
        gcov = covcor1$geno_cov,
        SI = 15,
        weights = rep(5, 1)
      )
      SH_selection <- SH$index
      #write.xlsx(SH_selection, "SH_selection.xlsx", sep = "/t")
      
      #Comparison of the genotype selected:
      coincidence <- coincidence_index(mgidi_index_noLP, FAI_noLP, SH, total = 13)
      coincidence
      
}

  #Linear model with interaction effect used to analyze data for multi-environment trials. Three seasons:
  {
    statFIELD <- desc_stat(Rall_field, stats = "all")
    #BLUPs
mixed_mod_Rall_field <- 
  gamem_met(Rall_field,
        ENV,
        GEN, 
        REP, 
        resp = c("GDD_F", "GDD_P", "Rust", "Biom_plant_g", "yield_plant_g"),
        block = BLOCK
        )
  BLUP_flowering <- mixed_mod_Rall_field$GDD_F$BLUPgen
  write.xlsx(BLUP_flowering, "BLUP_flowering.xlsx")
  BLUP_poding <- mixed_mod_Rall_field$GDD_P$BLUPgen
  write.xlsx(BLUP_poding, "BLUP_poding.xlsx")
  BLUP_rust <- mixed_mod_Rall_field$Rust$BLUPgen
  write.xlsx(BLUP_rust, "BLUP_rust.xlsx")
  BLUP_biomass <- mixed_mod_Rall_field$Biom_plant_g$BLUPgen
  write.xlsx(BLUP_biomass, "BLUP_biomass.xlsx")
  BLUP_yield <- mixed_mod_Rall_field$yield_plant_g$BLUPgen
  write.xlsx(BLUP_yield, "BLUP_yield.xlsx")
  
  get_model_data(mixed_mod_Rall_field)
  mixed_mod_Rall <- 
    gamem_met(Rall_T,
              ENV,
              GEN, 
              REP, 
              resp = "RustT",
              block = BLOCK
    )
  
  aksjdñ <- mixed_mod_Rall$RustT$BLUPgen %>%
    arrange(Predicted)
  Rust_GxE <- mixed_mod_Rall$RustT$MeansGxE
  mixed_mod_Rall$RustT$BLUPint
  mixed_mod_Rall$RustT$ESTIMATES
  mu <- ddply(Rust_GxE, "ENV", summarise, grp.mean=mean(Y))
  head(mu)
  ggplot(Rust_GxE, aes(x = Y))+
    geom_histogram(bins = 18,color="black", fill="white") +
    theme_classic() +
    facet_grid(~ ENV)+geom_vline(data=mu, aes(xintercept=grp.mean, color="red"),
                                 linetype="dashed") +
    theme(legend.position = "none",
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  
  mu1 <- ddply(BLUPcc1, "Trait", summarise, grp.mean=mean(Value))
  head(mu1)
  ggplot(BLUPcc1, aes(x = Value_nc))+
    geom_histogram(bins = 18,color="black", fill="white") +
    theme_classic() +
    facet_grid(~ Trait, scales = "free")+geom_vline(data=mu1, aes(xintercept=grp.mean, color="red"),
                                 linetype="dashed") +
    theme(legend.position = "none",
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  
  
  
  #write.xlsx(BLUP_all, "BLUP_all.xlsx", sep = "/t")

    #The WAASB object
waasb_model <- 
  waasb(Rall_field,
        env = ENV,
        gen = GEN,
        rep = REP,
        resp = c("GERM", "GDD_F", "GDD_P", "Rust", "Biom_plant_g", "yield_plant_g"),
        random = "gen", #Default
        verbose = TRUE) #Default

  data <- waasb_model$Rust$PCA
  print(data)
  
  data1 <- waasb_model$Rust$MeansGxE
  print(data1)
  
  WAASB_rust <- waasb_model$Rust$BLUPgen
  
  e <- plot_waasby(waasb_model)
  f <- plot_waasby(waasb_model, col.shape = c("gray20", "gray80"))
  arrange_ggplot(e, f, tag_levels = list(c("e", "f")))
  
  i <- plot_scores(waasb_model, type = 4)
  j <- plot_scores(waasb_model,
                   type = 4,
                   size.tex.gen = 1.5,
                   color = FALSE,
                   col.alpha.gen = 0,
                   col.alpha.env = 0,
                   plot_theme = theme_metan(color.background = "white"))
  arrange_ggplot(i, j, tag_levels = list(c("i", "j")))
  
  details_Rall <- get_model_data(mixed_mod_Rall_field, "genpar")
  get_model_data(mixed_mod_Rall_field, "lrt")
  get_model_data(mixed_mod_Rall_field, "details")
  get_model_data(mixed_mod_Rall_field, "estimates")
  #write.xlsx(BLUPgen_all, "BLUPgen_all.xlsx", sep = "/t")
  BLUPint_all <- mixed_mod_Rall_field$Rust$BLUPint
  mixed_mod_Rall_field$Rust$MeansGxE
  
  #Genotype-by-environment means
  make_mat(Rall_field, GEN, ENV, Rust) %>% print(rownames = TRUE)
  
  #Genotype-by-environment interaction effects
  ge_ef <- ge_effects(Rall_field, ENV, GEN, Rust)
  print(ge_ef$Rust)
  plot(ge_ef)
  
  #Genotype plus Genotype-by-environment effects
  gge_ef <- ge_effects(Rall_field, ENV, GEN, Rust, type = "gge")
  print(gge_ef$Rust)
  GGE_rust <- gge_ef$Rust
  #write.xlsx(GGE_rust, "GGE_rust.xlsx", sep = "/t")
  
  #superiority index
  super <- superiority(Rall_field, ENV, GEN, Rust)
  print(super$Rust$index)
  
  #environmental estratification
  fact <- ge_factanal(Rall_field, ENV, GEN, REP, Rust)
  plot(fact)
  
  #ranking genotypes
  gge_model <- gge(Rall_field, ENV, GEN, Rust)
  as_tibble(gge_model$Rust$ge_mat)
  #ploting
  o <- plot(gge_model, type = 8)
  p <- plot(gge_model,
            type = 8,
            col.gen = "black",
            col.env = "gray",
            size.text.gen = 6,
            plot_theme = theme_metan_minimal())
  arrange_ggplot(o, p, tag_levels = list(c("o", "p")))
  
  residualsR18_rust <- plot(gen_alphaR18,  which = c(2,5), labels = TRUE, var = 9)
  residualsR19_rust <-plot(gen_alphaR19,  which = c(2,5), labels = TRUE, var = 7)
  residualsR20_rust <-plot(gen_alphaR20,  which = c(2,5), labels = TRUE, var = 7)
  
  #Multi-trait stability index
    #Lin e Binns' superiority index
      out <- superiority(Rall_field, ENV, GEN, Rust)

  
    #Annicchiarico's genotypic confidence index
      Ann <- Annicchiarico(Rall_field,
                           env = ENV,
                           gen = GEN,
                           rep = REP,
                           resp = Rust)
      Ann_selection <- Ann$Rust$general
      #write.xlsx(Ann_selection, "Ann_selection.xlsx", sep = "/t")
  }
     #Same but with Rust % arc.sin transformed:
      {
      R18T <- filter(Rall_T, ENV == "R18") %>%
          gamem(GEN,
                REP,
                resp = "RustT",
                BLOCK)
      get_model_data(R18T)
      
      R19T <- filter(Rall_T, ENV == "R19") %>%
        gamem(GEN,
              REP,
              resp = "RustT",
              BLOCK)
      get_model_data(R19T)
      
      R20T <- filter(Rall_T, ENV == "R20") %>%
        gamem(GEN,
              REP,
              resp = "RustT",
              BLOCK)
      get_model_data(R20T)
      
      RCCT <- filter(Rall_T, ENV == "CC") %>%
        gamem(GEN,
              REP,
              resp = "RustT",
              BLOCK)
      get_model_data(RCCT)
      
      BLUP18T <- R18T$RustT$BLUPgen %>%
        select(GEN, BLUPg)
      BLUP19T <- R19T$RustT$BLUPgen %>%
        select(GEN, BLUPg)
      BLUP20T <- R20T$RustT$BLUPgen %>%
        select(GEN, BLUPg)
      BLUPCCT <- RCCT$RustT$BLUPgen %>%
        select(GEN, BLUPg)
      
      b18.19 <- inner_join(BLUP18T, BLUP19T, by = "GEN")
      b18.19.20 <- inner_join(b18.19, BLUP20T, by = "GEN")
      b18.19.20.CC <- inner_join(b18.19.20, BLUPCCT, by = "GEN")
      
      corr_plot(b18.19.20.CC)
      
      }

#GGE MODELS

gge_model <- gge(Rall_field, ENV, GEN, Biom_plant_g)

plot(gge_model, type = 8)

plot(gge_model,
     type = 9,
     sel_gen1 = "220",
     sel_gen2 = "246",
     col.gen = "black",
     title = FALSE)

plot(gge_model, type = 7, sel_gen = "59")

  #BLUPs correlations:
  {
BLUP18_19 <- left_join(BLUP18,BLUP19, by = "GEN", suffix = c(".18", ".19"))
BLUP18_19_20 <- left_join(BLUP18_19,BLUP20, by = "GEN")
BLUP18_19_20_all <- left_join(BLUP18_19_20,BLUP_all, by = "GEN")
BLUPALL <- left_join(BLUP18_19_20_all,BLUPcc, by = "GEN")
BLUPALL_info <- left_join(BLUPALL,Panel_info, by = "GEN")
#write.xlsx(BLUPALL_info, "BLUP_ALL_info.xlsx", sep = "/t")

summary(BLUPALL_info)
corr_plot(BLUPALL, c("Rust.18", "Rust.19", "Rust","Predicted","AUDPC_n", "PL50_n", "IF_n", "DS_n", "IT_n" ))
ggcorrmat(
  data     = BLUPALL,
  cor.vars = c("Rust.18", "Rust.19", "Rust","Predicted","AUDPC_n", "PL50_n", "IF_n", "DS_n", "IT_n" ),
  colors   = c("#cbac43", "white", "#550000"),
  title    = "Correlalogram",
  subtitle = "This is the subtitle",
  type = "robust",
  p.adjust.method = "bonferroni"
)

ggcorrmat(
  data     = BLUPALL,
  cor.vars = c("Rust.18", "Rust.19", "Rust","Predicted","AUDPC_n", "PL50_n", "IF_n", "DS_n", "IT_n" ),
  colors   = c("#cbac43", "white", "#550000"),
  title    = "Correlalogram",
  subtitle = "This is the subtitle",
  type = "np",
  p.adjust.method = "BY"
)

}



#Correlations
  #2018

corr_plot(BLUPALL_info, c("Rust.18", "Rust.19", "Rust","Predicted","AUDPC_n", "PL50_n", "IF_n", "DS_n", "IT_n" ),
          col.by = Material)

  #2019
material19 <- select(R19_nc, c(5,6,7))
material19 <- material19[1:322,1:3]

BLUP19 <- as.data.frame(get_model_data(gen_alphaR19, "blupg"))
BLUP19 <-cbind(BLUP19, material19)
corr_plot(BLUP19, GDD_F, GDD_P,GDD_M, Biom_plant_g, yield_plant_g, Rust, col.by = Material)

#2020
material20 <- select(R20_nc, c(5,6,7))
material20 <- material20[1:324,1:3]

BLUP20 <- as.data.frame(get_model_data(gen_alphaR20, "blupg"))
BLUP20 <-cbind(BLUP20, material20)
corr_plot(BLUP20, GDD_DTF, GDD_DTP,GDD_M, Biom_plant_g, yield_plant_g, Rust,Asco, Oidio, col.by = Material)

#Boxplots
  #2018
ggbetweenstats(x = Material, y = Rust.18, data = BLUPALL_info, 
               p.adjust.method = "holm", type = "p",
               bf.message = FALSE, var.equal = F,
               ggsignif.args = list(textsize = 1.5, tip_length = 0.01)) +
  theme(text = element_text(size = 8), plot.subtitle = element_text(size=8))

ggbetweenstats(x = Species, y = AUDPC_n, data = BLUPALL_info, 
               p.adjust.method = "none", type = "np",
               bf.message = FALSE, var.equal = F,
               ggsignif.args = list(textsize = 1.5, tip_length = 0.01)) +
  theme(text = element_text(size = 8), plot.subtitle = element_text(size=8))

#2019
ggbetweenstats(x = Material, y = Rust, data = BLUP19, 
               p.adjust.method = "none", type = "np",
               bf.message = FALSE, var.equal = F,
               ggsignif.args = list(textsize = 1.5, tip_length = 0.01)) +
  theme(text = element_text(size = 8), plot.subtitle = element_text(size=8))

ggbetweenstats(x = Structure, y = Rust, data = BLUP19, 
               p.adjust.method = "none", type = "np",
               bf.message = FALSE, var.equal = F,
               ggsignif.args = list(textsize = 1.5, tip_length = 0.01)) +
  theme(text = element_text(size = 8), plot.subtitle = element_text(size=8))

#2020
ggbetweenstats(x = Material, y = Rust, data = BLUP20, 
               p.adjust.method = "none", type = "np",
               bf.message = FALSE, var.equal = F,
               ggsignif.args = list(textsize = 1.5, tip_length = 0.01)) +
  theme(text = element_text(size = 8), plot.subtitle = element_text(size=8))

ggbetweenstats(x = Species, y = BLUP, data = BLUPs_anova, 
               p.adjust.method = "none", type = "np",
               bf.message = FALSE, var.equal = F,
               ggsignif.args = list(textsize = 1.5, tip_length = 0.01)) +
  theme(text = element_text(size = 8), plot.subtitle = element_text(size=8)) +
  coord_flip()


#intento de pca
blups.pca <- prcomp(BLUPALL_info[,c(40:42)], center = TRUE,scale. = TRUE)

summary(blups.pca)
ggbiplot(blups.pca, groups = BLUPALL_info$Species, ellipse = T)


ggplot(data = BLUPs_anova, aes(x = ENV, y = BLUP)) +
  geom_boxplot(aes(fill = Material), width = 0.8) + theme_classic()

ggplot(data = BLUPs_anova, aes(x = Species, y = BLUP)) +
  geom_boxplot(width = 0.8) + theme_classic()+ coord_flip()


#selected genotypes parameters (macro and microscopical)
colonies <- read.xlsx("U_pisi colonies.xlsx", sheet = "DEFINITIVO")
corr_plot(colonies)


#Efecto precocidad en roya:
BLUP <- read.xlsx(xlsxFile = "BLUP_field.xlsx", sep= "\t", rowNames = T, sheet = "BLUP_f+cc", )
head(BLUP)

corr_plot(BLUP, c("GDD_F", "GDD_P", "Rust", "Biomass","Yield", "AUDPC", "LP50", "IF", "IT", "DS"))

#genetic correlation:
DS_DS <- read_xlsx("data/CC_R19traits.xlsx", sheet = "Index_DS")

mixed_modsdfsdf <- 
  gamem_met(DS_DS,
            env = ENV,
            gen = GEN,
            resp = Rust,
            random = "gen", #Default
            verbose = TRUE) #Default
get_model_data(mixed_modsdfsdf, "genpar")

install.packages("quant_gen")
Pheno_rust <- data.frame(Pheno_rust)
Pheno_rust$R19 <- as.numeric(Pheno_rust$R19)
Pheno_rust$AUDPC <- as.numeric(Pheno_rust$AUDPC)
Pheno_rust$IF <- as.numeric(Pheno_rust$IF)
Pheno_rust$IT <- as.numeric(Pheno_rust$IT)
Pheno_rust$DS <- as.numeric(Pheno_rust$DS)
Pheno_rust$I_cc_FAI <- as.numeric(Pheno_rust$I_cc_FAI)
corr_plot(Pheno_rust, c("R19", "AUDPC", "IF", "IT", "DS", "I_cc_FAI"))

0.30/(sqrt(0.76)*sqrt(0.67))
