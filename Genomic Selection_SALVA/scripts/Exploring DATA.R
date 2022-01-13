##LOAD DATA##
library(readxl)
library(rAverage)

  #Season 2018 Salva with values already normalized with controls:
R18_nc <- read_excel("~/GitHub/Genomic-Selection/Genomic Selection_SALVA/data/R18_19_20.xlsx", 
                     sheet = "R18_nc", 
                     col_types = 
                       c("text", #ENV
                         "numeric", #REP
                         "numeric", #BLOCK
                         "numeric", #ROW
                         "text", #GEN
                         "numeric", #GERM
                         "numeric", #DtF
                         "numeric", #DtP
                         "numeric", #GDD_F	
                         "numeric", #GDD_P
                         "numeric", #biomass_pl
                         "numeric", #yield_pl
                         "numeric", #yield_KgHa
                         "numeric", #Asco_adj
                         "numeric", #Rust
                         "numeric", #Rust_adj
                         "numeric" #Oidio
                       ))
#asd
#Season 2019 Salva with values already normalized with controls:
R19_nc <- read_excel("~/GitHub/Genomic-Selection/Genomic Selection_SALVA/data/R18_19_20.xlsx", 
                     sheet = "R19_nc", 
                     col_types = 
                       c("text", #ENV
                         "numeric", #REP
                         "numeric", #BLOCK
                         "numeric", #ROW
                         "text", #GEN
                         "numeric", #GERM
                         "numeric", #GERMadj
                         "numeric", #DtFF
                         "numeric", #DtFFadj
                         "numeric", #GDD_FF
                         "numeric", #DtF
                         "numeric", #DtFadj
                         "numeric", #GDD_F
                         "numeric", #DtFP
                         "numeric", #DtFPadj
                         "numeric", #GDD_FP
                         "numeric", #DtP
                         "numeric", #DtPadj
                         "numeric", #GDD_P
                         "numeric", #DtM
                         "numeric", #DtMadj
                         "numeric", #GDD_M
                         "numeric", #rust
                         "numeric", #rustadj
                         "numeric", #biomass_pl
                         "numeric", #yield_pl
                         "numeric" #yield_KgHa
                       )) 
#Season 2020 Salva with values already normalized with controls:
R20_nc <- read_excel("~/GitHub/Genomic-Selection/Genomic Selection_SALVA/data/R18_19_20.xlsx", 
                  sheet = "R20_nc", 
                  col_types = 
                    c("text", #ENV
                      "numeric", #REP
                      "numeric", #BLOCK
                      "numeric", #ROW
                      "text", #GEN
                      "numeric", #GERM
                      "numeric", #DtFF
                      "numeric", #GDD_FF
                      "numeric", #DtF
                      "numeric", #GDD_F
                      "numeric", #DtFP
                      "numeric", #GDD_FP
                      "numeric", #DtP
                      "numeric", #GDD_P
                      "numeric", #DtM
                      "numeric", #GDD_M
                      "numeric", #rustadj
                      "numeric", #asco
                      "numeric", #asco_adj
                      "numeric", #oidio
                      "numeric", #Oidio_adj
                      "numeric", #Biom_plant_g
                      "numeric", #yield_plant_g	
                      "numeric"  #yield_kgHa
                    ))

#All season in a df with the common traits, values already normalized with controls:

Rall <- read_excel("~/GitHub/Genomic-Selection/Genomic Selection_SALVA/data/R18_19_20.xlsx", 
                   sheet = "R_all", 
                   col_types = 
                     c("text", #ENV
                       "numeric", #REP
                       "numeric", #BLOCK
                       "numeric", #ROW
                       "text", #GEN
                       "numeric", #GERM
                       "numeric", #DtF
                       "numeric", #GDD_F
                       "numeric", #DtP
                       "numeric", #GDD_P
                       "numeric", #Biom_pl
                       "numeric", #yield_pl
                       "numeric", #rustadj
                       "numeric", #asco_adj
                       "numeric" #Oidio_adj
                     )) 
df_cc <- read_excel("~/GitHub/Rust-collection/data/Rust_controlled_condition.xlsx", 
                     sheet = "CC_n", 
                     col_types = 
                       c("numeric", #REP
                         "numeric", #BLOCK
                         "text", #GEN
                         "numeric", #L_area
                         "numeric", #L_aphila SI = 1; NO = 0
                         "numeric", #L_border Lisa = 0; Dentada = 1; Semi = 0,5
                         "numeric", #AUDPC_n
                         "numeric", #PL50_n
                         "numeric", #IF_n
                         "numeric", #DS_n
                         "numeric" #IT_n  
                         ))


#functions#
library("metan")


#inspect
inspect(Rall, plot = T)
desc_stat(Rall)

gen_alpha <- 
  gamem(R18_nc, 
        GEN, 
        REP,
        c(Rust, Rust_adj), 
        block = BLOCK)

get_model_data(gen_alpha)
plot(gen_alpha, type = "re")


mixed_mod <- 
  gamem_met(Rall,
            env = ENV,
            gen = GEN,
            rep = REP,
            block = BLOCK,
            resp = c("GDD_F", "GDD_P"),
            random = "gen", #Default
            verbose = TRUE)
plot(mixed_mod)
plot(mixed_mod, type = "re")
get_model_data(mixed_mod, "lrt")
get_model_data(mixed_mod)
hist(mixed_mod$Rust_adj$BLUPgen$BLUPg)
