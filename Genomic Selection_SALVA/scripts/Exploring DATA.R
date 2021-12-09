    ##LOAD DATA##
library(readxl)
R18_19_20 <- read_excel("C:/Users/Salva/Desktop/DOCTORADO NUBE/ENSAYOS/GS/Genomic Selection_SALVA/data/R18_19_20.xlsx", 
                              sheet = "R_all", col_types = c("text", 
                                                                        "text", "text", "text", "text", "text", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "numeric"))


Rall <- data.frame(R18_19_20)
as_factor(Rall$ENV)
as_factor(Rall$REP)
as_factor(Rall$BLOCK)
as_factor(Rall$GEN)




#functions#

install.packages("metan"
)
library("metan")


#inspect

inspect(Rall, plot = T)
desc_stat(Rall)

gen_alpha <- 
  gamem(R18, 
        GEN, 
        REP,
        Biomass_pl_gr, 
        block = BLOCK)

plot(gen_alpha)
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
