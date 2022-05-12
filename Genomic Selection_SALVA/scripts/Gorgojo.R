##LOAD DATA##
library(readxl)
library(openxlsx)
library(metan)
library(ggstatsplot)
library(tidyverse)
library(GGEBiplots)
library(multcompView)
#####Gorgojo Data
Gorgojo_data <- read_excel("GitHub/Genomic-Selection/Genomic Selection_SALVA/data/Gorgojo_data.xlsx", 
                           col_types = c("numeric", "text", "numeric", "numeric"))
#(transformedvalue = 180/p · arcsine [sqroot(%/100)])
arcsin <- function(percentage) {
  Tvalue <- (180/pi)* asin(sqrt(percentage/100))
  return(Tvalue)
}
Gorgojo_data <- Gorgojo_data %>%
  select(GEN, ENV, REP, DS) %>%
  mutate(DS_t = arcsin(percentage = DS))
head(Gorgojo_data)
ggplot(data = Gorgojo_mean, aes(x = DS_t_mean)) +
  geom_histogram()+
  facet_wrap(~ENV, scales = "free")
ggplot(data = Gorgojo_data, aes(x = ENV, y = DS_t)) +
  geom_boxplot()
Gorgojo_mean <- Gorgojo_data %>%
  group_by(GEN, ENV) %>%
  summarise(DS_t_mean = mean(DS_t),
            DS_mean = mean(DS))
#some stats
desc_stat(Gorgojo_data, stats = "main", na.rm = T)

#AMMI model:
model <- performs_ammi(Gorgojo_data,
                       env = ENV,
                       gen = GEN,
                       rep = REP,
                       resp = DS, DS_t,
                       verbose = FALSE)
get_model_data(model)
    #Biplots:
a <- plot_scores(model,
                 highlight = c("93", "200"))
b <- plot_scores(model,
                 type = 2,
                 polygon = TRUE,
                 col.env = "gray70",
                 col.segm.env = "gray70",
                 col.gen = transparent_color(),
                 highlight = c("93", "200"),
                 axis.expand = 1.5)
c <- plot_scores(model, type = 4)
arrange_ggplot(a, b, c, tag_levels = "a", nrow = 1)

#GGE model and plots:
modelgge <- gge(Gorgojo_data, ENV, GEN, DS, max_iter = 100)
modelgge2 <- gge(Gorgojo_data, ENV, GEN, DS, svp = "genotype")
modelgge3 <- gge(Gorgojo_data, ENV, GEN, DS, svp = "symmetrical")
d <- plot(modelgge)
e <- plot(modelgge2, type = 8)
f <- plot(modelgge2,
          type = 2,
          col.gen = "black",
          col.env = "gray70",
          axis.expand = 1.5)
arrange_ggplot(d, e, f, tag_levels = list(c("d", "e", "f")), nrow = 1)

#BLUP model
model2 <- 
  gamem_met(Gorgojo_data,
            env = ENV,
            gen = GEN,
            rep = REP,
            resp = c(DS, DS_t))
g <- plot_blup(model2)
h <- plot_blup(model2,
               prob = 0.2,
               col.shape = c("gray20", "gray80"),
              invert = T)
arrange_ggplot(g, h, tag_levels = list(c("g", "h")))

#parametric and non-parametric indices
stats <- ge_stats(Gorgojo_data, ENV, GEN, REP, resp = c(DS, DS_t))
get_model_data(stats)


#mixed model:
mixed_mod <- 
  gamem_met(Gorgojo_data,
            env = ENV,
            gen = GEN,
            rep = REP,
            resp = c(DS, DS_t),
            random = "gen", #Default
            verbose = TRUE) #Default
plot(model2, type = "re", var = 2)
plot(model2, var = 1)
hist(model2$DS_t$BLUPgen$Predicted)
hist(model2$DS_t$BLUPint$Predicted)
hist(model2$DS_t$MeansGxE$Y)

table_var <- get_model_data(model2)
plot(model2$DS_t$BLUPgen$Predicted, model2$DS_t$BLUPgen$Y)
plot(model2$DS_t$BLUPint$Predicted, model2$DS_t$BLUPint$`BLUPg+ge`)
#variances of DS:
CVg = 37.2
CVr = 83.5
mu = 12.8
varg = ((CVg * mu)/100)^2
varr = ((CVr * mu)/100)^2
vari =(varg/0.333)-varg
CVi = (sqrt(vari)/mu)*100
#variances of DS_t:
CVg = 26.3
CVr = 51.5
mu = 12.8
varg = ((CVg * mu)/100)^2
varr = ((CVr * mu)/100)^2
vari =(varg/0.333)-varg
CVi = (sqrt(vari)/mu)*100



CVi_T <- tibble(Parameters = c("CVi", "vari", "varg", "varr"), 
                        DS = c(52.65,45.41, 22.67, 114.23), 
                        DS_t = c(37.22,22.7,11.33,43.45))
CVi_T
table_var <-rbind(table_var, CVi_T)
table_MeansGxE <- model2$DS_t$MeansGxE #No reps no NA's
table_BLUPint <- model2$DS_t$BLUPint #Predicted value with reps and no NA's
table_BLUPgen <- model2$DS_t$BLUPgen #Predicted value no ENV, no REP, no Na's

hist(table_BLUPgen$Y)
ggplot(data = table_MeansGxE, aes(x = ENV, y = Y))+
  geom_boxplot() +
  theme_classic() +
  coord_flip()

ggplot(data = table_BLUPint, aes(x = ENV, y = Predicted))+
  geom_boxplot(fill = REP)

ggplot(data = table_BLUPint , aes(x = Predicted)) +
  geom_histogram() +
  facet_wrap(~ENV)
  
table_BLUPint2 <- table_BLUPint %>%
  group_by(ENV, GEN) %>%
  summarise(Predicted_mean = mean(Predicted))
table_pivot <- table_BLUPint2 %>%
  pivot_wider(names_from = ENV, values_from = Predicted_mean)
corr_plot(table_pivot, `Co-18`,`Co-19`, `Co-19W`, `Co-20`, `Co-20W`,
          maxsize = 4,
          minsize = 2,
          smooth = T)
#write.csv(table_pivot, "table_gorgojo.txt")
table_info <- read.csv("table_gorgojo.txt", sep = "\t")
head(table_info)
ggplot(data = table_BLUPint2 , aes(x = Predicted_mean)) +
  geom_histogram() +
  facet_wrap(~ENV)

#Genotype-by-environment means:
make_mat(Gorgojo_data, GEN, ENV, DS_t)
#Genotype-by-environment interaction effects
ge_ef <- ge_effects(Gorgojo_data, ENV, GEN, DS_t)
GxE_inter <- ge_ef$DS_t
#Genotype plus Genotype-by-environment effects
gge_ef <- ge_effects(Gorgojo_data, ENV, GEN, DS_t, type = "gge")
GxEplusG <- gge_ef$DS_t
#ANOVA ind
ind <- anova_ind(Gorgojo_data, ENV, GEN, REP, DS_t)
ind$DS_t$individual
ind2 <- anova_ind(Gorgojo_data, ENV, GEN, REP, DS)
ind2$DS$individual

#pca
gorgojo.pca <- prcomp(table_info[6:10], center = T, scale. = T)
gorgojo.material <- table_info$Material
gorgojo.structure <- table_info$Structure
gorgojo.specie <- table_info$Specie
summary(gorgojo.pca)
ggbiplot::ggbiplot(gorgojo.pca, labels = table_info$GEN, ellipse = T, groups = gorgojo.specie)
corr_plot(table_info, Co.18, Co.19, Co.19W, Co.20, Co.20W, col.by = Material)

#histograms:
a <- ggplot(table_info, aes(x = Co.18))+
  geom_histogram(color="black", fill="white",
                 binwidth = 1.7) +
  theme_classic()
b <- ggplot(table_info, aes(x = Co.19))+
  geom_histogram(color="black", fill="white",
                 binwidth = 2.7) +
  theme_classic()
c <- ggplot(table_info, aes(x = Co.20))+
  geom_histogram(color="black", fill="white",
                 binwidth = 3.7) +
  theme_classic()
d <- ggplot(table_info, aes(x = Co.19W))+
  geom_histogram(color="black", fill="white",
                 binwidth = 3.7) +
  theme_classic()
e <- ggplot(table_info, aes(x = Co.20W))+
  geom_histogram(color="black", fill="white",
                 binwidth = 3.7) +
  theme_classic()
f <- ggplot(table_BLUPgen, aes(x = Predicted))+
  geom_histogram(color="black", fill="white",
                 binwidth = 1.7) +
  theme_classic()

arrange_ggplot(a,b,c,d,e,f, tag_levels = list(c("a", "b", "c", "d", "e", "f")))

ggbetweenstats(table_MeansGxE, x= ENV, y = Y, plot.type = "boxviolin")


#####Linear mixed model with other package:
install.packages("sommer")
library("sommer")
library("tidyverse")

head()
fit <-mmer(fixed = Y  ~ ENV, 
           random = ~ GEN,
           rcov = ~ units,
           data = table_MeansGxE, 
           verbose = FALSE)
summary(fit)
ratio <- fit$sigmaVector[1]/fit$sigmaVector[2]
fit$sigmaVector
12.26473/45.97501
data(DT_wheat)
DT <- DT_wheat
GT <- GT_wheat
hist(fit$U$GEN$Y)
blup_Y <- data.frame(fit$U$GEN$Y)
blup_Y %>%
  arrange(fit$U$GEN$Y)
Gorgojo_mean$GEN <- as.factor(Gorgojo_mean$GEN)
fit2 <-mmer(fixed = DS_t  ~ ENV + REP, 
           random = ~ GEN,
           rcov = ~ units,
           data = Gorgojo_data, 
           verbose = FALSE)
summary(fit2)
hist(fit2$U$GEN$DS_t_mean)
blup_Y2 <- data.frame(fit2$U$GEN$DS_t_mean)
blup_Y2 %>%
  arrange(fit2$U$GEN$DS_t_mean)
#LMM by ENV
Gorgojo_data$GEN <- as.factor(Gorgojo_data$GEN)
Gorgojo_data$ENV <- as.factor(Gorgojo_data$ENV)  
Gorgojo_data$REP <- as.factor(Gorgojo_data$REP)  

Co18 <- Gorgojo_data %>%
  filter(ENV == "Co-18")
Co19 <- Gorgojo_data %>%
  filter(ENV == "Co-19")
Co20 <- Gorgojo_data %>%
  filter(ENV == "Co-20")
Co19W <- Gorgojo_data %>%
  filter(ENV == "Co-19W")
Co20W <- Gorgojo_data %>%
  filter(ENV == "Co-20W")

  fit18 <- mmer(fixed = DS_t ~ REP,
               random = ~ GEN,
               rcov = ~units,
               data = Co18,
               verbose = F)
  summary(fit18)
  CVg18 <- fit18$sigmaVector[1]/fit18$sigmaVector[2]
  BLUP18 <- data.frame(fit18$U); hist(BLUP18$DS_t); #write.csv(BLUP18, "BLUP18_gorgojoLMM.txt")

  fit19 <- mmer(fixed = DS_t ~ REP,
                random = ~ GEN,
                rcov = ~units,
                data = Co19,
                verbose = F)
  summary(fit19)
  CVg19 <- fit19$sigmaVector[1]/fit19$sigmaVector[2]
  BLUP19 <- data.frame(fit19$U); hist(BLUP19$DS_t); #write.csv(BLUP19, "BLUP19_gorgojoLMM.txt")
  
  fit20 <- mmer(fixed = DS_t ~ REP,
                random = ~ GEN,
                rcov = ~units,
                data = Co20,
                verbose = F)
  summary(fit20)
  CVg20 <- fit20$sigmaVector[1]/fit20$sigmaVector[2]
  BLUP20 <- data.frame(fit20$U); hist(BLUP20$DS_t); #write.csv(BLUP20, "BLUP20_gorgojoLMM.txt")
  
  fit19W <- mmer(fixed = DS_t ~ REP,
                random = ~ GEN,
                rcov = ~units,
                data = Co19W,
                verbose = F)
  summary(fit19W)
  CVg19W <- fit19W$sigmaVector[1]/fit19W$sigmaVector[2]
  BLUP19W <- data.frame(fit19W$U); hist(BLUP19W$DS_t); #write.csv(BLUP19W, "BLUP19W_gorgojoLMM.txt")
  
  fit20W <- mmer(fixed = DS_t ~ REP,
                 random = ~ GEN,
                 rcov = ~units,
                 data = Co20W,
                 verbose = F)
  summary(fit20W)
  CVg20W <- fit20W$sigmaVector[1]/fit20W$sigmaVector[2]
  BLUP20W <- data.frame(fit20W$U); hist(BLUP20W$DS_t); #write.csv(BLUP20W, "BLUP20W_gorgojoLMM.txt")
BLUP_env <- data.frame(BLUP18, BLUP19, BLUP20, BLUP19W, BLUP20W)
row.names(BLUP_env) <- c(BLUP18, BLUP19, BLUP20, BLUP19W, BLUP20W)
#pca2
BLUP_env <- read.csv("BLUP_GorgojoLMM.txt", sep = "\t")
gorgojo.pca2 <- prcomp(BLUP_env[2:6], center = T, scale. = T)
summary(gorgojo.pca2)
ggbiplot::ggbiplot(gorgojo.pca2, ellipse = T, labels = BLUP_env$GEN)
corr_plot(BLUP_env, Co18, Co19, Co19W, Co20, Co20W)
#LMM total
BLUP_env2$GEN <- as.factor(BLUP_env2$GEN)
ans1 <- mmer(DS_t ~ ENV,
             random = ~ GEN + ENV:GEN,
             rcov = ~units,
             data = BLUP_env2)

summary(ans1)
#### calculate heritability
vpredict(fit19, h2 ~ V1/(V1+V2) )
CVg19

#### genetic correlation
vpredict(ans1, gen.cor ~ V2 / sqrt(V1*V3))

ansDG <- mmer(DS_t ~ ENV,
              random= ~ vs(ds(ENV),GEN, Gu=A),
              rcov= ~ units,
              data=BLUP_env2, verbose = FALSE)
summary(ansDG)

library(sommer)
data(DT_example)
DT <- DT_example
A <- as.matrix(read.csv("GENmatrix.csv", sep = ";", header = T))
colnames(A) <- c(1:320)
