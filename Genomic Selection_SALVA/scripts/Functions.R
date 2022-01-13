library(openxlsx)
library(data.table)
library(tidyverse)
library(ggstatsplot)
library(readxl)
library(hrbrthemes)

fn.dataframe(R18)

fn.media.dataframe <- function(df, argmivector, controls) {
  
  print("The mean of the vector selected in the DataFrame is")
  print(paste("Mean: ",round(mean(argmivector, ),2)))
  
  if (controls == T) {
    print(Control <- df %>%
            group_by(ROW, REP) %>%
            filter(df, GEN == "Control"))
  } else {
    print("No controls list shown")
  }
}

fn.media.dataframe(R18, R18$Germ, controls = T)

##-------------------------------------------------------------------------##
##  Primero hago un ejemplo manual creando la columna de FC para un trait  ##
##  Después, si funciona el manual, intentaré crear una función aplicable  ##
##                 PARA TODOS LOS TRAITS AL MISMO TIEMPO                   ##
##-------------------------------------------------------------------------##

#Cargar datos de prueba
R18 <- read_excel("~/GitHub/Genomic-Selection/Genomic Selection_SALVA/data/Ensayocampo2018.xlsx", sheet = "RawData", col_types = 
                    c("text", #ENV
                      "text", #REP
                      "numeric", #BLOCK
                      "numeric", #ROW
                      "text", #GEN
                      "numeric", #GERM
                      "numeric", #BIOM
                      "numeric", #BIOM_PL
                      "numeric", #YIELD
                      "numeric", #YIELD_PL
                      "numeric", #BIOM/YIELD
                      "numeric", #Rust in %
                      "numeric", #Dtf
                      "numeric", #DtFP
                      "numeric", #DtP
                      "numeric", #Ascochyta in %
                      "numeric"  #Oidio in %
                    )) 

#Grafico que visualiza si hay diferencias significativas entre réplicas para los controles
R18 %>%
  filter(GEN == "Control") %>%
  ggbetweenstats(
    x     = REP,
    y     = T_Rust,
    title = "Distribution of Germination across Replicates in the Controls"
  )

##Intentar hacer una función que haga solo esto:
F1 <- function(df, trait) {
  
  p <- df %>% filter(GEN == "Control") %>% ggbetweenstats(
    x     = REP,
    y     = {{trait}},
    title = c("Distribution of {{trait}} across Replicates in the Controls")
  )
  return(p)
}
F1(R18, "T_Germ")
F1(R18, "Asco")

#Tabla con las medias de los controles por ROW y REP
GERM_R <- R18 %>%
  group_by(REP, ROW) %>%
  filter(GEN == "Control") %>%
  summarise(T_Germ = mean(T_Germ))

GERM_R1 <- R18 %>%
  group_by(ROW) %>%
  filter(REP == "1") %>%
  filter(GEN == "Control") %>%
  summarise(T_Germ = mean(T_Germ)) %>%
  arrange(ROW)

GERM_R2 <- R18 %>%
  group_by(ROW) %>%
  filter(REP == "2") %>%
  filter(GEN == "Control") %>%
  summarise(T_Germ = mean(T_Germ)) %>%
  arrange(ROW)

GERM_R3 <- R18 %>%
  group_by(ROW) %>%
  filter(REP == "3") %>%
  filter(GEN == "Control") %>%
  summarise(T_Germ = mean(T_Germ)) %>%
  arrange(ROW)

head(GERM_R)

GERM_R %>% 
  ggplot(aes(x=REP,y=T_Germ, fill=REP)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2)

##Intentar hacer una función que haga solo esto:
F2 <- function(df) {
  list(
    l <- df %>%
      group_by(REP, ROW) %>%
      filter(GEN == "Control") %>%
      summarise(T_Germ = mean(T_Germ)),
    
    p <- df %>% 
      ggplot(aes(x=REP,y=T_Germ, fill=REP)) +
      geom_boxplot() + geom_jitter(width=0.1,alpha=0.2)
  )
}
F2(R18)


#Tabla con las medias de los controles por BLOCK y REP 
GERM_B <- R18 %>%
  group_by(REP, BLOCK) %>%
  filter(GEN == "Control") %>%
  summarise(T_Germ = mean(T_Germ))

GERM_B1 <- R18 %>%
  group_by(BLOCK) %>%
  filter(REP == "1") %>%
  filter(GEN == "Control") %>%
  summarise(T_Germ = mean(T_Germ)) %>%
  arrange(BLOCK)

GERM_B2 <- R18 %>%
  group_by(BLOCK) %>%
  filter(REP == "2") %>%
  filter(GEN == "Control") %>%
  summarise(T_Germ = mean(T_Germ)) %>%
  arrange(BLOCK)

GERM_B3 <- R18 %>%
  group_by(BLOCK) %>%
  filter(REP == "3") %>%
  filter(GEN == "Control") %>%
  summarise(T_Germ = mean(T_Germ)) %>%
  arrange(BLOCK)

head(GERM_B)

GERM_B %>% 
  ggplot(aes(x=REP,y=T_Germ, fill=REP)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2)

#Lista con el Factor de Corrección para cada R(i) y B(j)
{
as_tibble (FC_GERM <- c(
col1 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[1])+GERM_B$T_Germ[1:19])/2), #REP1
col2 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[2])+GERM_B$T_Germ[1:19])/2), #REP1
col3 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[3])+GERM_B$T_Germ[1:19])/2), #REP1
col4 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[4])+GERM_B$T_Germ[1:19])/2), #REP1
col5 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[5])+GERM_B$T_Germ[1:19])/2), #REP1
col6 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[6])+GERM_B$T_Germ[1:19])/2), #REP1
col7 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[7])+GERM_B$T_Germ[1:19])/2), #REP1
col8 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[8])+GERM_B$T_Germ[1:19])/2), #REP1
col9 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[9])+GERM_B$T_Germ[1:19])/2), #REP1
col10 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[10])+GERM_B$T_Germ[1:19])/2), #REP1
col11 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[11])+GERM_B$T_Germ[1:19])/2), #REP1
col12 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[12])+GERM_B$T_Germ[1:19])/2), #REP1
col13 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[13])+GERM_B$T_Germ[1:19])/2), #REP1
col14 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[14])+GERM_B$T_Germ[1:19])/2), #REP1
col15 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[15])+GERM_B$T_Germ[1:19])/2), #REP1
col16 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[16])+GERM_B$T_Germ[1:19])/2), #REP1
col17 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[17])+GERM_B$T_Germ[1:19])/2), #REP1
col18 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[18])+GERM_B$T_Germ[1:19])/2), #REP1
col19 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[19])+GERM_B$T_Germ[1:19])/2), #REP1

col20 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[20])+GERM_B$T_Germ[20:38])/2), #REP2
col21 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[21])+GERM_B$T_Germ[20:38])/2), #REP2
col22 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[22])+GERM_B$T_Germ[20:38])/2), #REP2
col23 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[23])+GERM_B$T_Germ[20:38])/2), #REP2
col24 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[24])+GERM_B$T_Germ[20:38])/2), #REP2
col25 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[25])+GERM_B$T_Germ[20:38])/2), #REP2
col26 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[26])+GERM_B$T_Germ[20:38])/2), #REP2
col27 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[27])+GERM_B$T_Germ[20:38])/2), #REP2
col28 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[28])+GERM_B$T_Germ[20:38])/2), #REP2
col29 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[29])+GERM_B$T_Germ[20:38])/2), #REP2
col30 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[30])+GERM_B$T_Germ[20:38])/2), #REP2
col31 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[31])+GERM_B$T_Germ[20:38])/2), #REP2
col32 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[32])+GERM_B$T_Germ[20:38])/2), #REP2
col33 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[33])+GERM_B$T_Germ[20:38])/2), #REP2
col34 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[34])+GERM_B$T_Germ[20:38])/2), #REP2
col35 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[35])+GERM_B$T_Germ[20:38])/2), #REP2
col36 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[36])+GERM_B$T_Germ[20:38])/2), #REP2
col37 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[37])+GERM_B$T_Germ[20:38])/2), #REP2
col38 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[38])+GERM_B$T_Germ[20:38])/2), #REP2

col39 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[39])+GERM_B$T_Germ[39:57])/2), #REP3
col40 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[40])+GERM_B$T_Germ[39:57])/2), #REP3
col41 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[41])+GERM_B$T_Germ[39:57])/2), #REP3
col42 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[42])+GERM_B$T_Germ[39:57])/2), #REP3
col43 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[43])+GERM_B$T_Germ[39:57])/2), #REP3
col44 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[44])+GERM_B$T_Germ[39:57])/2), #REP3
col45 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[45])+GERM_B$T_Germ[39:57])/2), #REP3
col46 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[46])+GERM_B$T_Germ[39:57])/2), #REP3
col47 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[47])+GERM_B$T_Germ[39:57])/2), #REP3
col48 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[48])+GERM_B$T_Germ[39:57])/2), #REP3
col49 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[49])+GERM_B$T_Germ[39:57])/2), #REP3
col50 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[50])+GERM_B$T_Germ[39:57])/2), #REP3
col51 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[51])+GERM_B$T_Germ[39:57])/2), #REP3
col52 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[52])+GERM_B$T_Germ[39:57])/2), #REP3
col53 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[53])+GERM_B$T_Germ[39:57])/2), #REP3
col54 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[54])+GERM_B$T_Germ[39:57])/2), #REP3
col55 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[55])+GERM_B$T_Germ[39:57])/2), #REP3
col56 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[56])+GERM_B$T_Germ[39:57])/2), #REP3
col57 <- mean(R18$T_Germ)/(((GERM_R$T_Germ[57])+GERM_B$T_Germ[39:57])/2) #REP3
)
)
  }

cbind(R18, FC_GERM)

R18_t <- R18 %>%
  mutate(GERM_t = (T_Germ * FC_GERM))



#Comprar los valores transformados con los controles y los origianes en un histograma:

comparison <- data.frame(type = c(rep("Germ", 1083),rep("Germ_T", 1083)), 
               value = c(R18_t$T_Germ, R18_t$GERM_t))

p <- comparison %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 10) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="")



##-------------------------------------------------------------------------##
##  Ya tengo el ejemplo manual que añade una nueva columna al df con las   ##
##  normalizaciones hechas con los controles según su posicion en el dis.  ##
##                 experimental. ahora intentar hacer FM                   ##
##-------------------------------------------------------------------------##





fun1 <- function(df, TRAIT, REP, BLOCK, ROW){
  print(TRAIT *((mean(TRAIT/(mean())))))
}


##cogida de internet##
myDF <- data.table(
  Grouping=rep(c("P1","P2"),each=6),
  type = rep(c(rep("samp",times=4),"CRTL","CRTL"),times=2),
  ID= rep(1:6, times=2),
  feat1 = rnorm(12),
  feat2 = rnorm(12)
)

cols_grouping=c('Grouping', 'type')
cols_features=c('feat1','feat2')

setkeyv(myDF,"Grouping")
myDF_norm=myDF[,lapply(.SD, median, rm.NA=TRUE), .SDcols=cols_features, by=cols_grouping]
setkeyv(myDF_norm,"Grouping")

crt_normalization = function(sub_table){
  for (col in cols_features) {
    i_col=paste0("i.",col)
    sub_table[[col]]=sub_table[[col]]/sub_table[[i_col]]
    sub_table[[i_col]]=NULL
  }
  return(sub_table)
}

myDF_norm = myDF_norm[
  myDF_norm[type == "CRTL",
            c("Grouping",cols_features),
            with=FALSE]
][,crt_normalization(.SD),by='Grouping']



#función para sustituir outliers por debajo y por encima del 5th y 95th percertil:

fun.percentile <- function(x){
  quantiles <- quantile( x, c(.05, .95 ) )
  x[ x < quantiles[1] ] <- quantiles[1]
  x[ x > quantiles[2] ] <- quantiles[2]
  x
}

F1 <- function(df, trait) {
  
#Is there differences between REPs in the control for the trait selected??
  p <- df %>% filter(GEN == "Control") %>% ggbetweenstats(
    x     = REP,
    y     = {{trait}},
    title = c("Distribution of {{trait}} across Replicates in the Controls")
  )
  return(p)
}

#Esta función crea una nueva columna con el valor corregido con los controles en el diseño experimental de SalvaGwas
F2 <- function(df, trait) {
  #Create the vectors for the BLOCK/ROW average by REP
  trait_B1 <- df %>%
    group_by(BLOCK) %>%
    filter(REP == "1") %>%
    filter(GEN == "Control") %>%
    summarise(trait = mean({{trait}})) %>%
    arrange(BLOCK)
  
  trait_B2 <- df %>%
    group_by(BLOCK) %>%
    filter(REP == "2") %>%
    filter(GEN == "Control") %>%
    summarise(trait = mean({{trait}})) %>%
    arrange(BLOCK)
  
  trait_B3 <- df %>%
    group_by(BLOCK) %>%
    filter(REP == "3") %>%
    filter(GEN == "Control") %>%
    summarise(trait = mean({{trait}})) %>%
    arrange(BLOCK)
  
  trait_R1 <- df %>%
    group_by(ROW) %>%
    filter(REP == "1") %>%
    filter(GEN == "Control") %>%
    summarise(trait = mean({{trait}})) %>%
    arrange(ROW)
  
  trait_R2 <- df %>%
    group_by(ROW) %>%
    filter(REP == "2") %>%
    filter(GEN == "Control") %>%
    summarise(trait = mean({{trait}})) %>%
    arrange(ROW)
  
  trait_R3 <- df %>%
    group_by(ROW) %>%
    filter(REP == "3") %>%
    filter(GEN == "Control") %>%
    summarise(trait = mean({{trait}})) %>%
    arrange(ROW)
  
  #Create the matrixes with the CF for the trait selected:
  RowsR1 <- c(trait_R1$trait)
  BlocksR1 <- c(trait_B1$trait)
  
  RowsR2 <- c(trait_R2$trait)
  BlocksR2 <- c(trait_B2$trait)
  
  RowsR3 <- c(trait_R3$trait)
  BlocksR3 <- c(trait_B3$trait)
  
  mean_trait <- mean(df[,trait])
  
  my_matR1 <- matrix(nrow=max(df$ROW), ncol=max(df$BLOCK))
  for (i in 1:length(RowsR1)){
    for (j in 1:length(BlocksR1)){
      my_matR1[i,j] <- mean_trait/((RowsR1[i]+BlocksR1[j])/2)     
    }
  }
  
  my_matR2 <- matrix(nrow=max(df$ROW), ncol=max(df$BLOCK))
  for (i in 1:length(RowsR2)){
    for (j in 1:length(BlocksR2)){
      my_matR2[i,j]<-mean_trait/((RowsR2[i]+BlocksR2[j])/2)       
    }
  }
  
  my_matR3 <- matrix(nrow=max(df$ROW), ncol=max(df$BLOCK))
  for (i in 1:length(RowsR3)){
    for (j in 1:length(BlocksR3)){
      my_matR3[i,j]<-mean_trait/((RowsR3[i]+BlocksR3[j])/2)       
    }
  }
  
  #Convert the matrixes to one column to add it in the df  
  CF1<- as.tibble(unlist(c(my_matR1)))
  CF2 <- as.tibble(unlist(c(my_matR2)))
  CF3 <- as.tibble(unlist(c(my_matR3)))
  FC_trait <- rbind(c(CF1, CF2, CF3))
  
  #Add to the data frame the correction factor column and multiply by its trait, creating a new column
  
  
  #df <- df %>% mutate(trait_t = as.numeric(unlist({{trait}} * FC_trait)))
  dfnew <- df %>% mutate(trait_t = as.numeric(unlist(df[,trait])) * as.numeric(FC_BIOM$value))
  # mutate(yield_plant_gT = as.numeric(unlist(yield_plant_g * FC_BIOM$value)))
  
  return(dfnew)  
  
}

