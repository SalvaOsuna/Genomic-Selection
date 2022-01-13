
                                                            ### 1. INSTALLING all required library in R ###

##install LDHeatMap
#install.packages("LDheatmap")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager:: install(c("snpStats","rtracklayer","GenomicRanges","GenomInfoDb","IRanges"))

# install.packages ("agricolae")
# install.packages ("ggplot2")
# install.packages ("qqperm")
# install.packages ("qvalue")
# BiocManager::install("qvalue")

#para instalar GAPIT3 packages:
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3",force=TRUE)



                                                          #### 2. Load all required library and GAPIT code 

library(LDheatmap)

 source("http://zzlab.net/GAPIT/GAPIT.library.R")
 source("http://zzlab.net/GAPIT/gapit_functions.txt")

library(GAPIT3) # alternativa a las dos lineas de codigo anterior para cargar las funciones y libreria necesaria para GAPIT
library(agricolae)
library (dplyr)
library(QQperm)
library (qvalue)
library(ggplot2)



                                                        ### 3. Import databases in Gapit ###

{

myY <- read.table("<nombre archivo fenotipo.txt>", sep="\t", head = TRUE) #import phenotype
myG <- read.table(file="GenPea_SilDArT_sort_def.hmp.txt",sep= "\t",head=FALSE) #import genotype in hapmap format
myKI <- read.table("GenPea_SilDArT_Kinship.txt", sep="\t",head = FALSE) #import Kinship matrix
myCV <- read.table("GenPea_SilDArT_PCA_23PC.txt", sep="\t", head = TRUE) # import Covariate Dataset - here PCAs up to 50% variance
}


                            
                                                    ### 4. GWAS ANALYSIS with GAPIT ###

## RUN GAPIT for the trait/environment of interest

myGAPIT <- GAPIT(
  Y=myY [c(1,2)],
  G=myG,
  CV=myCV,
  KI=myKI,
  model=c("GLM","MLM","MLMM","BLINK", "FarmCPU","SUPER"),
  Multiple_analysis=TRUE)


## Load GAPIT output for each model - modificar los nombre de los archivo cada vez


{
GLMga<-read.csv ("GAPIT.GLM.Mean_IR_PLHT_T2cor.GWAS.Results.csv", header=TRUE)  
MLMga<-read.csv ("GAPIT.MLM.Mean_IR_PLHT_T2cor.GWAS.Results.csv", header=TRUE)
MLMM<-read.csv ("GAPIT.MLMM.Mean_IR_PLHT_T2cor.GWAS.Results.csv", header=TRUE)
FarmCPU<-read.csv ("GAPIT.FarmCPU.Mean_IR_PLHT_T2cor.GWAS.Results.csv", header=TRUE)
Blink<-read.csv ("GAPIT.Blink.Mean_IR_PLHT_T2cor.GWAS.Results.csv", header=TRUE)
SUPER<-read.csv ("GAPIT.SUPER.Mean_IR_PLHT_T2cor.GWAS.Results.csv", header=TRUE)
distribTeor<-1:1955/1955#crea una lista numerica con la distribución de p teoricos de 1 al numero maximo de marcadores aqui: 1955
}

## Estimate and draw a graphical representation of the genomic inflation factor and save it to a single eps file

{setEPS()
  postscript(file="MultiLambdaplot_Mean_IR_PLHT_T2cor.eps",onefile=FALSE,width=24/2.54, height=16/2.54, family="Times")
  par(mfcol=c(2,3))
  estlambda2(GLMga$P.value,distribTeor, plot = TRUE, adjust.xy = TRUE)
  mtext("GLM",line=1, cex= 1.5, col=1)
  estlambda2(MLMga$P.value,distribTeor, plot = TRUE, adjust.xy = TRUE)
  mtext("MLM",line=1,cex=1.5,col=2)
  estlambda2(MLMM$P.value,distribTeor, plot = TRUE, adjust.xy = TRUE)
  mtext("MLMM",line=1,cex=1.5,col=3)
  estlambda2(FarmCPU$P.value,distribTeor, plot = TRUE, adjust.xy = TRUE)
  mtext("FarmCPU",line=1,cex=1.5,col=4)
  estlambda2(Blink$P.value,distribTeor, plot = TRUE, adjust.xy = TRUE)
  mtext("Blink",line=1,cex=1.5,col=5)
   estlambda2(SUPER$P.value,distribTeor, plot = TRUE, adjust.xy = TRUE)
  mtext("Super",line=1,cex=1.5,col=6)
  dev.off()
}

{
  GLM_lambda<- estlambda2(GLMga$P.value,distribTeor, plot = FALSE, adjust.xy = TRUE)
  MLM_lambda<- estlambda2(MLMga$P.value,distribTeor, plot = FALSE, adjust.xy = TRUE)
  MLMM_lambda<- estlambda2(MLMM$P.value,distribTeor, plot = FALSE, adjust.xy = TRUE)
  FarmCPU_lambda<- estlambda2(FarmCPU$P.value,distribTeor, plot = FALSE, adjust.xy = TRUE)
  Blink_lambda<- estlambda2(Blink$P.value,distribTeor, plot = FALSE, adjust.xy = TRUE)
  SUPER_lambda<- estlambda2(SUPER$P.value,distribTeor, plot = FALSE, adjust.xy = TRUE)
  }



## Draw the QQplots for each model and save them on a single eps file

{setEPS()
  postscript(file="MultiQQplot_Mean_IR_GYP_T2cor.eps",onefile=FALSE,width=24/2.54, height=15/2.54, family="Times")
  par(mfcol=c(2,3))
   qqplot(distribTeor,GLMga$P.value, adjust.xy=TRUE)
  qqplot(distribTeor,MLMga$P.value, adjust.xy=TRUE)
  mtext("GLM/MLM GAPIT",col=1)
  qqplot(distribTeor,MLMM$P.value, adjust.xy=TRUE)
  qqplot(distribTeor,FarmCPU$P.value, adjust.xy=TRUE)
  mtext("MLMM/FarmCPU",col=2)
  qqplot(distribTeor,Blink$P.value, adjust.xy=TRUE)
  qqplot(distribTeor,SUPER$P.value, adjust.xy=TRUE)
  mtext("Blink/Super",col=3)
  dev.off()
}

#### Q value analysis for FDR validation

{
GLMga_qval <- qvalue(p = GLMga$P.value) #estima los qvalues para cada lista de p-value
summary (GLMga_qval) #muestra el resumen de los datos de q-values
plot(GLMga_qval, rng = c(0,1)) #dibuja los graficos de Qvalues con el pio, y el numero de false discovery per valor de qvalue
}

plot(GLMga_qval, rng = c(0,0.96)) #dibuja los graficos por los valores de qvalue entre dos valores de q
GLMga_qlim<-0.7 #cutoff de qvalue


{
  MLMga_qval<- qvalue(p=MLMga$P.value)
  summary (MLMga_qval) #muestra el resumen de los datos de q-values
  plot(MLMga_qval, rng = c(0,1)) #dibuja los graficos de Qvalues con el pio, y el numero de false discovery per valor de qvalue
}

plot(MLMga_qval, rng = c(0,0.85)) #dibuja los graficos por los valores de qvalue entre dos valores de q
MLMga_qlim<-0.7 #cutoff de qvalue

{
  MLMM_qval<- qvalue(p=MLMM$P.value)
  summary (MLMM_qval) #muestra el resumen de los datos de q-values
  plot(MLMM_qval, rng = c(0,1)) #dibuja los graficos de Qvalues con el pio, y el numero de false discovery per valor de qvalue
}

plot(MLMM_qval, rng = c(0,0.7)) #dibuja los graficos por los valores de qvalue entre dos valores de q
MLMM_qlim<-0.3 #cutoff de qvalue

{
  FarmCPU_qval<- qvalue(p=FarmCPU$P.value)
  summary (FarmCPU_qval) #muestra el resumen de los datos de q-values
  plot(FarmCPU_qval, rng = c(0,1)) #dibuja los graficos de Qvalues con el pio, y el numero de false discovery per valor de qvalue
}

plot(FarmCPU_qval, rng = c(0,0.5)) #dibuja los graficos por los valores de qvalue entre dos valores de q
FarmCPU_qlim<-0.3 #cutoff de qvalue

{
  Blink_qval<- qvalue(p=Blink$P.value)
  summary (Blink_qval) #muestra el resumen de los datos de q-values
  plot(Blink_qval, rng = c(0,1)) #dibuja los graficos de Qvalues con el pio, y el numero de false discovery per valor de qvalue
}

plot(Blink_qval, rng = c(0,0.7)) #dibuja los graficos por los valores de qvalue entre dos valores de q
Blink_qlim<-0.3 #cutoff de qvalue

{
  SUPER_qval<- qvalue(p=SUPER$P.value)
  summary (SUPER_qval) #muestra el resumen de los datos de q-values
  plot(SUPER_qval, rng = c(0,1)) #dibuja los graficos de Qvalues con el pio, y el numero de false discovery per valor de qvalue
}

plot(SUPER_qval, rng = c(0,0.7)) #dibuja los graficos por los valores de qvalue entre dos valores de q
SUPER_qlim<-0.3 #cutoff de qvalue



## Create a table with significant marker-trait association for all models

{
  GLMga_Res<- GLMga %>% 
    mutate (Q_values=GLMga_qval$qvalues, Threshold_Q=GLMga_qlim, Model="GLM-GAPIT", Lambda=GLM_lambda$estimate) %>% 
    filter (Q_values<GLMga_qlim)
  
   MLMga_Res<- MLMga %>% 
    mutate (Q_values=MLMga_qval$qvalues, Threshold_Q=MLMga_qlim, Model="MLM-GAPIT",Lambda=MLM_lambda$estimate) %>%
    filter (Q_values<MLMga_qlim)
  
  MLMM_Res<- MLMM %>% 
   mutate (Q_values=MLMM_qval$qvalues, Threshold_Q=MLMM_qlim, Model="MLMM-GAPIT",Lambda=MLMM_lambda$estimate) %>%
   filter (Q_values<MLMM_qlim)
 
  FarmCPU_Res<- FarmCPU %>% 
    mutate (Q_values=FarmCPU_qval$qvalues, Threshold_Q=FarmCPU_qlim, Model="FarmCPU-GAPIT",Lambda=FarmCPU_lambda$estimate) %>%
    filter (Q_values<FarmCPU_qlim)
  
  Blink_Res<- Blink %>% 
    mutate (Q_values=Blink_qval$qvalues, Threshold_Q=Blink_qlim, Model="Blink-GAPIT",Lambda=Blink_lambda$estimate) %>%
    filter (Q_values<Blink_qlim)
  
  SUPER_Res<- SUPER %>% 
    mutate (Q_values=SUPER_qval$qvalues, Threshold_Q=SUPER_qlim, Model="SUPER-GAPIT",Lambda=SUPER_lambda$estimate) %>%
    filter (Q_values<SUPER_qlim)
  
  Sig<-GLMga_Res %>% rbind(MLMga_Res) %>% rbind(MLMM_Res) %>% rbind(FarmCPU_Res) %>% rbind(Blink_Res) %>% rbind(SUPER_Res) %>% arrange(SNP) %>% filter (Lambda<1.2 & Lambda>0.8)# Combina y ordena los dos archivos filtrado
  write.table (Sig, "GWAS_Sig_GAPIT_Mean_IR_PLHT_T2cor.txt", sep ="\t",row.names = FALSE) #escribe la nueva tabla en un archivo .txt
}


