# Q&A Alpha Lattice Design
library(agricolae)
library(tidyverse)

## see: http://subodhgreen.blogspot.com/2020/07/alpha-lattice-design-in-r-studio.html

data <-read.table(header = TRUE, text = "
Treatment	Replication	Block	Yield
1	1	1	6
2	1	1	7
3	1	1	5
4	1	1	8
5	1	1	6
6	1	2	16
7	1	2	12
8	1	2	12
9	1	2	13
10	1	2	8
11	1	3	17
12	1	3	7
13	1	3	7
14	1	3	9
15	1	3	14
16	1	4	18
17	1	4	16
18	1	4	13
19	1	4	13
20	1	4	14
21	1	5	14
22	1	5	15
23	1	5	11
24	1	5	14
25	1	5	14
1	2	1	24
6	2	1	13
11	2	1	24
16	2	1	11
21	2	1	8
2	2	2	21
7	2	2	11
12	2	2	14
17	2	2	11
22	2	2	23
3	2	3	16
8	2	3	4
13	2	3	12
18	2	3	12
23	2	3	12
4	2	4	17
9	2	4	10
14	2	4	30
19	2	4	9
24	2	4	23
5	2	5	15
10	2	5	15
15	2	5	22
20	2	5	16
25	2	5	19
")

data$Treatment=as.factor(data$Treatment)
data$Replication=as.factor(data$Replication)
data$Block=as.factor(data$Block)

str(data)
## 1st way
model<-aov(Yield~Replication+Treatment+Replication:Block, data=data)
summary(model)

out<-HSD.test(Yield, Treatment, DFerror = 16, MSerror = 13.66, alpha = 0.05,console=TRUE)

## 2nd way
modelLattice<-PBIB.test(Block,Treatment,Replication,Yield,k=5,console = TRUE,
                        method =c("VC"),test = "tukey",alpha = 0.05,group = TRUE)

bar.group(out$groups,ylim=c(0,20),xlab="treatment", ylab="Yield")
bar.err(out$means,ylim=c(0,20),xlab="treatment", ylab="Yield")

## Another example from
## https://www.rdocumentation.org/packages/agricolae/versions/1.3-3/topics/design.alpha

library(agricolae)
#Example one
trt<-1:30
t <- length(trt)
# size block k
k<-3
# Blocks s
s<-t/k
# replications r
r <- 2
outdesign<- design.alpha(trt,k,r,serie=2)
book<-outdesign$book
plots<-book[,1]
dim(plots)<-c(k,s,r)
for (i in 1:r) print(t(plots[,,i]))
outdesign$sketch
# Example two 
trt<-letters[1:12] 
t <- length(trt)
k<-3
r<-3
s<-t/k
outdesign<- design.alpha(trt,k,r,serie=2)
book<-outdesign$book
plots<-book[,1]
dim(plots)<-c(k,s,r)
for (i in 1:r) print(t(plots[,,i]))
outdesign$sketch
# }

## https://rdrr.io/cran/agridat/man/john.alpha.html
## uses: library(agridat)


library(agridat)
library(lattice)
library(lme4)
library(lucid)
# library(asreml)
library(asremlPlus)

data(john.alpha)
dat <- john.alpha

# RCB (no incomplete block)
(m0 <- lm(yield ~ 0 + gen + rep, data=dat))

# Block fixed (intra-block analysis) (bottom of table 7.4 in John)
m1 <- lm(yield ~ 0 + gen + rep + rep:block, dat)
anova(m1)

# Block random (combined inter-intra block analysis)
## libs(lme4, lucid)
(m2 <- lmer(yield ~ 0 + gen + rep + (1|rep:block), dat))

anova(m2)
## Analysis of Variance Table
##     Df Sum Sq Mean Sq  F value
## gen 24 380.43 15.8513 185.9942
## rep  2   1.57  0.7851   9.2123
vc(m2)
##        grp        var1 var2    vcov  sdcor
##  rep:block (Intercept) <NA> 0.06194 0.2489
##   Residual        <NA> <NA> 0.08523 0.2919


# Variety means.  John and Williams table 7.5.  Slight, constant
# difference for each method as compared to John and Williams.
means <- data.frame(rcb=coef(m0)[1:24],
                    ib=coef(m1)[1:24],
                    intra=fixef(m2)[1:24])
head(means)
##             rcb       ib    intra
## genG01 5.201233 5.268742 5.146433
## genG02 4.552933 4.665389 4.517265
## genG03 3.381800 3.803790 3.537934
## genG04 4.439400 4.728175 4.528828
## genG05 5.103100 5.225708 5.075944
## genG06 4.749067 4.618234 4.575394

##libs(lattice)
splom(means, main="john.alpha - means for RCB, IB, Intra-block")


# ----------
# asreml4

libs(asremlPus,lucid)

# Heritability calculation of Piepho & Mohring, Example 1

m3 <- asreml(yield ~ 1 + rep, data=dat, random=~ gen + rep:block)
sg2 <- summary(m3)$varcomp['gen','component'] # .142902

# Average variance of a difference of two adjusted means (BLUP)

p3 <- predict(m3, data=dat, classify="gen", sed=TRUE)
# Matrix of pair-wise SED values, squared
vdiff <- p3$sed^2
# Average variance of two DIFFERENT means (using lower triangular of vdiff)
vblup <- mean(vdiff[lower.tri(vdiff)]) # .05455038

# Note that without sed=TRUE, asreml reports square root of the average variance
# of a difference between the variety means, so the following gives the same value
# predict(m3, data=dat, classify="gen")$avsed ^ 2 # .05455038

# Average variance of a difference of two adjusted means (BLUE)
m4 <- asreml(yield ~ 1 + gen + rep, data=dat, random = ~ rep:block)
p4 <- predict(m4, data=dat, classify="gen", sed=TRUE)
vdiff <- p4$sed^2
vblue <- mean(vdiff[lower.tri(vdiff)]) # .07010875
# Again, could use predict(m4, data=dat, classify="gen")$avsed ^ 2

# H^2 Ad-hoc measure of heritability
sg2 / (sg2 + vblue/2) # .803

# H^2c Similar measure proposed by Cullis.
1-(vblup / 2 / sg2) # .809


# ----------

# Illustrate how to do the same calculations with lme4
# https://stackoverflow.com/questions/38697477

libs(lme4)

cov2sed <- function(x){
  # Convert var-cov matrix to SED matrix
  # sed[i,j] = sqrt( x[i,i] + x[j,j]- 2*x[i,j] )
  n <- nrow(x)
  vars <- diag(x)
  sed <- sqrt( matrix(vars, n, n, byrow=TRUE) +
                 matrix(vars, n, n, byrow=FALSE) - 2*x )
  diag(sed) <- 0
  return(sed)
}

# Same as asreml model m4. Note 'gen' must be first term
m5blue <- lmer(yield ~ 0 + gen + rep + (1|rep:block), dat)

libs(emmeans)
ls5blue <- emmeans(m5blue, "gen")
con <- ls5blue@linfct[,1:24] # contrast matrix for genotypes
# The 'con' matrix is identity diagonal, so we don't need to multiply,
# but do so for a generic approach
# sed5blue <- cov2sed(con 
tmp <- tcrossprod( crossprod(t(con), vcov(m5blue)[1:24,1:24]), con)
sed5blue <- cov2sed(tmp)


# vblue Average variance of difference between genotypes
vblue <- mean(sed5blue[upper.tri(sed5blue)]^2)
vblue # .07010875 matches 'vblue' from asreml

# Now blups
m5blup <- lmer(yield ~ 0 + (1|gen) + rep + (1|rep:block), dat)
# Need lme4::ranef in case ordinal is loaded
re5 <- lme4::ranef(m5blup,condVar=TRUE)
vv1 <- attr(re5$gen,"postVar")  
vblup <- 2*mean(vv1) # .0577 not exactly same as 'vblup' above
vblup

# H^2 Ad-hoc measure of heritability
sg2 <- c(lme4::VarCorr(m5blup)[["gen"]])  # 0.142902
sg2 / (sg2 + vblue/2) # .803 matches asreml

# H^2c Similar measure proposed by Cullis.
1-(vblup / 2 / sg2) # .809 from asreml, .800 from lme4



