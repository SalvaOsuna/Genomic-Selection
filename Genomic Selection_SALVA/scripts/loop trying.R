# Create a 19 x 19 matrix (of 19 rows and 19 blocks)
  #vector for every rep:
RowsR1 <- c(GERM_R1$T_Germ)
BlocksR1 <- c(GERM_B1$T_Germ)

RowsR2 <- c(GERM_R2$T_Germ)
BlocksR2 <- c(GERM_B2$T_Germ)

RowsR3 <- c(GERM_R3$T_Germ)
BlocksR3 <- c(GERM_B3$T_Germ)

  #matrixes for every rep:
my_matR1 <- matrix(nrow=max(R18$ROW), ncol=max(R18$BLOCK))
for (i in 1:length(RowsR1)){
  for (j in 1:length(BlocksR1)){
    my_matR1[i,j]<-(mean(R18$T_Germ))/((RowsR1[i]+BlocksR1[j])/2)       
  }
}

my_matR2 <- matrix(nrow=max(R18$ROW), ncol=max(R18$BLOCK))
for (i in 1:length(RowsR2)){
  for (j in 1:length(BlocksR2)){
    my_matR2[i,j]<-(mean(R18$T_Germ))/((RowsR2[i]+BlocksR2[j])/2)       
  }
}

my_matR3 <- matrix(nrow=max(R18$ROW), ncol=max(R18$BLOCK))
for (i in 1:length(RowsR3)){
  for (j in 1:length(BlocksR3)){
    my_matR3[i,j]<-(mean(R18$T_Germ))/((RowsR3[i]+BlocksR3[j])/2)       
  }
}

# Create a 19 x 19 x 3 array (of 19 Rows, 19 Blocks and 3 Reps)


Rows
Blocks

my_array <- array(NaN, dim=c(19, 19, 3))

for (i in 1:length(Rows)) {
  for (j in 1:length(Blocks)) {
    for (k in 1:max(Reps)) {
      if (Reps[k] == "1") {
      my_array[i,j]<-(mean(R18$T_Germ))/((Rows[i]+Blocks[j])/2)
    }
    }
  }
}
