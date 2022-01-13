######################################################################################
##                                                                                  ##
##                              PLOTS ARTICLE                                       ##
##                                                                                  ##
######################################################################################


## ENViRONMENT: 
  {
# What is the effect of the treatment on the value ?
model=lm( Rall_field$Rust ~ Rall_field$ENV )
ANOVA=aov(model)

# Tukey test to study each pair of treatment :
TUKEY <- TukeyHSD(x=ANOVA, 'Rall_field$ENV', conf.level=0.95)

# Tuckey test representation :
plot(TUKEY , las=1 , col="brown")


# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY, "Rall_field$ENV")


# A panel of colors to draw each group with the same color :
my_colors <- c( 
  rgb(143,199,74,maxColorValue = 255),
  rgb(242,104,34,maxColorValue = 255), 
  rgb(111,145,202,maxColorValue = 255)
)

# Draw the basic boxplot
a <- boxplot(Rall_field$Rust ~ Rall_field$ENV , ylim=c(min(Rall_field$Rust) , 1.1*max(Rall_field$Rust)) , col=my_colors[as.numeric(LABELS[,1])] , 
             ylab="Disease Severity (%)",
             xlab = "Environment",
             main="Distribution by environment",
             las = 1)

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )

#Add the labels
text( c(1:nlevels(Rall_field$ENV)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )
  }
  
## SPECIES
{
  # What is the effect of the treatment on the value ?
  BLUPs_anova$Species <- as.factor(BLUPs_anova$Species)
  model=lm( BLUPs_anova$BLUP ~ BLUPs_anova$Species)
  ANOVA=aov(model)
  
  # Tukey test to study each pair of treatment :
  TUKEY <- TukeyHSD(x=ANOVA, 'BLUPs_anova$Species', conf.level=0.95)
  
  # Tuckey test representation :
  plot(TUKEY , las=1 , col="brown")
  
  
  # I need to group the treatments that are not different each other together.
  generate_label_df <- function(TUKEY, variable){
    
    # Extract labels and factor levels from Tukey post-hoc 
    Tukey.levels <- TUKEY[[variable]][,4]
    Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
    
    #I need to put the labels in the same order as in the boxplot :
    Tukey.labels$treatment=rownames(Tukey.labels)
    Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
    return(Tukey.labels)
  }
  
  # Apply the function on my dataset
  LABELS <- generate_label_df(TUKEY, "BLUPs_anova$Species")
  
  
  # A panel of colors to draw each group with the same color :
  my_colors <- c( 
    rgb(143,199,74,maxColorValue = 255),
    rgb(242,104,34,maxColorValue = 255), 
    rgb(111,145,202,maxColorValue = 255)
  )
  
  # Draw the basic boxplot
  a <- boxplot(BLUPs_anova$BLUP ~ BLUPs_anova$Species , ylim=c(min(BLUPs_anova$BLUP) , 1.1*max(BLUPs_anova$BLUP)) , col=my_colors[as.numeric(LABELS[,1])] , 
               ylab="Disease Severity (%)",
               xlab = "Species",
               main="Distribution by Species",
               las = 1,
               horizontal = T)
  
  # I want to write the letter over each box. Over is how high I want to write it.
  over <- 0.1*max( a$stats[nrow(a$stats),] )
  
  #Add the labels
  text( c(1:nlevels(Rall_field$ENV)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )
}

