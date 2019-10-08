##########################################
##########################################
###RScript to make correlations between###
###Data of Primates from AnAge for their###
###lifespan and the computed ortholog ####
####lengths from ensembl#################
###########################################

#IMPORT libraries

library(ape)
library(geiger)
library(nlme)
library(phytools)
library(ggplot2)
library(qqman)
library(qualityTools)
library(Haplin)
library(gridBase)
library(tidyverse)
library(car)

###############
###############
#PHENOTIPIC DB#
###############

#Import Phenotypes from AnAge
Freqs_table <- read.csv("freqs_total_paired.tsv", header = TRUE, sep = "\t")
View(Freqs_table)

#Test of goodness of fit

scalar1 <- function(x) {x / sqrt(sum(x^2))}
vect_freqs_norm <- Freqs_table$Exp_prob/sum(Freqs_table$Exp_prob)
Freqs_table$Norm_exp_prob <- vect_freqs_norm

res <- chisq.test(Freqs_table$Total_observed, p = Freqs_table$Norm_exp_prob)
res$expected

Table_pairs <- data.frame(Combinations=Freqs_table$Combined_diseases,
                          Observed=Freqs_table$Total_observed, 
                          Expected=as.numeric((Freqs_table$Norm_exp_prob)*sum(Freqs_table$Total_observed)), 
                          stringsAsFactors=FALSE)
options(scipen=999)
View(Table_pairs)


#PLOT DATAFRAME
library(gridExtra)
library(grid)

par(mfrow=c(2,4))


#THEME
myt <- ttheme_default(base_size = 8,
  # Use hjust and x to left justify the text
  # Alternate the row fill colours
  core = list(fg_params=list(hjust = 1, x=1),
              bg_params=list(fill=c("white", "grey"))),
  
  # Change column header to white text and red background
  colhead = list(fg_params=list(col="white"),
                 bg_params=list(fill="black"))
)

  grid.newpage()
  grid.draw(tableGrob(format(Table_pairs[1:15, ], big.mark=","), theme=myt  ,rows=NULL))
  grid.draw(tableGrob(format(Table_pairs[16:30, ], big.mark=","), theme=myt  ,rows=NULL))
  grid.draw(tableGrob(format(Table_pairs[31:45, ], big.mark=","), theme=myt  ,rows=NULL))
  grid.draw(tableGrob(format(Table_pairs[46:60, ], big.mark=","), theme=myt  ,rows=NULL))
  grid.draw(tableGrob(format(Table_pairs[61:75, ], big.mark=","), theme=myt  ,rows=NULL))
  grid.draw(tableGrob(format(Table_pairs[76:90, ], big.mark=","), theme=myt  ,rows=NULL))
  grid.draw(tableGrob(format(Table_pairs[91:105, ], big.mark=","), theme=myt  ,rows=NULL))
  
#Post-Hoc test
  
#Test for Normality (Shapiro Test)
  
shapiro.test(Freqs_table$Total_observed)
  
pvalues <- c()
for (i in 1:nrow(Freqs_table)){
observed_posthoc <- c(Freqs_table$Total_observed[i], sum(Freqs_table$Total_observed[-i]))
expected_posthoc <- c(Freqs_table$Norm_exp_prob[i], sum(Freqs_table$Norm_exp_prob[-i]))
pv <- as.double(chisq.test(observed_posthoc, p=expected_posthoc)["p.value"])
pvalues<- append(pvalues, pv)
print(Freqs_table$Combined_diseases[i])
print(chisq.test(observed_posthoc, p=expected_posthoc))
}

#adjust FDR
names(pvalues) <- as.character(Freqs_table$Combined_diseases)
corrected_pvalues <- p.adjust(pvalues, method = "bonferroni", n = length(pvalues))
corrected_pvalues[corrected_pvalues<0.05]

ks.test(Freqs_table$Total_observed, Freqs_table$Norm_exp_prob)
ks.test(Freqs_table$Total_observed, as.numeric((Freqs_table$Norm_exp_prob)*sum(Freqs_table$Total_observed)))
