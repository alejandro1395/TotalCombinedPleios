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
Freqs_table <- read.csv("freqs_total.tsv", header = TRUE, sep = "\t")
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
sorted_corrected_pvalues <- sort(corrected_pvalues)
pleios_df <- as.data.frame(cbind(names(sorted_corrected_pvalues), as.numeric(sorted_corrected_pvalues)))
colnames(pleios_df) <- c("Combined_diseases", "Corrected_pval")
View(pleios_df) 

#histogram
par(mfrow=c(1,1))
ggplot(data=pleios_df[1:25,], aes(x= reorder(Combined_diseases, as.numeric(as.character(Corrected_pval))), 
                                  y=round(as.numeric(as.character(Corrected_pval)), 3), fill=Combined_diseases), 
       fill=Combined_diseases) + scale_fill_hue(l=10, c=45) +
  geom_bar(stat="identity", width=0.2) + guides(fill=FALSE) +
  theme_minimal() + theme(aspect.ratio = 0.9/1, 
        axis.text.x = element_text(angle = 40, hjust = 1)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red") + ylim(0, 1) + ylab("Corrected_p_values") +
  xlab("Combined diseases")

