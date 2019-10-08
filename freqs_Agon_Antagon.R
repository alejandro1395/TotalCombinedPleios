##########################################
##########################################
###RScript to make correlations between###
###Data of Primates from AnAge for their##
###lifespan and the computed ortholog ####
####lengths from ensembl##################
##########################################

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
library(readr)
library(dplyr)
library(reshape2)
library(tibble)

###############
###############
#PHENOTIPIC DB#
###############

#Import Phenotypes from AnAge
Counts_table <- read.csv("freqs_Agon_Antagon.tsv", header = TRUE, sep="\t")
View(Counts_table)

#Description plots

new_Counts_table <- melt(data = Counts_table, id.vars = "Combined_diseases", measure.vars = c("Agonistic", "Antagonistic"))

ggplot(data = new_Counts_table
                  , aes(x = Combined_diseases, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'dodge', alpha = 2/3) +  
  labs(x = 'Combination_of_disease', y = NULL, fill = 'Type'
       , title = 'Proportions of Pleiotropye type by Combination') +
  facet_wrap(~ variable) + scale_fill_manual(values=c("#997799", "#E69F00")) + theme(aspect.ratio = 0.9/1, 
      axis.text.x = element_text(angle = 40, hjust = 1, size = 5))

newsort_Counts_table <- new_Counts_table[order(new_Counts_table$Combined_diseases),]
ggplot(data = newsort_Counts_table,
        aes(x = 1, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = 'dodge', alpha = 2/3) +  
  labs(x = 'Combination_of_disease', y = NULL, fill = 'Type'
       , title = 'Proportions of Pleiotropye type by Combination') +
  facet_wrap(~ Combined_diseases) + scale_fill_manual(values=c("#997799", "#E69F00")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())# what's the reader looking at?

#Chi-square table and props

rownames(Counts_table) <- Counts_table[,1]
Counts_table <- Counts_table[,2:3]
chisq.test(Counts_table)


Counts_table<-Counts_table %>% rownames_to_column('new_column')
Counts_table_filtered = Counts_table[ rowSums(Counts_table)!=0, ]


#Post-Hoc test

#Test for Normality (Shapiro Test)

shapiro.test(Counts_table_filtered$Agonistic)
shapiro.test(Counts_table_filtered$Antagonistic)

pvalues <- c()
for (i in 1:nrow(Counts_table_filtered)){
  Agon_posthoc <- c(Counts_table_filtered$Agonistic[i], sum(Counts_table_filtered$Agonistic[-i]))
  Antagon_posthoc <- c(Counts_table_filtered$Antagonistic[i], sum(Counts_table_filtered$Antagonistic[-i]))
  pv <- as.double(chisq.test(cbind(Agon_posthoc, Antagon_posthoc))["p.value"])
  pvalues<- append(pvalues, pv)
  print(rownames(Counts_table_filtered)[i])
  print(chisq.test(cbind(Agon_posthoc, Antagon_posthoc)))
}

#adjust FDR
names(pvalues) <- as.character(rownames(Counts_table_filtered))
corrected_pvalues <- p.adjust(pvalues, method = "bonferroni", n = length(pvalues))
sorted_corrected_pvalues <- sort(corrected_pvalues)
pleios_df <- as.data.frame(cbind(names(sorted_corrected_pvalues), as.numeric(sorted_corrected_pvalues)))
colnames(pleios_df) <- c("Combined_diseases", "Corrected_pval")
View(pleios_df) 