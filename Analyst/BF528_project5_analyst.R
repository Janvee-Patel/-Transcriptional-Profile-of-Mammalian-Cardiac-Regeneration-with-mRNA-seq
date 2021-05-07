#----------------------------------------------
#Title: BF528 Individual Project 5 - Analyst
#Author: Janvee Patel
#Date: 05/04/2021
#----------------------------------------------

#Code for Analyst Part 6 - Identifying differentially expressed genes associated with myocyte differentiation

#load libraries
library(dplyr)
library(ggplot2)
library(gridExtra)

#Set the current working directory
setwd("/projectnb/bf528/users/frizzled/project_5_jpatel2/Analyst")

#read in the gene_exp.diff file from Cuffdiff output
#sort by q-value
gene_DE <- read.table("/projectnb/bf528/users/frizzled/project_5_jpatel2/Programmer/cuffdiff_out/gene_exp.diff", header=TRUE) %>% filter(status=="OK")
gene_DE <- gene_DE %>% arrange(q_value)

#generate histogram of log2 fold change for all genes
p1 <- ggplot(data=gene_DE, mapping=aes(x=log2.fold_change.)) + geom_histogram(bins=30, fill="azure2", color="azure4") + xlab("log2 Fold Change") + ylab("Frequency") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background = element_blank(), axis.line = element_line()) + labs(title="Histogram of log2 Fold Change\nfor All Genes", tag="A") + theme(plot.title = element_text(hjust = 0.5))

#subset by the significant term and determine total number of genes using this method
#generate histogram of log2 fold change for significant genes
sig_gene_DE <- subset(gene_DE, significant=="yes")
nrow(sig_gene_DE)
p2 <- ggplot(data=sig_gene_DE, mapping=aes(x=log2.fold_change.)) + geom_histogram(bins=35, fill="azure2", color="azure4") + xlab("log2 Fold Change") + ylab("Frequency") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background = element_blank(), axis.line = element_line()) + labs(title="Histogram of log2 Fold Change\nfor Significant Genes", tag="B") + theme(plot.title = element_text(hjust = 0.5))

#combine the 2 plots together
overall <- grid.arrange(p1, p2, ncol=2)

#Further subset the significant genes into up and down regulated genes
up_sig_gene_DE <- subset(sig_gene_DE, log2.fold_change. >= 0)
nrow(up_sig_gene_DE)

dwn_sig_gene_DE <- subset(sig_gene_DE, log2.fold_change. < 0)
nrow(dwn_sig_gene_DE)

#Subset the original gene_exp.diff file using the filter of p-value < 0.01
pval_gene_DE <- subset(gene_DE, p_value < 0.01)
nrow(pval_gene_DE)

#Further subset the significant genes from p-value < 0.01 filter into up and down-regulated genes
up_pval_gene_DE <- subset(pval_gene_DE, log2.fold_change. >= 0)
nrow(up_pval_gene_DE)

dwn_pval_gene_DE <- subset(pval_gene_DE, log2.fold_change. < 0)
nrow(dwn_pval_gene_DE)

#write out the up and down-regulated genes to CSV and txt files
#contains only the gene symbols
write(up_sig_gene_DE$gene, file="up_DEG.txt")
write(dwn_sig_gene_DE$gene, file="down_DEG.txt")

#contains all data including the gene symbols
write.csv(up_sig_gene_DE, file="up_DEG.csv")
write.csv(dwn_sig_gene_DE, file="down_DEG.csv")


#plot volcano plots of log2 FC and log10(p-value)
#remove infinite logFC values
gene_DE <- gene_DE[is.finite(gene_DE$log2.fold_change.),]

#split into not significant using the q-value < 0.05 filter or (significant label)
not_sig <- gene_DE[gene_DE$q_value > 0.05, ]
nsx = not_sig$log2.fold_change.
nsy = -log10(not_sig$p_value)

#split into significant genes
sig <- gene_DE[gene_DE$q_value < 0.05, ]
  
#split into positive fold change
positive_fc <- sig[sig$log2.fold_change. > 0, ]
px = positive_fc$log2.fold_change.
py = -log10(positive_fc$p_value)
  
#split into negative fold change
negative_fc <- sig[sig$log2.fold_change. < 0, ]
nx = negative_fc$log2.fold_change.
ny = -log10(negative_fc$p_value)
  
#plot the volcano plot
clrs <- c("Not Significant"="grey", "Increased"="red", "Decreased"="blue")
ggplot() + geom_point(data=not_sig, aes(x=nsx, y=nsy, color="Not Significant"), alpha=0.3) + geom_point(data=positive_fc, aes(x=px, y=py, color="Increased"), alpha=0.3) + geom_point(data=negative_fc, aes(x=nx, y=ny, color="Decreased"), alpha=0.3) + xlab("Log2 Fold Change") + theme_bw() + labs(title="Volcano Plot of Significant Differentially Expressed Genes", color="Differentially Expressed") + ylab("-log10(P-value)") + scale_color_manual(values=clrs) + theme(plot.title = element_text(size=17, hjust = 0.5))
