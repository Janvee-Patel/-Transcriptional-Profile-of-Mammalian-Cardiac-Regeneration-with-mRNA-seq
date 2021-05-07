#----------------------------------------------
#Title: BF528 Individual Project 5 - Programmer
#Author: Janvee Patel
#Date: 05/04/2021
#----------------------------------------------

#Code for Part 5 - Plotting histogram of FPKM values from the genes.fpkm_tracking table of Cufflinks output

#Set the current working directory
setwd("/projectnb/bf528/users/frizzled/project_5_jpatel2/Programmer/")

#load libraries
library(dplyr)
library(ggplot2)

#Read in the FPKM table from Cufflinks output
P01_FPKM <- read.table("P0_1_cufflinks/genes.fpkm_tracking", header=TRUE)
P01_FPKM <- P01_FPKM[!sapply(P01_FPKM, function(x) all(x == "-"))]
nrow(P01_FPKM)

#remove genes that have FPKM value < 1
P01_FPKM <- P01_FPKM %>% filter(FPKM > 1)
nrow(P01_FPKM)

#Plot histogram of FPKM values for all genes 
ggplot(data=P01_FPKM, mapping=aes(x=log10(FPKM)+1)) + geom_histogram(bins=25, fill="azure2", color="azure4") + xlab("log10(FPKM)+1") + ylab("Frequency") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background = element_blank(), axis.line = element_line()) + labs(title="Histogram of FPKM Values for All Genes\n(Filter of FPKM > 1)") + theme(plot.title = element_text(hjust = 0.5))
