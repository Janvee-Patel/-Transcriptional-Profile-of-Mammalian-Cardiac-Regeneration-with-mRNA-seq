#-------------------------------------
#title: BF528 Project5 Biologist Re-do
#author: Janvee Patel
#date: 05/04/2021
#-------------------------------------

#import libraries
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(gridExtra)
library(ggpubr)
library(dendextend)

setwd("/projectnb/bf528/users/frizzled/project_5_jpatel2/Analyst")

#import FPKM tracking tables for each of the replicates of the samples 
#P0
P01 <- read.table("/projectnb/bf528/users/frizzled/project_5_jpatel2/Programmer/P0_1_cufflinks/genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)
colnames(P01)[10] <- "P0_1_FPKM"
P02 <- read.table("FPKM_Tables/P02genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)
colnames(P02)[10] <- "P0_2_FPKM"

#P4
P41 <- read.table("FPKM_Tables/P41genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)
colnames(P41)[10] <- "P4_1_FPKM"
P42 <- read.table("FPKM_Tables/P42genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)
colnames(P42)[10] <- "P4_2_FPKM"

#P7
P71 <- read.table("FPKM_Tables/P71genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)
colnames(P71)[10] <- "P7_1_FPKM"
P72 <- read.table("FPKM_Tables/P72genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)
colnames(P72)[10] <- "P7_2_FPKM"

#Ad
Ad1 <- read.table("FPKM_Tables/Ad1genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)
colnames(Ad1)[10] <- "Ad_1_FPKM"
Ad2 <- read.table("FPKM_Tables/Ad2genes.fpkm_tracking", header=TRUE) %>% arrange(.,tracking_id)
colnames(Ad2)[10] <- "Ad_2_FPKM"

#not in function
`%notin%` <- Negate(`%in%`)


#7.1
#Generate plots of FPKM values for P0, P4, P7, Ad samples
#Extract all values for all genes from each of the time periods, then just merged across
#Plotted the averages of the replicates for each sample
clrs <- colorRampPalette(brewer.pal(8, "Dark2"))(16)
gen_fpkmplot <- function(gene_list, title_name, ltag) {
  P0 <- (P01[P01$gene_short_name %in% gene_list, "P0_1_FPKM"] + P02[P02$gene_short_name %in% gene_list, "P0_2_FPKM"]) / 2
  P4 <- (P41[P41$gene_short_name %in% gene_list, "P4_1_FPKM"] + P42[P42$gene_short_name %in% gene_list, "P4_2_FPKM"]) / 2
  P7 <- (P71[P71$gene_short_name %in% gene_list, "P7_1_FPKM"] + P72[P72$gene_short_name %in% gene_list, "P7_2_FPKM"]) / 2
  Ad <- (Ad1[Ad1$gene_short_name %in% gene_list, "Ad_1_FPKM"] + Ad2[Ad2$gene_short_name %in% gene_list, "Ad_2_FPKM"]) / 2
  fpkm_genes <- cbind(P0, P4, P7, Ad)  #merged across
  row.names(fpkm_genes) <- gene_list
  fpkm_genes <- data.frame(t(fpkm_genes))  #transpose
  fpkm_genes$timeperiod <- rownames(fpkm_genes)
  fpkm_genes <- melt(fpkm_genes, id.vars=c("timeperiod"))
  thm = theme(plot.title = element_text(hjust=0.5))
  ggplot(fpkm_genes, aes(x=factor(timeperiod, level=c("P0", "P4", "P7", "Ad")), y=value, colour=variable, group=variable)) + geom_line() + geom_point() + labs(title = title_name, tag=ltag, x="Maturation Time Stages", y="FPKM", colour="Genes") + theme_bw() + thm + scale_color_manual(values=clrs)
}


#Sarcomere
gene_sarc <- c("Tcap", "Des", "Myoz2", "Pdlim5", "Csrp3", "Cryab", "Pygm")
splot <- gen_fpkmplot(gene_sarc, "Sarcomere", "A")

#Mitochondria
#No Mpc1 in our FPKM tables
gene_mito <- c("Slc25a11", "Prdx3", "Echs1", "Phyh", "Acat1")
mplot <- gen_fpkmplot(gene_mito, "Mitochondria", "B")

#Cell Cycle
#No Bora in our FPKM tables
gene_cc <- c("Cdc45", "Cdc6", "Cdc27", "Aurkb", "Cdc23", "Cdc7", "Rad51", "E2f1", "E2f8", "Cdc26", "Cdk7")
cplot <- gen_fpkmplot(gene_cc, "Cell Cycle", "C")

#Arrange plots into one figure
overall <- grid.arrange(splot, mplot, cplot, nrow=3)
annotate_figure(overall, top= "FPKM Values of Representative Genes During Maturation Stages\n", fig.lab.size = 10)


#Created FPKM Matrix with the 8 samples above
#Merged the FPKM columns from each of the tracking tables from the FPKM matrix previously created into single dataframe
#Addressed situations with a duplicate tracking_id after merging the tables

#Get tables to be merged
P0_1_FPKM <- P01[c(1, 5, 10)]
fpkm_mat <- read.csv("/project/bf528/project_2/data/fpkm_matrix.csv", header=TRUE, sep="\t") %>% distinct()

#Merge tables
fpkm_combined <- merge(P0_1_FPKM, fpkm_mat, by="tracking_id")

#Contains all duplicated rows
duplct <- fpkm_combined %>% group_by(tracking_id) %>% dplyr::filter(n() > 1)

#Remove the duplicated rows
fpkm_combined = subset(fpkm_combined, fpkm_combined$tracking_id %notin% duplct$tracking_id)

#Contains Ensembl tracking_ids ENSMUSG00000074899 and ENSMUSG00000044083
duplct_add <- duplct %>% group_by(tracking_id) %>% dplyr::filter(n() > 4)
duplct <- subset(duplct, duplct$tracking_id %notin% duplct_add$tracking_id)

#Edit the duplicates to get appropriate rows
duplct <- duplct %>% group_by(tracking_id) %>% slice(c(1,4))

#Edit duplct_add to get appropriate rows
duplct_add <- duplct_add[c(3, 5, 7, 11, 14, 20, 25), ]

#Add duplicated rows back into fpkm_combined
fpkm_combined <- rbind(fpkm_combined, duplct, duplct_add)

#Find the most differentially expressed genes
#Read in gene_exp.diff file
#Filter out rows to get where the status is "OK"
gene_de <- read.table("/projectnb/bf528/users/frizzled/project_5_jpatel2/Programmer/cuffdiff_out/gene_exp.diff", header=TRUE) %>% filter(status=="OK")

#Get the top 1000 significantly differentially expressed genes
#Ordered by q_value
gene_de <- gene_de %>% dplyr::filter(significant=="yes")
top_genes <- head(gene_de[order(gene_de$q_value), ], n=1000)

#Find and remove which entries have multiple gene symbols listed
#Chose the first gene symbol for rows with multiple listed and added back into gene_symbol
s <- c(which(grepl(",", top_genes$gene)))
gene_symbol <- top_genes[-s, c("gene")]
v <- vector()

for (entry in s) {
  v <- append(v, str_split(top_genes[entry, "gene"], ",")[[1]][1])
}

#Contains 1000 gene symbols
gene_symbol <- c(gene_symbol, v)


#Additional code displaying how I attempted to implement the different Bioconductor packages to map gene symbol to Ensembl id
#These methods did not return all Ensembl ids for the gene symbols that were inputted for reasons such as gene symbol was gene synonym for given id, and different or deprecated ids
#Results from options 1, 2, 3 include rows from top_genes that had multiple gene symbols listed per row

#Option 1 (not taken)
#Returns 921 Ensembl tracking_ids for the 1000 gene symbols that were inputted
#library("biomaRt")
#searchDatasets(mart = ensembl, pattern = "musculus")
#ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#ensmusids <- getBM(attributes=c('ensembl_gene_id'), filters = 'mgi_symbol', values = gene_symbol, mart = ensembl)


#Option 2 (not taken)
#Returns 938 Ensembl tracking_ids for the 1000 gene symbols that were inputted
#library(EnsDb.Mmusculus.v79)
#ensmusids <- ensembldb::select(EnsDb.Mmusculus.v79, keys= gene_symbol, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))


#Option 3 (not taken)
#Returns 857 Ensembl tracking_ids for the 1000 gene symbols that were inputted
#library(org.Mm.eg.db)
#ensmusids <- as.data.frame(mapIds(org.Mm.eg.db, gene_symbol, keytype="SYMBOL", column = "ENSEMBL"))
#sum(!(is.na(ensmusids)))


#Option 4 (taken)
#As the above methods returned differing numbers of Ensembl ids, I utilized the gene symbol
#Using Ensembl tracking_ids would have been an ideal method within the above Bioconductor packages, however, there were inconsistencies with the tracking_ids that I observed

#Subset the FPKM combined matrix by the top 1000 genes found to be differentially expressed between P0 and Ad
fpkm_combined_sub <- subset(fpkm_combined, fpkm_combined$gene_short_name %in% gene_symbol)


#Remove the tracking_id and gene_short_name columns from the subsetted matrix
#Add gene names back as the row names and address duplicates by adding "_" to differentiate
f_ready <- fpkm_combined_sub[-1:-2]
rownames(f_ready) <- make.unique(fpkm_combined_sub[,2], sep="_")

#Addresses rows that have 0s for clustering setup
#Set as matrix
f_ready <- f_ready[apply(f_ready, 1, function(x) !all(x==0)), ]
f_ready <- as.matrix(f_ready)

#rearrange the order of the columns
f_ready<- f_ready[,c("P0_1_FPKM", "P0_2_FPKM", "P4_1_FPKM", "P4_2_FPKM", "P7_1_FPKM", "P7_2_FPKM", "Ad_1_FPKM", "Ad_2_FPKM")]

#Add labels for clustered heatmap annotations for columns
samples_col <- data.frame(Samples = rep(c("P0", "P4", "P7", "Ad"), c(2,2,2,2)))
row.names(samples_col) <- colnames(f_ready)

#Generated clustered heatmap with genes along rows and 8 samples along columns with dendrograms and labels
pal <- rev(colorRampPalette(c("red2", "white", "blue2"))(100))  #color palette
clrs <- list(Samples = c(P0="#37414A", P4="#727B84", P7="#9a9a9a", Ad="#D9E2E1"))
cold <- pheatmap::pheatmap(f_ready)[[2]]
cold <- dendextend::rotate(cold, order = c("P0_1_FPKM", "P0_2_FPKM", "P4_1_FPKM", "P4_2_FPKM", "P7_1_FPKM", "P7_2_FPKM", "Ad_1_FPKM", "Ad_2_FPKM"))  #order of columns
png(filename="heatmap_bio.png", width=930, height=1037, res=200)
pheatmap::pheatmap(f_ready, col=pal, scale="row", cluster_rows = TRUE, clustering_distance_rows = "euclidean", cluster_cols=as.hclust(cold), show_rownames = FALSE, annotation_col = samples_col, annotation_colors=clrs, main = "RNA-Seq Gene Expression\n", fontsize = 11, legend=TRUE, angle_col = 90, fontsize_col = 10)
dev.off()