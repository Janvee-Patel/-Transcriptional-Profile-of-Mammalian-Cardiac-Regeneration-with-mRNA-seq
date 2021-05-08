# BF528 Individual Project 5 Description
Transcriptional Profile of Mammalian Cardiac Regeneration with mRNA-Seq

The goals of this project were to align and normalize sequencing reads to the mm9 mouse genome (MGSCv37) using Tuxedo Suite software, assess quality control metrics, and perform differential expression analysis. 

Citation: O’Meara et al. Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration. Circ Res. Feb 2015. PMID: 25477501

# Repository Contents
Each folder indicates the role for the project. The folder contains the files used during that process of the analysis. 

Programmer:
   - run_tophat.qsub: script that aligns the RNA-Seq reads to the mm9 mouse reference genome
   - run_cufflinks.qsub: script that runs Cufflinks and quantifies how reads map to genomic regions defined by annotation
   - run_cuffdiff.qsub: script that runs Cuffdiff and quantifies differential gene expression
   - QCmet_RSeQC.qsub: script that runs quality control steps and read-mapping statistics
   - BF528_project5_programmer-5.R: R script that generates histograms of FPKM values

Analyst:
   - BF528_project5_analyst.R: R script that processes the differential gene expression output from Cuffdiff, generates histograms of log2 fold changes, and determines significant differentially expressed genes.

Biologist:
   - BF528_project5_biologist.R: R script that generates line plots of FPKM values for representative genes, and a clustered heatmap of the top 1,000 significant differentially expressed genes.
