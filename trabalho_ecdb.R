
#Install and load packages
install.packages(c("data.table", "tidyverse", "Bioconductor", "TCGAbiolinks"))


library(data.table)
library(tidyverse)
library(TCGAbiolinks)  # For TCGA data
library(SummarizedExperiment)


#This block of code is to gather data directly from The Cancer Genome Atlas

query <- GDCquery(
  project = "TCGA-UCEC", #specific code for: The Cancer Genome Atlas - Uterine Corpus Endometrial Carcinoma
  data.category = "Transcriptome Profiling", #Download data for genomic expression - RNA-Seq 
  data.type = "Gene Expression Quantification", #Genomic expression obtained from RNA-Seq
  workflow.type = "STAR - Counts"
)

GDCdownload(query)  #Function to download the query made

rna_data <- GDCprepare(query)  # Prepare the data for analysis 

#rna_data generates a class object named SummarizedExperiment containing assays and metadata (60660 elements):

  #=> unstranded
  #=> stranded_first
  #=> stranded_second
  #=> tpm_unstrand
  #=> fpkm_unstrand
  #=> fpkm_uq_unstrand


#rna_data contains assays as columns and rows as genes
rna_data
colnames(rna_data) # lists assays
rownames(rna_data) #lists all genes

#Convert rna_data to a dataframe 

rna_counts <- as.data.frame(assay(rna_data))  # matrix counting
metadata <- as.data.frame(colData(rna_data))  # Clinic metadata 


#Initial data exploration 
#Verify data structure

#check black spaces, inconsistencies, odd characters 
colnames(rna_counts) 
rownames(rna_counts)
colnames(metadata)

#check if the number of genes/assays are correct (genes x assays)
dim(rna_counts)

#check if there are essential columns missing on our metadata 
dim(metadata)

#check if there are NA values
colSums(is.na(rna_counts))

colSums(is.na(metadata)) 

#There are NA values in metadata, we need to verify how we can handle it 

#remove rows where all values are null
rna_counts <- rna_counts[rowSums(is.na(rna_counts)) < ncol(rna_counts), ]   
metadata <- metadata[rowSums(is.na(metadata)) < ncol(metadata), ] 
#check if the data type are correct 
str(rna_counts)
str(metadata)

#Summarize data 


