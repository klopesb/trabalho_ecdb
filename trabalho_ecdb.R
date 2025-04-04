
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
colSums(is.na(rna_counts)) #no NA

colSums(is.na(metadata)) 

#Considering there are NA values in our metadata, we are going to verify it's proportion 
metadata_na_percent <- colSums(is.na(metadata)) / nrow(metadata) * 100
metadata_na_percent

#remove columns with > 10% of NAs

metadata_clean <- metadata[, colMeans(is.na(metadata)) <= 0.1] #from 84, now we have 47 cols 

#check if rna_counts and metadata are corresponding 
all(colnames(rna_normalized) %in% rownames(metadata_clean))

#normalize data

#apply log to estabilize gene variation
#rna_counts_log <- log2(rna_counts + 1)  # +1 para evitar log(0)

rna_filtered <- rna_counts[apply(rna_counts,1 , sd) > 0, ] #remove lines where expression in is 0
rna_normalized <- t(scale(t(rna_filtered))) #this way we don't produce NaNs 
colSums(is.na(rna_normalized))


#data standardization - applied z-score to each row(gene)
#summary(rna_counts_scaled)
apply(rna_normalized, 1, mean)  # Média de cada gene (deve ser ~0)
apply(rna_normalized, 1, sd)    # Desvio padrão de cada gene (deve ser ~1)


#check if the data type is correct 
str(rna_counts)
str(metadata)

#Exploratory Analysis

  


