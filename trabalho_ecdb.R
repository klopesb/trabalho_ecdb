
#Install and load packages
install.packages(c("data.table", "tidyverse", "Bioconductor", "TCGAbiolinks", "edgeR", "DESeq2"))


library(data.table)
library(tidyverse)
library(TCGAbiolinks)  # For TCGA data
library(SummarizedExperiment)
library(ggplot2)
library(DESeq2)
library(graphics)
library(lattice)
library(edgeR)
library(dplyr)


#This block of code is to gather data directly from The Cancer Genome Atlas

if (file.exists("rna_data_TCGA_UCEC.rds")) {
  rna_data <- readRDS("rna_data_TCGA_UCEC.rds")
} else {
  query <- GDCquery(
    project = "TCGA-UCEC",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  GDCdownload(query)
  rna_data <- GDCprepare(query)
  saveRDS(rna_data, file = "rna_data_TCGA_UCEC.rds")
}

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


#geneExp <-  SummarizedExperiment::assay(rna_data, "unstranded")

#dim(geneExp)

#apply log to estabilize gene variation
#rna_counts_log <- log2(rna_counts + 1)  # +1 para evitar log(0)

rna_filtered <- rna_counts[apply(rna_counts,1 , sd) > 0, ] #remove lines where expression in is 0
rna_normalized <- t(scale(t(rna_filtered))) #this way we don't produce NaNs 
colSums(is.na(rna_normalized))


#data standardization - applied z-score to each row(gene)
#summary(rna_counts_scaled)
#apply(rna_normalized, 1, mean)  # Média de cada gene (deve ser ~0)
#apply(rna_normalized, 1, sd)    # Desvio padrão de cada gene (deve ser ~1)


#check if the data type is correct 
str(rna_counts)
str(metadata)

#Exploratory Analysis
#Vale a pena comparar essas variáveis 
#Tissue Type 
table(colData(rna_data)$tissue_type)
table(colData(rna_data)$classification_of_tumor)
table(colData(rna_data)$race)
table(colData(rna_data)$primary_diagnosis)
table(colData(rna_data)$vital_status)

#selecção das variáveis de interesse: vital_status, age_at_index, classification_of_tumor, tissue_type, primary_diagnosis 
  #Endometrioid adenocarcinoma, NOS (404)
  #Serous cystadenocarcinoma, NOS (124)

#Fazer nova tabela metadados só com essas variaveis 
metadata_subset <- dplyr::select(
  metadata_clean,
  vital_status,
  age_at_index,
  classification_of_tumor,
  tissue_type,
  primary_diagnosis
)


diagnosis_filtered <- metadata_subset %>%
  mutate(primary_diagnosis = ifelse(
    primary_diagnosis %in% c("Endometrioid adenocarcinoma, NOS", "Serous cystadenocarcinoma, NOS"),
    primary_diagnosis,
    "Others"
  ))


tumor_class <- metadata_subset %>%
  mutate(classification_of_tumor = ifelse(classification_of_tumor == "primary", "primary", "others"))


head(metadata_subset)

rm(metadata_clean)
rm(pca)
rm(dds)
rm(expr)
rm(geneExp)
rm(rna_counts)
rm(pca_df)

#Frequências absolutas e relativas:
# Tabelas de frequência
table(metadata_subset$vital_status)
prop.table(table(metadata_subset$vital_status)) * 100

table(tumor_class$classification_of_tumor)
prop.table(table(tumor_class$classification_of_tumor)) * 100

table(metadata_subset$tissue_type)
prop.table(table(metadata_subset$tissue_type)) * 100

table(metadata_subset$primary_diagnosis)
prop.table(table(metadata_subset$primary_diagnosis)) * 100

#Gráficos de barra 
# Exemplo: vital_status
ggplot(metadata_subset[!is.na(metadata_subset$vital_status), ], aes(x = vital_status, y = after_stat(prop), group = 1)) +
  geom_bar(fill = "lightcoral") +
  theme_minimal() +
  labs(title = "Distribuição: Vital Status", x = "", y = "Proporção") +
  theme(plot.title = element_text(hjust = 0.5))

# Exemplo: classification_of_tumor <- juntar primario vs outros 
ggplot(tumor_class[!is.na(tumor_class$classification_of_tumor), ], aes(x = classification_of_tumor)) +
  geom_bar(fill = "lightcoral") +
  theme_minimal() +
  labs(title = "Distribuição: Classification of Tumor", x = "", y = "Frequência") +
  theme(plot.title = element_text(hjust = 0.5))

# Exemplo: tissue_type
ggplot(metadata_subset, aes(x = tissue_type)) +
  geom_bar(fill = "lightcoral") +
  theme_minimal() +
  labs(title = "Distribuição: Tissue Type", x = "", y = "Frequência") +
  theme(plot.title = element_text(hjust = 0.5))

#primary_diagnosis
ggplot(diagnosis_filtered, aes(x = primary_diagnosis)) +
  geom_bar(fill = "lightcoral") +
  theme_minimal() +
  labs(title = "Distribuição: Primary Diagnosis", x = "", y = "Frequência") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))

#Análise descritiva de age_at_index 

summary(metadata_subset$age_at_index)
sd(metadata_subset$age_at_index, na.rm = TRUE)

# Histograma de Idade 
ggplot(metadata_subset[!is.na(metadata_subset$age_at_index), ], aes(x = age_at_index)) +
  geom_density(fill = "lightcoral") +
  theme_minimal() +
  labs(title = "Densidade da Idade dos Pacientes",
       x = "Idade",
       y = "Densidade")



# Análises cruzadas

#Idade média dos pacientes vivos vs. mortos:
ggplot(metadata_subset[!is.na(metadata_subset$age_at_index), ], aes(x = vital_status, y = age_at_index)) +
  geom_boxplot(fill = "lightcoral") +
  theme_minimal() +
  labs(title = "Idade por Status Vital", x = "Status Vital", y = "Idade") +
  theme(plot.title = element_text(hjust = 0.5))

tapply(metadata_subset$age_at_index, metadata_subset$vital_status, summary)


# Pacientes com tumor x tecido normal: idade média
ggplot(metadata_subset[!is.na(metadata_subset$age_at_index), ], aes(x = tissue_type, y = age_at_index)) +
  geom_boxplot(fill = "lightcoral") +
  theme_minimal() +
  labs(title = "Idade por Tipo de Tecido", x = "", y = "Idade") +
  theme(plot.title = element_text(hjust = 0.5))

tapply(metadata_subset$age_at_index, metadata_subset$tissue_type, summary)

#Frequência de vivos/mortos por tipo de tumor
ggplot(tumor_class[!is.na(tumor_class$age_at_index), ], aes(x = classification_of_tumor, fill = vital_status)) +
  geom_bar(position = "fill") +  # proporções
  theme_minimal() +
  labs(title = "Proporção de Status Vital por Classificação do Tumor",
       x = "Classificação do Tumor", y = "Proporção") +
  scale_y_continuous(labels = scales::percent)

tapply(tumor_class$age_at_index, tumor_class$vital_status, summary)

#testes de hipótese 
#verificar se a variável numérica segue uma distribuição normal 

#H0 (Hipótese nula): Os dados seguem uma distribuição normal.

#H1 (Hipótese alternativa): Os dados não seguem uma distribuição normal.

shapiro.test(metadata_subset$age_at_index)

#rejeito H0, dados não seguem distr normal 

#Vital Status (Categórica) vs. Age at Index (Numérica)
#Objetivo: Verificar se a idade média dos pacientes difere entre os grupos "vivos" e "mortos".
#Teste de Mann-Whitney (se a distribuição for não paramétrica).
#Hipóteses:

#H0 (Hipótese nula): A média da idade dos pacientes vivos é igual à média dos pacientes mortos.

#H1 (Hipótese alternativa): A média da idade dos pacientes vivos é diferente da média dos pacientes mortos.

wilcox.test(age_at_index ~ vital_status, data = metadata_subset)

#Age at Index (Numérica) vs. Tissue Type (Categórica)
#Objetivo: Verificar se a idade média dos pacientes difere entre os grupos de tecido normal e tumoral.

#Teste sugerido: Teste t de Student (se a distribuição for normal) ou Teste de Mann-Whitney (se não for normal).

#Hipóteses:
  
#H0 (Hipótese nula): A média da idade dos pacientes é igual entre os tecidos normais e tumorais.

#H1 (Hipótese alternativa): A média da idade dos pacientes difere entre os tecidos normais e tumorais.

wilcox.test(age_at_index ~ tissue_type, data = metadata_subset)


#Age at Index (Numérica) vs. Primary Diagnosis (Categórica)
#Objetivo: Verificar se a idade média dos pacientes difere entre os diagnósticos primários.

#Teste sugerido: ANOVA ou Teste de Kruskal-Wallis (se os dados não forem normalmente distribuídos).

  
#H0 (Hipótese nula): A média da idade dos pacientes é igual entre todos os diagnósticos primários.

#H1 (Hipótese alternativa): A média da idade dos pacientes difere entre pelo menos dois diagnósticos primários.
kruskal.test(age_at_index ~ primary_diagnosis, data = diagnosis_filtered)


#, comparar a idade entre os diferentes tipos de tumor (classification_of_tumor)
#H0: As idades entre os diferentes tipos de tumor são iguais (não há diferença significativa).

#H1: Pelo menos um tipo de tumor tem uma idade média diferente.

kruskal.test(age_at_index ~ classification_of_tumor, data = tumor_class)

#H0: Não há associação entre vital_status e classification_of_tumor (as duas variáveis são independentes).

#H1: Existe uma associação significativa entre vital_status e classification_of_tumor.

# Criar uma tabela de contingência
contingency_table <- table(tumor_class$vital_status, tumor_class$classification_of_tumor)

# Aplicar o teste de qui-quadrado
chisq.test(contingency_table)
