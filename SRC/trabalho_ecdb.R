install.packages(c("data.table", "tidyverse", "Bioconductor", "TCGAbiolinks", "edgeR", "DESeq2"))
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")

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
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyr)


#Buscar informações do The Cancer Genome Atlas

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


rna_data
colnames(rna_data) # lists assays
rownames(rna_data) #lists all genes

#Converter rna_data para um dataframe 

rna_counts <- as.data.frame(assay(rna_data))  # matrix counting
metadata <- as.data.frame(colData(rna_data))  # Clinic metadata 


#Exploração inicial de dados 
#Verificar estrutura dos dados 

colnames(rna_counts) 
rownames(rna_counts)
colnames(metadata)

#Verificar dimensões (genes x assays)

dim(rna_counts)
dim(metadata)

#Checar existência de valores omissos 

colSums(is.na(rna_counts)) 
colSums(is.na(metadata)) 

#Verificar porcentagem de valores omissos em cada coluna dos metadados

metadata_na_percent <- colSums(is.na(metadata)) / nrow(metadata) * 100
metadata_na_percent

#Remover colunas com valores omissos > 10%

metadata_clean <- metadata[, colMeans(is.na(metadata)) <= 0.1] 

#filtragem da contagem de dados em RNA 
rna_filtered <- rna_counts[apply(rna_counts,1 , sd) > 0, ] #remove lines where expression in is 0
rna_normalized <- t(scale(t(rna_filtered))) #this way we don't produce NaNs 
colSums(is.na(rna_normalized))

#Checar se colunas e linhas de rna_counts e metadata estão correspondentes
all(colnames(rna_normalized) %in% rownames(metadata_clean))
#data standardization - applied z-score to each row(gene)
#summary(rna_counts_scaled)
#apply(rna_normalized, 1, mean)  # Média de cada gene (deve ser ~0)
#apply(rna_normalized, 1, sd)    # Desvio padrão de cada gene (deve ser ~1)


#Checar se os tipos de dados estão corretos 
str(rna_counts)
str(metadata)

#Análise exploratória

#Verificar a quantidade de categorias em cada coluna
table(colData(rna_data)$tissue_type)
table(colData(rna_data)$classification_of_tumor)
table(colData(rna_data)$primary_diagnosis)
table(colData(rna_data)$vital_status)

#selecção das variáveis de interesse: 
  #vital_status, age_at_index, classification_of_tumor, tissue_type, primary_diagnosis (Endometrioid adenocarcinoma, NOS (404), Serous cystadenocarcinoma, NOS (124))

#Criar nova tabela metadados só com essas variaveis de interesse

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
ggplot(metadata_subset[!is.na(metadata_subset$vital_status), ], aes(x = vital_status)) +
  geom_bar(fill = "lightcoral") +
  theme_minimal() +
  labs(title = "Distribuição: Vital Status", x = "", y = "Frequência") +
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

#Testes de Hipótese 


#H0 (Hipótese nula): Os dados seguem uma distribuição normal.
#H1 (Hipótese alternativa): Os dados não seguem uma distribuição normal.

shapiro.test(metadata_subset$age_at_index)



#Vital Status (Categórica) vs. Age at Index (Numérica)

#H0 (Hipótese nula): A média da idade dos pacientes vivos é igual à média dos pacientes mortos.
#H1 (Hipótese alternativa): A média da idade dos pacientes vivos é diferente da média dos pacientes mortos.

wilcox.test(age_at_index ~ vital_status, data = metadata_subset)

#Age at Index (Numérica) vs. Tissue Type (Categórica)

#H0 (Hipótese nula): A média da idade dos pacientes é igual entre os tecidos normais e tumorais.
#H1 (Hipótese alternativa): A média da idade dos pacientes difere entre os tecidos normais e tumorais.

wilcox.test(age_at_index ~ tissue_type, data = metadata_subset)

#Age at Index (Numérica) vs. Primary Diagnosis (Categórica)

#H0 (Hipótese nula): A média da idade dos pacientes é igual entre todos os diagnósticos primários.
#H1 (Hipótese alternativa): A média da idade dos pacientes difere entre pelo menos dois diagnósticos primários.

kruskal.test(age_at_index ~ primary_diagnosis, data = diagnosis_filtered)


#Age at Index vs os classification of tumor)

#H0: As idades entre os diferentes tipos de tumor são iguais (não há diferença significativa).
#H1: Pelo menos um tipo de tumor tem uma idade média diferente.

kruskal.test(age_at_index ~ classification_of_tumor, data = tumor_class)

#Vital Status vs. Classification of Tumor

#H0: Não há associação entre vital_status e classification_of_tumor (as duas variáveis são independentes).
#H1: Existe uma associação significativa entre vital_status e classification_of_tumor.

# Criar uma tabela de contingência
contingency_table <- table(tumor_class$vital_status, tumor_class$classification_of_tumor)

# Aplicar o teste de qui-quadrado
chisq.test(contingency_table)



#Expressao diferencial

#vital_status x EA
# 1. Filtrar metadados: apenas Endometrioid adenocarcinoma, NOS e sem NAs
EA_meta <- metadata_subset %>%
  filter(primary_diagnosis == "Endometrioid adenocarcinoma, NOS" & !is.na(vital_status))

# 2. Garantir que os nomes das colunas em contagens correspondem às amostras
rna_EA <- rna_filtered[, rownames(EA_meta)]

# 3. Filtrar genes: manter genes com contagem ≥ 20 em pelo menos 4 amostras
rna_filtered_EA <- rowSums(rna_EA >= 20) >= 4
rna_counts_filtrado <- rna_EA[rna_filtered_EA, ]

# 4. Criar objeto DESeq2
dds_EA <- DESeqDataSetFromMatrix(countData = rna_counts_filtrado,
                                 colData = EA_meta,
                                 design = ~ vital_status)

# 5. Transformar caracteres em fatores
dds_EA$vital_status <- factor(dds_EA$vital_status)

# 6. Executar DESeq2
dds_EA <- DESeq(dds_EA)

# 7. Resultados: comparar Dead vs Alive (Dead - referência por padrão)
res_EA <- results(dds_EA, contrast = c("vital_status", "Dead", "Alive"))

# 8. Ver resumo e visualizar top genes
summary(res_EA)
head(res_EA[order(res_EA$pvalue), ])

sum(res_EA$padj < 0.1, na.rm = TRUE)


#plotar gráfico de DE
DESeq2::plotMA(res_EA, main="DESeq2", ylim=c(-6,6))

res_df <- as.data.frame(res_EA) 
res_df$gene <- rownames(res_EA)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significativo", "Não significativo")

# Transformar o objeto de resultados em um data frame
res_df_EA <- as.data.frame(res_EA)
res_df_EA$gene <- rownames(res_df_EA)

# Salvar a análise diferencial de expressão do EA
write.csv(res_df_EA, "DEG_results_EA.csv", row.names = FALSE)


# Criar a coluna 'regulation' com base em critérios de significância
res_df_EA <- res_df %>%
  mutate(regulation = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "Not Sig"
  ))


# Selecionar os genes mais significativos
top_genes_EA <- res_df_EA %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 20)

# Plotar gráfico de barras
ggplot(top_genes_EA, aes(x = reorder(gene, log2FoldChange), y = log2FoldChange, fill = regulation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "lightcoral", "Down" = "lightgreen")) +
  theme_minimal() +
  labs(title = "Top 20 genes diferencialmente expressos (Dead vs Alive) para Endometrioid adenocarcinoma",
       x = "Gene",
       y = "log2 Fold Change")

# Mapear os IDs ENSEMBL para o nome completo do gene para top_genes_EA
ensembl_ids_EA <- top_genes_EA$gene

# Limpar as versões dos IDs ENSEMBL
clean_ids_EA <- sub("\\..*", "", ensembl_ids_EA)

# Mapear os IDs ENSEMBL para o nome completo do gene
gene_names_EA <- mapIds(org.Hs.eg.db, keys = clean_ids_EA, column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")

# Adicionar os nomes dos genes ao dataframe
top_genes_EA$gene_name <- gene_names_EA

# Visualizar o resultado
head(top_genes_EA)

# Plotar o gráfico com os nomes dos genes, removendo os NAs diretamente
ggplot(na.omit(top_genes_EA), aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = regulation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "lightcoral", "Down" = "lightgreen")) +
  theme_minimal() +
  labs(title = "Top 20 genes diferencialmente expressos (Dead vs Alive) para Endometrial adenocarcinoma",
       x = "Gene",
       y = "log2 Fold Change")

#vital_status x SC 

# 1. Filtrar metadados: apenas Endometrioid adenocarcinoma, NOS e sem NAs
SC_meta <- metadata_subset %>%
  filter(primary_diagnosis == "Serous cystadenocarcinoma, NOS" & !is.na(vital_status))

# 2. Garantir que os nomes das colunas em contagens correspondem às amostras
rna_SC <- rna_filtered[, rownames(SC_meta)]

# 3. Filtrar genes: manter genes com contagem ≥ 20 em pelo menos 4 amostras
rna_filtered_SC <- rowSums(rna_SC >= 20) >= 4
rna_counts_filtrado_SC <- rna_SC[rna_filtered_SC, ]

# 4. Criar objeto DESeq2
dds_SC <- DESeqDataSetFromMatrix(countData = rna_counts_filtrado_SC,
                                 colData = SC_meta,
                                 design = ~ vital_status)

# 5. Transformar caracteres em fatores
dds_SC$vital_status <- factor(dds_SC$vital_status)

# 6. Executar DESeq2
dds_SC <- DESeq(dds_SC)

# 7. Resultados: comparar Dead vs Alive (Dead - referência por padrão)
res_SC <- results(dds_SC, contrast = c("vital_status", "Dead", "Alive"))

# 8. Ver resumo e visualizar top genes
summary(res_SC)
head(res_SC[order(res_SC$pvalue), ])

sum(res_SC$padj < 0.1, na.rm = TRUE)


#plotar gráfico de DE
DESeq2::plotMA(res_SC, main="DESeq2", ylim=c(-6,6))

res_df_SC <- as.data.frame(res_SC) 
res_df_SC$gene <- rownames(res_SC)
res_df_SC$significant <- ifelse(res_df_SC$padj < 0.05 & abs(res_df_SC$log2FoldChange) > 1, "Significativo", "Não significativo")

# Salvar a análise diferencial de expressão do SC
write.csv(res_df_SC, "DEG_results_SC.csv", row.names = FALSE)

# Transformar o objeto de resultados em um data frame
res_df_SC <- as.data.frame(res_SC)
res_df_SC$gene <- rownames(res_df_SC)

# Criar a coluna 'regulation' com base em critérios de significância
res_df_SC <- res_df_SC %>%
  mutate(regulation = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "Not Sig"
  ))


# Selecionar os genes mais significativos
top_genes_SC <- res_df_SC %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 20)

# Plotar gráfico de barras
ggplot(top_genes_SC, aes(x = reorder(gene, log2FoldChange), y = log2FoldChange, fill = regulation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "lightcoral", "Down" = "lightgreen")) +
  theme_minimal() +
  labs(title = "Top 20 genes diferencialmente expressos (Dead vs Alive) para Serous cystadenocarcinoma",
       x = "Gene",
       y = "log2 Fold Change")


#Plotar gráfico com os top 20 genes para SC

ensembl_ids_SC <- top_genes_SC$gene

#Limpar as versões dos IDs ENSEMBL
clean_ids_SC <- sub("\\..*", "", ensembl_ids_SC)

# Mapear os IDs ENSEMBL para o nome completo do gene (GENENAME)
gene_names_SC <- mapIds(org.Hs.eg.db, keys = clean_ids_SC, column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")

# Adicionar os nomes dos genes ao dataframe
top_genes_SC$gene_name <- gene_names_SC

# Visualizar o resultado
head(top_genes_SC)

# Plotar o gráfico com os nomes dos genes
ggplot(na.omit(top_genes_SC), aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = regulation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "lightcoral", "Down" = "lightgreen")) +
  theme_minimal() +
  labs(title = "Top 20 genes diferencialmente expressos (Dead vs Alive) para Serous cystadenocarcinoma",
       x = "Gene",
       y = "log2 Fold Change")

############################################################################################

##Enrequicimento

#Análise de Enriquecimento Funcional para EA

# Filtrar genes significativos (p < 0.05 e |log2FC| > 1)
EA_genes_sig <- res_df_EA %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
EA_genes_ids <- sub("\\..*", "", EA_genes_sig$gene)

# Converter ENSEMBL para ENTREZ ID
EA_entrez_ids <- mapIds(org.Hs.eg.db,
                        keys = EA_genes_ids,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")

# Remover NAs
EA_entrez_ids <- na.omit(EA_entrez_ids)

# Enriquecimento GO - Biological Process
EA_GO_enrich <- enrichGO(gene = EA_entrez_ids,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)


#Guardar EA
EA_go_df <- as.data.frame(EA_GO_enrich)
head(EA_go_df)
write.csv(EA_go_df, file = "GO_enrichment_EA.csv", row.names = FALSE)

#Enriquecimento para SC

SC_genes_sig <- res_df_SC %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
SC_genes_ids <- sub("\\..*", "", SC_genes_sig$gene)

SC_entrez_ids <- mapIds(org.Hs.eg.db,
                        keys = SC_genes_ids,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")

SC_entrez_ids <- na.omit(SC_entrez_ids)

SC_GO_enrich <- enrichGO(gene = SC_entrez_ids,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)

barplot(SC_GO_enrich, showCategory = 15, title = "GO BP - SC")

# Guardar em CSV
SC_go_df <- as.data.frame(SC_GO_enrich)
write.csv(SC_go_df, file = "GO_enrichment_SC.csv", row.names = FALSE)

##KEGG

EA_KEGG <- enrichKEGG(gene = EA_entrez_ids,
                      organism = 'hsa',
                      keyType = "ncbi-geneid",
                      pvalueCutoff = 0.05)

barplot(EA_KEGG, showCategory = 10, title = "KEGG - EA")

SC_KEGG <- enrichKEGG(gene = SC_entrez_ids,
                      organism = 'hsa',
                      keyType = "ncbi-geneid",
                      pvalueCutoff = 0.05)

barplot(SC_KEGG, showCategory = 10, title = "KEGG - SC")

#Guardar como csv EA
EA_KEGG_df <- as.data.frame(EA_KEGG)
write.csv(EA_KEGG_df, "KEGG_enrichment_EA.csv", row.names = FALSE)

#Guardar como csv SC
SC_KEGG_df <- as.data.frame(SC_KEGG)
write.csv(SC_KEGG_df, "KEGG_enrichment_SC.csv", row.names = FALSE)

head(SC_go_df)

# Tabela SC

# 1. Limpar IDs ENSEMBL no res_df_SC
res_df_SC$ensembl_clean <- sub("\\..*", "", res_df_SC$gene)

# 2. Mapear ENSEMBL → SYMBOL
res_df_SC$symbol <- mapIds(org.Hs.eg.db,
                           keys = res_df_SC$ensembl_clean,
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")

# 3. Separar os genes em cada termo GO
SC_GO_df_split <- as.data.frame(SC_GO_enrich) %>%
  filter(p.adjust < 0.05) %>%
  separate_rows(geneID, sep = "/") %>%
  rename(symbol = geneID)

# 4. Juntar com log2FoldChange dos genes
SC_GO_merged <- SC_GO_df_split %>%
  left_join(res_df_SC, by = "symbol")

# 5. Calcular a média do log2FoldChange por termo GO
SC_GO_summary <- SC_GO_merged %>%
  group_by(ID, Description) %>%
  summarise(
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    sd_log2FC = sd(log2FoldChange, na.rm = TRUE),
    Count = n(),
    p.adjust = first(p.adjust)
  ) %>%
  arrange(desc(abs(mean_log2FC)))  # ordenar pela magnitude da média

# 6. Visualizar resultado
print(SC_GO_summary)


# Tabela EA

# 1. Limpar os ENSEMBL IDs do res_df_EA
res_df_EA$ensembl_clean <- sub("\\..*", "", res_df_EA$gene)

# 2. Mapear ENSEMBL → símbolo de gene (SYMBOL)
res_df_EA$symbol <- mapIds(org.Hs.eg.db,
                           keys = res_df_EA$ensembl_clean,
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")

# 3. Separar os geneIDs de cada termo GO (enrichGO EA)
EA_GO_df_split <- as.data.frame(EA_GO_enrich) %>%
  filter(p.adjust < 0.05) %>%
  separate_rows(geneID, sep = "/") %>%
  rename(symbol = geneID)

# 4. Juntar os símbolos com os log2FoldChange
EA_GO_merged <- EA_GO_df_split %>%
  left_join(res_df_EA, by = "symbol")

# 5. Calcular a média e o desvio padrão de log2FC por termo GO
EA_GO_summary <- EA_GO_merged %>%
  group_by(ID, Description) %>%
  summarise(
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    sd_log2FC = sd(log2FoldChange, na.rm = TRUE),
    Count = n(),
    p.adjust = first(p.adjust)
  ) %>%
  arrange(desc(abs(mean_log2FC)))  # ordenar pela magnitude (positivo ou negativo)

# 6. Ver resultado no console
print(EA_GO_summary)

