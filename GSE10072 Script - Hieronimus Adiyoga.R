# Load Package 

library(GEOquery)
library(limma)
library(hgu133a.db)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)

# Download Dataset

options(download.file.method = "wininet")

gset <- getGEO("GSE10072", GSEMatrix = TRUE)
gset <- gset[[1]]

# Expression Matrix 

exprSet <- exprs(gset)
dim(exprSet)

# Define Group 

pheno <- pData(gset)

table(pheno$characteristics_ch1)

# Check Metadata Structure 
# 1st column consists of male and female 

colnames(pheno)
pheno[, grep("characteristics", colnames(pheno))]
table(pheno$characteristics_ch1.2)

# Check Sourcename

table(pheno$source_name_ch1)

# Define group 

group <- ifelse(grepl("normal",
                      pheno$source_name_ch1,
                      ignore.case = TRUE),
                "Normal",
                "Tumor")

group <- factor(group)

table(group)

# Matrix Model 

design <- model.matrix(~ group)
fit <- lmFit(exprSet, design)
fit <- eBayes(fit)

results <- topTable(fit,
                    coef = 2,
                    number = Inf,
                    adjust.method = "BH")

# Check 10 top genes

head(results)

# Check significant genes

deg <- subset(results,
              adj.P.Val < 0.05 & abs(logFC) > 1)

nrow(deg)

# Check expression types 

up   <- subset(deg, logFC > 1)
down <- subset(deg, logFC < -1)

nrow(up)
nrow(down)

colnames(design)

# Mapping Probe

library(hgu133a.db)
library(AnnotationDbi)

results$SYMBOL <- mapIds(hgu133a.db,
                         keys = rownames(results),
                         column = "SYMBOL",
                         keytype = "PROBEID",
                         multiVals = "first")

results <- results[!is.na(results$SYMBOL), ]

# Collapse Duplicate Gene 

results <- results[!duplicated(results$SYMBOL), ]

# Final DEG Table 
deg <- subset(results,
              adj.P.Val < 0.05 & abs(logFC) > 1)

nrow(deg)

# Save 
write.csv(deg, "DEG_GSE10072.csv")

# Volcano Plot 
library(ggplot2)

results$threshold <- "Not Significant"
results$threshold[results$adj.P.Val < 0.05 & results$logFC > 1]  <- "Up"
results$threshold[results$adj.P.Val < 0.05 & results$logFC < -1] <- "Down"

ggplot(results, aes(x = logFC,
                    y = -log10(adj.P.Val),
                    color = threshold)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  scale_color_manual(values = c("blue","grey","red")) +
  ggtitle("Volcano Plot - GSE10072")

# Heatmap Top 50 Genes 
install.packages("pheatmap")
library(pheatmap)

nrow(deg)
top50 <- rownames(deg)[1:50]
length(top50)

annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(exprSet)

all(rownames(annotation_col) == colnames(exprSet))

heatmap_data <- exprSet[top50, , drop = FALSE]

pheatmap(heatmap_data,
         scale = "row",
         annotation_col = annotation_col)

# GO Enrichment 
library(clusterProfiler)
library(org.Hs.eg.db)

# Convert Symbol 

gene_df <- bitr(deg$SYMBOL,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

# GO Analysis 

ego <- enrichGO(gene          = gene_df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05)

dotplot(ego, showCategory = 15)

# KEGG Pathway 

ekegg <- enrichKEGG(gene         = gene_df$ENTREZID,
                    organism     = "hsa",
                    pvalueCutoff = 0.05)

dotplot(ekegg, showCategory = 15)

# Data Retrieval

nrow(deg)
nrow(up)
nrow(down)

# GO KEGG Interpretation

head(ego)
head(ekegg)