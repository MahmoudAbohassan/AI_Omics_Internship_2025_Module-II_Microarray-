# Install and load required packages
if (!requireNamespace("GEOquery", quietly = TRUE)) install.packages("GEOquery")
if (!requireNamespace("affy", quietly = TRUE)) install.packages("affy")
if (!requireNamespace("arrayQualityMetrics", quietly = TRUE)) install.packages("arrayQualityMetrics")
if (!requireNamespace("limma", quietly = TRUE)) install.packages("limma")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(limma)
library(ggplot2)


# 1. Load dataset from GEO

gse <- getGEO("GSE13911", GSEMatrix = TRUE)
gse <- gse[[1]]


# 2. Perform Quality Control (before normalization)

exprs_raw <- exprs(gse)
pheno <- pData(gse)

# Generate quality control report (before normalization)
arrayQualityMetrics(expressionset = gse,
                    outdir = "QC_before_norm",
                    force = TRUE)


# 3. Normalize the data

exprs_norm <- normalizeBetweenArrays(exprs_raw, method = "quantile")

# Replace raw expression values with normalized values
exprs(gse) <- exprs_norm

# Generate quality control report (after normalization)
arrayQualityMetrics(expressionset = gse,
                    outdir = "QC_after_norm",
                    force = TRUE)

# 4. Filtering low-intensity probes

cutoff <- apply(exprs_norm, 1, function(x) mean(x) > 5)
exprs_filtered <- exprs_norm[cutoff, ]

cat("Number of transcripts remaining after filtering:", nrow(exprs_filtered), "\n")


# 5. Define phenotype groups (Normal vs Cancer)

pheno$Group <- ifelse(grepl("normal", pheno$source_name_ch1, ignore.case = TRUE),
                      "Normal", "Cancer")

table(pheno$Group)


# 6. Principal Component Analysis (PCA)

pca <- prcomp(t(exprs_filtered), scale. = TRUE)

pca_data <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = pheno$Group
)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA plot (GSE13911)", x = "PC1", y = "PC2") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank()
  )