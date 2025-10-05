# Microarray Preprocessing Workflow – GSE13911

This repository contains an R script for **preprocessing microarray gene expression data** as part of the **AI and Omics Research Internship (2025)**.  
The workflow uses the publicly available dataset **GSE13911**, which includes samples from **gastric mucosa (normal)** and **gastric adenocarcinoma (cancer)** tissues.

---

# Objectives

1. Perform **Quality Control (QC)** before and after normalization  
2. Apply **RMA / Quantile Normalization** to correct technical variations  
3. **Filter** low-intensity probes to remove background noise  
4. Define **phenotype groups** (Normal vs Cancer)  
5. Visualize data variation using **Principal Component Analysis (PCA)**

---

# Workflow Overview

| Step | Description |
|------|--------------|
| **1. Data Download** | Retrieve dataset `GSE13911` from NCBI GEO using `GEOquery` |
| **2. QC (Before Normalization)** | Detect outlier arrays using `arrayQualityMetrics` |
| **3. Normalization** | Apply quantile normalization (`normalizeBetweenArrays`) |
| **4. QC (After Normalization)** | Ensure intensity distributions are uniform |
| **5. Filtering** | Remove low-expressed probes (mean expression < 5) |
| **6. Grouping** | Label samples as *Normal* or *Cancer* |
| **7. PCA** | Visualize clustering between groups |

---

# Required R Packages

```r
GEOquery
affy
arrayQualityMetrics
limma
ggplot2


# Author

Mahmoud Ahmed Abohassan
AI and Omics Research Internship (2025)
Contact: mahmoudabohassan03@gmail.com


