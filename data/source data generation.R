## create data source for deconvolution
library(Seurat)
library(ggplot2)
library(usethis)
library(dplyr)
library(readr)
library(tidyverse)

###Signature matrix profile are generated from our full-length snRNA-Seq data set WAT aging data set
###https://onlinelibrary.wiley.com/doi/10.1111/acel.14287

##Load functions
Convert_To_ENSG2 <- function(df) {
  gene_info <- read.csv("~/OneDrive - AdventHealth/Single cell seq/BIOage/CogentAP output/gene_info_incl_introns.csv")
  gene_info <- gene_info %>%
    dplyr::select(2, 3)  ##subset to only relevant columns
  df <- dplyr::left_join(df, gene_info, by = c("Gene_Name"))
  Ensemble_ID <- df$Ensembl_ID  ###Remove columns with Gene_ID and Gene_name
  return(Ensemble_ID)
}


##Load finalized seurat object
BA <- readRDS("~/Library/CloudStorage/OneDrive-AdventHealth/Single cell seq/BIOage/decontx harmony integration/BioAge.harmony.sct.V3.RDS")


###Create signature matrix for all cell types using filtered assay and normalized data
Sig_matrix <- AggregateExpression(BA, return.seurat = FALSE, assay = 'filtered', slot = 'data', group.by = "CellType")
Sig_matrix <- Sig_matrix$filtered
dim(Sig_matrix)
write.csv(Sig_matrix, 'BA.signature.matrix.normalized.csv')


###Calculate HVG of different iterations
##start with merged but unintegrated object to define which HVG are used for clustering and therefore can be used to optimize signature matrix
Merge <- readRDS("~/OneDrive - AdventHealth/Single cell seq/BIOage/Batch correction comparison/decontx/seurat_Unadjusted_BioAgedecont.RDS")


merge.list <- SplitObject(Merge, split.by="orig.ident")
merge.list <- lapply(X = merge.list, 
                     FUN = SCTransform, 
                     method = "glmGamPoi", 
                     return.only.var.genes = FALSE)

###The results showed that 6000 HVG produced the best results, but different iterations can be formed here
Gene_Name <- SelectIntegrationFeatures(object.list = merge.list, nfeatures = 6000)
df <- as.data.frame(Gene_Name)
df <- Convert_To_ENSG2(df)
write.csv(df, "6000HVG.csv")




