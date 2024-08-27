library(dplyr)
library(ggplot2)
library(WGCNA)
library(tidyverse)
library(granulator)

###Load required functions
remove_duplicate_genes<-function(eset, column_of_symbol, method = "mean"){
  eset<-as.data.frame(eset)
  rownames(eset)<-NULL
  dups <- dim(eset)[1] - length(unique(eset[,column_of_symbol]))
  if(dups==0){
    eset<-tibble:: column_to_rownames(eset,var = column_of_symbol)
    return(eset)
  }else{
    if(method=="mean"){
      order_index=apply(eset[,setdiff(colnames(eset),column_of_symbol)],1,function(x) mean(x,na.rm=T))
      eset<-eset[order(order_index,decreasing=T),]
      eset<-eset %>%dplyr:: distinct(!!sym(column_of_symbol),.keep_all = TRUE) %>%
        tibble:: column_to_rownames(.,var = column_of_symbol)
      return(eset)
    }else if(method == "sd"){
      order_index = apply(eset[,setdiff(colnames(eset),column_of_symbol)],1,function(x) sd(x,na.rm=T))
      eset<-eset[order(order_index,decreasing=T),]
      eset<-eset %>% distinct(!!sym(column_of_symbol),.keep_all = TRUE) %>%
        tibble:: column_to_rownames(.,var = column_of_symbol)
      return(eset)
    }
  }
}

Convert_To_ENSG <- function(df) {
  gene_info <- read.csv("~/OneDrive - AdventHealth/Single cell seq/BIOage/CogentAP output/gene_info_incl_introns.csv")
  gene_info <- gene_info %>%
    dplyr::select(2, 3)  ##subset to only relevant columns
  df$Gene_Name <- rownames(df)
  df <- dplyr::left_join(df, gene_info, by = c("Gene_Name"))
  df <- remove_duplicate_genes(df, column_of_symbol = 'Gene_Name', method = 'mean') ##remove duplicate genes that have the same symbol based on which one has the greatest mean 
  rownames(df) <- df$Ensembl_ID  ##make ENSG rownames
  sum(duplicated(df$Gene_ID))
  df$Gene_ID[duplicated(sig$Gene_ID)]
  df <- subset(df, select = -c(Ensembl_ID))  ###Remove columns with Gene_ID and Gene_name
}


####To run the deconvolution pipeline you will need to download the data files from the data folder and save it in your local drive
##To see how the data files were generated please view the source data generation file.
##The gene_info_incl_introns.csv is a direct output from TakaraBio's CogentAP pipeline.

##In this example we will deconvolute the Petersen et al (2024) data which is presented in the manuscript and can be downloaded from GSE 
##Paper https://www.cell.com/cell-metabolism/fulltext/S1550-4131(24)00081-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1550413124000810%3Fshowall%3Dtrue
##GSE https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244118

bulk <- read.delim("~/OneDrive - AdventHealth/Single cell seq/Deconvolution/Petersen 2024/GSE244118_abdominal.fat_all.gene_CPM.txt", 
                  sep = '\t')
dim(bulk)
bulk2 <- bulk[,c(1, 8:60)]
bulk2[1:6, 1:6]
bulk2 <- remove_duplicate_genes(bulk2, column_of_symbol = 'ensembl_gene_id', method = 'mean')
mat <- as.matrix(bulk2)


##Load the signature matrix, this is currently unfiltered
sig <- read.csv("~/OneDrive - AdventHealth/Single cell seq/Deconvolution/BioAge harmony/BA.signature.matrix.normalized.csv", row.names = 1)
head(sig)
sig[1:6, 1:5]

sig <- Convert_To_ENSG(sig)

##In our paper we had the best results from 6000 HVG list so we will subset our signature matrix to these 6000 HVG
hvg6000 <- read.csv("~/Library/CloudStorage/OneDrive-AdventHealth/Single cell seq/Deconvolution/BioAge harmony/6000HVG.csv", row.names = 1)
hvg6000s <- hvg6000$x 

sig2 <- sig[rownames(sig) %in% hvg6000s, ]
sig2 <-  as.matrix(sig2)
dim(sig2)
sig2[1:6, 1:6]


###Run deconvolution
decon <- deconvolute(m = mat, sigMatrix = sig2, methods = c("dtangle"))
str(decon)

##extract proportions
props <- decon$proportions$dtangle_sig1
props
##plot proportions
plot_proportions(deconvoluted = decon, method = 'dtangle', signature = 'sig1') + 
  theme(axis.text.x=element_blank(),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(colour = "black", size = 12)) +
  labs(fill = "Cell Type", title = '6000 HVG') +scale_fill_brewer(palette = "Paired")

##Save sample type proportions
write.csv(props, 'cell.type.proportions.csv')
