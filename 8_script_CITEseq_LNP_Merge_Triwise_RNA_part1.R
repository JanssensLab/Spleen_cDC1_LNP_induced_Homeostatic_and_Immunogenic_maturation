## Script for analyzing LNP cDC1 CITE-seq project data
## Triwise analysis RNA part 1

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library('ggplot2')
library('openxlsx')
library('cowplot')
library('muscat')
library('purrr')
library('tidyverse')
library('Matrix')
library('DESeq2')
library('triwise')
library('ggrepel')
library("htmlwidgets")
library('openxlsx')

################################################################################
## LOAD DATASET
################################################################################

setwd('/home/clintdn/VIB/DATA/Sophie/scRNA-seq_Victor/LNP_CITEseq_experiment/')

sampleFolder <- "VBO_LNP_merge/"
sampleName <- "VBO_4-12"
experiment <- sampleName

########################################
##### Some variables
########################################

listLabels<-list('VBO004','VBO005','VBO006','VBO007',
                 'VBO008','VBO009','VBO010','VBO011','VBO012')

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')
source('~/VIB/DATA/Roos/Daan 1/script_featurePlots.R')

### Modified function for robust edgeR
source('~/VIB/DATA/Roos/Cara/pbDSrob.R')
source('~/VIB/DATA/Roos/Cara/pbDS_DESeq2.R')

###First letter upper case
firstup <- function(x) {
  x<-tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

################################################################################
## TRIWISE COMPARING BETWEEN SAMPLES
################################################################################

## Read in muscat data
res<-readRDS("VBO_LNP_merge/results/Robjects/res_Muscat_DeSeq2_new_DS_Simple_RNA_Harmony_clustering_VBO_4-12.rds")
expression_2h <- assays(varianceStabilizingTransformation(res[["data"]][["Early mature cDC1s"]]))[[1]] #actually contains 2+8hr
expression_8h <- assays(varianceStabilizingTransformation(res[["data"]][["Late mature cDC1s"]]))[[1]] #actually contains 2+8hr

saveRDS(expression_2h,"VBO_LNP_merge/results/Robjects/Expr_2h_Early_mat_cDC1s_Muscat_DeSeq2_new_DS_Simple_RNA_Harmony_clustering_VBO_4-12.rds")
saveRDS(expression_8h,"VBO_LNP_merge/results/Robjects/Expr_8h_Late_mat_cDC1s_Muscat_DeSeq2_new_DS_Simple_RNA_Harmony_clustering_VBO_4-12.rds")

## Extra: (not used here)
res_ADT<-readRDS("VBO_LNP_merge/results/Robjects_ADT_muscat/res_Muscat_DeSeq2_new_ADT_DS_Simple_RNA_Harmony_clustering_VBO_4-12.rds") ##DESeq2
expression_2h_ADT <- assays(varianceStabilizingTransformation(res_ADT[["data"]][["Early mature cDC1s"]]))[[1]] #actually contains 2+8hr
expression_8h_ADT <- assays(varianceStabilizingTransformation(res_ADT[["data"]][["Late mature cDC1s"]]))[[1]] #actually contains 2+8hr

saveRDS(expression_2h_ADT,"VBO_LNP_merge/results/Robjects/Expr_2h_ADT_Early_mat_cDC1s_Muscat_DeSeq2_new_DS_Simple_RNA_Harmony_clustering_VBO_4-12.rds")
saveRDS(expression_8h_ADT,"VBO_LNP_merge/results/Robjects/Expr_8h__ADT_Late_mat_cDC1s_Muscat_DeSeq2_new_DS_Simple_RNA_Harmony_clustering_VBO_4-12.rds")

##### Get columns
cols_eLNP_8h<-grep("VBO005",colnames(expression_8h))
cols_pIC_LNP_8h<-grep("VBO006",colnames(expression_8h))
cols_CpG_LNP_8h<-grep("VBO007",colnames(expression_8h))
cols_pIC_alone_8h<-grep("VBO008",colnames(expression_8h))

cols_eLNP_2h<-grep("VBO009",colnames(expression_2h))
cols_pIC_LNP_2h<-grep("VBO010",colnames(expression_2h))
cols_CpG_LNP_2h<-grep("VBO011",colnames(expression_2h))
cols_pIC_alone_2h<-grep("VBO012",colnames(expression_2h))

##### Create expTableMean -> filter by sample (hr)
eLNP_8h_mean<-apply(expression_8h[,cols_eLNP_8h],1,mean)
pIC_LNP_8h_mean<-apply(expression_8h[,cols_pIC_LNP_8h],1,mean)
CpG_LNP_8h_mean<-apply(expression_8h[,cols_CpG_LNP_8h],1,mean)
pIC_alone_8h_mean<-apply(expression_8h[,cols_pIC_alone_8h],1,mean)

expTable_8h_mean_v2<-cbind(eLNP_8h_mean,pIC_LNP_8h_mean,CpG_LNP_8h_mean)
colnames(expTable_8h_mean_v2)<-c('eLNP_8h','pIC_LNP_8h','CpG_LNP_8h')

expTable_8h_mean_v3<-cbind(pIC_LNP_8h_mean,CpG_LNP_8h_mean,pIC_alone_8h_mean)
colnames(expTable_8h_mean_v3)<-c('pIC_LNP_8h','CpG_LNP_8h','pIC_alone_8h')

expTable_8h_mean_v4<-cbind(eLNP_8h_mean,pIC_LNP_8h_mean,pIC_alone_8h_mean)
colnames(expTable_8h_mean_v4)<-c('eLNP_8h','pIC_LNP_8h','pIC_alone_8h')

eLNP_2h_mean<-apply(expression_2h[,cols_eLNP_2h],1,mean)
pIC_LNP_2h_mean<-apply(expression_2h[,cols_pIC_LNP_2h],1,mean)
CpG_LNP_2h_mean<-apply(expression_2h[,cols_CpG_LNP_2h],1,mean)
pIC_alone_2h_mean<-apply(expression_2h[,cols_pIC_alone_2h],1,mean)

expTable_2h_mean_v2<-cbind(eLNP_2h_mean,pIC_LNP_2h_mean,CpG_LNP_2h_mean)
colnames(expTable_2h_mean_v2)<-c('eLNP_2h','pIC_LNP_2h','CpG_LNP_2h')

expTable_2h_mean_v3<-cbind(pIC_LNP_2h_mean,CpG_LNP_2h_mean,pIC_alone_2h_mean)
colnames(expTable_2h_mean_v3)<-c('pIC_LNP_2h','CpG_LNP_2h','pIC_alone_2h')

expTable_2h_mean_v4<-cbind(eLNP_2h_mean,pIC_LNP_2h_mean,pIC_alone_2h_mean)
colnames(expTable_2h_mean_v4)<-c('eLNP_2h','pIC_LNP_2h','pIC_alone_2h')

## Save exptables
write.xlsx(as.data.frame(expTable_2h_mean_v2), "VBO_LNP_merge/results/Triwise/ExpTable_mean_2h_EM_v2.xlsx", row.names = T)
write.xlsx(as.data.frame(expTable_2h_mean_v3), "VBO_LNP_merge/results/Triwise/ExpTable_mean_2h_EM_v3.xlsx", row.names = T)
write.xlsx(as.data.frame(expTable_2h_mean_v4), "VBO_LNP_merge/results/Triwise/ExpTable_mean_2h_EM_v4.xlsx", row.names = T)
write.xlsx(as.data.frame(expTable_8h_mean_v2), "VBO_LNP_merge/results/Triwise/ExpTable_mean_8h_LM_v2.xlsx", row.names = T)
write.xlsx(as.data.frame(expTable_8h_mean_v3), "VBO_LNP_merge/results/Triwise/ExpTable_mean_8h_LM_v3.xlsx", row.names = T)
write.xlsx(as.data.frame(expTable_8h_mean_v4), "VBO_LNP_merge/results/Triwise/ExpTable_mean_8h_LM_v4.xlsx", row.names = T)

######## CREATE BARYCOORDINATES
barycoords_8h_v2 <- transformBarycentric(expTable_8h_mean_v2)
barycoords_8h_v3 <- transformBarycentric(expTable_8h_mean_v3)
barycoords_8h_v4 <- transformBarycentric(expTable_8h_mean_v4)

barycoords_2h_v2 <- transformBarycentric(expTable_2h_mean_v2)
barycoords_2h_v3 <- transformBarycentric(expTable_2h_mean_v3)
barycoords_2h_v4 <- transformBarycentric(expTable_2h_mean_v4)

################################################################################
## TRIWISE CONTINUED FOR EITHER (first batch: only v2 and v3)
################################################################################

################
### DEG
################

## Update February 2024!! Switch to stricter list (clint suffix: p_adj.glb < 0.01, abs(logFC) > 1, baseMean > 50)!!!!!!!
CpG_LNP_8h_vs_eLNP_8h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_8h-eLNPs_8h_VBO_4-12_clint.rds")
eLNP_8h_vs_pIC_LNP_8h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_eLNPs_8h-pIC_LNPs_8h_VBO_4-12_clint.rds")

CpG_LNP_8h_vs_pIC_LNP_8h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_8h-pIC_LNPs_8h_VBO_4-12_clint.rds")

CpG_LNP_8h_vs_pIC_alone_8h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_8h-pIC_8h_VBO_4-12_clint.rds")
pIC_alone_8h_vs_pIC_LNP_8h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_pIC_8h-pIC_LNPs_8h_VBO_4-12_clint.rds")

Gdiffexp_CpG_LNP_8h_vs_eLNP_8h <- c(CpG_LNP_8h_vs_eLNP_8h$`Late mature cDC1s`$gene)
Gdiffexp_eLNP_8h_vs_pIC_LNP_8h <- c(eLNP_8h_vs_pIC_LNP_8h$`Late mature cDC1s`$gene)

Gdiffexp_CpG_LNP_8h_vs_pIC_LNP_8h <- c(CpG_LNP_8h_vs_pIC_LNP_8h$`Late mature cDC1s`$gene)

Gdiffexp_CpG_LNP_8h_vs_pIC_alone_8h <- c(CpG_LNP_8h_vs_pIC_alone_8h$`Late mature cDC1s`$gene)
Gdiffexp_pIC_alone_8h_vs_pIC_LNP_8h <- c(pIC_alone_8h_vs_pIC_LNP_8h$`Late mature cDC1s`$gene)

## Update February 2024!! Switch to stricter list (clint suffix: p_adj.glb < 0.01, abs(logFC) > 1, baseMean > 50)!!!!!!!
CpG_LNP_2h_vs_eLNP_2h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_2h-eLNPs_2h_VBO_4-12_clint.rds")
eLNP_2h_vs_pIC_LNP_2h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_eLNPs_2h-pIC_LNPs_2h_VBO_4-12_clint.rds")

CpG_LNP_2h_vs_pIC_LNP_2h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_2h-pIC_LNPs_2h_VBO_4-12_clint.rds")

CpG_LNP_2h_vs_pIC_alone_2h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_2h-pIC_2h_VBO_4-12_clint.rds")
pIC_alone_2h_vs_pIC_LNP_2h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_pIC_2h-pIC_LNPs_2h_VBO_4-12_clint.rds")

Gdiffexp_CpG_LNP_2h_vs_eLNP_2h <- c(CpG_LNP_2h_vs_eLNP_2h$`Early mature cDC1s`$gene)
Gdiffexp_eLNP_2h_vs_pIC_LNP_2h <- c(eLNP_2h_vs_pIC_LNP_2h$`Early mature cDC1s`$gene)

Gdiffexp_CpG_LNP_2h_vs_pIC_LNP_2h <- c(CpG_LNP_2h_vs_pIC_LNP_2h$`Early mature cDC1s`$gene)

Gdiffexp_CpG_LNP_2h_vs_pIC_alone_2h <- c(CpG_LNP_2h_vs_pIC_alone_2h$`Early mature cDC1s`$gene)
Gdiffexp_pIC_alone_2h_vs_pIC_LNP_2h <- c(pIC_alone_2h_vs_pIC_LNP_2h$`Early mature cDC1s`$gene)

##########################
### GENES TO COLOR ON PLOT
##########################

#NONE
genelist = NULL

#TOP UPREGULATED GENES
NRTOPGENES <- 25 #MARK TOP 25 GENES
genelist <- barycoords %>% tibble::rownames_to_column("gene_id") %>% top_n(NRTOPGENES, r) %>% pull(gene_id)

#OWN GENES OF INTEREST
genelist <- c("Ccl5")
C103_list<-read.xlsx("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Victor/2020_analysis/Immunity 2016 Ardouin-2_Sophie.xlsx",
                    sheet = "C103")
genelist<-firstup(C103_list$Gene.Symbol)
C49_list<-read.xlsx("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Victor/2020_analysis/Immunity 2016 Ardouin-2_Sophie.xlsx",
                     sheet = "C49")
genelist<-firstup(C49_list$Gene.Symbol)
C61_list<-read.xlsx("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Victor/2020_analysis/Immunity 2016 Ardouin-2_Sophie.xlsx",
                    sheet = "C61")
genelist<-firstup(C61_list$Gene.Symbol)

##########################
### GENES TO MARK ON PLOT
##########################

# GOI (GENES OF INTEREST) ARE TOP EXPRESSED OR OWN LIST

##########################
### GENERATE PLOT
##########################
dir.create("VBO_LNP_merge/results/Triwise")

## 8h
p1 <- plotDotplot(barycoords_8h_v2, Gdiffexp = Gdiffexp_CpG_LNP_8h_vs_eLNP_8h, Goi= genelist, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v2_Late_mature_cDC1s_CpG_LNP_8h_vs_eLNP_8h_paper_2024.pdf", width = 15, height = 15)
print(p1)
dev.off()

p2 <- plotDotplot(barycoords_8h_v2, Gdiffexp = Gdiffexp_eLNP_8h_vs_pIC_LNP_8h, Goi= genelist, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v2_Late_mature_cDC1s_eLNP_8h_vs_pIC_LNP_8h_paper_2024.pdf", width = 15, height = 15)
print(p2)
dev.off()

p3 <- plotDotplot(barycoords_8h_v2, Gdiffexp = Gdiffexp_CpG_LNP_8h_vs_pIC_LNP_8h, Goi= genelist, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v2_Late_mature_cDC1s_CpG_LNP_8h_vs_pIC_LNP_8h_paper_2024.pdf", width = 15, height = 15)
print(p3)
dev.off()

p4 <- plotDotplot(barycoords_8h_v3, Gdiffexp = Gdiffexp_CpG_LNP_8h_vs_pIC_LNP_8h, Goi= genelist, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v3_Late_mature_cDC1s_CpG_LNP_8h_vs_pIC_LNP_8h_paper_2024.pdf", width = 15, height = 15)
print(p4)
dev.off()

p5 <- plotDotplot(barycoords_8h_v3, Gdiffexp = Gdiffexp_CpG_LNP_8h_vs_pIC_alone_8h, Goi= genelist, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v3_Late_mature_cDC1s_CpG_LNP_8h_vs_pIC_alone_8h_paper_2024.pdf", width = 15, height = 15)
print(p5)
dev.off()

p6 <- plotDotplot(barycoords_8h_v3, Gdiffexp = Gdiffexp_pIC_alone_8h_vs_pIC_LNP_8h, Goi= genelist, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v3_Late_mature_cDC1s_pIC_alone_8h_vs_pIC_LNP_8h_paper_2024.pdf", width = 15, height = 15)
print(p6)
dev.off()

####

## 2h
p1_2h <- plotDotplot(barycoords_2h_v2, Gdiffexp = Gdiffexp_CpG_LNP_2h_vs_eLNP_2h, Goi= genelist, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v2_Early_mature_cDC1s_CpG_LNP_2h_vs_eLNP_2h_paper_2024.pdf", width = 15, height = 15)
print(p1_2h)
dev.off()

p2_2h <- plotDotplot(barycoords_2h_v2, Gdiffexp = Gdiffexp_eLNP_2h_vs_pIC_LNP_2h, Goi= genelist, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v2_Early_mature_cDC1s_eLNP_2h_vs_pIC_LNP_2h_paper_2024.pdf", width = 15, height = 15)
print(p2_2h)
dev.off()

p3_2h <- plotDotplot(barycoords_2h_v2, Gdiffexp = Gdiffexp_CpG_LNP_2h_vs_pIC_LNP_2h, Goi= genelist, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v2_Early_mature_cDC1s_CpG_LNP_2h_vs_pIC_LNP_2h_paper_2024.pdf", width = 15, height = 15)
print(p3_2h)
dev.off()

p4_2h <- plotDotplot(barycoords_2h_v3, Gdiffexp = Gdiffexp_CpG_LNP_2h_vs_pIC_LNP_2h, Goi= genelist, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/pdfPlot_v3_Early_mature_cDC1s_CpG_LNP_2h_vs_pIC_LNP_2h_paper_2024.pdf", width = 15, height = 15)
print(p4_2h)
dev.off()

p5_2h <- plotDotplot(barycoords_2h_v3, Gdiffexp = Gdiffexp_CpG_LNP_2h_vs_pIC_alone_2h, Goi= genelist, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v3_Early_mature_cDC1s_CpG_LNP_2h_vs_pIC_alone_2h_paper_2024.pdf", width = 15, height = 15)
print(p5_2h)
dev.off()

p6_2h <- plotDotplot(barycoords_2h_v3, Gdiffexp = Gdiffexp_pIC_alone_2h_vs_pIC_LNP_2h, Goi= genelist, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v3_Early_mature_cDC1s_pIC_alone_2h_vs_pIC_LNP_2h_paper_2024.pdf", width = 15, height = 15)
print(p6_2h)
dev.off()

# Interactive
## 8h
p1<-interactiveDotplot(expTable_8h_mean_v2, Gdiffexp=Gdiffexp_CpG_LNP_8h_vs_eLNP_8h, Goi= genelist, plotLocalEnrichment=FALSE, rmax=5)
print(p1)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v2_Late_mature_cDC1s_CpG_LNP_8h_vs_eLNP_8h_paper_2024.html")
saveWidget(p1,file=fileName) ##needs full path!

p2<-interactiveDotplot(expTable_8h_mean_v2, Gdiffexp=Gdiffexp_eLNP_8h_vs_pIC_LNP_8h, Goi= genelist, plotLocalEnrichment=FALSE, rmax=5)
print(p2)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v2_Late_mature_cDC1s_eLNP_8h_vs_pIC_LNP_8h_paper_2024.html")
saveWidget(p2,file=fileName) ##needs full path!

p3<-interactiveDotplot(expTable_8h_mean_v2, Gdiffexp=Gdiffexp_CpG_LNP_8h_vs_pIC_LNP_8h, Goi= genelist, plotLocalEnrichment=FALSE, rmax=5)
print(p3)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v2_Late_mature_cDC1s_CpG_LNP_8h_vs_pIC_LNP_8h_paper_2024.html")
saveWidget(p3,file=fileName) ##needs full path!

p4<-interactiveDotplot(expTable_8h_mean_v3, Gdiffexp=Gdiffexp_CpG_LNP_8h_vs_pIC_LNP_8h, Goi= genelist, plotLocalEnrichment=FALSE, rmax=5)
print(p4)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v3_Late_mature_cDC1s_CpG_LNP_8h_vs_pIC_LNP_8h_paper_2024.html")
saveWidget(p4,file=fileName) ##needs full path!

p5<-interactiveDotplot(expTable_8h_mean_v3, Gdiffexp=Gdiffexp_CpG_LNP_8h_vs_pIC_alone_8h, Goi= genelist, plotLocalEnrichment=FALSE, rmax=5)
print(p5)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v3_Late_mature_cDC1s_CpG_LNP_8h_vs_pIC_alone_8h_paper_2024.html")
saveWidget(p5,file=fileName) ##needs full path!

p6<-interactiveDotplot(expTable_8h_mean_v3, Gdiffexp=Gdiffexp_pIC_alone_8h_vs_pIC_LNP_8h, Goi= genelist, plotLocalEnrichment=FALSE, rmax=5)
print(p6)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v3_Late_mature_cDC1s_pIC_alone_8h_vs_pIC_LNP_8h_paper_2024.html")
saveWidget(p6,file=fileName) ##needs full path!


## 2h
p1_2h<-interactiveDotplot(expTable_2h_mean_v2, Gdiffexp=Gdiffexp_CpG_LNP_2h_vs_eLNP_2h, Goi= genelist, plotLocalEnrichment=FALSE, rmax=5)
print(p1_2h)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v2_Early_mature_cDC1s_CpG_LNP_2h_vs_eLNP_2h_paper_2024.html")
saveWidget(p1_2h,file=fileName) ##needs full path!

p2_2h<-interactiveDotplot(expTable_2h_mean_v2, Gdiffexp=Gdiffexp_eLNP_2h_vs_pIC_LNP_2h, Goi= genelist, plotLocalEnrichment=FALSE, rmax=5)
print(p2_2h)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v2_Early_mature_cDC1s_eLNP_2h_vs_pIC_LNP_2h_paper_2024.html")
saveWidget(p2_2h,file=fileName) ##needs full path!

p3_2h<-interactiveDotplot(expTable_2h_mean_v2, Gdiffexp=Gdiffexp_CpG_LNP_2h_vs_pIC_LNP_2h, Goi= genelist, plotLocalEnrichment=FALSE, rmax=5)
print(p3_2h)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v2_Early_mature_cDC1s_CpG_LNP_2h_vs_pIC_LNP_2h_paper_2024.html")
saveWidget(p3_2h,file=fileName) ##needs full path!

p4_2h<-interactiveDotplot(expTable_2h_mean_v3, Gdiffexp=Gdiffexp_CpG_LNP_2h_vs_pIC_LNP_2h, Goi= genelist, plotLocalEnrichment=FALSE, rmax=5)
print(p4_2h)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v3_Early_mature_cDC1s_CpG_LNP_2h_vs_pIC_LNP_2h_paper_2024.html")
saveWidget(p4_2h,file=fileName) ##needs full path!

p5_2h<-interactiveDotplot(expTable_2h_mean_v3, Gdiffexp=Gdiffexp_CpG_LNP_2h_vs_pIC_alone_2h, Goi= genelist, plotLocalEnrichment=FALSE, rmax=5)
print(p5_2h)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v3_Early_mature_cDC1s_CpG_LNP_2h_vs_pIC_alone_2h_paper_2024.html")
saveWidget(p5_2h,file=fileName) ##needs full path!

p6_2h<-interactiveDotplot(expTable_2h_mean_v3, Gdiffexp=Gdiffexp_pIC_alone_2h_vs_pIC_LNP_2h, Goi= genelist, plotLocalEnrichment=FALSE, rmax=5)
print(p6_2h)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v3_Early_mature_cDC1s_pIC_alone_2h_vs_pIC_LNP_2h_paper_2024.html")
saveWidget(p6_2h,file=fileName) ##needs full path!

###########################################################################################################
###########################################################################################################

################################################################################
######### TRIWISE PLOTS (second batch for v2/v3/v4 -> colored lists)
################################################################################

########################################
##### Preparation
########################################

wantedColors7<-c(nodiffall="gray55",diffall="black",nodiffSet1="indianred1",diffSet1="red", 
                 nodiffSet2="limegreen",diffSet2="darkgreen", nodiffSet3="lightblue", diffSet3="darkblue", 
                 nodiffSet4="darkorchid1", diffSet4="darkmagenta", nodiffSet5="lightsalmon", diffSet5="darkorange", 
                 nodiffSet6="lightskyblue", diffSet6="turquoise4", nodiffSet7="bisque", diffSet7="burlywood4")

## Update February 2024!! Switch to stricter list (clint suffix: p_adj.glb < 0.01, abs(logFC) > 1, baseMean > 50)!!!!!!!
CpG_LNP_8h_vs_eLNP_8h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_8h-eLNPs_8h_VBO_4-12_clint.rds")
eLNP_8h_vs_pIC_LNP_8h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_eLNPs_8h-pIC_LNPs_8h_VBO_4-12_clint.rds")

CpG_LNP_8h_vs_pIC_LNP_8h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_8h-pIC_LNPs_8h_VBO_4-12_clint.rds")

CpG_LNP_8h_vs_pIC_alone_8h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_8h-pIC_8h_VBO_4-12_clint.rds")
pIC_alone_8h_vs_pIC_LNP_8h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_pIC_8h-pIC_LNPs_8h_VBO_4-12_clint.rds")

eLNP_8h_vs_pIC_alone_8h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_eLNPs_8h-pIC_8h_VBO_4-12_clint.rds")


Gdiffexp_CpG_LNP_8h_vs_eLNP_8h <- c(CpG_LNP_8h_vs_eLNP_8h$`Late mature cDC1s`$gene)
Gdiffexp_eLNP_8h_vs_pIC_LNP_8h <- c(eLNP_8h_vs_pIC_LNP_8h$`Late mature cDC1s`$gene)

Gdiffexp_CpG_LNP_8h_vs_pIC_LNP_8h <- c(CpG_LNP_8h_vs_pIC_LNP_8h$`Late mature cDC1s`$gene)

Gdiffexp_CpG_LNP_8h_vs_pIC_alone_8h <- c(CpG_LNP_8h_vs_pIC_alone_8h$`Late mature cDC1s`$gene)
Gdiffexp_pIC_alone_8h_vs_pIC_LNP_8h <- c(pIC_alone_8h_vs_pIC_LNP_8h$`Late mature cDC1s`$gene)

Gdiffexp_eLNP_8h_vs_pIC_alone_8h <- c(eLNP_8h_vs_pIC_alone_8h$`Late mature cDC1s`$gene)

## Update February 2024!! Switch to stricter list (clint suffix: p_adj.glb < 0.01, abs(logFC) > 1, baseMean > 50)!!!!!!!
CpG_LNP_2h_vs_eLNP_2h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_2h-eLNPs_2h_VBO_4-12_clint.rds")
eLNP_2h_vs_pIC_LNP_2h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_eLNPs_2h-pIC_LNPs_2h_VBO_4-12_clint.rds")

CpG_LNP_2h_vs_pIC_LNP_2h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_2h-pIC_LNPs_2h_VBO_4-12_clint.rds")

CpG_LNP_2h_vs_pIC_alone_2h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_2h-pIC_2h_VBO_4-12_clint.rds")
pIC_alone_2h_vs_pIC_LNP_2h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_pIC_2h-pIC_LNPs_2h_VBO_4-12_clint.rds")

eLNP_2h_vs_pIC_alone_2h<-readRDS("VBO_LNP_merge/results/Robjects/tbl_fil_muscat_DESeq2_new_eLNPs_2h-pIC_2h_VBO_4-12_clint.rds")

Gdiffexp_CpG_LNP_2h_vs_eLNP_2h <- c(CpG_LNP_2h_vs_eLNP_2h$`Early mature cDC1s`$gene)
Gdiffexp_eLNP_2h_vs_pIC_LNP_2h <- c(eLNP_2h_vs_pIC_LNP_2h$`Early mature cDC1s`$gene)

Gdiffexp_CpG_LNP_2h_vs_pIC_LNP_2h <- c(CpG_LNP_2h_vs_pIC_LNP_2h$`Early mature cDC1s`$gene)

Gdiffexp_CpG_LNP_2h_vs_pIC_alone_2h <- c(CpG_LNP_2h_vs_pIC_alone_2h$`Early mature cDC1s`$gene)
Gdiffexp_pIC_alone_2h_vs_pIC_LNP_2h <- c(pIC_alone_2h_vs_pIC_LNP_2h$`Early mature cDC1s`$gene)

Gdiffexp_eLNP_2h_vs_pIC_alone_2h <- c(eLNP_2h_vs_pIC_alone_2h$`Early mature cDC1s`$gene)


allDEgenes_8h_v2<-unique(c(Gdiffexp_CpG_LNP_8h_vs_eLNP_8h,Gdiffexp_eLNP_8h_vs_pIC_LNP_8h,Gdiffexp_CpG_LNP_8h_vs_pIC_LNP_8h))
length(allDEgenes_8h_v2)

allDEgenes_8h_v3<-unique(c(Gdiffexp_CpG_LNP_8h_vs_pIC_alone_8h,Gdiffexp_pIC_alone_8h_vs_pIC_LNP_8h,Gdiffexp_CpG_LNP_8h_vs_pIC_LNP_8h))
length(allDEgenes_8h_v3)

allDEgenes_8h_v4<-unique(c(Gdiffexp_eLNP_8h_vs_pIC_alone_8h,Gdiffexp_pIC_alone_8h_vs_pIC_LNP_8h,Gdiffexp_eLNP_8h_vs_pIC_LNP_8h))
length(allDEgenes_8h_v4)


allDEgenes_2h_v2<-unique(c(Gdiffexp_CpG_LNP_2h_vs_eLNP_2h,Gdiffexp_eLNP_2h_vs_pIC_LNP_2h,Gdiffexp_CpG_LNP_2h_vs_pIC_LNP_2h))
length(allDEgenes_2h_v2)

allDEgenes_2h_v3<-unique(c(Gdiffexp_CpG_LNP_2h_vs_pIC_alone_2h,Gdiffexp_pIC_alone_2h_vs_pIC_LNP_2h,Gdiffexp_CpG_LNP_2h_vs_pIC_LNP_2h))
length(allDEgenes_2h_v3)

allDEgenes_2h_v4<-unique(c(Gdiffexp_eLNP_2h_vs_pIC_alone_2h,Gdiffexp_pIC_alone_2h_vs_pIC_LNP_2h,Gdiffexp_eLNP_2h_vs_pIC_LNP_2h))
length(allDEgenes_2h_v4)

###Triwise plot combined DEGs
p_v2_8h <- plotDotplot(barycoords_8h_v2, Gdiffexp = allDEgenes_8h_v2, showlabels = F, sizevalues = 50, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v2_8h_Late_mature_cDC1s_all_DEGs_paper_2024.pdf", width = 15, height = 15)
print(p_v2_8h)
dev.off()

p_v2_8h <- plotDotplot(barycoords_8h_v2, Gdiffexp = allDEgenes_8h_v2, showlabels = F, sizevalues = 50, rmax = 8) #higher rmax 04/2024
pdf("VBO_LNP_merge/results/Triwise/Plot_v2_8h_Late_mature_cDC1s_all_DEGs_paper_2024_high_rmax.pdf", width = 15, height = 15)
print(p_v2_8h)
dev.off()

p_v3_8h <- plotDotplot(barycoords_8h_v3, Gdiffexp = allDEgenes_8h_v3, showlabels = F, sizevalues = 50, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v3_8h_Late_mature_cDC1s_all_DEGs_paper_2024.pdf", width = 15, height = 15)
print(p_v3_8h)
dev.off()

p_v4_8h <- plotDotplot(barycoords_8h_v4, Gdiffexp = allDEgenes_8h_v4, showlabels = F, sizevalues = 50, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v4_8h_Late_mature_cDC1s_all_DEGs_paper_2024.pdf", width = 15, height = 15)
print(p_v4_8h)
dev.off()

p_v2_2h <- plotDotplot(barycoords_2h_v2, Gdiffexp = allDEgenes_2h_v2, showlabels = F, sizevalues = 50, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v2_2h_Early_mature_cDC1s_all_DEGs_paper_2024.pdf", width = 15, height = 15)
print(p_v2_2h)
dev.off()

p_v2_2h <- plotDotplot(barycoords_2h_v2, Gdiffexp = allDEgenes_2h_v2, showlabels = F, sizevalues = 50, rmax = 8)
pdf("VBO_LNP_merge/results/Triwise/Plot_v2_2h_Early_mature_cDC1s_all_DEGs_paper_2024_high_rmax.pdf", width = 15, height = 15)
print(p_v2_2h)
dev.off()

p_v3_2h <- plotDotplot(barycoords_2h_v3, Gdiffexp = allDEgenes_2h_v3, showlabels = F, sizevalues = 50, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v3_2h_Early_mature_cDC1s_all_DEGs_paper_2024.pdf", width = 15, height = 15)
print(p_v3_2h)
dev.off()

p_v4_2h <- plotDotplot(barycoords_2h_v4, Gdiffexp = allDEgenes_2h_v4, showlabels = F, sizevalues = 50, rmax = 5)
pdf("VBO_LNP_merge/results/Triwise/Plot_v4_2h_Early_mature_cDC1s_all_DEGs_paper_2024.pdf", width = 15, height = 15)
print(p_v4_2h)
dev.off()

###Rose plot
p_v2_8h <- plotRoseplot(barycoords_8h_v2, Gdiffexp = allDEgenes_8h_v2, showlabels = F)
pdf("VBO_LNP_merge/results/Triwise/Roseplot_v2_8h_Late_mature_cDC1s_all_DEGs_paper_2024.pdf", width = 15, height = 15)
print(p_v2_8h)
dev.off()

p_v3_8h <- plotRoseplot(barycoords_8h_v3, Gdiffexp = allDEgenes_8h_v3, showlabels = F)
pdf("VBO_LNP_merge/results/Triwise/Roseplot_v3_8h_Late_mature_cDC1s_all_DEGs_paper_2024.pdf", width = 15, height = 15)
print(p_v3_8h)
dev.off()

p_v4_8h <- plotRoseplot(barycoords_8h_v4, Gdiffexp = allDEgenes_8h_v4, showlabels = F)
pdf("VBO_LNP_merge/results/Triwise/Roseplot_v4_8h_Late_mature_cDC1s_all_DEGs_paper_2024.pdf", width = 15, height = 15)
print(p_v4_8h)
dev.off()

p_v2_2h <- plotRoseplot(barycoords_2h_v2, Gdiffexp = allDEgenes_2h_v2, showlabels = F)
pdf("VBO_LNP_merge/results/Triwise/Roseplot_v2_2h_Early_mature_cDC1s_all_DEGs_paper_2024.pdf", width = 15, height = 15)
print(p_v2_2h)
dev.off()

p_v3_2h <- plotRoseplot(barycoords_2h_v3, Gdiffexp = allDEgenes_2h_v3, showlabels = F)
pdf("VBO_LNP_merge/results/Triwise/Roseplot_v3_2h_Early_mature_cDC1s_all_DEGs_paper_2024.pdf", width = 15, height = 15)
print(p_v3_2h)
dev.off()

p_v4_2h <- plotRoseplot(barycoords_2h_v4, Gdiffexp = allDEgenes_2h_v4, showlabels = F)
pdf("VBO_LNP_merge/results/Triwise/Roseplot_v4_2h_Early_mature_cDC1s_all_DEGs_paper_2024.pdf", width = 15, height = 15)
print(p_v4_2h)
dev.off()

###Interactive plot
p_v2_8h<-interactiveDotplot(expTable_8h_mean_v2, Gdiffexp=allDEgenes_8h_v2, plotLocalEnrichment=FALSE, rmax=5)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v2_8h_Late_mature_cDC1s_all_DEGs_paper_2024.html")
saveWidget(p_v2_8h,file=fileName) ##needs full path!

p_v3_8h<-interactiveDotplot(expTable_8h_mean_v3, Gdiffexp=allDEgenes_8h_v3, plotLocalEnrichment=FALSE, rmax=5)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v3_8h_Late_mature_cDC1s_all_DEGs_paper_2024.html")
saveWidget(p_v3_8h,file=fileName) ##needs full path!

p_v4_8h<-interactiveDotplot(expTable_8h_mean_v4, Gdiffexp=allDEgenes_8h_v4, plotLocalEnrichment=FALSE, rmax=5)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v4_8h_Late_mature_cDC1s_all_DEGs_paper_2024.html")
saveWidget(p_v4_8h,file=fileName) ##needs full path!


p_v2_2h<-interactiveDotplot(expTable_2h_mean_v2, Gdiffexp=allDEgenes_2h_v2, plotLocalEnrichment=FALSE, rmax=5)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v2_2h_Early_mature_cDC1s_all_DEGs_paper_2024.html")
saveWidget(p_v2_2h,file=fileName) ##needs full path!

p_v3_2h<-interactiveDotplot(expTable_2h_mean_v3, Gdiffexp=allDEgenes_2h_v3, plotLocalEnrichment=FALSE, rmax=5)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v3_2h_Early_mature_cDC1s_all_DEGs_paper_2024.html")
saveWidget(p_v3_2h,file=fileName) ##needs full path!

p_v4_2h<-interactiveDotplot(expTable_2h_mean_v4, Gdiffexp=allDEgenes_2h_v4, plotLocalEnrichment=FALSE, rmax=5)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/interactivePlot_v4_2h_Early_mature_cDC1s_all_DEGs_paper_2024.html")
saveWidget(p_v4_2h,file=fileName) ##needs full path!


###Extra: create colored plots with gene lists! (only eBayes)
###Seperate the DE genes of the different comparisons
#set1 (red)=DE genes only present in CpG_LNP_8h_vs_eLNP_8h
#set2 (green)=DE genes only present in eLNP_8h_vs_pIC_LNP_8h
#set3 (blue)=DE genes only present in CpG_LNP_8h_vs_pIC_LNP_8h
#set4 (purple)=DE genes present in all 3 comparisons
#set5 (orange)=DE genes present in CpG_LNP_8h_vs_eLNP_8h and eLNP_8h_vs_pIC_LNP_8h
#set6 (turquoise)=DE genes present in eLNP_8h_vs_pIC_LNP_8h and CpG_LNP_8h_vs_pIC_LNP_8h
#set7 (brown)=DE genes present in CpG_LNP_8h_vs_eLNP_8h and CpG_LNP_8h_vs_pIC_LNP_8h

DEgenesGroup1 <- CpG_LNP_8h_vs_eLNP_8h$`Late mature cDC1s`
rownames(DEgenesGroup1) <- DEgenesGroup1$gene
DEgenesGroup2 <- eLNP_8h_vs_pIC_LNP_8h$`Late mature cDC1s`
rownames(DEgenesGroup2) <- DEgenesGroup2$gene
DEgenesGroup3 <- CpG_LNP_8h_vs_pIC_LNP_8h$`Late mature cDC1s`
rownames(DEgenesGroup3) <- DEgenesGroup3$gene

genesSet5tmp<-intersect(rownames(DEgenesGroup1),rownames(DEgenesGroup2))
genesSet6tmp<-intersect(rownames(DEgenesGroup2),rownames(DEgenesGroup3))
genesSet7tmp<-intersect(rownames(DEgenesGroup1),rownames(DEgenesGroup3))
genesSet4<-intersect(genesSet5tmp, rownames(DEgenesGroup3))
genesSet5<-setdiff(genesSet5tmp,genesSet4)
genesSet6<-setdiff(genesSet6tmp,genesSet4)
genesSet7<-setdiff(genesSet7tmp,genesSet4)
genesSet1<-setdiff(rownames(DEgenesGroup1),unique(c(rownames(DEgenesGroup2),rownames(DEgenesGroup3))))
genesSet2<-setdiff(rownames(DEgenesGroup2),unique(c(rownames(DEgenesGroup1),rownames(DEgenesGroup3))))
genesSet3<-setdiff(rownames(DEgenesGroup3),unique(c(rownames(DEgenesGroup1),rownames(DEgenesGroup2))))

# 8 changed to 9: p adj loc -> p adj glb!!
genesSet1matrix<-DEgenesGroup1[genesSet1,c(1,3,4,9)]
genesSet2matrix<-DEgenesGroup2[genesSet2,c(1,3,4,9)]
genesSet3matrix<-DEgenesGroup3[genesSet3,c(1,3,4,9)]
genesSet4matrix<-cbind(DEgenesGroup1[genesSet4,c(1,3,4,9)],DEgenesGroup2[genesSet4,c(1,3,4,9)],DEgenesGroup3[genesSet4,c(1,3,4,9)])
genesSet5matrix<-cbind(DEgenesGroup1[genesSet5,c(1,3,4,9)],DEgenesGroup2[genesSet5,c(1,3,4,9)])
genesSet6matrix<-cbind(DEgenesGroup2[genesSet6,c(1,3,4,9)],DEgenesGroup3[genesSet6,c(1,3,4,9)])
genesSet7matrix<-cbind(DEgenesGroup1[genesSet7,c(1,3,4,9)],DEgenesGroup3[genesSet7,c(1,3,4,9)]) 

geneset_list<-tibble::lst(genesSet1matrix,genesSet2matrix,genesSet3matrix,genesSet4matrix,
                          genesSet5matrix,genesSet6matrix,genesSet7matrix)
names(geneset_list)<-c("red","green","blue", "purple","orange","turquoise","brown")

### Write results
write.xlsx(geneset_list, "VBO_LNP_merge/results/Triwise/triwisePlot_withColors_8h_v2_paper_2024.xlsx")

p_v2_8h<-plotDotplot(barycoords_8h_v2, Gdiffexp=allDEgenes_8h_v2, showlabels = F, rmax = 8, Goi=list(Set1=genesSet1,Set2=genesSet2,Set3=genesSet3,Set4=genesSet4,
                                                                                                     Set5=genesSet5,Set6=genesSet6,Set7=genesSet7), colorvalues=wantedColors7) +
  theme(legend.position = "none")
pdf("VBO_LNP_merge/results/Triwise/Plot_v2_8h_Late_mature_cDC1s_all_DEGs_colored_paper_2024_high_rmax.pdf", width = 15, height = 15)
print(p_v2_8h)
dev.off()

myList<-list("set1"=genesSet1,"set2"=genesSet2,"set3"=genesSet3,"set4"=genesSet4,"set5"=genesSet5,"set6"=genesSet6,"set7"=genesSet7)
saveRDS(myList,file="VBO_LNP_merge/results/Triwise/listGenes_coloredTriwise_v2_8h_paper_2024.rds")

## Extra paper Feb 2024: annotate certain genes on this colored plot
p_v2_8h <- p_v2_8h + scale_y_continuous(expand = expansion(mult = 0.2)) + scale_x_continuous(expand = expansion(mult = 0.2))   #change scaling

#MARK GENES ON PLOT
C103_list<-read.xlsx("/home/clintdn/VIB/DATA/Sophie/RNA-seq_Victor/2020_analysis/Immunity 2016 Ardouin-2_Sophie.xlsx",
                     sheet = "C103")
genelist<-firstup(C103_list$Gene.Symbol)

# FROM FIND MAKERS
NegXBarycoords <- subset(barycoords_8h_v2[intersect(unlist(genelist),allDEgenes_8h_v2),],barycoords_8h_v2[intersect(unlist(genelist),allDEgenes_8h_v2),1]<0)
PosXBarycoords <- subset(barycoords_8h_v2[intersect(unlist(genelist),allDEgenes_8h_v2),],barycoords_8h_v2[intersect(unlist(genelist),allDEgenes_8h_v2),1]>0)
p_v2_8h <- p_v2_8h +
  geom_text_repel(NegXBarycoords[],
                  mapping=aes(label=rownames(NegXBarycoords), x=NegXBarycoords[,1], y=NegXBarycoords[,2]),
                  force = 1, size=2, segment.alpha=0.5, segment.colour = '#999999',
                  nudge_x = (0-NegXBarycoords[,1]-7),
                  box.padding =0.25) +
  geom_text_repel(PosXBarycoords[,],
                  mapping=aes(label=rownames(PosXBarycoords), x=PosXBarycoords[,1], y=PosXBarycoords[,2]),
                  force = 1, size=2, segment.alpha=0.5, segment.colour = '#999999',
                  nudge_x = (0-PosXBarycoords[,1]+7),
                  box.padding =0.25)

p_v2_8h

######

###Extra: create colored plots with gene lists! (only eBayes)
###Seperate the DE genes of the different comparisons
#set1 (red)=DE genes only present in CpG_LNP_2h_vs_eLNP_2h
#set2 (green)=DE genes only present in eLNP_2h_vs_pIC_LNP_2h
#set3 (blue)=DE genes only present in CpG_LNP_2h_vs_pIC_LNP_2h
#set4 (purple)=DE genes present in all 3 comparisons
#set5 (orange)=DE genes present in CpG_LNP_2h_vs_eLNP_2h and eLNP_2h_vs_pIC_LNP_2h
#set6 (turquoise)=DE genes present in eLNP_2h_vs_pIC_LNP_2h and CpG_LNP_2h_vs_pIC_LNP_2h
#set7 (brown)=DE genes present in CpG_LNP_2h_vs_eLNP_2h and CpG_LNP_2h_vs_pIC_LNP_2h

DEgenesGroup1 <- CpG_LNP_2h_vs_eLNP_2h$`Early mature cDC1s`
rownames(DEgenesGroup1) <- DEgenesGroup1$gene
DEgenesGroup2 <- eLNP_2h_vs_pIC_LNP_2h$`Early mature cDC1s`
rownames(DEgenesGroup2) <- DEgenesGroup2$gene
DEgenesGroup3 <- CpG_LNP_2h_vs_pIC_LNP_2h$`Early mature cDC1s`
rownames(DEgenesGroup3) <- DEgenesGroup3$gene

genesSet5tmp<-intersect(rownames(DEgenesGroup1),rownames(DEgenesGroup2))
genesSet6tmp<-intersect(rownames(DEgenesGroup2),rownames(DEgenesGroup3))
genesSet7tmp<-intersect(rownames(DEgenesGroup1),rownames(DEgenesGroup3))
genesSet4<-intersect(genesSet5tmp, rownames(DEgenesGroup3))
genesSet5<-setdiff(genesSet5tmp,genesSet4)
genesSet6<-setdiff(genesSet6tmp,genesSet4)
genesSet7<-setdiff(genesSet7tmp,genesSet4)
genesSet1<-setdiff(rownames(DEgenesGroup1),unique(c(rownames(DEgenesGroup2),rownames(DEgenesGroup3))))
genesSet2<-setdiff(rownames(DEgenesGroup2),unique(c(rownames(DEgenesGroup1),rownames(DEgenesGroup3))))
genesSet3<-setdiff(rownames(DEgenesGroup3),unique(c(rownames(DEgenesGroup1),rownames(DEgenesGroup2))))

genesSet1matrix<-DEgenesGroup1[genesSet1,c(1,3,4,9)]
genesSet2matrix<-DEgenesGroup2[genesSet2,c(1,3,4,9)]
genesSet3matrix<-DEgenesGroup3[genesSet3,c(1,3,4,9)]
genesSet4matrix<-cbind(DEgenesGroup1[genesSet4,c(1,3,4,9)],DEgenesGroup2[genesSet4,c(1,3,4,9)],DEgenesGroup3[genesSet4,c(1,3,4,9)])
genesSet5matrix<-cbind(DEgenesGroup1[genesSet5,c(1,3,4,9)],DEgenesGroup2[genesSet5,c(1,3,4,9)])
genesSet6matrix<-cbind(DEgenesGroup2[genesSet6,c(1,3,4,9)],DEgenesGroup3[genesSet6,c(1,3,4,9)])
genesSet7matrix<-cbind(DEgenesGroup1[genesSet7,c(1,3,4,9)],DEgenesGroup3[genesSet7,c(1,3,4,9)]) 

geneset_list<-tibble::lst(genesSet1matrix,genesSet2matrix,genesSet3matrix,genesSet4matrix,
                          genesSet5matrix,genesSet6matrix,genesSet7matrix)
names(geneset_list)<-c("red","green","blue", "purple","orange","turquoise","brown")

### Write results
write.xlsx(geneset_list, "VBO_LNP_merge/results/Triwise/triwisePlot_withColors_2h_v2_paper_2024.xlsx")

p_v2_2h<-plotDotplot(barycoords_2h_v2, Gdiffexp=allDEgenes_2h_v2, showlabels = F, rmax = 8, Goi=list(Set1=genesSet1,Set2=genesSet2,Set3=genesSet3,Set4=genesSet4,
                                                                                                     Set5=genesSet5,Set6=genesSet6,Set7=genesSet7), colorvalues=wantedColors7) +
  theme(legend.position = "none")
pdf("VBO_LNP_merge/results/Triwise/Plot_v2_2h_Early_mature_cDC1s_all_DEGs_colored_paper_2024_high_rmax.pdf", width = 15, height = 15)
print(p_v2_2h)
dev.off()

myList<-list("set1"=genesSet1,"set2"=genesSet2,"set3"=genesSet3,"set4"=genesSet4,"set5"=genesSet5,"set6"=genesSet6,"set7"=genesSet7)
saveRDS(myList,file="VBO_LNP_merge/results/Triwise/listGenes_coloredTriwise_v2_2h_paper_2024.rds")

## Extra paper Feb 2024: annotate certain genes on this colored plot

##################################################

## Extra: plot genelists meta analysis on the triwise of second batch
Genelist_tibble<-readRDS("../../Meta_analysis/Tibble_all_genelists_heatmap_v7_with_top_200_LNP_lists.rds")

#OWN GENES OF INTEREST
for (i in 1:length(Genelist_tibble)){
  genelist<-Genelist_tibble[[i]]
  name<-names(Genelist_tibble[i])
  
  #Triwise
  p_v2_2h <- plotDotplot(barycoords_2h_v2, Gdiffexp = allDEgenes_2h_v2, Goi= genelist, showlabels = F, sizevalues = 50, rmax = 5)
  pdf(paste0("VBO_LNP_merge/results/Triwise/Plot_v2_2h_Early_mature_cDC1s_",name,"_paper_2024.pdf"), width = 15, height = 15)
  print(p_v2_2h)
  dev.off()
  
  p_v2_8h <- plotDotplot(barycoords_8h_v2, Gdiffexp = allDEgenes_8h_v2, Goi= genelist, showlabels = F, sizevalues = 50, rmax = 5)
  pdf(paste0("VBO_LNP_merge/results/Triwise/Plot_v2_8h_Late_mature_cDC1s_",name,"_paper_2024.pdf"), width = 15, height = 15)
  print(p_v2_8h)
  dev.off()
  
  #Roseplot
  p2_v2_2h <- plotRoseplot(barycoords_2h_v2, Gdiffexp = allDEgenes_2h_v2, Goi= genelist, showlabels = F)
  pdf(paste0("VBO_LNP_merge/results/Triwise/Roseplot_v2_2h_Early_mature_cDC1s_",name,"_paper_2024.pdf"), width = 15, height = 15)
  print(p2_v2_2h)
  dev.off()
  
  p2_v2_8h <- plotRoseplot(barycoords_8h_v2, Gdiffexp = allDEgenes_8h_v2, Goi = genelist, showlabels = F)
  pdf(paste0("VBO_LNP_merge/results/Triwise/Roseplot_v2_8h_Late_mature_cDC1s_",name,"_paper_2024.pdf"), width = 15, height = 15)
  print(p2_v2_8h)
  dev.off()
}
