## Script for analyzing LNP cDC1 CITE-seq project data
## Muscat script for comparing between conditions part 1 ADT
## LM one condition vs LM other condition
## EM one condition vs EM other condition

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

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd('/home/clintdn/VIB/DATA/Sophie/scRNA-seq_Victor/LNP_CITEseq_experiment/')

sampleName<-"VBO_4-12"
sampleFolder<-"VBO_LNP_merge/"

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

## Read objects
seuratObj <- readRDS(file = paste0(sampleFolder,"results/Robjects/seuratObj_",sampleName,".rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"results/Robjects/diagnostics_",sampleName,".rds"))

########################################################################################################################

## Simplify annotation
seuratObj$annotated_clusters_Muscat<-seuratObj$annotated_clusters
levels(seuratObj$annotated_clusters_Muscat)<-c('Early mature cDC1s',"Late mature cDC1s","Early mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Late mature cDC1s",
                                               "Proliferating cDC1s","Proliferating cDC1s","Late mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Other cDC1s","Late mature cDC1s",
                                               "Early mature cDC1s","Late mature cDC1s","Doublets","Late mature cDC1s","Early mature cDC1s")

DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat")

## Combine MULTI_ID with Orig.ident or Genotype?? Need to put genotype in correct order: otherwise WT after DKO
seuratObj@meta.data[["MULTI_ID_merge"]]
seuratObj@meta.data$MULTI_ID_merge<-as.factor(seuratObj@meta.data$MULTI_ID_merge)

########################################
##### Load data (already filtered) 
########################################

## Perform on ADT slot!!!
seuratObj[['ADT']]@counts <- as.matrix(seuratObj[['ADT']]@counts)
seuratObj[['ADT']]@data <- as.matrix(seuratObj[['ADT']]@data) 

### Convert to sce ### ADT here
sce <- as.SingleCellExperiment(seuratObj, assay = "ADT")
########################################
##### Prepare sce object
########################################

### Create test metadata table
sample_id<-levels(as.factor(seuratObj@meta.data$MULTI_ID_merge))
group_id<-c(rep("WT",4),rep("eLNPs_8h",4),rep("pIC_LNPs_8h",4),rep("CpG_LNPs_8h",3),rep("pIC_8h",4),
            rep("eLNPs_2h",4),rep("pIC_LNPs_2h",4),rep("CpG_LNPs_2h",4),rep("pIC_2h",4)) 
n_cells<-rep(1,35)

### Create metaData matrix
metaData<-data.frame(sample_id,group_id,n_cells)

### Look through samples (1->8) for amount of cells after filtering
for(i in 1:length(sample_id)){
  toSearch<-sample_id[i]
  metaData[i, "n_cells"]<-length(grep(paste0("\\<",toSearch,"\\>"),seuratObj@meta.data$MULTI_ID_merge)) ##Need to do specific grep. Otherwise fault between Hashtag 1 and Hashtag 10!!!
}

### Retrieve the sizes of the two groups
WT_size<-sum(metaData$n_cells[1:4]) #14539
eLNP8h_size<-sum(metaData$n_cells[5:8]) #12207
pICLNP8h_size<-sum(metaData$n_cells[9:12]) #13483
CpGLNP8h_size<-sum(metaData$n_cells[13:15]) #12840
pIC8h_size<-sum(metaData$n_cells[16:19]) #13414
eLNP2h_size<-sum(metaData$n_cells[20:23]) #9187
pICLNP2h_size<-sum(metaData$n_cells[24:27]) #14315
CpGLNP2h_size<-sum(metaData$n_cells[28:31]) #12766
pIC2h_size<-sum(metaData$n_cells[32:35]) #13837

### Add group_id to sce object
sce$group_id<-c(rep("WT",WT_size),rep("eLNPs_8h",eLNP8h_size),rep("pIC_LNPs_8h",pICLNP8h_size),
                rep("CpG_LNPs_8h",CpGLNP8h_size),rep("pIC_8h",pIC8h_size),
                rep("eLNPs_2h",eLNP2h_size),rep("pIC_LNPs_2h",pICLNP2h_size),
                rep("CpG_LNPs_2h",CpGLNP2h_size),rep("pIC_2h",pIC2h_size))

### Add new sample_id to sce object (combination of orig.ident and group_id)
sce$sample_id<-sce$MULTI_ID_merge

### Finish prep object
(sce <- prepSCE(sce, 
                kid = "annotated_clusters_Muscat", # cell population assignments
                gid = "group_id",   # group IDs (resp/Non_Resp)
                sid = "sample_id",    # sample IDs (group+orig.ident)
                drop = TRUE))        # drop all other colData columns

### Store ids separately for easy access
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

############################
##### Visualization
############################

### Calculate number of cells per cluster-sample (can already exclude here or later in aggregate step)
t(table(sce$cluster_id, sce$sample_id))

##############################################################################
##### Differential state analysis (only on sample level for this experiment)
##############################################################################
# Aggregation-based method that acts on *pseudobulk* data. Each gene is tested for state changes in each cluster.
# Thus, a total of #(genes)\times\#(clusters) tests will be performed per comparison of interest.

# First aggregate the measurements for each sample (in each cluster) to obtain pseudobulk data.
# In general, `aggregateData()` will aggregate the data by the `colData` variables specified with argument `by`, 
# and return a `SingleCellExperiment` containing pseudobulk data.  

# For DS analysis, measurements must be aggregated at the cluster-sample level (default `by = c("cluster_id", "sample_id"`). 
# In this case, the returned `SingleCellExperiment` will contain one assay per cluster, where rows = genes and columns = samples. 
# Arguments `assay` and `fun` specify the input data and summary statistic, respectively, to use for aggregation.  
# Default it's applied to the sum of the raw counts, but there are also other possiblities (input data (raw/(log-)normalized counts, 
# CPM ect.) and summary statistics (sum, mean, median))

### Aggregation
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)
# pseudobulks for 1st subpopulation
t(head(assay(pb)))

### Pseudobulk-level MDS plot
# Multi-dimensional scaling (MDS) plot of aggregated signal to explore overall sample similarities.
# Ideally, such a representation of the data should separate both clusters and groups from one another. -> Not the case!!
(pb_mds <- pbMDS(pb) + scale_shape_manual(values=1:nlevels(pb$group_id))) #Add function for more than 6 shapes (9 groups)

# Removing cluster-sample instance(s) ‘Other cDC1s’-‘VBO004_Hashtag1’, ‘Other cDC1s’-‘VBO004_Hashtag2’, ‘Other cDC1s’-‘VBO004_Hashtag3’, ‘Other cDC1s’-‘VBO004_Hashtag4’, 
# ‘Other cDC1s’-‘VBO005_Hashtag1’, ‘Other cDC1s’-‘VBO005_Hashtag2’, ‘Other cDC1s’-‘VBO005_Hashtag3’, ‘Other cDC1s’-‘VBO005_Hashtag4’, ‘Other cDC1s’-‘VBO008_Hashtag4’, 
# ‘Other cDC1s’-‘VBO009_Hashtag1’, ‘Other cDC1s’-‘VBO009_Hashtag2’, ‘Other cDC1s’-‘VBO009_Hashtag3’, ‘Other cDC1s’-‘VBO010_Hashtag1’, ‘Other cDC1s’-‘VBO010_Hashtag2’, 
# ‘Other cDC1s’-‘VBO010_Hashtag3’, ‘Other cDC1s’-‘VBO010_Hashtag4’, ‘Other cDC1s’-‘VBO011_Hashtag1’, ‘Other cDC1s’-‘VBO011_Hashtag2’, ‘Other cDC1s’-‘VBO011_Hashtag3’, 
# ‘Other cDC1s’-‘VBO011_Hashtag4’, ‘Other cDC1s’-‘VBO012_Hashtag1’, ‘Other cDC1s’-‘VBO012_Hashtag2’, ‘Other cDC1s’-‘VBO012_Hashtag3’, ‘Other cDC1s’-‘VBO012_Hashtag4’, 
# ‘Doublets’-‘VBO006_Hashtag1’, ‘Doublets’-‘VBO006_Hashtag2’, ‘Doublets’-‘VBO006_Hashtag3’, ‘Doublets’-‘VBO006_Hashtag4’, ‘Doublets’-‘VBO007_Hashtag1’, 
# ‘Doublets’-‘VBO009_Hashtag1’, ‘Doublets’-‘VBO009_Hashtag2’, ‘Doublets’-‘VBO009_Hashtag3’, ‘Doublets’-‘VBO009_Hashtag4’, ‘Doublets’-‘VBO010_Hashtag1’, 
# ‘Doublets’-‘VBO010_Hashtag2’, ‘Doublets’-‘VBO010_Hashtag3’, ‘Doublets’-‘VBO010_Hashtag4’, ‘Doublets’-‘VBO011_Hashtag1’, ‘Doublets’-‘VBO011_Hashtag3’, 
# ‘Doublets’-‘VBO012_Hashtag2’

dir.create(paste0(sampleFolder,"results/Muscat_ADT/"))
dir.create(paste0(sampleFolder,"results/Robjects_ADT_muscat/"))

Comparison<-"ADT_DS_Simple_RNA_Harmony_clustering"

png(file=paste0(sampleFolder,"results/Muscat_ADT/meanVariancePlot_subset_",Comparison,".png"),width = 3500, height = 2000, res = 300)
pb_mds + scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set2"))
dev.off()

### Experimental design
# We can provide `pbDS` with a design matrix capturing the experimental design using `model.matrix` (package `r Rpackage("stats")`), 
# and a contrast matrix that specifies our comparison of interesting using `makeContrasts` from the `r Biocpkg("limma")` package. 
# Alternatively, the comparison(s) of interest (or a list thereof) can be specified with via `coefs` (see `?glmQLFTest` for details). 

# Here, we want to carry out a single comparison of DKO against WT samples, 
# thus placing `"WT"` on the right-hand side as the reference condition.
library("limma")
# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id) #model.matrix(~ 0 + ei$group_id + ei$gender) Include covariates here!!
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("CpG_LNPs_2h-WT", "CpG_LNPs_8h-WT", "eLNPs_2h-WT", "eLNPs_8h-WT",
                            "pIC_2h-WT", "pIC_8h-WT", "pIC_LNPs_2h-WT", "pIC_LNPs_8h-WT",
                            "CpG_LNPs_2h-eLNPs_2h", "CpG_LNPs_2h-pIC_2h", "CpG_LNPs_2h-pIC_LNPs_2h", 
                            "eLNPs_2h-pIC_2h","eLNPs_2h-pIC_LNPs_2h", "pIC_2h-pIC_LNPs_2h", 
                            "CpG_LNPs_8h-eLNPs_8h", "CpG_LNPs_8h-pIC_8h", "CpG_LNPs_8h-pIC_LNPs_8h", 
                            "eLNPs_8h-pIC_8h","eLNPs_8h-pIC_LNPs_8h", "pIC_8h-pIC_LNPs_8h",
                          levels = mm)

# run DS analysis
### DeSeq2 method ###
res_DESeq2<-pbDS(pb, design = mm, contrast = contrast,method = "DESeq2")

## Save res object
saveRDS(res_DESeq2,file=paste0(sampleFolder,"results/Robjects_ADT_muscat/res_Muscat_DeSeq2_new_",Comparison,"_",sampleName,".rds"))

## Choose Res object!
res<-readRDS(file=paste0(sampleFolder,"results/Robjects_ADT_muscat/res_Muscat_DeSeq2_new_",Comparison,"_",sampleName,".rds")) ##DESeq2

# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)
# view results for 1st cluster
k1 <- tbl[[1]]
head(format(k1[, -ncol(k1)], digits = 2))

## Create loop for all comparisons
for (i in 1:length(res$table)){
  tbl <- res$table[[i]]
  Comp<-names(res$table)[i]
  write.xlsx(tbl, paste0(sampleFolder,"results/Muscat_ADT/tbl_full_muscat_DESeq2_new_",Comp,"_",sampleName,".xlsx"))
  
  # filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
  tbl_fil <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
    dplyr::arrange(u, p_adj.loc)
  })
  
  # filter clint -> here more strict because too many DE genes!!!!!!!!!!!!!
  tbl_fil_clint <- lapply(tbl, function(u) { ## EMDS stronger filtering for heatmaps!! (12/10/2023)
    u <- dplyr::filter(u, p_adj.glb < 0.01, abs(logFC) > 1, baseMean > 50)
    dplyr::arrange(u, p_adj.glb)
  })
  
  saveRDS(tbl_fil,file=paste0(sampleFolder,"results/Robjects_ADT_muscat/tbl_fil_muscat_DESeq2_new_",Comp,"_",sampleName,".rds"))
  write.xlsx(tbl_fil, paste0(sampleFolder,"results/Muscat_ADT/tbl_fil_muscat_DESeq2_new_",Comp,"_",sampleName,".xlsx"))
  
  saveRDS(tbl_fil_clint,file=paste0(sampleFolder,"results/Robjects_ADT_muscat/tbl_fil_muscat_DESeq2_new_",Comp,"_",sampleName,"_clint.rds"))
  write.xlsx(tbl_fil_clint, paste0(sampleFolder,"results/Muscat_ADT/tbl_fil_muscat_DESeq2_new_",Comp,"_",sampleName,"_clint.xlsx"))
  
}

### Extra ###
tbl_freq<- resDS(sce, res, frq = TRUE)
write.xlsx(tbl_freq, paste0(sampleFolder,"results/Muscat_ADT/tbl_full_freq_muscat_DESeq2_new_",Comparison,"_",sampleName,".xlsx"))

## Get list of cell populations
listCells<-levels(sce$cluster_id)[c(1,2)] #Empty lists!!

names(tbl_fil_clint)

pbHeatmap(sce, res, top_n = 20)

for (i in c(1,3,5,7,9,10,12,13,14)){ 
  H2 <- pbHeatmap(sce, res, top_n = 25, c = names(res$table[i]), k = listCells[1], normalize = T) 
  
  Comp<-names(res$table)[i]
  
  pdf(file=paste0(sampleFolder,"results/Muscat_ADT/Heatmap_Overview_DESeq2_new_normalized_",Comp,".pdf"), width = 8, height = 8)
  print(H2)
  dev.off()
  
}

for (i in c(2,4,6,8,15,16,17,18,19,20)){ 
  H2 <- pbHeatmap(sce, res, top_n = 25, c = names(res$table[i]), k = listCells[2], normalize = T) 
  
  Comp<-names(res$table)[i]
  
  pdf(file=paste0(sampleFolder,"results/Muscat_ADT/Heatmap_Overview_DESeq2_new_normalized_",Comp,".pdf"), width = 8, height = 8)
  print(H2)
  dev.off()
  
}

### Calculating expression frequencies
# It is often worthwhile to also consider the expression frequencies of each gene,
# i.e., the fraction of cells that express a given gene in each sample and/or group.
# `calcExprFreqs` computes cluster-sample/-group wise expression frequencies.
# Here, a gene is considered to be expressed when the specified measurement value (argument `assay`)
# falls above a certain threshold (argument `th`).
# Note that, `assay = "counts"` and `th = 0` (default) amounts to the fraction of cells for which a respective gene has been detected.
# `calcExprFreqs` will return a `r Biocpkg("SingleCellExperiment")` object, where sheets (assays) = clusters, rows = genes,
# and columns = samples (and groups, if `group_id`s are present in the `colData` of the input SCE).

frq <- calcExprFreqs(sce, assay = "counts", th = 0)
# one sheet per cluster
assayNames(frq)
# expression frequencies in each
# sample & group; 1st cluster
t(head(assay(frq), 5))

#Interested in groups: 
gids <- levels(sce$group_id)

#Check if larger than 10% expressed in one of those groups -> logical (T<->F)
frq10 <- vapply(as.list(assays(frq)),
                function(u) apply(u[, gids] > 0.1, 1, any),
                logical(nrow(sce)))
t(head(frq10))

#Filter table from before (already filtered on logFC and FDR) further on expr -> loop through subpops and only retain genes which are
#in the frq10 object for that subpop
## Create loop for all comparisons
for (i in 1:length(res$table)){ 
  tbl <- res$table[[i]]
  
  Comp<-names(res$table)[i]
  
  # filter clint -> here more strict because too many DE genes!!!!!!!!!!!!!
  tbl_fil_clint <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.glb < 0.01, abs(logFC) > 1, baseMean > 50)
    dplyr::arrange(u, p_adj.glb)
  })
  
  kids_v2<-names(tbl)
  tbl_fil2 <- lapply(kids_v2, function(k) dplyr::filter(tbl_fil_clint[[k]], gene %in% names(which(frq10[, k])))) ## EMDS stronger filtering for heatmaps!! (12/10/2023)
  names(tbl_fil2)<-kids_v2
  
  ## Save table:
  saveRDS(tbl_fil2,file=paste0(sampleFolder,"results/Robjects_ADT_muscat/tbl_fil2_muscat_DESeq2_new_",Comp,"_",sampleName,".rds"))
  write.xlsx(tbl_fil2, paste0(sampleFolder,"results/Muscat_ADT/tbl_fil2_muscat_DESeq2_new_",Comp,"_",sampleName,".xlsx"))
  
}

################################################################################################
################################################################################################
################################################################################################

## EMDS stronger filtering for heatmaps!! (12/10/2023)

## Get gene lists filtered stronger (global adj PV and higher basemean!!) 
genes_2h_cpgLNP_vs_eLNP<-readRDS("VBO_LNP_merge/results/Robjects_ADT_muscat/tbl_fil2_muscat_DESeq2_new_CpG_LNPs_2h-eLNPs_2h_VBO_4-12.rds")
genelist_2h_cpgLNP_vs_eLNP<-genes_2h_cpgLNP_vs_eLNP$`Early mature cDC1s`$gene
genes_2h_picLNP_vs_eLNP<-readRDS("VBO_LNP_merge/results/Robjects_ADT_muscat/tbl_fil2_muscat_DESeq2_new_eLNPs_2h-pIC_LNPs_2h_VBO_4-12.rds")
genelist_2h_picLNP_vs_eLNP<-genes_2h_picLNP_vs_eLNP$`Early mature cDC1s`$gene

genelist_2h<-c(setdiff(genelist_2h_cpgLNP_vs_eLNP,genelist_2h_picLNP_vs_eLNP),intersect(genelist_2h_cpgLNP_vs_eLNP,genelist_2h_picLNP_vs_eLNP),
               setdiff(genelist_2h_picLNP_vs_eLNP,genelist_2h_cpgLNP_vs_eLNP))

genes_8h_cpgLNP_vs_eLNP<-readRDS("VBO_LNP_merge/results/Robjects_ADT_muscat/tbl_fil2_muscat_DESeq2_new_CpG_LNPs_8h-eLNPs_8h_VBO_4-12.rds")
genelist_8h_cpgLNP_vs_eLNP<-genes_8h_cpgLNP_vs_eLNP$`Late mature cDC1s`$gene
genes_8h_picLNP_vs_eLNP<-readRDS("VBO_LNP_merge/results/Robjects_ADT_muscat/tbl_fil2_muscat_DESeq2_new_eLNPs_8h-pIC_LNPs_8h_VBO_4-12.rds")
genelist_8h_picLNP_vs_eLNP<-genes_8h_picLNP_vs_eLNP$`Late mature cDC1s`$gene

genelist_8h<-c(setdiff(genelist_8h_cpgLNP_vs_eLNP,genelist_8h_picLNP_vs_eLNP),intersect(genelist_8h_cpgLNP_vs_eLNP,genelist_8h_picLNP_vs_eLNP),
               setdiff(genelist_8h_picLNP_vs_eLNP,genelist_8h_cpgLNP_vs_eLNP))

## Extra test with seurat average heatmaps (05/04/23)
## Subset to only keep pops to compare
seuratObj$annotated_clusters_Muscat_v2<-paste0(seuratObj$orig.ident,"_",seuratObj$annotated_clusters_Muscat)

DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2")

Idents(seuratObj)<-seuratObj$annotated_clusters_Muscat_v2
seuratObj_subset_EM<-subset(seuratObj, idents = levels(Idents(seuratObj))[c(2,29,36,37,42)]) 
seuratObj_subset_LM<-subset(seuratObj, idents = levels(Idents(seuratObj))[c(4,6,11,17,22)]) 

DimPlot(seuratObj_subset_EM, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2", split.by = "orig.ident", ncol = 3)
DimPlot(seuratObj_subset_LM, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2", split.by = "orig.ident", ncol = 3)

seuratObj_subset_EM$Condition_v2<-as.character(seuratObj_subset_EM$Condition)
seuratObj_subset_LM$Condition_v2<-as.character(seuratObj_subset_LM$Condition)

seuratObj_subset_EM$Condition_v2<-gsub("-","_",seuratObj_subset_EM$Condition_v2)
seuratObj_subset_EM$Condition_v2<-gsub(" ","_",seuratObj_subset_EM$Condition_v2)
seuratObj_subset_LM$Condition_v2<-gsub("-","_",seuratObj_subset_LM$Condition_v2)
seuratObj_subset_LM$Condition_v2<-gsub(" ","_",seuratObj_subset_LM$Condition_v2)

seuratObj_subset_EM$Condition_v2<-factor(as.character(seuratObj_subset_EM$Condition_v2), levels = c("Steady_state","eLNPs_2h","pIC_LNPs_2h","CpG_LNPs_2h","pIC_alone_2h"))
seuratObj_subset_LM$Condition_v2<-factor(as.character(seuratObj_subset_LM$Condition_v2), levels = c("Steady_state","eLNPs_8h","pIC_LNPs_8h","CpG_LNPs_8h","pIC_alone_8h"))

## Combine MULTI_ID with Orig.ident or Genotype?? Need to put genotype in correct order: otherwise WT after DKO
seuratObj_subset_EM@meta.data[["MULTI_ID_merge"]]
seuratObj_subset_EM@meta.data$MULTI_ID_merge<-as.factor(as.character(seuratObj_subset_EM@meta.data$MULTI_ID_merge))
seuratObj_subset_LM@meta.data[["MULTI_ID_merge"]]
seuratObj_subset_LM@meta.data$MULTI_ID_merge<-as.factor(as.character(seuratObj_subset_LM@meta.data$MULTI_ID_merge))

## Update muscat annotation
seuratObj_subset_EM$annotated_clusters_Muscat<-gsub(" ","_",as.character(seuratObj_subset_EM$annotated_clusters_Muscat))
seuratObj_subset_LM$annotated_clusters_Muscat<-gsub(" ","_",as.character(seuratObj_subset_LM$annotated_clusters_Muscat))

seuratObj_subset_EM$annotated_clusters_Muscat<-factor(as.character(seuratObj_subset_EM$annotated_clusters_Muscat), 
                                                      levels = c("Early_mature_cDC1s"))
seuratObj_subset_LM$annotated_clusters_Muscat<-factor(as.character(seuratObj_subset_LM$annotated_clusters_Muscat), 
                                                      levels = c("Late_mature_cDC1s"))

## Try new sample_ID
seuratObj_subset_EM$sample_ID<-paste0(as.character(seuratObj_subset_EM$MULTI_ID),"_",as.character(seuratObj_subset_EM$Condition_v2),
                                      "_",as.character(seuratObj_subset_EM$annotated_clusters_Muscat))
seuratObj_subset_LM$sample_ID<-paste0(as.character(seuratObj_subset_LM$MULTI_ID),"_",as.character(seuratObj_subset_LM$Condition_v2),
                                      "_",as.character(seuratObj_subset_LM$annotated_clusters_Muscat))

seuratObj_subset_EM$sample_ID<-factor(seuratObj_subset_EM$sample_ID, #Added Imm replicates
                                      levels = c("Hashtag1_Steady_state_Early_mature_cDC1s","Hashtag2_Steady_state_Early_mature_cDC1s","Hashtag3_Steady_state_Early_mature_cDC1s","Hashtag4_Steady_state_Early_mature_cDC1s",
                                                 "Hashtag1_eLNPs_2h_Early_mature_cDC1s","Hashtag2_eLNPs_2h_Early_mature_cDC1s","Hashtag3_eLNPs_2h_Early_mature_cDC1s","Hashtag4_eLNPs_2h_Early_mature_cDC1s",       
                                                 "Hashtag1_pIC_LNPs_2h_Early_mature_cDC1s","Hashtag2_pIC_LNPs_2h_Early_mature_cDC1s","Hashtag3_pIC_LNPs_2h_Early_mature_cDC1s","Hashtag4_pIC_LNPs_2h_Early_mature_cDC1s",
                                                 "Hashtag1_CpG_LNPs_2h_Early_mature_cDC1s","Hashtag2_CpG_LNPs_2h_Early_mature_cDC1s","Hashtag3_CpG_LNPs_2h_Early_mature_cDC1s","Hashtag4_CpG_LNPs_2h_Early_mature_cDC1s",
                                                 "Hashtag1_pIC_alone_2h_Early_mature_cDC1s","Hashtag2_pIC_alone_2h_Early_mature_cDC1s","Hashtag3_pIC_alone_2h_Early_mature_cDC1s","Hashtag4_pIC_alone_2h_Early_mature_cDC1s"))
levels(seuratObj_subset_EM$sample_ID)

seuratObj_subset_LM$sample_ID<-factor(seuratObj_subset_LM$sample_ID, #Added Imm replicates
                                      levels = c("Hashtag1_Steady_state_Late_mature_cDC1s","Hashtag2_Steady_state_Late_mature_cDC1s","Hashtag3_Steady_state_Late_mature_cDC1s","Hashtag4_Steady_state_Late_mature_cDC1s",
                                                 "Hashtag1_eLNPs_8h_Late_mature_cDC1s","Hashtag2_eLNPs_8h_Late_mature_cDC1s","Hashtag3_eLNPs_8h_Late_mature_cDC1s","Hashtag4_eLNPs_8h_Late_mature_cDC1s",
                                                 "Hashtag1_pIC_LNPs_8h_Late_mature_cDC1s","Hashtag2_pIC_LNPs_8h_Late_mature_cDC1s","Hashtag3_pIC_LNPs_8h_Late_mature_cDC1s","Hashtag4_pIC_LNPs_8h_Late_mature_cDC1s",
                                                 "Hashtag1_CpG_LNPs_8h_Late_mature_cDC1s","Hashtag2_CpG_LNPs_8h_Late_mature_cDC1s","Hashtag3_CpG_LNPs_8h_Late_mature_cDC1s",
                                                 "Hashtag1_pIC_alone_8h_Late_mature_cDC1s","Hashtag2_pIC_alone_8h_Late_mature_cDC1s","Hashtag3_pIC_alone_8h_Late_mature_cDC1s","Hashtag4_pIC_alone_8h_Late_mature_cDC1s"))
levels(seuratObj_subset_LM$sample_ID)

## Update clusters
seuratObj_subset_EM$annotated_clusters_Muscat_v2<-factor(paste0(seuratObj_subset_EM$Condition_v2,"_",seuratObj_subset_EM$annotated_clusters_Muscat),
                                                         levels = c("Steady_state_Early_mature_cDC1s",
                                                                    "eLNPs_2h_Early_mature_cDC1s","pIC_LNPs_2h_Early_mature_cDC1s",
                                                                    "CpG_LNPs_2h_Early_mature_cDC1s","pIC_alone_2h_Early_mature_cDC1s"))
DimPlot(seuratObj_subset_EM, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2", split.by = "Condition_v2", ncol = 3)

seuratObj_subset_LM$annotated_clusters_Muscat_v2<-factor(paste0(seuratObj_subset_LM$Condition_v2,"_",seuratObj_subset_LM$annotated_clusters_Muscat),
                                                         levels = c("Steady_state_Late_mature_cDC1s",
                                                                    "eLNPs_8h_Late_mature_cDC1s","pIC_LNPs_8h_Late_mature_cDC1s",
                                                                    "CpG_LNPs_8h_Late_mature_cDC1s","pIC_alone_8h_Late_mature_cDC1s"))
DimPlot(seuratObj_subset_LM, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2",split.by = "Condition_v2", ncol = 3)

## Basic annotation
seuratObj_subset_EM$Basic_annotation<-"cDC1s"
seuratObj_subset_LM$Basic_annotation<-"cDC1s"


## Look at average expression per sample
Idents(seuratObj_subset_EM)<-seuratObj_subset_EM$sample_ID
seuratObj_subset_EM_average <- AverageExpression(seuratObj_subset_EM, return.seurat = T)

Idents(seuratObj_subset_LM)<-seuratObj_subset_LM$sample_ID
seuratObj_subset_LM_average <- AverageExpression(seuratObj_subset_LM, return.seurat = T)

seuratObj_subset_EM_average <- ScaleData(seuratObj_subset_EM_average, assay = "ADT")
seuratObj_subset_LM_average <- ScaleData(seuratObj_subset_LM_average, assay = "ADT")


## Create heatmaps (same code as above -> same variables for genelist, but now filtered stricter. New filenames to save!)
library("RColorBrewer")
Colset_EM<-rep(brewer.pal(n = 5, name = "Set3"), each = 4)
Colset_LM<-Colset_EM[-13] #One sample less for CpG_LNPs_8h
Colset_EM_3samples<-Colset_EM[-c(1:4,17:20)] #Only 3 groups
Colset_LM_3samples<-Colset_LM[-c(1:4,16:19)] #Only 3 groups

## Second batch
H2.2h <- DoHeatmap(seuratObj_subset_EM_average, assay = "ADT", cells = levels(Idents(seuratObj_subset_EM_average))[c(5:16)],
                   features = intersect(genelist_2h_cpgLNP_vs_eLNP,genelist_2h_picLNP_vs_eLNP), size = 3,
                   draw.lines = FALSE, group.colors = Colset_EM) + 
  scale_color_manual(values = Colset_EM_3samples) +
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill")

H1.8h <- DoHeatmap(seuratObj_subset_LM_average, assay = "ADT", cells = levels(Idents(seuratObj_subset_LM_average))[c(5:15)],
                   features = setdiff(genelist_8h_cpgLNP_vs_eLNP,genelist_8h_picLNP_vs_eLNP), size = 3,
                   draw.lines = FALSE, group.colors = Colset_LM) + 
  scale_color_manual(values = Colset_LM_3samples) +
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill")
H2.8h <- DoHeatmap(seuratObj_subset_LM_average, assay = "ADT", cells = levels(Idents(seuratObj_subset_LM_average))[c(5:15)],
                   features = intersect(genelist_8h_cpgLNP_vs_eLNP,genelist_8h_picLNP_vs_eLNP), size = 3,
                   draw.lines = FALSE, group.colors = Colset_LM) + 
  scale_color_manual(values = Colset_LM_3samples) +
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill")
H3.8h <- DoHeatmap(seuratObj_subset_LM_average, assay = "ADT", cells = levels(Idents(seuratObj_subset_LM_average))[c(5:15)],
                   features = setdiff(genelist_8h_picLNP_vs_eLNP,genelist_8h_cpgLNP_vs_eLNP), size = 3,
                   draw.lines = FALSE, group.colors = Colset_LM) + 
  scale_color_manual(values = Colset_LM_3samples) +
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                        mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill")

pdf(file=paste0(sampleFolder,"results/Muscat_ADT/Heatmap_seurat_2h_Sofie_common_CpG_and_pIC_LNP_vs_eLNP_3samples_EMDS_strict.pdf"), width = 20, height = 15)
print(H2.2h)
dev.off()
pdf(file=paste0(sampleFolder,"results/Muscat_ADT/Heatmap_seurat_8h_Sofie_specific_CpG_LNP_vs_eLNP_3samples_EMDS_strict.pdf"), width = 20, height = 15)
print(H1.8h)
dev.off()
pdf(file=paste0(sampleFolder,"results/Muscat_ADT/Heatmap_seurat_8h_Sofie_common_CpG_andpIC_LNP_vs_eLNP_3samples_EMDS_strict.pdf"), width = 20, height = 15)
print(H2.8h)
dev.off()
pdf(file=paste0(sampleFolder,"results/Muscat_ADT/Heatmap_seurat_8h_Sofie_specific_pIC_LNP_vs_eLNP_3samples_EMDS_strict.pdf"), width = 20, height = 15)
print(H3.8h)
dev.off()