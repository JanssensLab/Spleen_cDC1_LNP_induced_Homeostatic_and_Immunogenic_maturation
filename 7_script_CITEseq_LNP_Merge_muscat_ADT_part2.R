## Script for analyzing LNP cDC1 CITE-seq project data
## Muscat script for comparing between conditions and maturation stages part 2 ADT
## See good diff between common/homeo/immuno genes: LM eLNP with Imm SS and LM pIC-LNP

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

### Modified function for robust edgeR
source('~/VIB/DATA/Roos/Cara/pbDSrob.R')
source('~/VIB/DATA/Roos/Cara/pbDS_DESeq2.R')

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

## Update annotation
seuratObj$annotated_clusters_Muscat_v2<-gsub(" ","_",seuratObj$annotated_clusters_Muscat)

## Update condition
seuratObj$Condition_v2<-as.factor(as.character(seuratObj$Condition))
seuratObj$Condition_v2<-gsub("-","_",seuratObj$Condition_v2)
seuratObj$Condition_v2<-gsub(" ","_",seuratObj$Condition_v2)

## Try new sample_ID
seuratObj$sample_ID<-paste0(as.character(seuratObj$MULTI_ID),"_",as.character(seuratObj$Condition_v2),
                            "_",as.character(seuratObj$annotated_clusters_Muscat_v2))
seuratObj$sample_ID<-as.factor(seuratObj$sample_ID)
levels(seuratObj$sample_ID)

## Update clusters
seuratObj$annotated_clusters_Muscat_v3<-as.factor(paste0(seuratObj$Condition_v2,"_",seuratObj$annotated_clusters_Muscat_v2))
DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v3",split.by = "Condition_v2", ncol = 3)

## Basic annotation
seuratObj$Basic_annotation<-"cDC1s"

## Subset the seuratObj to only the samples we want to compare
Idents(seuratObj)<-seuratObj$sample_ID
seuratObj_subset<-subset(seuratObj, idents = c(levels(seuratObj$sample_ID)[grep("CpG_LNPs_2h_Immature_cDC1s",levels(seuratObj$sample_ID))], 
                                               levels(seuratObj$sample_ID)[grep("CpG_LNPs_2h_Early_mature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("CpG_LNPs_8h_Late_mature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("eLNPs_2h_Immature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("eLNPs_2h_Early_mature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("eLNPs_8h_Late_mature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("pIC_LNPs_2h_Immature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("pIC_LNPs_2h_Early_mature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("pIC_LNPs_8h_Late_mature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("pIC_alone_2h_Immature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("pIC_alone_2h_Early_mature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("pIC_alone_8h_Late_mature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("Steady_state_Immature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("Steady_state_Early_mature_cDC1s",levels(seuratObj$sample_ID))],
                                               levels(seuratObj$sample_ID)[grep("Steady_state_Late_mature_cDC1s",levels(seuratObj$sample_ID))]))

seuratObj_subset$sample_ID<-as.factor(as.character(seuratObj_subset$sample_ID))
Idents(seuratObj_subset)<-seuratObj_subset$sample_ID
seuratObj_subset$annotated_clusters_Muscat_v3<-factor(as.character(seuratObj_subset$annotated_clusters_Muscat_v3), 
                                                      levels = c("Steady_state_Immature_cDC1s","Steady_state_Early_mature_cDC1s","Steady_state_Late_mature_cDC1s",
                                                                 "eLNPs_2h_Immature_cDC1s","eLNPs_2h_Early_mature_cDC1s","eLNPs_8h_Late_mature_cDC1s", 
                                                                 "pIC_alone_2h_Immature_cDC1s","pIC_alone_2h_Early_mature_cDC1s","pIC_alone_8h_Late_mature_cDC1s",
                                                                 "CpG_LNPs_2h_Immature_cDC1s","CpG_LNPs_2h_Early_mature_cDC1s","CpG_LNPs_8h_Late_mature_cDC1s",
                                                                 "pIC_LNPs_2h_Immature_cDC1s","pIC_LNPs_2h_Early_mature_cDC1s","pIC_LNPs_8h_Late_mature_cDC1s"))
DimPlot(seuratObj_subset, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v3", split.by = "Condition")

########################################
##### Load data (already filtered) 
########################################

### Convert to sce ### RNA or SCT!!!!!!
sce <- as.SingleCellExperiment(seuratObj_subset, assay = "ADT")

########################################
##### Prepare sce object
########################################

### Create test metadata table
sample_id<-levels(as.factor(as.character(seuratObj_subset$sample_ID)))
group_id<-as.character()
for (i in 1:length(sample_id)){
  group_id[i]<-gsub(substr(sample_id[i],1,9), "", sample_id[i])
}
n_cells<-rep(1,59)

### Create metaData matrix
metaData<-data.frame(sample_id,group_id,n_cells)

### Look through samples (1->8) for amount of cells after filtering
for(i in 1:length(sample_id)){
  toSearch<-sample_id[i]
  metaData[i, "n_cells"]<-length(grep(paste0("\\<",toSearch,"\\>"),seuratObj_subset$sample_ID)) ##Need to do specific grep. Otherwise fault between Hashtag 1 and Hashtag 10!!!
}

## Use update annotation for group_id
sce$group_id<-sce$annotated_clusters_Muscat_v3

### Add new sample_id to sce object (combination of orig.ident and group_id)
sce$sample_id<-sce$sample_ID

### Finish prep object
(sce <- prepSCE(sce, 
                kid = "Basic_annotation", # cell population assignments
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


dir.create(paste0(sampleFolder,"results/Muscat_ADT_Thesis_Victor/"))

Comparison<-"DS_ADT_thesis_Victor"

png(file=paste0(sampleFolder,"results/Muscat_ADT_Thesis_Victor/meanVariancePlot_subset_",Comparison,".png"),width = 3500, height = 2000, res = 300)
pb_mds #+ scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set2"))
dev.off()

## Extra August 2023 for paper:
## Create subset of sce without other Imm cDC1 clusters (only keep SS)
sce_subset<-sce[,which(!sce$group_id %in% c("eLNPs_2h_Immature_cDC1s","pIC_alone_2h_Immature_cDC1s","CpG_LNPs_2h_Immature_cDC1s","pIC_LNPs_2h_Immature_cDC1s"))]
pb_subset <- aggregateData(sce_subset,
                           assay = "counts", fun = "sum",
                           by = c("cluster_id", "sample_id"))
(pb_mds_subset <- pbMDS(pb_subset) + scale_shape_manual(values=1:nlevels(pb_subset$group_id)))

pdf(file=paste0(sampleFolder,"results/Muscat_ADT_Thesis_Victor/meanVariancePlot_subset_paper_",Comparison,".pdf"),width = 15, height = 8)
pb_mds_subset + scale_color_manual(values = c("cadetblue1",brewer.pal(10, "Paired")))
dev.off()

####################################################################

#################################################
###### First analysis: basic design matrix ######
#################################################

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
contrast <- makeContrasts("CpG_LNPs_8h_Late_mature_cDC1s-Steady_state_Immature_cDC1s",
                          "eLNPs_8h_Late_mature_cDC1s-Steady_state_Immature_cDC1s",
                          "pIC_LNPs_8h_Late_mature_cDC1s-Steady_state_Immature_cDC1s",
                          "Steady_state_Late_mature_cDC1s-Steady_state_Immature_cDC1s",
                          "CpG_LNPs_8h_Late_mature_cDC1s-eLNPs_8h_Late_mature_cDC1s",
                          "pIC_LNPs_8h_Late_mature_cDC1s-eLNPs_8h_Late_mature_cDC1s",
                          "CpG_LNPs_8h_Late_mature_cDC1s-Steady_state_Late_mature_cDC1s",
                          "pIC_LNPs_8h_Late_mature_cDC1s-Steady_state_Late_mature_cDC1s",
                          levels = mm)

# run DS analysis
### DeSeq2 method ###
res_DESeq2<-pbDS(pb, design = mm, contrast = contrast,method = "DESeq2")

## Save res object
saveRDS(res_DESeq2,file=paste0(sampleFolder,"results/Robjects_Thesis_Victor/res_Muscat_DeSeq2_new_",Comparison,"_",sampleName,".rds"))

## Choose Res object!
res<-readRDS(file=paste0(sampleFolder,"results/Robjects_Thesis_Victor/res_Muscat_DeSeq2_new_",Comparison,"_",sampleName,".rds")) ##DESeq2

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
  write.xlsx(tbl, paste0(sampleFolder,"results/Muscat_ADT_Thesis_Victor/tbl_full_muscat_DESeq2_new_",Comparison,"_",Comp,"_",sampleName,".xlsx"))
  
  # filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
  tbl_fil <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
    dplyr::arrange(u, p_adj.loc)
  })
  
  # filter clint -> here more strict because too many DE genes!!!!!!!!!!!!!
  tbl_fil_clint <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.glb < 0.01, abs(logFC) > 1, baseMean > 50)
    dplyr::arrange(u, p_adj.glb)
  })
  
  saveRDS(tbl_fil,file=paste0(sampleFolder,"results/Robjects_Thesis_Victor/tbl_fil_muscat_DESeq2_new_",Comparison,"_",Comp,"_",sampleName,".rds"))
  write.xlsx(tbl_fil, paste0(sampleFolder,"results/Muscat_ADT_Thesis_Victor/tbl_fil_muscat_DESeq2_new_",Comparison,"_",Comp,"_",sampleName,".xlsx"))
  
  saveRDS(tbl_fil_clint,file=paste0(sampleFolder,"results/Robjects_Thesis_Victor/tbl_fil_muscat_DESeq2_new_",Comparison,"_",Comp,"_",sampleName,"_clint.rds"))
  write.xlsx(tbl_fil_clint, paste0(sampleFolder,"results/Muscat_ADT_Thesis_Victor/tbl_fil_muscat_DESeq2_new_",Comparison,"_",Comp,"_",sampleName,"_clint.xlsx"))
  
}

### Extra ###
tbl_freq<- resDS(sce, res, frq = TRUE)
write.xlsx(tbl_freq, paste0(sampleFolder,"results/Muscat_ADT_Thesis_Victor/tbl_full_freq_muscat_DESeq2_new_",Comparison,"_",sampleName,".xlsx"))

## Get list of cell populations
names(tbl_fil_clint)
listCells<-levels(sce$cluster_id)

pbHeatmap(sce, res, top_n = 20)

for (i in 1:length(res$table)){
  H2 <- pbHeatmap(sce, res, top_n = 25, c = names(res$table[i]), k = listCells, normalize = T) 
  
  Comp<-names(res$table)[i]
  
  pdf(file=paste0(sampleFolder,"results/Muscat_ADT_Thesis_Victor/Heatmap_Overview_DESeq2_new_normalized_",Comparison,"_",Comp,".pdf"), width = 10, height = 8)
  print(H2)
  dev.off()
  
}

