## Script for analyzing LNP cDC1 CITE-seq project data
## Muscat script for comparing between conditions part 1 RNA
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
library("limma")
library("clusterProfiler")
library("org.Mm.eg.db")

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

## Combine MULTI_ID with Orig.ident or Genotype
seuratObj@meta.data[["MULTI_ID_merge"]]
seuratObj@meta.data$MULTI_ID_merge<-as.factor(seuratObj@meta.data$MULTI_ID_merge)

########################################
##### Load data (already filtered) 
########################################

#Convert sparse matrix to matrix 
seuratObj[['RNA']]@counts <- as.matrix(seuratObj[['RNA']]@counts) #Necessary to avoid issue with calcExprFreqs step in DS analysis
seuratObj[['RNA']]@data <- as.matrix(seuratObj[['RNA']]@data) #Necessary to avoid issue with calcExprFreqs step in DS analysis

### Convert to sce ### 
sce <- as.SingleCellExperiment(seuratObj, assay = "RNA")

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

### Look through samples for amount of cells after filtering
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

# Removing cluster-sample instance(s) ‘Ly6c2+ Mature NK cells’-‘DKO_Hashtag3’, ‘Ifitm+ Mature NK cells’-‘WT_Hashtag3’, 
# ‘Ifitm+ Mature NK cells’-‘DKO_Hashtag3’, ‘NK cells’-‘DKO_Hashtag3’, ‘Spp1+Ifng+Ccl3+Ccl4+ Immature NK cells’-‘WT_Hashtag3’, 
# ‘Spp1+Ifng+Ccl3+Ccl4+ Immature NK cells’-‘DKO_Hashtag3’, ‘Gzmc+Ifitm+ Mature NK cells’-‘WT_Hashtag3’, 
# ‘Gzmc+Ifitm+ Mature NK cells’-‘DKO_Hashtag3’, ‘Klra5+ Mature NK cells’-‘DKO_Hashtag3’

dir.create(paste0(sampleFolder,"results/Muscat/"))

Comparison<-"DS_Simple_RNA_Harmony_clustering"

png(file=paste0(sampleFolder,"results/Muscat/meanVariancePlot_subset_",Comparison,".png"),width = 3500, height = 2000, res = 300)
pb_mds + scale_color_manual(values = RColorBrewer::brewer.pal(6, "Set2"))
dev.off()

### Experimental design
# We can provide `pbDS` with a design matrix capturing the experimental design using `model.matrix` (package `r Rpackage("stats")`), 
# and a contrast matrix that specifies our comparison of interesting using `makeContrasts` from the `r Biocpkg("limma")` package. 
# Alternatively, the comparison(s) of interest (or a list thereof) can be specified with via `coefs` (see `?glmQLFTest` for details). 

# Here, we want to carry out a single comparison of DKO against WT samples, 
# thus placing `"WT"` on the right-hand side as the reference condition.
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
saveRDS(res_DESeq2,file=paste0(sampleFolder,"results/Robjects/res_Muscat_DeSeq2_new_",Comparison,"_",sampleName,".rds"))

## Choose Res object!
res<-readRDS(file=paste0(sampleFolder,"results/Robjects/res_Muscat_DeSeq2_new_",Comparison,"_",sampleName,".rds")) ##DESeq2

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
  write.xlsx(tbl, paste0(sampleFolder,"results/Muscat/tbl_full_muscat_DESeq2_new_",Comp,"_",sampleName,".xlsx"))
  
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
  
  saveRDS(tbl_fil,file=paste0(sampleFolder,"results/Robjects/tbl_fil_muscat_DESeq2_new_",Comp,"_",sampleName,".rds"))
  write.xlsx(tbl_fil, paste0(sampleFolder,"results/Muscat/tbl_fil_muscat_DESeq2_new_",Comp,"_",sampleName,".xlsx"))
  
  saveRDS(tbl_fil_clint,file=paste0(sampleFolder,"results/Robjects/tbl_fil_muscat_DESeq2_new_",Comp,"_",sampleName,"_clint.rds"))
  write.xlsx(tbl_fil_clint, paste0(sampleFolder,"results/Muscat/tbl_fil_muscat_DESeq2_new_",Comp,"_",sampleName,"_clint.xlsx"))
  
}

### Extra ###
tbl_freq<- resDS(sce, res, frq = TRUE)
write.xlsx(tbl_freq, paste0(sampleFolder,"results/Muscat/tbl_full_freq_muscat_DESeq2_new_",Comparison,"_",sampleName,".xlsx"))


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
  saveRDS(tbl_fil2,file=paste0(sampleFolder,"results/Robjects/tbl_fil2_muscat_DESeq2_new_",Comp,"_",sampleName,".rds"))
  write.xlsx(tbl_fil2, paste0(sampleFolder,"results/Muscat/tbl_fil2_muscat_DESeq2_new_",Comp,"_",sampleName,".xlsx"))
  
}


## Get list of cell populations
listCells<-levels(sce$cluster_id)[c(1,2)] #Empty lists!!

names(tbl_fil_clint)

pbHeatmap(sce, res, top_n = 20)

for (i in 1:length(res$table)){
  H2 <- pbHeatmap(sce, res, top_n = 25, c = names(res$table[i]), k = listCells, normalize = T) 
  
  Comp<-names(res$table)[i]
  
  pdf(file=paste0(sampleFolder,"results/Muscat/Heatmap_Overview_DESeq2_new_normalized_",Comp,".pdf"), width = 8, height = 15)
  print(H2)
  dev.off()
  
}

####################################################

## GSEA in R:

#Create dir
dir.create(paste0(sampleFolder,"results/GSEA_Muscat/"))

#Load markers
## Choose Res object!
res<-readRDS(file=paste0(sampleFolder,"results/Robjects/res_Muscat_DeSeq2_new_",Comparison,"_",sampleName,".rds")) ##DESeq2

for (i in 1:length(res$table)){
  
  tbl <- res$table[[i]]
  
  # filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
  tbl_fil <- lapply(tbl, function(u) {
    u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
    dplyr::arrange(u, p_adj.loc)
  })
  
  Comp<-names(res$table)[i]
  
  ## All cell clusters
  EnrichGO<-tibble::lst()
  for (cell_pop in 1:length(tbl_fil)) {
    #Background:
    Background_scRNAseq<-as.character(tbl[[cell_pop]]$gene) #full table!!
    #EnrichGO:
    EnrichGO[cell_pop]<-enrichGO(
      as.character(tbl_fil[[cell_pop]]$gene),
      'org.Mm.eg.db',
      keyType = "SYMBOL",
      ont = "ALL",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe = Background_scRNAseq,
      qvalueCutoff = 0.2,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE,
      pool = T
    )
  }
  
  ## Save results
  saveRDS(EnrichGO, paste0(sampleFolder,"results/GSEA_Muscat/EnrichGO_results_background_scRNAseq_",Comp,"_",sampleName,".rds"))

  for (cell_pop in (1:length(tbl_fil))) { #1:length(tbl_fil) -> skip 10 and 13 because 0 results!!
    ## Create dotplot top 10 each category
    D_background<-dotplot(EnrichGO[[cell_pop]], split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    
    ## Save files
    pdf(file=paste0(sampleFolder,"results/GSEA_Muscat/Dotplot_background_scRNAseq_",names(tbl_fil)[cell_pop],"_",Comp,"_",sampleName,".pdf"), width = 15, height = 10)
    print(D_background)
    dev.off()
    
    write.xlsx(EnrichGO[[cell_pop]]@result,file=paste0(sampleFolder,"results/GSEA_Muscat/EnrichGO_background_scRNAseq_",names(tbl_fil)[cell_pop],"_",Comp,"_",sampleName,".xlsx"))
  }
  
}

#####################

## Extra check polyIC vs polyIC_LNP (13/12/22)
tbl_pIC <- res$table[[20]]$`Late mature cDC1s`

# filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
tbl_pIC_fil <- dplyr::filter(tbl_pIC, p_adj.loc < 0.05, abs(logFC) > 1)
tbl_pIC_fil<-dplyr::arrange(tbl_pIC_fil, p_adj.loc)

tbl_pIC_fil_UP <- dplyr::filter(tbl_pIC, p_adj.loc < 0.05, logFC > 1)
tbl_pIC_fil_DOWN <- dplyr::filter(tbl_pIC, p_adj.loc < 0.05, logFC < -1)

Comp<-names(res$table)[20]

#Background:
Background_scRNAseq<-as.character(tbl_pIC$gene) #full table!!
#EnrichGO:
EnrichGO_pIC_UP<-enrichGO(
  as.character(tbl_pIC_fil_UP$gene),
  'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = Background_scRNAseq,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T
)

EnrichGO_pIC_DOWN<-enrichGO(
  as.character(tbl_pIC_fil_DOWN$gene),
  'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = Background_scRNAseq,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T
)

D_background_pIC_UP<-dotplot(EnrichGO_pIC_UP, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
D_background_pIC_DOWN<-dotplot(EnrichGO_pIC_DOWN, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Save files
pdf(file=paste0(sampleFolder,"results/GSEA_Muscat/Dotplot_background_scRNAseq_Late_mature_cDC1s_UP_",Comp,"_",sampleName,".pdf"), width = 15, height = 10)
print(D_background_pIC_UP)
dev.off()

pdf(file=paste0(sampleFolder,"results/GSEA_Muscat/Dotplot_background_scRNAseq_Late_mature_cDC1s_DOWN_",Comp,"_",sampleName,".pdf"), width = 15, height = 10)
print(D_background_pIC_DOWN)
dev.off()

write.xlsx(EnrichGO_pIC_UP@result,file=paste0(sampleFolder,"results/GSEA_Muscat/EnrichGO_background_scRNAseq_Late_mature_cDC1s_UP_",Comp,"_",sampleName,".xlsx"))
write.xlsx(EnrichGO_pIC_DOWN@result,file=paste0(sampleFolder,"results/GSEA_Muscat/EnrichGO_background_scRNAseq_Late_mature_cDC1s_DOWN_",Comp,"_",sampleName,".xlsx"))
