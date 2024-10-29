## Script for analyzing LNP cDC1 CITE-seq project data
## Detailed pipeline run on the merge of the 9 samples

######################################################################
################ SPECIFY INPUT  AND OUTPUT ###########################
######################################################################

setwd('/home/clintdn/VIB/DATA/Sophie/scRNA-seq_Victor/LNP_CITEseq_experiment/VBO_LNP_merge/')

output.dir <- "results/"

samplename <- "VBO_4-12"

experiment <- samplename

source('/home/clintdn/VIB/DATA/Sophie/RNA-seq_Sandra/CITEseq_Test/RAW_DATA/script_functions_COVID.R') 

################################################################################
############################## RUN SCRIPT  ##################################### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ START HTO
################################################################################

library("tidyverse")
library("Seurat")
library("SingleCellExperiment")
library("scater")
library("data.table")
library("ggpubr")
library("gridExtra")
library("fitdistrplus") 
library("dplyr")
library("plyr")
library("caret")
library("openxlsx")
library("ggplot2")
library("tibble")
library("scran")
library("limma")
library("ggrepel")
library("intrinsicDimension")
library("tidyr")
library("purrr")
library("ggrepel")
library("clustree")
library("limma")
library("patchwork")
library("RColorBrewer")
library('cowplot')
library("harmony")
library("DisneyTools")
library("gplots")

library(future)
plan("multiprocess", workers = 6)

dir.create(output.dir)
dir.create(paste0(output.dir,"Robjects"))
dir.create(paste0(output.dir,"Plots"))
dir.create(paste0(output.dir,"Plots/RNA"))

diagnostics<-list()

## Load in seuratobjects
seuratObj_VBO004 <- readRDS(file="../VBO004/results/Robjects/seuratObj_VBO004_clean.rds")
seuratObj_VBO005 <- readRDS(file="../VBO005/results/Robjects/seuratObj_VBO005_clean.rds")
seuratObj_VBO006 <- readRDS(file="../VBO006/results/Robjects/seuratObj_VBO006_clean.rds")
seuratObj_VBO007 <- readRDS(file="../VBO007/results/Robjects/seuratObj_VBO007_clean.rds")
seuratObj_VBO008 <- readRDS(file="../VBO008/results/Robjects/seuratObj_VBO008_clean.rds")
seuratObj_VBO009 <- readRDS(file="../VBO009/results/Robjects/seuratObj_VBO009_clean.rds")
seuratObj_VBO010 <- readRDS(file="../VBO010/results/Robjects/seuratObj_VBO010_clean.rds")
seuratObj_VBO011 <- readRDS(file="../VBO011/results/Robjects/seuratObj_VBO011_clean.rds")
seuratObj_VBO012 <- readRDS(file="../VBO012/results/Robjects/seuratObj_VBO012_clean.rds")

seuratObj_VBO004<-DietSeurat(seuratObj_VBO004)
seuratObj_VBO005<-DietSeurat(seuratObj_VBO005)
seuratObj_VBO006<-DietSeurat(seuratObj_VBO006)
seuratObj_VBO007<-DietSeurat(seuratObj_VBO007)
seuratObj_VBO008<-DietSeurat(seuratObj_VBO008)
seuratObj_VBO009<-DietSeurat(seuratObj_VBO009)
seuratObj_VBO010<-DietSeurat(seuratObj_VBO010)
seuratObj_VBO011<-DietSeurat(seuratObj_VBO011)
seuratObj_VBO012<-DietSeurat(seuratObj_VBO012)

gc()

## Perform merge
seuratObj <- merge(seuratObj_VBO004, y = c(seuratObj_VBO005,seuratObj_VBO006,seuratObj_VBO007,seuratObj_VBO008,
                                           seuratObj_VBO009,seuratObj_VBO010,seuratObj_VBO011,seuratObj_VBO012),
                   project = samplename, merge.data = F) # 116588 cells

## Clean up
rm(seuratObj_VBO004)
rm(seuratObj_VBO005)
rm(seuratObj_VBO006)
rm(seuratObj_VBO007)
rm(seuratObj_VBO008)
rm(seuratObj_VBO009)
rm(seuratObj_VBO010)
rm(seuratObj_VBO011)
rm(seuratObj_VBO012)

gc()

### List samples
listLabels<-list('VBO004','VBO005','VBO006','VBO007',
                 'VBO008','VBO009','VBO010','VBO011','VBO012')

## Important Note: Skip all SCT steps!! Make object smaller
## Remove SCT slot!!
seuratObj[["SCT"]]<-NULL
DefaultAssay(seuratObj)<-"RNA"

##################################################
########## Normalize according to Seurat 
##################################################
# Basic Seurat log normalization: only cDC1s included, so no benefit to scran normalization (+faster and no sce conversion!)
seuratObj <- NormalizeData(object = seuratObj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(seuratObj))

### Get more info about HVGs (mean, dispersion and dispersion scaled)
head(HVFInfo(seuratObj))

### Plot variable features with and without labels
top10 <- head(x = VariableFeatures(object = seuratObj), 10)
plot1 <- VariableFeaturePlot(object = seuratObj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

png(file=paste0(output.dir,"Plots/RNA/05_hvg.png"), width = 850, height = 642)
CombinePlots(plots = list(plot1, plot2))
dev.off()

##### Add to diagnostics #####
diagnostics[['varGenes']]<-length(VariableFeatures(seuratObj))

###########################
########## SCALING & PCA
###########################
# Scaling the RNA assay data, not the SCT!
seuratObj <- ScaleData(seuratObj, assay = "RNA")

# Run PCA on rna normalized through scran/scater
seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj), 
                    npcs = 150, ndims.print = 1:5, nfeatures.print = 10, assay = "RNA",
                    reduction.name = "RNA_pca",reduction.key = "rnaPC_")

########################################
########## PCA PLOT
########################################
pdf(file=paste0(output.dir,"Plots/RNA/07_PCA.pdf"), width = 10)
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "RNA_pca", dims = c(1,3))
dev.off()

####################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
####################################################

# based on Pipecomp paper
int.dim <- maxLikGlobalDimEst(seuratObj@reductions$RNA_pca@cell.embeddings,
                              k=20,unbiased = TRUE, neighborhood.aggregation = 'robust')

est.PC <-round(int.dim[[1]])
est.PC
diagnostics[['est.dimsPC.maxlik']] <- est.PC

### Create PCElbowplot
png(file=paste0(output.dir,"Plots/RNA/08b_selectPC_RNA.png"), width = 850, height = 642)
ElbowPlot(object = seuratObj, ndims = 50, reduction = "RNA_pca") + geom_vline(aes(xintercept=est.PC))
dev.off()


dimsToTry<-c(est.PC,seq(20,40,by=5))

resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "RNA_pca", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj,  graph.name = "RNA_snn", resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "RNA", reduction ="RNA_pca",
                       reduction.name = "RNA_umap", reduction.key = "rnaUMAP_")
  umapPlot<-DimPlot(seuratObj, reduction = "RNA_umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "RNA_umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

################################################################################
################################################################################
### MANUAL PART
################################################################################
################################################################################

### Final
dimsToTry<-c(32)
resToUse<-0.8
diagnostics[['dimsPC']]<-dimsToTry
diagnostics[['res']]<-resToUse

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "RNA_pca", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj,  graph.name = "RNA_snn", resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "RNA", reduction ="RNA_pca",
                       reduction.name = "RNA_umap", reduction.key = "rnaUMAP_")
  umapPlot<-DimPlot(seuratObj, reduction = "RNA_umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "RNA_umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/10b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

names(seuratObj)

### Clustering: trying out clusTree
Perplexity<-32
Resolution<-0.8
Perplexity_UMAP<-32

seuratObj <- FindNeighbors(object = seuratObj, reduction = "RNA_pca", dims = 1:Perplexity)
resolutions <- seq(0,1,by=0.1)

for(res in resolutions){
  seuratObj <- FindClusters(object = seuratObj,  graph.name = "RNA_snn", resolution = res)
}

pdf(file=paste0(output.dir,"Plots/RNA/10c_Clustree.pdf"))
clustree(seuratObj, prefix = "RNA_snn_res.")
dev.off()

# Final Resolution and final clusters
res <- 0.8
diagnostics[['res']]<-res
seuratObj$RNA_clusters <- seuratObj$RNA_snn_res.0.8
Idents(seuratObj) <- seuratObj$RNA_snn_res.0.8 

################################################################################
################################################################################
### AUTOMATIC PART
################################################################################
################################################################################

umapPlot<-DimPlot(seuratObj, reduction = "RNA_umap", label = T, group.by= "RNA_clusters", label.size = 6)
seuratObj@active.assay

pdf(file=paste0(output.dir,"Plots/RNA/11_tSNE_UMAP.pdf"), width = 17*0.45, height = 12.4*0.45)
umapPlot
# tsnePlot
dev.off()

# Diagnostic plots
# a series of plots to assess how well the normalization has worked
# for this we will plot the library size and calculate the correlation with the libsize

# diagnostic plots
# first calculate fraction of zeros
seuratObj$fzero<- Matrix::colSums(seuratObj@assays$RNA@counts==0)/seuratObj@assays$RNA@counts@Dim[1]

diagplot <-function(object,reduction,metric){
  data <- Embeddings(object = object[[reduction]])
  data <- as.data.frame(x = data)
  assay <- "RNA"
  method <- "lognormalize" #"sctransform"
  log10UMI <- log10(FetchData(object,vars="sum"))[[1]]
  zero.fraction <- FetchData(object,vars="fzero")[[1]]
  PC_1 <- FetchData(object,vars="rnaPC_1")[[1]]
  if(metric=="umi"){
    if(reduction=="RNA_pca"){
      p <- ggplot() +
        geom_point(aes(x=rnaPC_1,y=rnaPC_2, colour=log10UMI), data=data, size=2, shape=20)
    }
    if(reduction=="RNA_umap"){
      p <- ggplot() +
        geom_point(aes(x=rnaUMAP_1,y=rnaUMAP_2, colour=log10UMI), data=data, size=2, shape=20)
    }
    p <- p + scale_colour_gradientn(colours = c("darkblue","darkgreen","green","yellow")) +
      ggtitle(paste0("UMI (",reduction,")\n",method)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  }
  if(metric=="zeros"){
    if(reduction=="RNA_pca"){
      p <- ggplot() +
        geom_point(aes(x=rnaPC_1,y=rnaPC_2, colour=zero.fraction), data=data, size=2, shape=20)
    }
    
    if(reduction=="RNA_umap"){
      p <- ggplot() +
        geom_point(aes(x=rnaUMAP_1,y=rnaUMAP_2, colour=zero.fraction), data=data, size=2, shape=20)
    }
    p <- p + scale_colour_gradientn(colours = c("red","blue")) +
      ggtitle(paste0("zero fraction  (",reduction,")\n",method)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  }
  if(metric=="both"){
    p <- ggplot() +
      geom_point(aes(x=log10UMI,y=rnaPC_1, colour=zero.fraction), data=data, size=2, shape=20)
    as <- cor(log10UMI,PC_1)
    p <- p + scale_colour_gradientn(colours = c("red","blue")) +
      ggtitle(paste0("UMI & zero fraction r=",round(as,2)," \n",method)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  }
  return(p +  theme_classic() )
}

pca.umi.RNA <- diagplot(seuratObj,'RNA_pca','umi')
pca.zeros.RNA <- diagplot(seuratObj,'RNA_pca','zeros')
pca.both.RNA <- diagplot(seuratObj,'RNA_pca','both')

umap.umi.RNA <- diagplot(seuratObj,'RNA_umap','umi')
umap.zeros.RNA <- diagplot(seuratObj,'RNA_umap','zeros')



pdf(file=paste0(output.dir,"Plots/RNA/12a_norm_diagnostics.pdf"), width = 17*0.45, height = 12.4*0.45)
VlnPlot(seuratObj, features =c("sum")) + ggtitle("UMI count")
VlnPlot(seuratObj, features =c("fzero")) + ggtitle("fraction of zero genes")
VlnPlot(seuratObj, features =c("detected")) + ggtitle("Detected gene count")
pca.umi.RNA
pca.zeros.RNA
pca.both.RNA
umap.umi.RNA
umap.zeros.RNA
dev.off()  

pdf(file=paste0(output.dir,"Plots/RNA/12b_mito_contamination.pdf"), width = 17*0.45, height = 12.4*0.45)
FeaturePlot(seuratObj, features = "subsets_Mito_percent")
VlnPlot(seuratObj, features = "subsets_Mito_percent")
dev.off() 

# adding the pca.drop to the metadata of the seuratobject
seuratObj$pca.drop <- metaData$pca.drop[!metaData$final.drop]
pdf(file=paste0(output.dir,"Plots/RNA/12c_PCA_outliers.pdf"), width = 17*0.45, height = 12.4*0.45)
FeaturePlot(seuratObj, features = "pca.drop")
VlnPlot(seuratObj, features = "pca.drop")
dev.off() 

pdf(file=paste0(output.dir,"Plots/RNA/13_RBC_contamination.pdf"), width = 17*0.45, height = 12.4*0.45)
FeaturePlot(seuratObj, features = "subsets_RBC_percent")
VlnPlot(seuratObj, features = "subsets_RBC_percent")
dev.off()  

pdf(file=paste0(output.dir,"Plots/RNA/14_COVID_infection.pdf"), width = 17*0.45, height = 12.4*0.45)
FeaturePlot(seuratObj, features = "subsets_COVID_percent")
VlnPlot(seuratObj, features = "subsets_COVID_percent")
dev.off()  

###########################################################################################

##### Save object
saveRDS(seuratObj, file=paste0(output.dir,"Robjects/seuratObj_",experiment,"_SCT.rds"))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

##### Read object
seuratObj <- readRDS(file=paste0(output.dir,"Robjects/seuratObj_",experiment,"_SCT.rds"))
diagnostics <- readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

################################################################################
############################## START ADT SCRIPT  ############################### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ START RNA
################################################################################

dir.create(paste0(output.dir,"Plots/ADT"))

###################################################
########## CLR transformation
###################################################
# we are doing the built-in clr normalization
# clr depends on the number of cells in a sample so possibly better to do it for each sample separately
# or try for a arcsinh(5) or logicle transform

seuratObj <- NormalizeData(seuratObj, assay = "ADT", normalization.method = "CLR", verbose=T)
Key(seuratObj)
names(seuratObj)

###################################################
########## Informative Antibodies
###################################################
# we are going to calculate the HVABs to see what ABs contain the most information
# this will create an ordering of the ABs
seuratObj <- FindVariableFeatures(seuratObj,assay="ADT",selection.method = "vst", nfeatures=nrow(seuratObj@assays$ADT)) #rawDataADT
HVABs <- VariableFeatures(seuratObj, assay = "ADT")
length(VariableFeatures(seuratObj,assay = "ADT")) # should be the number of ADTs added
head(HVFInfo(seuratObj,assay = "ADT"))

### Plot variable features with and without labels
top10 <- head(x = VariableFeatures(object = seuratObj, assay = "ADT"), 10)
plot1 <- VariableFeaturePlot(object = seuratObj, assay = "ADT")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

png(file=paste0(output.dir,"Plots/ADT/01_hvg.png"), width = 850, height = 642)
plot1 + plot2
dev.off()

##### Add to diagnostics #####
diagnostics[['varGenes']]<-length(VariableFeatures(seuratObj))

################################################################################
########## SCALING & PCA
################################################################################

# Scaling the ADT assay data
seuratObj <- ScaleData(seuratObj, assay = "ADT")

# PCA only meaningfull for large panels 
seuratObj <- RunPCA(object = seuratObj, assay = "ADT", features = VariableFeatures(seuratObj, assay = "ADT"), 
                    npcs = nrow(seuratObj@assays$ADT), ndims.print = 1:5, nfeatures.print = 10, #rawDataADT
                    reduction.name = "ADT_pca",reduction.key = "adtPC_")

names(seuratObj)

########################################
########## PCA PLOT
########################################
pdf(file=paste0(output.dir,"Plots/ADT/02_PCA.pdf"), width = 10)
DimPlot(object = seuratObj, reduction = "ADT_pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "ADT_pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "ADT_pca", dims = c(1,3))
dev.off()

########################################
########## HEATMAP OF PCs
########################################

# plotting residual variance ADT (HVG) vs PC-loadings
# build the data-frame
loadings <- rownames_to_column(as.data.frame(seuratObj@reductions$ADT_pca@feature.loadings),var="AB")
var.load <- seuratObj@assays$ADT@meta.features %>%
  rownames_to_column(var="AB") %>%
  filter(AB %in% HVABs) %>%
  # select(gene,sct.residual_variance) %>%
  left_join(loadings,by="AB")

var.load.plot <- ggplot(var.load,aes(vst.variance.standardized,adtPC_2,color=adtPC_2)) +
  geom_point() +
  geom_abline(aes(slope = 0, intercept = 0)) +
  geom_label_repel(aes(label=ifelse(abs(adtPC_2)>0.1,as.character(AB),'')),color="black") +
  theme_classic() +
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red" )
var.load.plot

# Make PC plots for inspection
pc_list <- colnames(seuratObj@reductions$ADT_pca@feature.loadings)
pc_plots <- list()
for (i in 1:dim(seuratObj@reductions$ADT_pca@feature.loadings)[2]) {
  labels <-ifelse(abs(var.load[,2+i])>0.1,as.character(var.load$AB),'')
  p <- ggplot(var.load,aes_string("vst.variance.standardized",pc_list[[i]][1],color=pc_list[[i]][1])) +
    geom_point() +
    geom_abline(aes(slope = 0, intercept = 0)) +
    geom_label_repel(aes(label=ifelse(abs(.data[[pc_list[[i]]]])>0.1,as.character(AB),'')),color="black") +
    theme_classic() +
    scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red" )
  pc_plots[[i]] <- p
}

pdf(file=paste0(output.dir,"Plots/ADT/03a_PCloadings.pdf"))
for (i in 1:dim(seuratObj@reductions$ADT_pca@feature.loadings)[2]) {
  print(pc_plots[[i]])
}
dev.off()

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

### Create PCElbowplot
png(file=paste0(output.dir,"Plots/ADT/03b_selectPC_ADT.png"), width = 850, height = 642)
ElbowPlot(object = seuratObj, ndims = dim(seuratObj@reductions$ADT_pca@feature.loadings)[2], reduction = "ADT_pca") #+ geom_vline(aes(xintercept=est.PC))
dev.off()

################################################################################
########## CLUSTER THE CELLS
################################################################################

dimsToTry<-c(seq(10,20,by=5)) 
resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "ADT_pca", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj,  graph.name = "ADT_snn", resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "ADT", reduction ="ADT_pca",
                       reduction.name = "ADT_umap", reduction.key = "adtUMAP_")
  umapPlot<-DimPlot(seuratObj, reduction = "ADT_umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "ADT_umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/ADT/04b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

################################################################################
################################################################################
### MANUAL PART
################################################################################
################################################################################

dimsToTry<-10 #est.PC
resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "ADT_pca", dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj,  graph.name = "ADT_snn", resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "ADT", reduction ="ADT_pca",
                       reduction.name = "ADT_umap", reduction.key = "adtUMAP_")
  umapPlot<-DimPlot(seuratObj, reduction = "ADT_umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "ADT_umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/ADT/04b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 20, height=10)
  
}

names(seuratObj)

### Clustering: trying out clusTree
Perplexity <- 10 #est.PC
Resolution<-0.8
resolutions <- seq(0,1,by=0.1)
for(res in resolutions){
  seuratObj <- FindClusters(object = seuratObj,  graph.name = "ADT_snn", resolution = res)
}

pdf(file=paste0(output.dir,"Plots/ADT/4c_Clustree.pdf"))
clustree(seuratObj, prefix = "ADT_snn_res.")
dev.off()

# Final Resolution and final clusters
res <- 0.8
diagnostics[['res']]<-res
seuratObj$ADT_clusters <- seuratObj$ADT_snn_res.0.8

################################################################################
################################################################################
### AUTOMATIC PART
################################################################################
################################################################################

umapPlot<-DimPlot(seuratObj, reduction = "ADT_umap", label = T, group.by= "ADT_clusters", label.size = 6)

pdf(file=paste0(output.dir,"Plots/ADT/05_tSNE_UMAP_ADT.pdf"), width = 17*0.45, height = 12.4*0.45)
umapPlot
dev.off()

###################################
#### Cross-Modality UMAPS and tSNE
##################################
umap_rna_rna <-  DimPlot(seuratObj, reduction = "RNA_umap",group.by="RNA_clusters",label=T)  + ggtitle("RNA clusters")
umap_rna_adt <-  DimPlot(seuratObj, reduction = "RNA_umap",group.by="ADT_clusters",label=T)  + ggtitle("ADT clusters")
umap_adt_adt <-  DimPlot(seuratObj, reduction = "ADT_umap",group.by="ADT_clusters",label=T)  + ggtitle("ADT clusters")
umap_adt_rna <-  DimPlot(seuratObj, reduction = "ADT_umap",group.by="RNA_clusters",label=T)  + ggtitle("RNA clusters")

pdf(file=paste0(output.dir,"Plots/ADT/06_overview_UMAP_tSNE.pdf"), width = 17, height = 12.4)
CombinePlots(plots = list(umap_rna_rna, umap_adt_rna, umap_rna_adt, umap_adt_adt), ncol = 2)
dev.off()

####################################
#### Clustering correspondence table
####################################
clustering.table <- table(seuratObj$RNA_clusters, seuratObj$ADT_clusters)
clus.heatmap <- pheatmap::pheatmap(clustering.table, scale = "column", display_numbers = clustering.table, main="RNA (rows) by ADT (columns) clusters") 

pdf(file = paste0(output.dir,"Plots/ADT/09_ADT_RNA_heatmap.pdf"), width = 11.69, height = 8.27)
clus.heatmap
dev.off()

#################################
#### Thresholded counttables
#################################

filtered.ADT<-data.frame("AB" = rownames(seuratObj@assays$ADT@counts),
                         "ADT_counts" = rowSums(seuratObj@assays$ADT@counts)) %>%
  dplyr::filter(ADT_counts>10) 

RNA.counts <- data.frame("gene" = rownames(seuratObj@assays$RNA@counts),
                         "RNA_counts" = rowSums(seuratObj@assays$RNA@counts))

write.xlsx(filtered.ADT, file = paste0(output.dir,"Plots/ADT/10_Countsummary_ADT_",experiment, ".xlsx"))
write.xlsx(RNA.counts, file = paste0(output.dir,"Plots/ADT/10_Countsummary_RNA_",experiment, ".xlsx"))

saveRDS(seuratObj, file = paste0(output.dir,"Robjects/seuratObj_",experiment,"_ADT.rds"))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment, ".rds"))

############################
dir.create(paste0(output.dir,"Plots/HTO"))

seuratObj$MULTI_ID_merge<-paste0(seuratObj$orig.ident,"_",seuratObj$MULTI_ID)
pdf(file = paste0(output.dir,"Plots/HTO/Overview_",experiment,".pdf"), width = 10, height = 12)
DimPlot(seuratObj,group.by = "MULTI_ID_merge", label = F, pt.size = 0.1) + ggtitle("MULTIseqDemux")
DimPlot(seuratObj,split.by = "MULTI_ID_merge", group.by = "MULTI_ID_merge", label = F, pt.size = 0.1, ncol = 4) + ggtitle("MULTIseqDemux")
dev.off()

################################################################################
########## RUN HARMONY
################################################################################

## Perform harmony so samples overlap better
DefaultAssay(seuratObj) <- "RNA"

########## Create vlnPlot before running Harmony ##########
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = seuratObj, reduction = "RNA_pca", pt.size = 0.2, group.by = "orig.ident")
p2 <- VlnPlot(object = seuratObj, features = "rnaPC_1", pt.size = 0.2, group.by = "orig.ident") 
plot_grid(p1,p2)
ggsave(plot_grid(p1, p2), file=paste0(output.dir,"Plots/RNA/1a_vlnPlot_beforeAlignment.png"))

########## Run Harmony ##########
### Increase theta parameter in case of bad overlap!
options(repr.plot.height = 3, repr.plot.width = 6)
seuratObj<-RunHarmony(seuratObj, reduction = "RNA_pca", assay.use = "RNA", group.by.vars = "orig.ident", theta = 2, 
                      plot_convergence = TRUE, nclust = 50, max.iter.cluster = 100, max.iter.harmony = 20, dims.use=1:40,
                      reduction.save = "RNA_harmony")

### Get embeddings
harmony_embeddings <- Embeddings(seuratObj, 'RNA_harmony')
harmony_embeddings[1:5, 1:5]

########## Create vlnPlot after running Harmony ##########
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = seuratObj, reduction = "RNA_harmony", pt.size = 0.2, group.by = "orig.ident") 
p2 <- VlnPlot(object = seuratObj, features = "RNAharmony_1", pt.size = 0.2, group.by = "orig.ident") 
plot_grid(p1,p2)
ggsave(plot_grid(p1, p2), file=paste0(output.dir,"Plots/RNA/1b_vlnPlot_afterAlignment.png"))

########################################
########## Choose dims
########################################

########## Via PCelbowplot ##########
ElbowPlot(object = seuratObj, reduction = "RNA_harmony", ndims = 40)

dimsToTry<-c(seq(15,30,by=5))

resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "RNA_harmony", dims = dimsToUse, graph.name = "RNA_harmony_snn")
  seuratObj <- FindClusters(object = seuratObj, graph.name = "RNA_harmony_snn", resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "RNA", reduction ="RNA_harmony",
                       reduction.name = "RNA_harmony_umap", reduction.key = "rnaharmonyUMAP_")
  umapPlot<-DimPlot(seuratObj, reduction = "RNA_harmony_umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "RNA_harmony_umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/16b_UMAP_",min(dimsToUse),"-",max(dimsToUse),"_RNA.png"), width = 20, height=10)
  
}

################################################################################
################################################################################
### MANUAL PART
################################################################################
################################################################################

### Final
dimsToTry<-c(20)
resToUse<-0.8
diagnostics[['dimsPC']]<-dimsToTry
diagnostics[['res']]<-resToUse

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "RNA_harmony", dims = dimsToUse, graph.name = "RNA_harmony_snn")
  seuratObj <- FindClusters(object = seuratObj, graph.name = "RNA_harmony_snn", resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "RNA", reduction ="RNA_harmony",
                       reduction.name = "RNA_harmony_umap", reduction.key = "rnaharmonyUMAP_")
  umapPlot<-DimPlot(seuratObj, reduction = "RNA_harmony_umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "RNA_harmony_umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/RNA/16b_UMAP_",min(dimsToUse),"-",max(dimsToUse),"_RNA.png"), width = 20, height=10)
  
}

names(seuratObj)

### Clustering: trying out clusTree
Perplexity<-20
Resolution<-0.8
Perplexity_UMAP<-20

seuratObj <- FindNeighbors(object = seuratObj, reduction = "RNA_harmony", dims = 1:Perplexity, graph.name = "RNA_harmony_snn")
resolutions <- seq(0,1,by=0.1)

for(res in resolutions){
  seuratObj <- FindClusters(object = seuratObj, graph.name = "RNA_harmony_snn",  resolution = res)
}

pdf(file=paste0(output.dir,"Plots/RNA/16c_Clustree.pdf"))
clustree(seuratObj, prefix = "RNA_harmony_snn_res.")
dev.off()

# Final Resolution and final clusters
res <- 0.8
diagnostics[['res']]<-res
seuratObj$harmony_clusters <- seuratObj$RNA_harmony_snn_res.0.8

####################################################################

################################################################################
########## RUN HARMONY
################################################################################

## Perform harmony son ADT too
DefaultAssay(seuratObj) <- "ADT"

########## Create vlnPlot before running Harmony ##########
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = seuratObj, reduction = "ADT_pca", pt.size = 0.2, group.by = "orig.ident")
p2 <- VlnPlot(object = seuratObj, features = "adtPC_1", pt.size = 0.2, group.by = "orig.ident") 
plot_grid(p1,p2)
ggsave(plot_grid(p1, p2), file=paste0(output.dir,"Plots/ADT/11a_vlnPlot_beforeAlignment.png"))

########## Run Harmony ##########
### Increase theta parameter in case of bad overlap!
options(repr.plot.height = 3, repr.plot.width = 6)
seuratObj<-RunHarmony(seuratObj, reduction = "ADT_pca", assay.use = "ADT", group.by.vars = "orig.ident", theta = 2, 
                      plot_convergence = TRUE, nclust = 50, max.iter.cluster = 100, max.iter.harmony = 20, dims.use=1:10,
                      reduction.save = "ADT_harmony") #Only use first 10 dims for ADT!!

### Get embeddings
harmony_embeddings_ADT <- Embeddings(seuratObj, 'ADT_harmony')
harmony_embeddings_ADT[1:5, 1:5]

########## Create vlnPlot after running Harmony ##########
options(repr.plot.height = 6, repr.plot.width = 12)
p1 <- DimPlot(object = seuratObj, reduction = "ADT_harmony", pt.size = 0.2, group.by = "orig.ident") 
p2 <- VlnPlot(object = seuratObj, features = "ADTharmony_1", pt.size = 0.2, group.by = "orig.ident") 
plot_grid(p1,p2)
ggsave(plot_grid(p1, p2), file=paste0(output.dir,"Plots/ADT/11b_vlnPlot_afterAlignment.png"))


########################################
########## Choose dims
########################################

########## Via PCelbowplot ##########
ElbowPlot(object = seuratObj, reduction = "ADT_harmony", ndims = 40)

dimsToTry<-c(seq(10,30,by=5))

resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "ADT_harmony", dims = dimsToUse, graph.name = "ADT_harmony_snn")
  seuratObj <- FindClusters(object = seuratObj, graph.name = "ADT_harmony_snn", resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "ADT", reduction ="ADT_harmony",
                       reduction.name = "ADT_harmony_umap", reduction.key = "adtharmonyUMAP_")
  umapPlot<-DimPlot(seuratObj, reduction = "ADT_harmony_umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "ADT_harmony_umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/ADT/12b_UMAP_",min(dimsToUse),"-",max(dimsToUse),"_ADT.png"), width = 20, height=10)
  
}

################################################################################
################################################################################
### MANUAL PART
################################################################################
################################################################################

### Final
dimsToTry<-c(15)
resToUse<-0.8
diagnostics[['dimsPC_ADT_harmony']]<-dimsToTry
diagnostics[['res_ADT_harmony']]<-resToUse

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, reduction = "ADT_harmony", dims = dimsToUse, graph.name = "ADT_harmony_snn")
  seuratObj <- FindClusters(object = seuratObj, graph.name = "ADT_harmony_snn", resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30, assay = "ADT", reduction ="ADT_harmony",
                       reduction.name = "ADT_harmony_umap", reduction.key = "adtharmonyUMAP_")
  umapPlot<-DimPlot(seuratObj, reduction = "ADT_harmony_umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "ADT_harmony_umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(output.dir,"Plots/ADT/16b_UMAP_",min(dimsToUse),"-",max(dimsToUse),"_ADT.png"), width = 20, height=10)
  
}

names(seuratObj)

### Clustering: trying out clusTree
Perplexity<-15
Resolution<-0.8
Perplexity_UMAP<-15

seuratObj <- FindNeighbors(object = seuratObj, reduction = "ADT_harmony", dims = 1:Perplexity, graph.name = "ADT_harmony_snn")
resolutions <- seq(0,1,by=0.1)

for(res in resolutions){
  seuratObj <- FindClusters(object = seuratObj, graph.name = "ADT_harmony_snn",  resolution = res)
}

pdf(file=paste0(output.dir,"Plots/ADT/16c_Clustree.pdf"))
clustree(seuratObj, prefix = "ADT_harmony_snn_res.")
dev.off()

# Final Resolution and final clusters
res <- 0.8
diagnostics[['res_ADT_harmony']]<-res
seuratObj$harmony_clusters_ADT <- seuratObj$ADT_harmony_snn_res.0.8

####################################################################

## Switch back
DefaultAssay(seuratObj) <- "RNA"

####################################################################
dir.create(paste0(output.dir,"Annotation/"))

## Check annotation spleen atlas
D4<-DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "predicted.id", label = TRUE,
            label.size = 2, repel = TRUE)  + ggtitle("Query transferred labels")

pdf(file = paste0(output.dir,"Annotation/Dimplot_spleen_query_predicted_annotation.pdf"),
    height = 10, width = 20)
D4
dev.off()

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')
source('~/VIB/DATA/Roos/Daan 1/script_featurePlots.R')

# DimPlot(seuratObj, reduction = "SCT_umap", label = T, repel = T, label.size = 4, group.by = "seurat_clusters")
U1<-DimPlot(seuratObj, reduction = "RNA_harmony_umap", label = T, repel = T, label.size = 4, group.by = "harmony_clusters") +
  labs(title = "Harmony clusters on Harmony UMAP")
U2<-DimPlot(seuratObj, reduction = "ADT_harmony_umap", label = T, repel = T, label.size = 4, group.by = "harmony_clusters_ADT") +
  labs(title = "ADT clusters on ADT UMAP")

#ADT clustering on SCT UMAP and vice versa
U3<-DimPlot(seuratObj, reduction = "RNA_harmony_umap", label = T, repel = T, label.size = 4, group.by = "harmony_clusters_ADT") +
  labs(title = "ADT clusters on Harmony UMAP")
U4<-DimPlot(seuratObj, reduction = "ADT_harmony_umap", label = T, repel = T, label.size = 4, group.by = "harmony_clusters") +
  labs(title = "Harmony clusters on ADT UMAP")

################################################################################
########## ANNOTATION
################################################################################

pdf(file=paste0("results/Annotation/1_annotation_",experiment,".pdf"), width = 15)
grid.arrange(U1, U3, ncol=2)
grid.arrange(U2, U4, ncol=2)
dev.off()


################################################################################
########## COLOR CELLS ACCORDING TO ADT on SCT + vice versa
################################################################################
DimPlot(seuratObj, reduction = "RNA_harmony_umap", label = F, group.by="harmony_clusters")
DimPlot(seuratObj, reduction = "RNA_harmony_umap", label = F, group.by="harmony_clusters_ADT")

clusterMatrix<-seuratObj@meta.data
umapTable<-as.data.frame(seuratObj[['RNA_harmony_umap']]@cell.embeddings, stringsAsFactors = F)

#ADT clusters
pdf(file=paste0("results/Annotation/1_annotation_Color_ADT_clusters_on_UMAP_",experiment,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$harmony_clusters_ADT))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$harmony_clusters_ADT==i),])))
  C1<-C1+ggtitle(paste0("ADT_cluster_",i))
  print(C1)
}
dev.off()

#SCT clusters
pdf(file=paste0("results/Annotation/1_annotation_Color_Harmony_clusters_on_UMAP_",experiment,".pdf"), width = 15)
for (i in 0:(length(levels(seuratObj@meta.data$harmony_clusters))-1)) {
  C1<-colorSomeCells(clusterMatrix, umapTable, WhichCells(seuratObj, cells = rownames(seuratObj@meta.data[which(seuratObj@meta.data$harmony_clusters==i),])))
  C1<-C1+ggtitle(paste0("RNA_cluster_",i))
  print(C1)
}
dev.off()

####################################################################################

### Find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(seuratObj) <- "ADT"
Idents(seuratObj)<-seuratObj$harmony_clusters

allADTMarkers <- FindAllMarkers(seuratObj, min.pct = 0.10,logfc.threshold = 0.30, return.thresh = 0.01, only.pos = TRUE)
table(allADTMarkers$cluster)
saveRDS(allADTMarkers, file=paste0(output.dir,"Robjects/allADTMarkers_",experiment,".rds"))

### Add to diagnostics
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))
diagnostics[['ADTmarkersPerCluster']]<-paste0(table(allADTMarkers$cluster)," markers for cluster ",rownames(table(allADTMarkers$cluster)))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

### Create list with markers
totalNrClusters<-max(as.numeric(names(table(allADTMarkers$cluster))))
totalNrClustersPlusOne<-totalNrClusters+1
ADTmarkersList<-list()

for(i in 1:totalNrClustersPlusOne){
  clusterNr<-i-1
  
  tmp<-allADTMarkers[allADTMarkers$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  ADTmarkersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(ADTmarkersList)<-paste0("cluster",0:totalNrClusters)

### Write to Excel

write.xlsx(ADTmarkersList, file =paste0(output.dir, "/allClusters_ADT_markers_",experiment,".xlsx"))

####################################################################################

### Find markers for every cluster compared to all remaining cells, report only the positive ones
Idents(seuratObj)<-seuratObj$harmony_clusters
DefaultAssay(seuratObj) <- "RNA"

allMarkers <- FindAllMarkers(seuratObj, min.pct = 0.10,logfc.threshold = 0.30, return.thresh = 0.01, only.pos = TRUE)
table(allMarkers$cluster)
saveRDS(allMarkers, file=paste0(output.dir,"Robjects/allRNAMarkers_",experiment,".rds"))

### Add to diagnostics
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))
diagnostics[['RNAmarkersPerCluster']]<-paste0(table(allMarkers$cluster)," markers for cluster ",rownames(table(allMarkers$cluster)))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

### Create list with markers
totalNrClusters<-max(as.numeric(names(table(allMarkers$cluster))))
totalNrClustersPlusOne<-totalNrClusters+1
markersList<-list()

for(i in 1:totalNrClustersPlusOne){
  clusterNr<-i-1
  
  tmp<-allMarkers[allMarkers$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  markersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(markersList)<-paste0("cluster",0:totalNrClusters)

### Write to Excel

write.xlsx(markersList, file =paste0(output.dir, "/allClusters_RNA_markers_",experiment,".xlsx"))

###############################################

U0<-DimPlot(seuratObj, reduction = "RNA_harmony_umap", label = T, repel = T, label.size = 4, group.by = "Original_annotation")
ggsave(U0, file=paste0(output.dir,"Annotation/2_Original_annotation_UMAP_",experiment,".png"), width =20, height= 12, dpi = "retina")

## Check clusters 
Test1<-FindMarkers(seuratObj, ident.1 = 12, ident.2 = 14)

### Create annotated UMAP
seuratObj@meta.data$annotated_clusters <- seuratObj$harmony_clusters
seuratObj@meta.data$annotated_clusters<- factor(seuratObj@meta.data$annotated_clusters,levels(seuratObj@meta.data$annotated_clusters)[c(1,2,11:18,3:10)]) #reorder levels
levels(seuratObj@meta.data$annotated_clusters) <- c('Il12b+ cDC1s 1',"Tnfrsf4 hi late mature cDC1s","Cxcl9-11 hi cDC1s",
                                                    "Ms4a4c+ late mature cDC1s","Late immature cDC1s","Tmem176+ late mature cDC1s",
                                                    "Proliferating cDC1s 1","Proliferating cDC1s 2","Hamp+ late mature cDC1s",
                                                    "Ccl5 late mature cDC1s","Serpina3+ cDC1s","Other cDC1s","H2-M2+ late mature cDC1s",
                                                    "IFN+ cDC1s","Zmynd15+ late mature cDC1s","Doublets","Small cluster 1","Small cluster 2")

U1 <- DimPlot(seuratObj, reduction = "RNA_harmony_umap", label = T, repel = T, label.size = 4, group.by = "annotated_clusters")
ggsave(U1, file=paste0(output.dir,"Annotation/2_Annotated_UMAP_",experiment,".png"), width =20, height= 12, dpi = "retina")

U2 <- DimPlot(seuratObj, reduction = "RNA_harmony_umap", label = T, repel = T, split.by = "orig.ident", ncol = 3, label.size = 2, group.by = "annotated_clusters")
ggsave(U2, file=paste0(output.dir,"Annotation/4_Annotated_UMAP_split_",experiment,".png"), width =20, height= 15, dpi = "retina")

U3 <- DimPlot(seuratObj, reduction = "RNA_umap", label = F, repel = T, split.by = "orig.ident", ncol = 3, label.size = 4, group.by = "annotated_clusters")
ggsave(U3, file=paste0(output.dir,"Annotation/4_Annotated_UMAP_split_non_harmony_",experiment,".png"), width =20, height= 15, dpi = "retina")

Idents(seuratObj)<-seuratObj$Condition

U4 <- DimPlot(seuratObj, reduction = "RNA_umap", label = F, repel = T, order = rev(c("Steady state","eLNPs 2h","eLNPs 8h","pIC-LNPs 2h","pIC-LNPs 8h",
                                                                                     "CpG-LNPs 2h","CpG-LNPs 8h","pIC alone 2h","pIC alone 8h")),
              ncol = 3, label.size = 4, cols = brewer.pal(10, "Paired")[2:10])
ggsave(U4, file=paste0(output.dir,"Annotation/4_Annotated_UMAP_split_non_harmony_v2_",experiment,".png"), width =20, height= 15, dpi = "retina")

pdf(file=paste0(output.dir,"Annotation/4_Annotated_UMAP_split_non_harmony_v2_",experiment,".pdf"), width =20, height= 15)
U4
dev.off()

seuratObj$Condition_v3<-factor(seuratObj$Condition,levels = c("Steady state","eLNPs 2h","eLNPs 8h","pIC-LNPs 2h","pIC-LNPs 8h",
                                                              "CpG-LNPs 2h","CpG-LNPs 8h","pIC alone 2h","pIC alone 8h"))
Idents(seuratObj)<-seuratObj$Condition_v3

U5 <- DimPlot(seuratObj, reduction = "RNA_umap", label = F, repel = T, split.by = "Condition_v3", ncol = 3, label.size = 4,
              cols = brewer.pal(10, "Paired")[2:10])
ggsave(U5, file=paste0(output.dir,"Annotation/4_Annotated_UMAP_split_non_harmony_v3_",experiment,".png"), width =20, height= 15, dpi = "retina")

pdf(file=paste0(output.dir,"Annotation/4_Annotated_UMAP_split_non_harmony_v3_",experiment,".pdf"), width =20, height= 15)
U5
dev.off()

##############################################

## Annotated markers
Idents(seuratObj)<-seuratObj@meta.data$annotated_clusters

names(markersList)<-levels(Idents(seuratObj))

### Write to Excel
write.xlsx(markersList, file =paste0(output.dir, "allClusters_annotated_",experiment,".xlsx"), overwrite = T)

################################################################################

## Sample distribution across clusters
# create a dataset
Sample <- seuratObj@meta.data$orig.ident
cluster <- seuratObj@active.ident

data <- data.frame(table(Sample, cluster))

# Stacked
png(file=paste0(output.dir,"Annotation/3_SampleDistribution_ggplot2_annot_1_",experiment,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(output.dir,"Annotation/3_SampleDistribution_ggplot2_annot_2_",experiment,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
dev.off()

##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
data_split1<-data[which(data$Sample=="VBO004"),]
data_split2<-data[which(data$Sample=="VBO005"),]
data_split3<-data[which(data$Sample=="VBO006"),]
data_split4<-data[which(data$Sample=="VBO007"),]
data_split5<-data[which(data$Sample=="VBO008"),]
data_split6<-data[which(data$Sample=="VBO009"),]
data_split7<-data[which(data$Sample=="VBO010"),]
data_split8<-data[which(data$Sample=="VBO011"),]
data_split9<-data[which(data$Sample=="VBO012"),]
data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-round((x/sum(data_split1$Freq))*100,2)})
data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-round((x/sum(data_split2$Freq))*100,2)})
data_split3$Freq<-lapply(data_split3$Freq,function(x){x<-round((x/sum(data_split3$Freq))*100,2)})
data_split4$Freq<-lapply(data_split4$Freq,function(x){x<-round((x/sum(data_split4$Freq))*100,2)})
data_split5$Freq<-lapply(data_split5$Freq,function(x){x<-round((x/sum(data_split5$Freq))*100,2)})
data_split6$Freq<-lapply(data_split6$Freq,function(x){x<-round((x/sum(data_split6$Freq))*100,2)})
data_split7$Freq<-lapply(data_split7$Freq,function(x){x<-round((x/sum(data_split7$Freq))*100,2)})
data_split8$Freq<-lapply(data_split8$Freq,function(x){x<-round((x/sum(data_split8$Freq))*100,2)})
data_split9$Freq<-lapply(data_split9$Freq,function(x){x<-round((x/sum(data_split9$Freq))*100,2)})

data_new<-rbind(data_split1,data_split2,data_split3,data_split4,
                data_split5,data_split6,data_split7,data_split8,
                data_split9)

png(file=paste0(output.dir,"Annotation/3_SampleDistribution_ggplot2_annot_1_",experiment,"_adjusted.png"), width = 2000, height = 1500, res = 300)
ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

################################################################################

# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using seuratObj[['weighted.nn']]
# The WNN graph can be accessed at seuratObj[["wknn"]], 
# and the SNN graph used for clustering at seuratObj[["wsnn"]]
# Cell-specific modality weights can be accessed at seuratObj$RNA.weight
seuratObj <- FindMultiModalNeighbors(
  seuratObj, reduction.list = list("RNA_harmony", "ADT_harmony"),  #Use harmony for both RNA and ADT!!
  dims.list = list(1:20, 1:15), modality.weight.name = "RNA.weight",
  k.nn = 30,  knn.range = 500, prune.SNN = 1/20 #Too many isolated clusters https://github.com/satijalab/seurat/issues/4793
)

# We can now use these results for downstream analysis, such as visualization and clustering. 
# For example, we can create a UMAP visualization of the data based on a weighted combination of RNA and protein data. 
# We can also perform graph-based clustering and visualize these results on the UMAP, alongside a set of cell annotations.

seuratObj <- RunUMAP(seuratObj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seuratObj <- FindClusters(seuratObj, graph.name = "wsnn", algorithm = 3, resolution = 0.8, future.seed = TRUE, verbose = FALSE) #Avoid warning set future.seed = T
seuratObj <- FindClusters(seuratObj, graph.name = "wsnn", algorithm = 3, resolution = 0.5, future.seed = TRUE, verbose = FALSE) #Avoid warning set future.seed = T
seuratObj <- FindClusters(seuratObj, graph.name = "wsnn", algorithm = 3, resolution = 0.3, future.seed = TRUE, verbose = FALSE) #Avoid warning set future.seed = T


p1 <- DimPlot(seuratObj, reduction = 'wnn.umap', group.by = 'wsnn_res.0.8', label = TRUE, repel = TRUE, label.size = 4) + NoLegend() + ggtitle("Multimodal UMAP with new numbered clusters")
p2 <- DimPlot(seuratObj, reduction = 'wnn.umap', group.by = 'annotated_clusters', label = TRUE, repel = TRUE, label.size = 4) + NoLegend() + ggtitle("Multimodal UMAP with harmony annotation")
p3 <- DimPlot(seuratObj, reduction = 'RNA_harmony_umap', group.by = 'annotated_clusters', label = TRUE, repel = TRUE, label.size = 4) + NoLegend() + ggtitle("Harmony UMAP with harmony annotation")
p4 <- DimPlot(seuratObj, reduction = 'ADT_harmony_umap', group.by = 'annotated_clusters', label = TRUE, repel = TRUE, label.size = 4) + NoLegend() + ggtitle("ADT UMAP with harmony annotation")

pdf(file = paste0(output.dir,"Annotation/Multimodal_plot_UMAP_",experiment,"_annotation.pdf"),
    height = 20, width = 30)
p1+p2+p3+p4
dev.off()

# Finally, we can visualize the modality weights that were learned for each cell. 
# Each of the populations with the highest RNA weights represent progenitor cells, 
# while the populations with the highest protein weights represent T cells. 
# This is in line with our biological expectations, 
# as the antibody panel does not contain markers that can distinguish between different progenitor populations.
pdf(file = paste0(output.dir,"Annotation/Multimodal_vlnplot_RNA_weights_",experiment,"_full_annotation.pdf"),
    height = 10, width = 15)
VlnPlot(seuratObj, features = "RNA.weight", group.by = 'annotated_clusters', sort = TRUE, pt.size = 0.1) +
  NoLegend() + ggtitle("RNA weights annotation clusters")
VlnPlot(seuratObj, features = "ADT.weight", group.by = 'annotated_clusters', sort = TRUE, pt.size = 0.1) +
  NoLegend() + ggtitle("ADT weights annotation clusters")
dev.off()

### Find ADT markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(seuratObj)<-"ADT"
Idents(seuratObj)<-seuratObj$wsnn_res.0.8

allMarkers_ADT <- FindAllMarkers(seuratObj, min.pct = 0.10,logfc.threshold = 0.30, return.thresh = 0.01, only.pos = TRUE)
table(allMarkers_ADT$cluster)
saveRDS(allMarkers_ADT, file=paste0(output.dir,"Robjects/allADTMarkers_Multimodal_",experiment,".rds"))

### Add to diagnostics
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))
diagnostics[['ADTmarkersPerMultimodalCluster']]<-paste0(table(allMarkers_ADT$cluster)," markers for Multimodal cluster ",rownames(table(allMarkers_ADT$cluster)))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

### Create list with markers
totalNrClusters_ADT<-max(as.numeric(names(table(allMarkers_ADT$cluster))))
totalNrClustersADTPlusOne<-totalNrClusters_ADT+1
ADTmarkersList<-list()

for(i in 1:totalNrClustersADTPlusOne){
  clusterNr<-i-1
  
  tmp<-allMarkers_ADT[allMarkers_ADT$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  ADTmarkersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(ADTmarkersList)<-paste0("cluster",0:totalNrClusters_ADT)

### Write to Excel
write.xlsx(ADTmarkersList, file =paste0(output.dir, "/allMultimodalClusters_ADT_markers_",experiment,".xlsx"))

### Find RNA markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(seuratObj)<-"RNA"
Idents(seuratObj)<-seuratObj$wsnn_res.0.8

allMarkers <- FindAllMarkers(seuratObj, min.pct = 0.10,logfc.threshold = 0.30, return.thresh = 0.01, only.pos = TRUE)
table(allMarkers$cluster)
saveRDS(allMarkers, file=paste0(output.dir,"Robjects/allRNAMarkers_Multimodal_",experiment,".rds"))

### Add to diagnostics
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))
diagnostics[['RNAmarkersPerMultimodalCluster']]<-paste0(table(allMarkers$cluster)," markers for Multimodal cluster ",rownames(table(allMarkers$cluster)))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

### Create list with markers
totalNrClusters<-max(as.numeric(names(table(allMarkers$cluster))))
totalNrClustersPlusOne<-totalNrClusters+1
markersList<-list()

for(i in 1:totalNrClustersPlusOne){
  clusterNr<-i-1
  
  tmp<-allMarkers[allMarkers$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  markersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(markersList)<-paste0("cluster",0:totalNrClusters)

### Write to Excel
write.xlsx(markersList, file =paste0(output.dir, "/allMultimodalClusters_RNA_markers_",experiment,".xlsx"))

##############################################

### Create annotated UMAP
seuratObj@meta.data$annotated_clusters_multimodal <- seuratObj$wsnn_res.0.8
seuratObj@meta.data$annotated_clusters_multimodal<- factor(seuratObj@meta.data$annotated_clusters_multimodal,levels(seuratObj@meta.data$annotated_clusters_multimodal)[c(1,2,8:15,3:7)]) #reorder levels
levels(seuratObj@meta.data$annotated_clusters_multimodal) <- c('Il12b+ cDC1s 1',"Tmem176+ late mature cDC1s 1","Tmem176- late mature cDC1s",
                                                               "CD43+ cDC1s","Late immature cDC1s","Other cDC1s 1",
                                                               "Tmem176+ late mature cDC1s 2","Proliferating cDC1s 1","Proliferating cDC1s 2",
                                                               "Other cDC1s 2","Proliferating cDC1s 3","Hamp+ late mature cDC1s",
                                                               'Il12b+ cDC1s 2',"Serpina3+ cDC1s","Tmem176+ late mature cDC1s 3")

U1 <- DimPlot(seuratObj, reduction = "wnn.umap", label = T, repel = T, label.size = 4, group.by = "annotated_clusters_multimodal")
ggsave(U1, file=paste0(output.dir,"Annotation/2_Annotated_UMAP_multimodal_",experiment,".png"), width =20, height= 12, dpi = "retina")

##############################################

## Annotated markers
Idents(seuratObj)<-seuratObj$annotated_clusters_multimodal

names(markersList)<-levels(Idents(seuratObj))

### Write to Excel
write.xlsx(markersList, file =paste0(output.dir, "allMultimodalClusters_annotated_",experiment,".xlsx"), overwrite = T)

##############################################

### Find ADT markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(seuratObj)<-"ADT"
Idents(seuratObj)<-seuratObj$wsnn_res.0.3

allMarkers_ADT <- FindAllMarkers(seuratObj, min.pct = 0.10,logfc.threshold = 0.30, return.thresh = 0.01, only.pos = TRUE)
table(allMarkers_ADT$cluster)
saveRDS(allMarkers_ADT, file=paste0(output.dir,"Robjects/allADTMarkers_Multimodal_simple_",experiment,".rds"))

### Add to diagnostics
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))
diagnostics[['ADTmarkersPerMultimodalClusterSimple']]<-paste0(table(allMarkers_ADT$cluster)," markers for Multimodal cluster ",rownames(table(allMarkers_ADT$cluster)))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

### Create list with markers
totalNrClusters_ADT<-max(as.numeric(names(table(allMarkers_ADT$cluster))))
totalNrClustersADTPlusOne<-totalNrClusters_ADT+1
ADTmarkersList<-list()

for(i in 1:totalNrClustersADTPlusOne){
  clusterNr<-i-1
  
  tmp<-allMarkers_ADT[allMarkers_ADT$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  ADTmarkersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(ADTmarkersList)<-paste0("cluster",0:totalNrClusters_ADT)

### Write to Excel
write.xlsx(ADTmarkersList, file =paste0(output.dir, "/allMultimodalClusters_ADT_markers_simple_",experiment,".xlsx"))

### Find RNA markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(seuratObj)<-"RNA"
Idents(seuratObj)<-seuratObj$wsnn_res.0.3

allMarkers <- FindAllMarkers(seuratObj, min.pct = 0.10,logfc.threshold = 0.30, return.thresh = 0.01, only.pos = TRUE)
table(allMarkers$cluster)
saveRDS(allMarkers, file=paste0(output.dir,"Robjects/allRNAMarkers_Multimodal_simple_",experiment,".rds"))

### Add to diagnostics
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))
diagnostics[['RNAmarkersPerMultimodalClusterSimple']]<-paste0(table(allMarkers$cluster)," markers for Multimodal cluster ",rownames(table(allMarkers$cluster)))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

### Create list with markers
totalNrClusters<-max(as.numeric(names(table(allMarkers$cluster))))
totalNrClustersPlusOne<-totalNrClusters+1
markersList<-list()

for(i in 1:totalNrClustersPlusOne){
  clusterNr<-i-1
  
  tmp<-allMarkers[allMarkers$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  markersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(markersList)<-paste0("cluster",0:totalNrClusters)

### Write to Excel
write.xlsx(markersList, file =paste0(output.dir, "/allMultimodalClusters_RNA_markers_simple_",experiment,".xlsx"))

### Create simple annotated UMAP
seuratObj@meta.data$annotated_clusters_multimodal_simple <- seuratObj$wsnn_res.0.3
levels(seuratObj@meta.data$annotated_clusters_multimodal_simple) <- c("Late mature cDC1s","Immature cDC1s",'Il12b+ cDC1s',
                                                                      "Intermediate cDC1s","Proliferating cDC1s 2","Proliferating cDC1s 1",
                                                                      "Other cDC1s","Hamp+ late mature cDC1s")

U1 <- DimPlot(seuratObj, reduction = "wnn.umap", label = T, repel = T, label.size = 6, group.by = "annotated_clusters_multimodal_simple", raster = F)
ggsave(U1, file=paste0(output.dir,"Annotation/2_Annotated_UMAP_multimodal_simple_",experiment,".png"), width =20, height= 12, dpi = "retina")

pdf(file=paste0(output.dir,"Annotation/2_Annotated_UMAP_multimodal_simple_",experiment,".pdf"), width =20, height= 12)
print(U1)
dev.off()

##############################################

## Annotated markers
Idents(seuratObj)<-seuratObj$annotated_clusters_multimodal_simple

names(markersList)<-levels(Idents(seuratObj))

### Write to Excel
write.xlsx(markersList, file =paste0(output.dir, "allMultimodalClusters_annotated_simple_",experiment,".xlsx"), overwrite = T)

################################################################################

## Sample distribution across clusters

# create a dataset
Sample <- seuratObj@meta.data$orig.ident
cluster <- seuratObj@active.ident

data <- data.frame(table(Sample, cluster))

# Stacked
png(file=paste0(output.dir,"Annotation/7_SampleDistribution_ggplot2_MMannot_1_",experiment,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

png(file=paste0(output.dir,"Annotation/7_SampleDistribution_ggplot2_MMannot_2_",experiment,".png"), width = 2000, height = 1500, res = 300)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
dev.off()

##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
data_split1<-data[which(data$Sample=="VBO004"),]
data_split2<-data[which(data$Sample=="VBO005"),]
data_split3<-data[which(data$Sample=="VBO006"),]
data_split4<-data[which(data$Sample=="VBO007"),]
data_split5<-data[which(data$Sample=="VBO008"),]
data_split6<-data[which(data$Sample=="VBO009"),]
data_split7<-data[which(data$Sample=="VBO010"),]
data_split8<-data[which(data$Sample=="VBO011"),]
data_split9<-data[which(data$Sample=="VBO012"),]
data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-round((x/sum(data_split1$Freq))*100,2)})
data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-round((x/sum(data_split2$Freq))*100,2)})
data_split3$Freq<-lapply(data_split3$Freq,function(x){x<-round((x/sum(data_split3$Freq))*100,2)})
data_split4$Freq<-lapply(data_split4$Freq,function(x){x<-round((x/sum(data_split4$Freq))*100,2)})
data_split5$Freq<-lapply(data_split5$Freq,function(x){x<-round((x/sum(data_split5$Freq))*100,2)})
data_split6$Freq<-lapply(data_split6$Freq,function(x){x<-round((x/sum(data_split6$Freq))*100,2)})
data_split7$Freq<-lapply(data_split7$Freq,function(x){x<-round((x/sum(data_split7$Freq))*100,2)})
data_split8$Freq<-lapply(data_split8$Freq,function(x){x<-round((x/sum(data_split8$Freq))*100,2)})
data_split9$Freq<-lapply(data_split9$Freq,function(x){x<-round((x/sum(data_split9$Freq))*100,2)})

data_new<-rbind(data_split1,data_split2,data_split3,data_split4,
                data_split5,data_split6,data_split7,data_split8,
                data_split9)

png(file=paste0(output.dir,"Annotation/7_SampleDistribution_ggplot2_MMannot_1_",experiment,"_adjusted.png"), width = 2000, height = 1500, res = 300)
ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + #scale_fill_manual(values=c('Blue','Red')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

#############################

## Add conditions to object!
seuratObj$Condition<-as.factor(seuratObj$orig.ident)
levels(seuratObj$Condition)<-c("Steady state","eLNPs 8h","pIC-LNPs 8h","CpG-LNPs 8h","pIC alone 8h",
                               "eLNPs 2h","pIC-LNPs 2h","CpG-LNPs 2h","pIC alone 2h")

#############################
dir.create(paste0(output.dir,"Feature_plots/"))

## Featureplots of interesting marker genes and proteins
features<-c("CD63","CD274","CD278","CD43","CD117","CD62L","ESAM","CD103","CD207-mh",
            "CD197","CD80","CD86","XCR1")

features<-c("Cd63","Ccr7","Fscn1","Sell","Itgae","Esam","Cd207","Cd80","Cd86","Xcr1",
            "Cd69","Mki67","Tmem176a","Tmem176b","Il12b","Tnfrsf4","Cxcl9","Cxcl11",
            "Ms4a4c","Egr3","Hamp","Ccl5","Serpina3f","Serpina3g","H2-M2","Ifna2","Zmynd15")

features<-c("Fscn1","Cd80","CD80","Ccl22","CD63","Cd63","CD197","Ccr7",
            "Icosl","Scube3","Slco5a1","Tmem176b",
            "Ifit2","Isg20","Oas3","Ifi47",
            "Dnase1l3","Ch25h","Nlrp3","Ido1","Ern1","Josd1",
            "Srebf1","Srebf2","Nr1h2","Nr1h3",
            "Abca1","Abcg1","Apol7c","Apoe")

features<-c("Havcr2","Bag6","Lgals9")

features<-c("Nlrp3","Casp4","Casp8")

features<-c("Smpdl3a")

V2 <- VlnPlot(seuratObj, features = "Smpdl3a", split.by = "Condition", group.by = "annotated_clusters_Muscat") 

pdf(file=paste0(output.dir,"Feature_plots/Violin_plot_split_Smpdl3a_",samplename,".pdf"), height = 10, width = 20)
V2
dev.off()

SCT_markers_WT<-read.xlsx("../../../RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Marker_lists/SCTmarkersList_SCTclus_SAM2and3_WT_subset_final_2021.xlsx",
                          sheet = "Early mature cDC1s")
features<-c(SCT_markers_WT$gene[c(1,4,5,6,7,8,9,19,20,21,22)],"Ifit1","Ifit2","Ifit3")

pdf(file=paste0(output.dir,"Feature_plots/Feature_plot_markers_Sophie_Caro_ISG_cluster_RNA_harmony_UMAP_",samplename,"_blue_grey.pdf"), height = 7, width = 10)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features =feature, cols = c("grey", "blue"), 
                  reduction = "RNA_harmony_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T) 
  print(F1)
}
dev.off()

pdf(file=paste0(output.dir,"Feature_plots/Feature_plot_markers_Sophie_Caro_ISG_cluster_split_RNA_harmony_UMAP_",samplename,"_blue_grey.pdf"), height = 5, width = 40)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features =feature, cols = c("grey", "blue"), 
                  reduction = "RNA_harmony_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T,
                  split.by = "Condition", raster = T) & theme(legend.position = "right") 
  print(F1)
}
dev.off()

ISG_signature <- list(c(SCT_markers_WT$gene[c(1,4,5,6,7,8,9,19,20,21,22)],"Ifit1","Ifit2","Ifit3"))

seuratObj <- AddModuleScore(object = seuratObj, features = ISG_signature, name = "ISG_signature_score")

## Simplify annotation
seuratObj$annotated_clusters_Muscat<-seuratObj$annotated_clusters
levels(seuratObj$annotated_clusters_Muscat)<-c('Early mature cDC1s',"Late mature cDC1s","Early mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Late mature cDC1s",
                                               "Proliferating cDC1s","Proliferating cDC1s","Late mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Other cDC1s","Late mature cDC1s",
                                               "Early mature cDC1s","Late mature cDC1s","Doublets","Late mature cDC1s","Early mature cDC1s")

DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat")
Idents(seuratObj)<-seuratObj$annotated_clusters_Muscat

F1 <- FeaturePlot(object = seuratObj, features = "ISG_signature_score1", label = T, repel =T, order = T, cols = c("Yellow","Red"),
                  reduction = "RNA_harmony_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, raster = T)  #+ scale_color_viridis(option = "C") #wnn.umap
F2<-FeaturePlot(object = seuratObj, features = "ISG_signature_score1", label = F, repel =T, order = T, cols = c("Yellow","Red"),
                reduction = "RNA_harmony_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, #wnn.umap
                split.by = "Condition", raster = T) & theme(legend.position = "right") 

pdf(file = paste0(output.dir,"Feature_plots/Feature_plot_markers_Sophie_Caro_ISG_modulescore_RNA_harmony_UMAP_",samplename,".pdf"), width = 7, height = 10)
F1
dev.off()

pdf(file = paste0(output.dir,"Feature_plots/Feature_plot_markers_Sophie_Caro_ISG_modulescore_split_RNA_harmony_UMAP_",samplename,".pdf"), width = 40, height = 5)
F2
dev.off()

###########################################################################

## Simplified muscat annotation
pdf(file=paste0(output.dir,"Annotation/8_Annotated_UMAP_RNA_harmony_detailed_",experiment,".pdf"), width =20, height= 12)
DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters")
dev.off()

pdf(file=paste0(output.dir,"Annotation/8_Annotated_UMAP_RNA_harmony_muscat_annotation_",experiment,".pdf"), width =19, height= 12)
DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat")
dev.off()

pdf(file=paste0(output.dir,"Annotation/8_Annotated_UMAP_RNA_harmony_muscat_annotation_split_",experiment,".pdf"), width =20, height= 20)
DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat", split.by = "Condition", ncol = 3)
dev.off()

## Repeat sample distribution across clusters

# create a dataset
Sample <- seuratObj@meta.data$Condition
cluster <- seuratObj$annotated_clusters_Muscat

data <- data.frame(table(Sample, cluster))

# Stacked
pdf(paste0(output.dir,"Annotation/8_SampleDistribution_ggplot2_RNA_harmony_muscat_annot_1_",experiment,".pdf"), width = 20, height = 15)
ggplot(data, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + scale_fill_manual(values=brewer.pal(10, "Paired")[2:10]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

pdf(paste0(output.dir,"Annotation/8_SampleDistribution_ggplot2_RNA_harmony_muscat_annot_2_",experiment,".pdf"), width = 20, height = 15)
ggplot(data, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
dev.off()

##Controlled for amount of cells per sample (split by sample and divide counts by total cell count sample)
data_split1<-data[which(data$Sample=="Steady state"),]
data_split2<-data[which(data$Sample=="eLNPs 8h"),]
data_split3<-data[which(data$Sample=="pIC-LNPs 8h"),]
data_split4<-data[which(data$Sample=="CpG-LNPs 8h"),]
data_split5<-data[which(data$Sample=="pIC alone 8h"),]
data_split6<-data[which(data$Sample=="eLNPs 2h"),]
data_split7<-data[which(data$Sample=="pIC-LNPs 2h"),]
data_split8<-data[which(data$Sample=="CpG-LNPs 2h"),]
data_split9<-data[which(data$Sample=="pIC alone 2h"),]
data_split1$Freq<-lapply(data_split1$Freq,function(x){x<-round((x/sum(data_split1$Freq))*100,2)})
data_split2$Freq<-lapply(data_split2$Freq,function(x){x<-round((x/sum(data_split2$Freq))*100,2)})
data_split3$Freq<-lapply(data_split3$Freq,function(x){x<-round((x/sum(data_split3$Freq))*100,2)})
data_split4$Freq<-lapply(data_split4$Freq,function(x){x<-round((x/sum(data_split4$Freq))*100,2)})
data_split5$Freq<-lapply(data_split5$Freq,function(x){x<-round((x/sum(data_split5$Freq))*100,2)})
data_split6$Freq<-lapply(data_split6$Freq,function(x){x<-round((x/sum(data_split6$Freq))*100,2)})
data_split7$Freq<-lapply(data_split7$Freq,function(x){x<-round((x/sum(data_split7$Freq))*100,2)})
data_split8$Freq<-lapply(data_split8$Freq,function(x){x<-round((x/sum(data_split8$Freq))*100,2)})
data_split9$Freq<-lapply(data_split9$Freq,function(x){x<-round((x/sum(data_split9$Freq))*100,2)})

data_new<-rbind(data_split1,data_split2,data_split3,data_split4,
                data_split5,data_split6,data_split7,data_split8,
                data_split9)

pdf(paste0(output.dir,"Annotation/8_SampleDistribution_ggplot2_RNA_harmony_muscat_annot_1_",experiment,"_adjusted.pdf"), width = 20, height = 15)
ggplot(data_new, aes(fill=Sample, y=Freq, x=cluster)) + theme_bw() + scale_fill_manual(values=brewer.pal(10, "Paired")[2:10]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6)) +
  geom_bar(position="fill", stat="identity", colour="white")
dev.off()

pdf(paste0(output.dir,"Annotation/8_SampleDistribution_ggplot2_RNA_harmony_muscat_annot_2_",experiment,"_adjusted.pdf"), width = 20, height = 15)
ggplot(data_new, aes(fill=cluster, y=Freq, x=Sample)) + theme_bw() +
  geom_bar(position="fill", stat="identity", colour="white")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.6))
dev.off()

###########################################################################

## Simplify clustering
seuratObj$annotated_clusters_Muscat<-seuratObj$annotated_clusters
levels(seuratObj$annotated_clusters_Muscat)<-c('Early mature cDC1s',"Late mature cDC1s","Early mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Late mature cDC1s",
                                               "Proliferating cDC1s","Proliferating cDC1s","Late mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Other cDC1s","Late mature cDC1s",
                                               "Early mature cDC1s","Late mature cDC1s","Doublets","Late mature cDC1s","Early mature cDC1s")

DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat")

seuratObj$Condition_v2<-as.factor(as.character(seuratObj$Condition))
seuratObj$Condition_v2<-gsub("-","_",seuratObj$Condition_v2)
seuratObj$Condition_v2<-gsub(" ","_",seuratObj$Condition_v2)

### Create new clusters: split on source
seuratObj@meta.data$newClustersTmp<-seuratObj$annotated_clusters_Muscat
seuratObj@meta.data$newClusters<-paste0(seuratObj@meta.data$newClustersTmp,"_",seuratObj@meta.data$Condition_v2)
head(seuratObj@meta.data)

### Use the new clusters
seuratObj@meta.data$newClusters<- as.factor(seuratObj@meta.data$newClusters) #reorder levels
Idents(seuratObj)<-seuratObj@meta.data$newClusters

DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "newClusters")

#################################################################
########## GET DIFF ADT MARKERS between conditions ##########
#################################################################

getDE_ADT<-function(ident1, ident2){
  markersDiff <- FindMarkers(seuratObj, ident.1 = ident1, ident.2 = ident2, assay = "ADT",
                             min.pct = 0.10) 
  markersDiff<-markersDiff[markersDiff$p_val_adj < 0.01,]
  markersDiff<-markersDiff[order(markersDiff$avg_log2FC, decreasing = T),]
  
  markersDiff$geneSymbol<-rownames(markersDiff)
  markersDiff$pct.1<-markersDiff$pct.1+0.001
  markersDiff$pct.2<-markersDiff$pct.2+0.001
  
  markersDiff<-rbind(markersDiff[markersDiff$avg_log2FC > 0,] %>% dplyr::mutate(.,score=pct.1/pct.2*avg_log2FC),
                     markersDiff[markersDiff$avg_log2FC < 0,] %>% dplyr::mutate(.,score=pct.2/pct.1*avg_log2FC))
  markersDiff<-markersDiff[order(markersDiff$score, decreasing = T),]
  return(markersDiff)
}

#### Get diff markers #####
DE_ADT_eLNP_vs_WT_EM_2h<-getDE_ADT("Early mature cDC1s_eLNPs_2h","Early mature cDC1s_Steady_state")
DE_ADT_eLNP_vs_WT_EM_2h<-DE_ADT_eLNP_vs_WT_EM_2h[order(DE_ADT_eLNP_vs_WT_EM_2h$avg_log2FC,decreasing = T),]
head(DE_ADT_eLNP_vs_WT_EM_2h)
dim(DE_ADT_eLNP_vs_WT_EM_2h)

DE_ADT_eLNP_vs_WT_LM_8h<-getDE_ADT("Late mature cDC1s_eLNPs_8h","Late mature cDC1s_Steady_state")
DE_ADT_eLNP_vs_WT_LM_8h<-DE_ADT_eLNP_vs_WT_LM_8h[order(DE_ADT_eLNP_vs_WT_LM_8h$avg_log2FC,decreasing = T),]
head(DE_ADT_eLNP_vs_WT_LM_8h)
dim(DE_ADT_eLNP_vs_WT_LM_8h)

DE_ADT_CpG_LNP_vs_WT_EM_2h<-getDE_ADT("Early mature cDC1s_CpG_LNPs_2h","Early mature cDC1s_Steady_state")
DE_ADT_CpG_LNP_vs_WT_EM_2h<-DE_ADT_CpG_LNP_vs_WT_EM_2h[order(DE_ADT_CpG_LNP_vs_WT_EM_2h$avg_log2FC,decreasing = T),]
head(DE_ADT_CpG_LNP_vs_WT_EM_2h)
dim(DE_ADT_CpG_LNP_vs_WT_EM_2h)

DE_ADT_CpG_LNP_vs_WT_LM_8h<-getDE_ADT("Late mature cDC1s_CpG_LNPs_8h","Late mature cDC1s_Steady_state")
DE_ADT_CpG_LNP_vs_WT_LM_8h<-DE_ADT_CpG_LNP_vs_WT_LM_8h[order(DE_ADT_CpG_LNP_vs_WT_LM_8h$avg_log2FC,decreasing = T),]
head(DE_ADT_CpG_LNP_vs_WT_LM_8h)
dim(DE_ADT_CpG_LNP_vs_WT_LM_8h)

DE_ADT_pIC_alone_vs_WT_EM_2h<-getDE_ADT("Early mature cDC1s_pIC_alone_2h","Early mature cDC1s_Steady_state")
DE_ADT_pIC_alone_vs_WT_EM_2h<-DE_ADT_pIC_alone_vs_WT_EM_2h[order(DE_ADT_pIC_alone_vs_WT_EM_2h$avg_log2FC,decreasing = T),]
head(DE_ADT_pIC_alone_vs_WT_EM_2h)
dim(DE_ADT_pIC_alone_vs_WT_EM_2h)

DE_ADT_pIC_alone_vs_WT_LM_8h<-getDE_ADT("Late mature cDC1s_pIC_alone_8h","Late mature cDC1s_Steady_state")
DE_ADT_pIC_alone_vs_WT_LM_8h<-DE_ADT_pIC_alone_vs_WT_LM_8h[order(DE_ADT_pIC_alone_vs_WT_LM_8h$avg_log2FC,decreasing = T),]
head(DE_ADT_pIC_alone_vs_WT_LM_8h)
dim(DE_ADT_pIC_alone_vs_WT_LM_8h)

DE_ADT_pIC_LNP_vs_WT_EM_2h<-getDE_ADT("Early mature cDC1s_pIC_LNPs_2h","Early mature cDC1s_Steady_state")
DE_ADT_pIC_LNP_vs_WT_EM_2h<-DE_ADT_pIC_LNP_vs_WT_EM_2h[order(DE_ADT_pIC_LNP_vs_WT_EM_2h$avg_log2FC,decreasing = T),]
head(DE_ADT_pIC_LNP_vs_WT_EM_2h)
dim(DE_ADT_pIC_LNP_vs_WT_EM_2h)

DE_ADT_pIC_LNP_vs_WT_LM_8h<-getDE_ADT("Late mature cDC1s_pIC_LNPs_8h","Late mature cDC1s_Steady_state")
DE_ADT_pIC_LNP_vs_WT_LM_8h<-DE_ADT_pIC_LNP_vs_WT_LM_8h[order(DE_ADT_pIC_LNP_vs_WT_LM_8h$avg_log2FC,decreasing = T),]
head(DE_ADT_pIC_LNP_vs_WT_LM_8h)
dim(DE_ADT_pIC_LNP_vs_WT_LM_8h)

##add to list
listDiffMarkers_ADT<-tibble::lst(DE_ADT_eLNP_vs_WT_EM_2h,DE_ADT_eLNP_vs_WT_LM_8h,DE_ADT_CpG_LNP_vs_WT_EM_2h,
                                 DE_ADT_CpG_LNP_vs_WT_LM_8h,DE_ADT_pIC_alone_vs_WT_EM_2h,DE_ADT_pIC_alone_vs_WT_LM_8h,
                                 DE_ADT_pIC_LNP_vs_WT_EM_2h,DE_ADT_pIC_LNP_vs_WT_LM_8h)

lapply(listDiffMarkers_ADT, dim)
listDiffMarkers_ADT<-lapply(listDiffMarkers_ADT,function(x){x<-x[order(x$score, decreasing=T),]})

###############

#Save RDS 
saveRDS(listDiffMarkers_ADT,file=paste0(output.dir,"Robjects/markers_ADT_DiffSamples_Full_",samplename,".rds"))

### Write to Excel
write.xlsx(listDiffMarkers_ADT, file = paste0(output.dir,"summaryDiffMarkers_ADT_Full_",samplename,".xlsx"))


#### Get diff markers v2 (maturation stage) #####
DE_ADT_eLNP_LM_8h_vs_EM_2h<-getDE_ADT("Late mature cDC1s_eLNPs_8h","Early mature cDC1s_eLNPs_2h")
DE_ADT_eLNP_LM_8h_vs_EM_2h<-DE_ADT_eLNP_LM_8h_vs_EM_2h[order(DE_ADT_eLNP_LM_8h_vs_EM_2h$avg_log2FC,decreasing = T),]
head(DE_ADT_eLNP_LM_8h_vs_EM_2h)
dim(DE_ADT_eLNP_LM_8h_vs_EM_2h)

DE_ADT_CpG_LNP_LM_8h_vs_EM_2h<-getDE_ADT("Late mature cDC1s_CpG_LNPs_8h","Early mature cDC1s_CpG_LNPs_2h")
DE_ADT_CpG_LNP_LM_8h_vs_EM_2h<-DE_ADT_CpG_LNP_LM_8h_vs_EM_2h[order(DE_ADT_CpG_LNP_LM_8h_vs_EM_2h$avg_log2FC,decreasing = T),]
head(DE_ADT_CpG_LNP_LM_8h_vs_EM_2h)
dim(DE_ADT_CpG_LNP_LM_8h_vs_EM_2h)

DE_ADT_WT_LM_8h_vs_EM_2h<-getDE_ADT("Late mature cDC1s_Steady_state","Early mature cDC1s_Steady_state")
DE_ADT_WT_LM_8h_vs_EM_2h<-DE_ADT_WT_LM_8h_vs_EM_2h[order(DE_ADT_WT_LM_8h_vs_EM_2h$avg_log2FC,decreasing = T),]
head(DE_ADT_WT_LM_8h_vs_EM_2h)
dim(DE_ADT_WT_LM_8h_vs_EM_2h)

DE_ADT_pIC_alone_LM_8h_vs_EM_2h<-getDE_ADT("Late mature cDC1s_pIC_alone_8h","Early mature cDC1s_pIC_alone_2h")
DE_ADT_pIC_alone_LM_8h_vs_EM_2h<-DE_ADT_pIC_alone_LM_8h_vs_EM_2h[order(DE_ADT_pIC_alone_LM_8h_vs_EM_2h$avg_log2FC,decreasing = T),]
head(DE_ADT_pIC_alone_LM_8h_vs_EM_2h)
dim(DE_ADT_pIC_alone_LM_8h_vs_EM_2h)

DE_ADT_pIC_LNP_LM_8h_vs_EM_2h<-getDE_ADT("Late mature cDC1s_pIC_LNPs_8h","Early mature cDC1s_pIC_LNPs_2h")
DE_ADT_pIC_LNP_LM_8h_vs_EM_2h<-DE_ADT_pIC_LNP_LM_8h_vs_EM_2h[order(DE_ADT_pIC_LNP_LM_8h_vs_EM_2h$avg_log2FC,decreasing = T),]
head(DE_ADT_pIC_LNP_LM_8h_vs_EM_2h)
dim(DE_ADT_pIC_LNP_LM_8h_vs_EM_2h)

##add to list
listDiffMarkers_ADT_maturation<-tibble::lst(DE_ADT_eLNP_LM_8h_vs_EM_2h,DE_ADT_CpG_LNP_LM_8h_vs_EM_2h,DE_ADT_WT_LM_8h_vs_EM_2h,
                                 DE_ADT_pIC_alone_LM_8h_vs_EM_2h,DE_ADT_pIC_LNP_LM_8h_vs_EM_2h)

lapply(listDiffMarkers_ADT_maturation, dim)
listDiffMarkers_ADT_maturation<-lapply(listDiffMarkers_ADT_maturation,function(x){x<-x[order(x$score, decreasing=T),]})

###############

#Save RDS 
saveRDS(listDiffMarkers_ADT_maturation,file=paste0(output.dir,"Robjects/markers_ADT_Diff_maturation_Full_",samplename,".rds"))

### Write to Excel
write.xlsx(listDiffMarkers_ADT_maturation, file = paste0(output.dir,"summaryDiffMarkers_maturation_ADT_Full_",samplename,".xlsx"))

###########################################################################

## Check for Malissen visualization
## Simplify annotation
seuratObj$annotated_clusters_Muscat<-seuratObj$annotated_clusters
levels(seuratObj$annotated_clusters_Muscat)<-c('Early mature cDC1s',"Late mature cDC1s","Early mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Late mature cDC1s",
                                               "Proliferating cDC1s","Proliferating cDC1s","Late mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Other cDC1s","Late mature cDC1s",
                                               "Early mature cDC1s","Late mature cDC1s","Doublets","Late mature cDC1s","Early mature cDC1s")

seuratObj$annotated_clusters_Muscat_v2<-paste0(seuratObj$orig.ident,"_",seuratObj$annotated_clusters_Muscat)

DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2")

## Subset to only keep pops to compare
Idents(seuratObj)<-seuratObj$annotated_clusters_Muscat_v2
seuratObj_subset<-subset(seuratObj, idents = levels(Idents(seuratObj))[c(1,2,4,6,11,17,22,29,36,37,42)]) #Also included 1 (Imm SS)

DimPlot(seuratObj_subset, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2", split.by = "orig.ident", ncol = 3)

seuratObj_subset$Condition_v2<-seuratObj_subset$Condition
levels(seuratObj_subset$Condition_v2)<-c("Steady state","eLNPs","pIC-LNPs","CpG-LNPs","pIC alone",
                                  "eLNPs","pIC-LNPs","CpG-LNPs","pIC alone")
seuratObj_subset$Condition_v2<-gsub("-","_",seuratObj_subset$Condition_v2)
seuratObj_subset$Condition_v2<-gsub(" ","_",seuratObj_subset$Condition_v2)

## Combine MULTI_ID with Orig.ident or Genotype?? Need to put genotype in correct order: otherwise WT after DKO
seuratObj_subset@meta.data[["MULTI_ID_merge"]]
seuratObj_subset@meta.data$MULTI_ID_merge<-as.factor(seuratObj_subset@meta.data$MULTI_ID_merge)

## Update muscat annotation
seuratObj_subset$annotated_clusters_Muscat<-factor(as.character(seuratObj_subset$annotated_clusters_Muscat), 
                                                   levels = c("Immature cDC1s", "Early mature cDC1s", "Late mature cDC1s")) #Adapted for Imm
levels(seuratObj_subset$annotated_clusters_Muscat)<-c("Immature","Early_mature","Late_mature")

## Try new sample_ID
seuratObj_subset$sample_ID<-paste0(as.character(seuratObj_subset$MULTI_ID),"_",as.character(seuratObj_subset$Condition_v2),
                            "_",as.character(seuratObj_subset$annotated_clusters_Muscat))
seuratObj_subset$sample_ID<-factor(seuratObj_subset$sample_ID, #Added Imm replicates
                                   levels = c("Hashtag1_Steady_state_Immature","Hashtag2_Steady_state_Immature","Hashtag3_Steady_state_Immature","Hashtag4_Steady_state_Immature",
                                              "Hashtag1_Steady_state_Early_mature","Hashtag2_Steady_state_Early_mature","Hashtag3_Steady_state_Early_mature","Hashtag4_Steady_state_Early_mature",
                                              "Hashtag1_Steady_state_Late_mature","Hashtag2_Steady_state_Late_mature","Hashtag3_Steady_state_Late_mature","Hashtag4_Steady_state_Late_mature",
                                              "Hashtag1_eLNPs_Early_mature","Hashtag2_eLNPs_Early_mature","Hashtag3_eLNPs_Early_mature","Hashtag4_eLNPs_Early_mature",       
                                              "Hashtag1_eLNPs_Late_mature","Hashtag2_eLNPs_Late_mature","Hashtag3_eLNPs_Late_mature","Hashtag4_eLNPs_Late_mature",
                                              "Hashtag1_pIC_LNPs_Early_mature","Hashtag2_pIC_LNPs_Early_mature","Hashtag3_pIC_LNPs_Early_mature","Hashtag4_pIC_LNPs_Early_mature",
                                              "Hashtag1_pIC_LNPs_Late_mature","Hashtag2_pIC_LNPs_Late_mature","Hashtag3_pIC_LNPs_Late_mature","Hashtag4_pIC_LNPs_Late_mature",
                                              "Hashtag1_CpG_LNPs_Early_mature","Hashtag2_CpG_LNPs_Early_mature","Hashtag3_CpG_LNPs_Early_mature","Hashtag4_CpG_LNPs_Early_mature",
                                              "Hashtag1_CpG_LNPs_Late_mature","Hashtag2_CpG_LNPs_Late_mature","Hashtag3_CpG_LNPs_Late_mature",
                                              "Hashtag1_pIC_alone_Early_mature","Hashtag2_pIC_alone_Early_mature","Hashtag3_pIC_alone_Early_mature","Hashtag4_pIC_alone_Early_mature",
                                              "Hashtag1_pIC_alone_Late_mature","Hashtag2_pIC_alone_Late_mature","Hashtag3_pIC_alone_Late_mature","Hashtag4_pIC_alone_Late_mature"))
levels(seuratObj_subset$sample_ID)

## Update clusters
seuratObj_subset$annotated_clusters_Muscat_v2<-factor(paste0(seuratObj_subset$Condition_v2,"_",seuratObj_subset$annotated_clusters_Muscat),
                                                      levels = c("Steady_state_Immature","Steady_state_Early_mature","Steady_state_Late_mature",
                                                                 "eLNPs_Early_mature","eLNPs_Late_mature","pIC_LNPs_Early_mature","pIC_LNPs_Late_mature",
                                                                 "CpG_LNPs_Early_mature","CpG_LNPs_Late_mature" ,"pIC_alone_Early_mature","pIC_alone_Late_mature"))
DimPlot(seuratObj_subset, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2",split.by = "Condition_v2", ncol = 3)

## Basic annotation
seuratObj_subset$Basic_annotation<-"cDC1s"


## Look at average expression per sample
Idents(seuratObj_subset)<-seuratObj_subset$sample_ID
seuratObj_subset_average <- AverageExpression(seuratObj_subset, return.seurat = T)

## Look at modulescores for various genelists
tbl_fil_CpG_2h<-readRDS(file=paste0("results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_2h-WT_VBO_4-12.rds"))
tbl_fil_CpG_8h<-readRDS(file=paste0("results/Robjects/tbl_fil_muscat_DESeq2_new_CpG_LNPs_8h-WT_VBO_4-12.rds"))
tbl_fil_eLNP_2h<-readRDS(file=paste0("results/Robjects/tbl_fil_muscat_DESeq2_new_eLNPs_2h-WT_VBO_4-12.rds"))
tbl_fil_eLNP_8h<-readRDS(file=paste0("results/Robjects/tbl_fil_muscat_DESeq2_new_eLNPs_8h-WT_VBO_4-12.rds"))
tbl_fil_pIC_LNP_2h<-readRDS(file=paste0("results/Robjects/tbl_fil_muscat_DESeq2_new_pIC_LNPs_2h-WT_VBO_4-12.rds"))
tbl_fil_pIC_LNP_8h<-readRDS(file=paste0("results/Robjects/tbl_fil_muscat_DESeq2_new_pIC_LNPs_8h-WT_VBO_4-12.rds"))
tbl_fil_pIC_2h<-readRDS(file=paste0("results/Robjects/tbl_fil_muscat_DESeq2_new_pIC_2h-WT_VBO_4-12.rds"))
tbl_fil_pIC_8h<-readRDS(file=paste0("results/Robjects/tbl_fil_muscat_DESeq2_new_pIC_8h-WT_VBO_4-12.rds"))

Genelist_CpG_2h<-list(head(tbl_fil_CpG_2h[["Early mature cDC1s"]][which(tbl_fil_CpG_2h$`Early mature cDC1s`$logFC > 0),"gene"],25))
Genelist_CpG_8h<-list(head(tbl_fil_CpG_8h[["Late mature cDC1s"]][which(tbl_fil_CpG_8h$`Late mature cDC1s`$logFC > 0),"gene"],25))
Genelist_eLNP_2h<-list(head(tbl_fil_eLNP_2h[["Early mature cDC1s"]][which(tbl_fil_eLNP_2h$`Early mature cDC1s`$logFC > 0),"gene"],25))
Genelist_eLNP_8h<-list(head(tbl_fil_eLNP_8h[["Late mature cDC1s"]][which(tbl_fil_eLNP_8h$`Late mature cDC1s`$logFC > 0),"gene"],25))
Genelist_pIC_LNP_2h<-list(head(tbl_fil_pIC_LNP_2h[["Early mature cDC1s"]][which(tbl_fil_pIC_LNP_2h$`Early mature cDC1s`$logFC > 0),"gene"],25))
Genelist_pIC_LNP_8h<-list(head(tbl_fil_pIC_LNP_8h[["Late mature cDC1s"]][which(tbl_fil_pIC_LNP_8h$`Late mature cDC1s`$logFC > 0),"gene"],25))
Genelist_pIC_2h<-list(head(tbl_fil_pIC_2h[["Early mature cDC1s"]][which(tbl_fil_pIC_2h$`Early mature cDC1s`$logFC > 0),"gene"],25))
Genelist_pIC_8h<-list(head(tbl_fil_pIC_8h[["Late mature cDC1s"]][which(tbl_fil_pIC_8h$`Late mature cDC1s`$logFC > 0),"gene"],25))
Genelist_WT_2h<-list(head(tbl_fil_eLNP_2h[["Early mature cDC1s"]][which(tbl_fil_eLNP_2h$`Early mature cDC1s`$logFC < 0),"gene"],25))
Genelist_WT_8h<-list(head(tbl_fil_eLNP_8h[["Late mature cDC1s"]][which(tbl_fil_eLNP_8h$`Late mature cDC1s`$logFC < 0),"gene"],25))

seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = Genelist_CpG_2h, name = "CpG_2h_score")
seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = Genelist_CpG_8h, name = "CpG_8h_score")
seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = Genelist_eLNP_2h, name = "eLNP_2h_score")
seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = Genelist_eLNP_8h, name = "eLNP_8h_score")
seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = Genelist_pIC_LNP_2h, name = "pIC_LNP_2h_score")
seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = Genelist_pIC_LNP_8h, name = "pIC_LNP_8h_score")
seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = Genelist_pIC_2h, name = "pIC_2h_score")
seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = Genelist_pIC_8h, name = "pIC_8h_score")
seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = Genelist_WT_2h, name = "WT_2h_score")
seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = Genelist_WT_8h, name = "WT_8h_score")

Groups<-c(rep("SS_Imm",4),rep("SS_EM",4),rep("SS_LM",4),
          rep("eLNPs_EM",4),rep("eLNPs_LM",4),
          rep("pIC_LNPs_EM",4),rep("pIC_LNPs_LM",4),
          rep("CpG_LNPs_EM",4),rep("CpG_LNPs_LM",3),
          rep("pIC_EM",4),rep("pIC_LM",4))
Groups_v2<-seq(1:43)

Test<-tibble::lst(as.data.frame(cbind(seuratObj_subset_average$CpG_2h_score1,Groups,Groups_v2)),
                  as.data.frame(cbind(seuratObj_subset_average$CpG_8h_score1,Groups,Groups_v2)),
                  as.data.frame(cbind(seuratObj_subset_average$eLNP_2h_score1,Groups,Groups_v2)),
                  as.data.frame(cbind(seuratObj_subset_average$eLNP_8h_score1,Groups,Groups_v2)),
                  as.data.frame(cbind(seuratObj_subset_average$pIC_LNP_2h_score1,Groups,Groups_v2)),
                  as.data.frame(cbind(seuratObj_subset_average$pIC_LNP_8h_score1,Groups,Groups_v2)),
                  as.data.frame(cbind(seuratObj_subset_average$pIC_2h_score1,Groups,Groups_v2)),
                  as.data.frame(cbind(seuratObj_subset_average$pIC_8h_score1,Groups,Groups_v2)),
                  as.data.frame(cbind(seuratObj_subset_average$WT_2h_score1,Groups,Groups_v2)),
                  as.data.frame(cbind(seuratObj_subset_average$WT_8h_score1,Groups,Groups_v2)))
names(Test)<-c("CpG_2h_score","CpG_8h_score","eLNP_2h_score","eLNP_8h_score","pIC_LNP_2h_score",
               "pIC_LNP_8h_score","pIC_2h_score","pIC_8h_score","WT_2h_score","WT_8h_score")
for (j in 1:length(Test)){
  Test[[j]][["Groups_v2"]]<-factor(Test[[j]][["Groups_v2"]], levels = seq(1:43))
  Test[[j]][["Groups"]]<-factor(Test[[j]][["Groups"]], levels = c("SS_Imm","SS_EM","SS_LM","eLNPs_EM","eLNPs_LM",
                                                                  "pIC_LNPs_EM","pIC_LNPs_LM","CpG_LNPs_EM","CpG_LNPs_LM",
                                                                  "pIC_EM","pIC_LM"))
  Test[[j]][["V1"]]<-as.numeric(Test[[j]][["V1"]])
}


dir.create(paste0(output.dir,"Genelist_expression_plots/"))

for (j in 1:length(Test)){
  p1<- ggplot(Test[[j]], aes(x=as.factor(Groups), y = V1, fill=as.factor(Groups) )) + 
    geom_bar(stat = "summary", fun = mean) +
    scale_fill_manual(values =c("cadetblue1",brewer.pal(10, "Paired"))) +
    theme(legend.position="none") +
    # theme(axis.text.x=element_blank(),
    #       axis.ticks.x=element_blank()) +
    scale_y_continuous(name = "Relative mRNA expression",
                       limits=c(-0.6,1.6),
                       breaks=c(-0.5,0,0.5,1,1.5)) +
    # scale_shape_manual(values=c(seq(0,9))) +
    labs(title=paste0(names(Test)[j]," genelist expression across various conditions"),
         x ="Samples", y = "Relative mRNA expression")
  
  p2 <- ggplot(Test[[j]], aes(x = as.factor(Groups_v2), y =V1)) +
    geom_point(aes(colour = factor(Groups), shape = factor(Groups)), size = 5) +
    # facet_wrap(~ Groups) +
    theme_light()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_y_continuous(name = "Relative mRNA expression",
                       limits=c(-0.6,1.6),
                       breaks=c(-0.5,0,0.5,1,1.5)) +
    scale_shape_manual(values=c(seq(0,10))) +
    # scale_size_manual(values = rep(10,43)) +
    scale_color_manual(values = c("cadetblue1",brewer.pal(10, "Paired"))) +
    labs(title=paste0(names(Test)[j]," genelist expression various conditions"),
         x ="Samples", y = "Relative mRNA expression",
         colour = "Cell populations", shape = "Cell populations")
  
  pdf(file=paste0(output.dir,"Genelist_expression_plots/Expression_plot_",names(Test)[j],"_average.pdf"), height = 10, width = 15)
  print(p1)
  dev.off()
  
  pdf(file=paste0(output.dir,"Genelist_expression_plots/Expression_plot_",names(Test)[j],"_individual.pdf"), height = 10, width = 15)
  print(p2)
  dev.off()
}

#########################################################################

## Extra 05/2023: meta analysis gene signatures
## Simplify annotation
seuratObj$annotated_clusters_Muscat<-seuratObj$annotated_clusters
levels(seuratObj$annotated_clusters_Muscat)<-c('Early mature cDC1s',"Late mature cDC1s","Early mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Late mature cDC1s",
                                               "Proliferating cDC1s","Proliferating cDC1s","Late mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Other cDC1s","Late mature cDC1s",
                                               "Early mature cDC1s","Late mature cDC1s","Doublets","Late mature cDC1s","Early mature cDC1s")

seuratObj$annotated_clusters_Muscat_v2<-paste0(seuratObj$orig.ident,"_",seuratObj$annotated_clusters_Muscat)

DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2")

## Subset to only keep pops to compare
Idents(seuratObj)<-seuratObj$annotated_clusters_Muscat_v2
seuratObj_subset<-subset(seuratObj, idents = levels(Idents(seuratObj))[c(1,2,4,6,11,17,22,29,36,37,42)]) #Also included 1 (Imm SS)

DimPlot(seuratObj_subset, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2", split.by = "orig.ident", ncol = 3)

seuratObj_subset$Condition_v2<-seuratObj_subset$Condition
levels(seuratObj_subset$Condition_v2)<-c("Steady state","eLNPs","pIC-LNPs","CpG-LNPs","pIC alone",
                                         "eLNPs","pIC-LNPs","CpG-LNPs","pIC alone")
seuratObj_subset$Condition_v2<-gsub("-","_",seuratObj_subset$Condition_v2)
seuratObj_subset$Condition_v2<-gsub(" ","_",seuratObj_subset$Condition_v2)

## Combine MULTI_ID with Orig.ident or Genotype?? Need to put genotype in correct order: otherwise WT after DKO
seuratObj_subset@meta.data[["MULTI_ID_merge"]]
seuratObj_subset@meta.data$MULTI_ID_merge<-as.factor(seuratObj_subset@meta.data$MULTI_ID_merge)

## Update muscat annotation
seuratObj_subset$annotated_clusters_Muscat<-factor(as.character(seuratObj_subset$annotated_clusters_Muscat), 
                                                   levels = c("Immature cDC1s", "Early mature cDC1s", "Late mature cDC1s")) #Adapted for Imm
levels(seuratObj_subset$annotated_clusters_Muscat)<-c("Immature","Early_mature","Late_mature")

## Try new sample_ID
seuratObj_subset$sample_ID<-paste0(as.character(seuratObj_subset$MULTI_ID),"_",as.character(seuratObj_subset$Condition_v2),
                                   "_",as.character(seuratObj_subset$annotated_clusters_Muscat))
seuratObj_subset$sample_ID<-factor(seuratObj_subset$sample_ID, #Added Imm replicates
                                   levels = c("Hashtag1_Steady_state_Immature","Hashtag2_Steady_state_Immature","Hashtag3_Steady_state_Immature","Hashtag4_Steady_state_Immature",
                                              "Hashtag1_Steady_state_Early_mature","Hashtag2_Steady_state_Early_mature","Hashtag3_Steady_state_Early_mature","Hashtag4_Steady_state_Early_mature",
                                              "Hashtag1_Steady_state_Late_mature","Hashtag2_Steady_state_Late_mature","Hashtag3_Steady_state_Late_mature","Hashtag4_Steady_state_Late_mature",
                                              "Hashtag1_eLNPs_Early_mature","Hashtag2_eLNPs_Early_mature","Hashtag3_eLNPs_Early_mature","Hashtag4_eLNPs_Early_mature",       
                                              "Hashtag1_eLNPs_Late_mature","Hashtag2_eLNPs_Late_mature","Hashtag3_eLNPs_Late_mature","Hashtag4_eLNPs_Late_mature",
                                              "Hashtag1_pIC_LNPs_Early_mature","Hashtag2_pIC_LNPs_Early_mature","Hashtag3_pIC_LNPs_Early_mature","Hashtag4_pIC_LNPs_Early_mature",
                                              "Hashtag1_pIC_LNPs_Late_mature","Hashtag2_pIC_LNPs_Late_mature","Hashtag3_pIC_LNPs_Late_mature","Hashtag4_pIC_LNPs_Late_mature",
                                              "Hashtag1_CpG_LNPs_Early_mature","Hashtag2_CpG_LNPs_Early_mature","Hashtag3_CpG_LNPs_Early_mature","Hashtag4_CpG_LNPs_Early_mature",
                                              "Hashtag1_CpG_LNPs_Late_mature","Hashtag2_CpG_LNPs_Late_mature","Hashtag3_CpG_LNPs_Late_mature",
                                              "Hashtag1_pIC_alone_Early_mature","Hashtag2_pIC_alone_Early_mature","Hashtag3_pIC_alone_Early_mature","Hashtag4_pIC_alone_Early_mature",
                                              "Hashtag1_pIC_alone_Late_mature","Hashtag2_pIC_alone_Late_mature","Hashtag3_pIC_alone_Late_mature","Hashtag4_pIC_alone_Late_mature"))
levels(seuratObj_subset$sample_ID)

## Update clusters
seuratObj_subset$annotated_clusters_Muscat_v2<-factor(paste0(seuratObj_subset$Condition_v2,"_",seuratObj_subset$annotated_clusters_Muscat),
                                                      levels = c("Steady_state_Immature","Steady_state_Early_mature","Steady_state_Late_mature",
                                                                 "eLNPs_Early_mature","eLNPs_Late_mature","pIC_LNPs_Early_mature","pIC_LNPs_Late_mature",
                                                                 "CpG_LNPs_Early_mature","CpG_LNPs_Late_mature" ,"pIC_alone_Early_mature","pIC_alone_Late_mature"))
DimPlot(seuratObj_subset, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2",split.by = "Condition_v2", ncol = 3)

## Basic annotation
seuratObj_subset$Basic_annotation<-"cDC1s"

## Look at average expression per sample
Idents(seuratObj_subset)<-seuratObj_subset$sample_ID
seuratObj_subset_average <- AverageExpression(seuratObj_subset, return.seurat = T)

## Gene lists
Genelist_tibble<-readRDS("../../../Meta_analysis/Tibble_all_genelists.rds")

for (i in 1:length(Genelist_tibble)){
  genelist<-Genelist_tibble[[i]]
  name<-names(Genelist_tibble[i])
  
  seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = list(intersect(genelist,rownames(seuratObj_subset_average))), name = name)
}

Groups<-c(rep("SS_Imm",4),rep("SS_EM",4),rep("SS_LM",4),
          rep("eLNPs_EM",4),rep("eLNPs_LM",4),
          rep("pIC_LNPs_EM",4),rep("pIC_LNPs_LM",4),
          rep("CpG_LNPs_EM",4),rep("CpG_LNPs_LM",3),
          rep("pIC_EM",4),rep("pIC_LM",4))
Groups_v2<-seq(1:43)

Test<-tibble::lst()
counter = 1
for (j in c(8:34)){
  Test[[counter]]<-as.data.frame(cbind(seuratObj_subset_average@meta.data[,j],Groups,Groups_v2), row.names = rownames(seuratObj_subset_average@meta.data))
  counter=counter+1
}

names(Test)<-names(Genelist_tibble)

for (j in 1:length(Test)){
  Test[[j]][["Groups_v2"]]<-factor(Test[[j]][["Groups_v2"]], levels = seq(1:43))
  Test[[j]][["Groups"]]<-factor(Test[[j]][["Groups"]], levels = c("SS_Imm","SS_EM","SS_LM","eLNPs_EM","eLNPs_LM",
                                                                  "pIC_LNPs_EM","pIC_LNPs_LM","CpG_LNPs_EM","CpG_LNPs_LM",
                                                                  "pIC_EM","pIC_LM"))
  Test[[j]][["V1"]]<-as.numeric(Test[[j]][["V1"]])
}


dir.create(paste0(output.dir,"Genelist_expression_plots/"))

for (j in 1:length(Test)){
  p1<- ggplot(Test[[j]], aes(x=as.factor(Groups), y = V1, fill=as.factor(Groups) )) + 
    geom_bar(stat = "summary", fun = mean) +
    scale_fill_manual(values =c("cadetblue1",brewer.pal(10, "Paired"))) +
    theme(legend.position="none") +
    # theme(axis.text.x=element_blank(),
    #       axis.ticks.x=element_blank()) +
    scale_y_continuous(name = "Relative mRNA expression",
                       limits=c(-0.6,1.6),
                       breaks=c(-0.5,0,0.5,1,1.5)) +
    # scale_shape_manual(values=c(seq(0,9))) +
    labs(title=paste0(names(Test)[j]," genelist expression across various conditions"),
         x ="Samples", y = "Relative mRNA expression")
  
  p2 <- ggplot(Test[[j]], aes(x = as.factor(Groups_v2), y =V1)) +
    geom_point(aes(colour = factor(Groups), shape = factor(Groups)), size = 5) +
    # facet_wrap(~ Groups) +
    theme_light()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_y_continuous(name = "Relative mRNA expression",
                       limits=c(-0.6,1.6),
                       breaks=c(-0.5,0,0.5,1,1.5)) +
    scale_shape_manual(values=c(seq(0,10))) +
    # scale_size_manual(values = rep(10,43)) +
    scale_color_manual(values = c("cadetblue1",brewer.pal(10, "Paired"))) +
    labs(title=paste0(names(Test)[j]," genelist expression various conditions"),
         x ="Samples", y = "Relative mRNA expression",
         colour = "Cell populations", shape = "Cell populations")
  
  pdf(file=paste0(output.dir,"Genelist_expression_plots/Expression_plot_",names(Test)[j],"_average.pdf"), height = 10, width = 15)
  print(p1)
  dev.off()
  
  pdf(file=paste0(output.dir,"Genelist_expression_plots/Expression_plot_",names(Test)[j],"_individual.pdf"), height = 10, width = 15)
  print(p2)
  dev.off()
}

#########################################################################

## Extra 08/2023: IFN gene signature
## Simplify annotation
seuratObj$annotated_clusters_Muscat<-seuratObj$annotated_clusters
levels(seuratObj$annotated_clusters_Muscat)<-c('Early mature cDC1s',"Late mature cDC1s","Early mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Late mature cDC1s",
                                               "Proliferating cDC1s","Proliferating cDC1s","Late mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Other cDC1s","Late mature cDC1s",
                                               "Early mature cDC1s","Late mature cDC1s","Doublets","Late mature cDC1s","Early mature cDC1s")

seuratObj$annotated_clusters_Muscat_v2<-paste0(seuratObj$orig.ident,"_",seuratObj$annotated_clusters_Muscat)

DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2")

## Subset to only keep pops to compare
Idents(seuratObj)<-seuratObj$annotated_clusters_Muscat_v2
seuratObj_subset<-subset(seuratObj, idents = levels(Idents(seuratObj))[c(1,2,4,6,11,17,22,29,36,37,42)]) #Also included 1 (Imm SS)

DimPlot(seuratObj_subset, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2", split.by = "orig.ident", ncol = 3)

seuratObj_subset$Condition_v2<-seuratObj_subset$Condition
levels(seuratObj_subset$Condition_v2)<-c("Steady state","eLNPs","pIC-LNPs","CpG-LNPs","pIC alone",
                                         "eLNPs","pIC-LNPs","CpG-LNPs","pIC alone")
seuratObj_subset$Condition_v2<-gsub("-","_",seuratObj_subset$Condition_v2)
seuratObj_subset$Condition_v2<-gsub(" ","_",seuratObj_subset$Condition_v2)

## Combine MULTI_ID with Orig.ident or Genotype?? Need to put genotype in correct order: otherwise WT after DKO
seuratObj_subset@meta.data[["MULTI_ID_merge"]]
seuratObj_subset@meta.data$MULTI_ID_merge<-as.factor(seuratObj_subset@meta.data$MULTI_ID_merge)

## Update muscat annotation
seuratObj_subset$annotated_clusters_Muscat<-factor(as.character(seuratObj_subset$annotated_clusters_Muscat), 
                                                   levels = c("Immature cDC1s", "Early mature cDC1s", "Late mature cDC1s")) #Adapted for Imm
levels(seuratObj_subset$annotated_clusters_Muscat)<-c("Immature","Early_mature","Late_mature")

## Try new sample_ID
seuratObj_subset$sample_ID<-paste0(as.character(seuratObj_subset$MULTI_ID),"_",as.character(seuratObj_subset$Condition_v2),
                                   "_",as.character(seuratObj_subset$annotated_clusters_Muscat))
seuratObj_subset$sample_ID<-factor(seuratObj_subset$sample_ID, #Added Imm replicates
                                   levels = c("Hashtag1_Steady_state_Immature","Hashtag2_Steady_state_Immature","Hashtag3_Steady_state_Immature","Hashtag4_Steady_state_Immature",
                                              "Hashtag1_Steady_state_Early_mature","Hashtag2_Steady_state_Early_mature","Hashtag3_Steady_state_Early_mature","Hashtag4_Steady_state_Early_mature",
                                              "Hashtag1_Steady_state_Late_mature","Hashtag2_Steady_state_Late_mature","Hashtag3_Steady_state_Late_mature","Hashtag4_Steady_state_Late_mature",
                                              "Hashtag1_eLNPs_Early_mature","Hashtag2_eLNPs_Early_mature","Hashtag3_eLNPs_Early_mature","Hashtag4_eLNPs_Early_mature",       
                                              "Hashtag1_eLNPs_Late_mature","Hashtag2_eLNPs_Late_mature","Hashtag3_eLNPs_Late_mature","Hashtag4_eLNPs_Late_mature",
                                              "Hashtag1_pIC_LNPs_Early_mature","Hashtag2_pIC_LNPs_Early_mature","Hashtag3_pIC_LNPs_Early_mature","Hashtag4_pIC_LNPs_Early_mature",
                                              "Hashtag1_pIC_LNPs_Late_mature","Hashtag2_pIC_LNPs_Late_mature","Hashtag3_pIC_LNPs_Late_mature","Hashtag4_pIC_LNPs_Late_mature",
                                              "Hashtag1_CpG_LNPs_Early_mature","Hashtag2_CpG_LNPs_Early_mature","Hashtag3_CpG_LNPs_Early_mature","Hashtag4_CpG_LNPs_Early_mature",
                                              "Hashtag1_CpG_LNPs_Late_mature","Hashtag2_CpG_LNPs_Late_mature","Hashtag3_CpG_LNPs_Late_mature",
                                              "Hashtag1_pIC_alone_Early_mature","Hashtag2_pIC_alone_Early_mature","Hashtag3_pIC_alone_Early_mature","Hashtag4_pIC_alone_Early_mature",
                                              "Hashtag1_pIC_alone_Late_mature","Hashtag2_pIC_alone_Late_mature","Hashtag3_pIC_alone_Late_mature","Hashtag4_pIC_alone_Late_mature"))
levels(seuratObj_subset$sample_ID)

## Update clusters
seuratObj_subset$annotated_clusters_Muscat_v2<-factor(paste0(seuratObj_subset$Condition_v2,"_",seuratObj_subset$annotated_clusters_Muscat),
                                                      levels = c("Steady_state_Immature","Steady_state_Early_mature","Steady_state_Late_mature",
                                                                 "eLNPs_Early_mature","eLNPs_Late_mature","pIC_LNPs_Early_mature","pIC_LNPs_Late_mature",
                                                                 "CpG_LNPs_Early_mature","CpG_LNPs_Late_mature" ,"pIC_alone_Early_mature","pIC_alone_Late_mature"))
DimPlot(seuratObj_subset, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2",split.by = "Condition_v2", ncol = 3)

## Basic annotation
seuratObj_subset$Basic_annotation<-"cDC1s"

## Look at average expression per sample
Idents(seuratObj_subset)<-seuratObj_subset$sample_ID
seuratObj_subset_average <- AverageExpression(seuratObj_subset, return.seurat = T)

## Gene lists
Genelist_IFN<-read.xlsx("../Documentation/Amigo2_GO_0009615_Mus_musculus.xlsx")
Genes_IFN<-unique(Genelist_IFN$Gene.symbol)

SCT_markers_WT<-read.xlsx("../../../RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Marker_lists/SCTmarkersList_SCTclus_SAM2and3_WT_subset_final_2021.xlsx",
                          sheet = "Early mature cDC1s")
Genes_IFN_v2<-c(SCT_markers_WT$gene[c(1,4,5,6,7,8,9,19,20,21,22)],"Ifit1","Ifit2","Ifit3")

seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = list(intersect(Genes_IFN,rownames(seuratObj_subset_average))), name = "IFN_signature")
seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = list(intersect(Genes_IFN_v2,rownames(seuratObj_subset_average))), name = "IFN_signature_v2")

Groups<-c(rep("SS_Imm",4),rep("SS_EM",4),rep("SS_LM",4),
          rep("eLNPs_EM",4),rep("eLNPs_LM",4),
          rep("pIC_LNPs_EM",4),rep("pIC_LNPs_LM",4),
          rep("CpG_LNPs_EM",4),rep("CpG_LNPs_LM",3),
          rep("pIC_EM",4),rep("pIC_LM",4))
Groups_v2<-seq(1:43)

Test<-tibble::lst()
counter = 1
for (j in c(8,9)){
  Test[[counter]]<-as.data.frame(cbind(seuratObj_subset_average@meta.data[,j],Groups,Groups_v2), row.names = rownames(seuratObj_subset_average@meta.data))
  counter=counter+1
}

names(Test)<-c("IFN_signature","IFN_signature_v2")

for (j in 1:length(Test)){
  Test[[j]][["Groups_v2"]]<-factor(Test[[j]][["Groups_v2"]], levels = seq(1:43))
  Test[[j]][["Groups"]]<-factor(Test[[j]][["Groups"]], levels = c("SS_Imm","SS_EM","SS_LM","eLNPs_EM","eLNPs_LM",
                                                                  "pIC_LNPs_EM","pIC_LNPs_LM","CpG_LNPs_EM","CpG_LNPs_LM",
                                                                  "pIC_EM","pIC_LM"))
  Test[[j]][["V1"]]<-as.numeric(Test[[j]][["V1"]])
}


dir.create(paste0(output.dir,"Genelist_expression_plots/"))

for (j in 1:length(Test)){
  p1<- ggplot(Test[[j]], aes(x=as.factor(Groups), y = V1, fill=as.factor(Groups) )) + 
    geom_bar(stat = "summary", fun = mean) +
    scale_fill_manual(values =c("cadetblue1",brewer.pal(10, "Paired"))) +
    theme(legend.position="none") +
    # theme(axis.text.x=element_blank(),
    #       axis.ticks.x=element_blank()) +
    scale_y_continuous(name = "Relative mRNA expression",
                       limits=c(-0.6,1.6),
                       breaks=c(-0.5,0,0.5,1,1.5)) +
    # scale_shape_manual(values=c(seq(0,9))) +
    labs(title=paste0(names(Test)[j]," genelist expression across various conditions"),
         x ="Samples", y = "Relative mRNA expression")
  
  p2 <- ggplot(Test[[j]], aes(x = as.factor(Groups_v2), y =V1)) +
    geom_point(aes(colour = factor(Groups), shape = factor(Groups)), size = 5) +
    # facet_wrap(~ Groups) +
    theme_light()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_y_continuous(name = "Relative mRNA expression",
                       limits=c(-0.6,1.6),
                       breaks=c(-0.5,0,0.5,1,1.5)) +
    scale_shape_manual(values=c(seq(0,10))) +
    # scale_size_manual(values = rep(10,43)) +
    scale_color_manual(values = c("cadetblue1",brewer.pal(10, "Paired"))) +
    labs(title=paste0(names(Test)[j]," genelist expression various conditions"),
         x ="Samples", y = "Relative mRNA expression",
         colour = "Cell populations", shape = "Cell populations")
  
  pdf(file=paste0(output.dir,"Genelist_expression_plots/Expression_plot_",names(Test)[j],"_average.pdf"), height = 10, width = 15)
  print(p1)
  dev.off()
  
  pdf(file=paste0(output.dir,"Genelist_expression_plots/Expression_plot_",names(Test)[j],"_individual.pdf"), height = 10, width = 15)
  print(p2)
  dev.off()
}

###########################################################################

## Finalization for paper 2024
table(seuratObj$annotated_clusters) #Only 2 cells in last 2 clusters each?! Strange clustering artifact? No markers for these clusters of course...

seuratObj@meta.data$annotated_clusters_paper_2024<-seuratObj@meta.data$annotated_clusters
levels(seuratObj@meta.data$annotated_clusters_paper_2024) <- c('Il12b+ cDC1s',"Tnfrsf4 hi late mature cDC1s","Cxcl9-11 hi cDC1s",
                                                              "Ms4a4c+ late mature cDC1s","Immature cDC1s","Tmem176+ late mature cDC1s",
                                                              "Proliferating cDC1s 1","Proliferating cDC1s 2","Hamp+ late mature cDC1s",
                                                              "Ccl5 hi late mature cDC1s","Serpina3+ cDC1s","Other cDC1s 1","H2-M2+ late mature cDC1s",
                                                              "IFN+ cDC1s","Zmynd15+ late mature cDC1s","Other cDC1s 2","Other late mature cDC1s","Other early mature cDC1s")

seuratObj$annotated_clusters_Muscat_paper_2024<-seuratObj$annotated_clusters_Muscat
levels(seuratObj$annotated_clusters_Muscat_paper_2024)<-c('Early mature cDC1s',"Late mature cDC1s","Immature cDC1s","Proliferating cDC1s","Other cDC1s 1","Other cDC1s 2")

seuratObj$annotated_clusters_Muscat_v2_paper_2024<-factor(paste0(seuratObj$Condition_v2,"_",seuratObj$annotated_clusters_Muscat_paper_2024), levels = c("Steady_state_Immature cDC1s","Steady_state_Proliferating cDC1s",
                                                                                                                                                        "Steady_state_Other cDC1s 2","Steady_state_Early mature cDC1s",
                                                                                                                                                        "Steady_state_Late mature cDC1s","eLNPs_2h_Immature cDC1s",
                                                                                                                                                        "eLNPs_2h_Proliferating cDC1s","eLNPs_2h_Other cDC1s 1" ,
                                                                                                                                                        "eLNPs_2h_Early mature cDC1s","eLNPs_2h_Late mature cDC1s",
                                                                                                                                                        "eLNPs_8h_Immature cDC1s","eLNPs_8h_Proliferating cDC1s",
                                                                                                                                                        "eLNPs_8h_Other cDC1s 2","eLNPs_8h_Early mature cDC1s",
                                                                                                                                                        "eLNPs_8h_Late mature cDC1s","pIC_alone_2h_Immature cDC1s",
                                                                                                                                                        "pIC_alone_2h_Proliferating cDC1s","pIC_alone_2h_Other cDC1s 2",
                                                                                                                                                        "pIC_alone_2h_Early mature cDC1s","pIC_alone_2h_Late mature cDC1s",
                                                                                                                                                        "pIC_alone_8h_Immature cDC1s","pIC_alone_8h_Proliferating cDC1s",
                                                                                                                                                        "pIC_alone_8h_Other cDC1s 1","pIC_alone_8h_Other cDC1s 2",
                                                                                                                                                        "pIC_alone_8h_Early mature cDC1s","pIC_alone_8h_Late mature cDC1s",
                                                                                                                                                        "pIC_LNPs_2h_Immature cDC1s", "pIC_LNPs_2h_Proliferating cDC1s",
                                                                                                                                                        "pIC_LNPs_2h_Early mature cDC1s", "pIC_LNPs_2h_Late mature cDC1s",
                                                                                                                                                        "pIC_LNPs_8h_Immature cDC1s","pIC_LNPs_8h_Proliferating cDC1s",
                                                                                                                                                        "pIC_LNPs_8h_Other cDC1s 1",
                                                                                                                                                        "pIC_LNPs_8h_Early mature cDC1s", "pIC_LNPs_8h_Late mature cDC1s", 
                                                                                                                                                        "CpG_LNPs_2h_Immature cDC1s", "CpG_LNPs_2h_Proliferating cDC1s",
                                                                                                                                                        "CpG_LNPs_2h_Other cDC1s 2",
                                                                                                                                                        "CpG_LNPs_2h_Early mature cDC1s","CpG_LNPs_2h_Late mature cDC1s",
                                                                                                                                                        "CpG_LNPs_8h_Immature cDC1s","CpG_LNPs_8h_Proliferating cDC1s",
                                                                                                                                                        "CpG_LNPs_8h_Other cDC1s 1","CpG_LNPs_8h_Other cDC1s 2",
                                                                                                                                                        "CpG_LNPs_8h_Early mature cDC1s", "CpG_LNPs_8h_Late mature cDC1s"))
                                                      
seuratObj@meta.data$newClusters_paper_2024<-paste0(seuratObj@meta.data$annotated_clusters_Muscat_paper_2024,"_",seuratObj@meta.data$Condition_v2)

seuratObj$Condition<-factor(seuratObj$Condition, levels = c("Steady state","eLNPs 2h", "eLNPs 8h","pIC alone 2h","pIC alone 8h",
                                                            "pIC-LNPs 2h","pIC-LNPs 8h","CpG-LNPs 2h","CpG-LNPs 8h"))

pdf(file=paste0(output.dir,"Annotation/9_Paper_2024_annotated_UMAP_RNA_harmony_detailed_",experiment,".pdf"), width =20, height= 12)
DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_paper_2024", label = T, label.size = 6, repel = T)
dev.off()

pdf(file=paste0(output.dir,"Annotation/9_Paper_2024_annotated_UMAP_RNA_harmony_detailed_no_label_",experiment,".pdf"), width =20, height= 12)
DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_paper_2024", label = F)
dev.off()

pdf(file=paste0(output.dir,"Annotation/9_Paper_2024_annotated_UMAP_RNA_harmony_muscat_annotation_",experiment,".pdf"), width =19, height= 12)
DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_paper_2024", label = T, label.size = 6, repel = T)
dev.off()

pdf(file=paste0(output.dir,"Annotation/9_Paper_2024_annotated_UMAP_RNA_harmony_muscat_annotation_no_label_",experiment,".pdf"), width =19, height= 12)
DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_paper_2024", label = F)
dev.off()

pdf(file=paste0(output.dir,"Annotation/9_Paper_2024_annotated_UMAP_RNA_harmony_muscat_annotation_split_",experiment,".pdf"), width =40, height= 15)
DimPlot(seuratObj, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_paper_2024", split.by = "Condition", ncol = 5)
dev.off()

pdf(file=paste0(output.dir,"Annotation/9_Paper_2024_annotated_UMAP_RNA_no_label_",experiment,".pdf"), width =15, height= 12)
DimPlot(seuratObj, reduction = "RNA_umap", label = F, group.by = "Condition",  order = rev(c("Steady state","eLNPs 2h","eLNPs 8h","pIC alone 2h","pIC alone 8h",
                                                                                             "pIC-LNPs 2h","pIC-LNPs 8h","CpG-LNPs 2h","CpG-LNPs 8h")),
       cols = brewer.pal(10, "Paired")[c(2:4,9:10,5:8)])
dev.off()

pdf(file=paste0(output.dir,"Annotation/9_Paper_2024_annotated_UMAP_RNA_split_per_sample_no_label_",experiment,".pdf"), width =40, height= 15)
DimPlot(seuratObj, reduction = "RNA_umap", label = F, split.by = "Condition", group.by = "Condition", ncol = 5, cols = brewer.pal(10, "Paired")[c(2:4,9:10,5:8)])
dev.off()

## Extra plot SS with annotation on RNA UMAP
Idents(seuratObj)<-seuratObj$Condition
seuratObj_SS<-subset(seuratObj, idents = "Steady state")

pdf(file=paste0(output.dir,"Annotation/9_Paper_2024_annotated_RNA_UMAP_muscat_annotation_split_",experiment,".pdf"), width =40, height= 15)
DimPlot(seuratObj, reduction = "RNA_umap", group.by = "annotated_clusters_Muscat_paper_2024", split.by = "Condition", ncol = 5)
dev.off()

pdf(file=paste0(output.dir,"Annotation/9_Paper_2024_annotated_RNA_UMAP_muscat_annotation_SS_",experiment,".pdf"), width =13, height= 12)
DimPlot(seuratObj_SS, reduction = "RNA_umap", group.by = "annotated_clusters_Muscat_paper_2024")
dev.off()

###############

### Find ADT markers for detailed annotated Harmony RNA clusters paper
DefaultAssay(seuratObj) <- "ADT"
Idents(seuratObj)<-seuratObj$annotated_clusters_paper_2024

allADTMarkers <- FindAllMarkers(seuratObj, min.pct = 0.10,logfc.threshold = 0.30, return.thresh = 0.01, only.pos = TRUE)
table(allADTMarkers$cluster)
saveRDS(allADTMarkers, file=paste0(output.dir,"Robjects/allADTMarkers_detailed_annotation_RNA_harmony_",experiment,".rds"))

### Add to diagnostics
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))
diagnostics[['ADTmarkersPerClusterDetailedHarmonyAnnotation']]<-paste0(table(allADTMarkers$cluster)," markers for cluster detailed Harmony annotation ",rownames(table(allADTMarkers$cluster)))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

### Create list with markers
totalNrClusters<-names(table(allADTMarkers$cluster))
ADTmarkersList<-list()

for(i in totalNrClusters){
  tmp<-allADTMarkers[allADTMarkers$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  ADTmarkersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(ADTmarkersList)<-levels(Idents(seuratObj))

### Write to Excel
write.xlsx(ADTmarkersList, file =paste0(output.dir, "/allClusters_ADT_markers_detailed_annotation_RNA_harmony_",experiment,".xlsx"))

############

### Find RNA markers for detailed annotated Harmony RNA clusters paper
DefaultAssay(seuratObj) <- "RNA"

allRNAMarkers <- FindAllMarkers(seuratObj, min.pct = 0.10,logfc.threshold = 0.30, return.thresh = 0.01, only.pos = TRUE)
table(allRNAMarkers$cluster)
saveRDS(allRNAMarkers, file=paste0(output.dir,"Robjects/allRNAMarkers_detailed_annotation_RNA_harmony_",experiment,".rds"))

### Add to diagnostics
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))
diagnostics[['RNAmarkersPerClusterDetailedHarmonyAnnotation']]<-paste0(table(allRNAMarkers$cluster)," markers for cluster detailed Harmony annotation ",rownames(table(allRNAMarkers$cluster)))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

### Create list with markers
totalNrClusters<-names(table(allRNAMarkers$cluster))
RNAmarkersList<-list()

for(i in totalNrClusters){
  tmp<-allRNAMarkers[allRNAMarkers$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  RNAmarkersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(RNAmarkersList)<-levels(Idents(seuratObj))

### Write to Excel
write.xlsx(RNAmarkersList, file =paste0(output.dir, "/allClusters_RNA_markers_detailed_annotation_RNA_harmony_",experiment,".xlsx"))

#################

### Find ADT markers for simple muscat annotated Harmony RNA clusters paper
DefaultAssay(seuratObj) <- "ADT"
Idents(seuratObj)<-seuratObj$annotated_clusters_Muscat_paper_2024

allADTMarkers <- FindAllMarkers(seuratObj, min.pct = 0.10,logfc.threshold = 0.30, return.thresh = 0.01, only.pos = TRUE)
table(allADTMarkers$cluster)
saveRDS(allADTMarkers, file=paste0(output.dir,"Robjects/allADTMarkers_simple_muscat_annotation_RNA_harmony_",experiment,".rds"))

### Add to diagnostics
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))
diagnostics[['ADTmarkersPerClusterSimpleHarmonyAnnotation']]<-paste0(table(allADTMarkers$cluster)," markers for cluster simple Harmony annotation ",rownames(table(allADTMarkers$cluster)))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

### Create list with markers
totalNrClusters<-names(table(allADTMarkers$cluster))
ADTmarkersList<-list()

for(i in totalNrClusters){
  tmp<-allADTMarkers[allADTMarkers$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  ADTmarkersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(ADTmarkersList)<-levels(Idents(seuratObj))

### Write to Excel
write.xlsx(ADTmarkersList, file =paste0(output.dir, "/allClusters_ADT_markers_simple_muscat_annotation_RNA_harmony_",experiment,".xlsx"))

############

### Find RNA markers for simple muscat annotated Harmony RNA clusters paper
DefaultAssay(seuratObj) <- "RNA"

allRNAMarkers <- FindAllMarkers(seuratObj, min.pct = 0.10,logfc.threshold = 0.30, return.thresh = 0.01, only.pos = TRUE)
table(allRNAMarkers$cluster)
saveRDS(allRNAMarkers, file=paste0(output.dir,"Robjects/allRNAMarkers_simple_muscat_annotation_RNA_harmony_",experiment,".rds"))

### Add to diagnostics
diagnostics<-readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))
diagnostics[['RNAmarkersPerClusterSimpleHarmonyAnnotation']]<-paste0(table(allRNAMarkers$cluster)," markers for cluster simple Harmony annotation ",rownames(table(allRNAMarkers$cluster)))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

### Create list with markers
totalNrClusters<-names(table(allRNAMarkers$cluster))
RNAmarkersList<-list()

for(i in totalNrClusters){
  tmp<-allRNAMarkers[allRNAMarkers$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  RNAmarkersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(RNAmarkersList)<-levels(Idents(seuratObj))

### Write to Excel
write.xlsx(RNAmarkersList, file =paste0(output.dir, "/allClusters_RNA_markers_simple_muscat_annotation_RNA_harmony_",experiment,".xlsx"))


#########################################################################

## Paper gene signatures with final colors and symbols
DimPlot(seuratObj, reduction = "RNA_harmony_umap", label = T, group.by = "annotated_clusters_Muscat_v2_paper_2024") + NoLegend()

## Subset to only keep pops to compare
Idents(seuratObj)<-seuratObj$annotated_clusters_Muscat_v2_paper_2024
seuratObj_subset<-subset(seuratObj, idents = levels(Idents(seuratObj))[c(1,4,5,9,15,19,26,29,35,39,46)]) #Also included 1 (Imm SS)

DimPlot(seuratObj_subset, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2_paper_2024", split.by = "Condition", ncol = 3)

## Update muscat annotation
seuratObj_subset$annotated_clusters_Muscat_paper_2024 <-factor(as.character(seuratObj_subset$annotated_clusters_Muscat_paper_2024),
                                                               levels = c("Immature cDC1s", "Early mature cDC1s", "Late mature cDC1s")) #Adapted for Imm
levels(seuratObj_subset$annotated_clusters_Muscat_paper_2024)<-c("Immature","Early_mature","Late_mature")

## Try new sample_ID
seuratObj_subset$sample_ID<-paste0(as.character(seuratObj_subset$MULTI_ID),"_",as.character(seuratObj_subset$Condition_v2),
                                   "_",as.character(seuratObj_subset$annotated_clusters_Muscat_paper_2024))
seuratObj_subset$sample_ID<-factor(seuratObj_subset$sample_ID, #Added Imm replicates
                                   levels = c("Hashtag1_Steady_state_Immature","Hashtag2_Steady_state_Immature","Hashtag3_Steady_state_Immature","Hashtag4_Steady_state_Immature",
                                              "Hashtag1_Steady_state_Early_mature","Hashtag2_Steady_state_Early_mature","Hashtag3_Steady_state_Early_mature","Hashtag4_Steady_state_Early_mature",
                                              "Hashtag1_Steady_state_Late_mature","Hashtag2_Steady_state_Late_mature","Hashtag3_Steady_state_Late_mature","Hashtag4_Steady_state_Late_mature",
                                              "Hashtag1_eLNPs_2h_Early_mature","Hashtag2_eLNPs_2h_Early_mature","Hashtag3_eLNPs_2h_Early_mature","Hashtag4_eLNPs_2h_Early_mature",       
                                              "Hashtag1_eLNPs_8h_Late_mature","Hashtag2_eLNPs_8h_Late_mature","Hashtag3_eLNPs_8h_Late_mature","Hashtag4_eLNPs_8h_Late_mature",
                                              "Hashtag1_pIC_alone_2h_Early_mature","Hashtag2_pIC_alone_2h_Early_mature","Hashtag3_pIC_alone_2h_Early_mature","Hashtag4_pIC_alone_2h_Early_mature",
                                              "Hashtag1_pIC_alone_8h_Late_mature","Hashtag2_pIC_alone_8h_Late_mature","Hashtag3_pIC_alone_8h_Late_mature","Hashtag4_pIC_alone_8h_Late_mature",
                                              "Hashtag1_pIC_LNPs_2h_Early_mature","Hashtag2_pIC_LNPs_2h_Early_mature","Hashtag3_pIC_LNPs_2h_Early_mature","Hashtag4_pIC_LNPs_2h_Early_mature",
                                              "Hashtag1_pIC_LNPs_8h_Late_mature","Hashtag2_pIC_LNPs_8h_Late_mature","Hashtag3_pIC_LNPs_8h_Late_mature","Hashtag4_pIC_LNPs_8h_Late_mature",
                                              "Hashtag1_CpG_LNPs_2h_Early_mature","Hashtag2_CpG_LNPs_2h_Early_mature","Hashtag3_CpG_LNPs_2h_Early_mature","Hashtag4_CpG_LNPs_2h_Early_mature",
                                              "Hashtag1_CpG_LNPs_8h_Late_mature","Hashtag2_CpG_LNPs_8h_Late_mature","Hashtag3_CpG_LNPs_8h_Late_mature"))
levels(seuratObj_subset$sample_ID)
seuratObj_subset$sample_ID[is.na(seuratObj_subset$sample_ID)]

## Look at average expression per sample
Idents(seuratObj_subset)<-seuratObj_subset$sample_ID
seuratObj_subset_average <- AverageExpression(seuratObj_subset, return.seurat = T)

## Gene lists Ardouin
Genelist_tibble<-readRDS("../../../Meta_analysis/Tibble_all_genelists_v3.rds") #Still has the Ardouin lists from their paper!! Homeo, Common, Immuno (Checked nr of genes!)

for (i in 11:13){ #Ardouin lists
  genelist<-Genelist_tibble[[i]]
  name<-names(Genelist_tibble[i])
  
  seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = list(intersect(genelist,rownames(seuratObj_subset_average))), name = name)
}

## IFN gene lists
Genelist_IFN<-read.xlsx("../Documentation/Amigo2_GO_0009615_Mus_musculus.xlsx")
Genes_IFN<-unique(Genelist_IFN$Gene.symbol) #330 genes!!

SCT_markers_WT<-read.xlsx("../../../RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Marker_lists/SCTmarkersList_SCTclus_SAM2and3_WT_subset_final_2021.xlsx",
                          sheet = "Early mature cDC1s")
Genes_IFN_v2<-c(SCT_markers_WT$gene[c(1,4,5,6,7,8,9,19,20,21,22)],"Ifit1","Ifit2","Ifit3") #14 genes
# [1] "Cxcl10" "Cxcl9"  "Ifi47"  "Cd40"   "Gbp2"   "Gbp5"   "Isg15"  "Nfkbia" "Ifi204"
# [10] "Ifi211" "Irf7"   "Ifit1"  "Ifit2"  "Ifit3" 

seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = list(intersect(Genes_IFN,rownames(seuratObj_subset_average))), name = "IFN_signature")
seuratObj_subset_average <- AddModuleScore(seuratObj_subset_average, features = list(intersect(Genes_IFN_v2,rownames(seuratObj_subset_average))), name = "IFN_signature_v2")

## Prep plot
Groups<-c(rep("SS_Imm",4),rep("SS_EM",4),rep("SS_LM",4),
          rep("eLNPs_2h_EM",4),rep("eLNPs_8h_LM",4),
          rep("pIC_2h_EM",4),rep("pIC_8h_LM",4),
          rep("pIC_LNPs_2h_EM",4),rep("pIC_LNPs_8h_LM",4),
          rep("CpG_LNPs_2h_EM",4),rep("CpG_LNPs_8h_LM",3))
Groups_v2<-seq(1:43)

## Check metadata seuratobj_subset_average
Test<-tibble::lst()
counter = 1
for (j in c(8:12)){ 
  Test[[counter]]<-as.data.frame(cbind(seuratObj_subset_average@meta.data[,j],Groups,Groups_v2), row.names = rownames(seuratObj_subset_average@meta.data))
  counter=counter+1
}

names(Test)<-c(names(Genelist_tibble)[11:13],"IFN_signature","IFN_signature_v2")

for (j in 1:length(Test)){
  Test[[j]][["Groups_v2"]]<-factor(Test[[j]][["Groups_v2"]], levels = seq(1:43))
  Test[[j]][["Groups"]]<-factor(Test[[j]][["Groups"]], levels = c("SS_Imm","SS_EM","SS_LM","eLNPs_2h_EM","eLNPs_8h_LM","pIC_2h_EM","pIC_8h_LM",
                                                                  "pIC_LNPs_2h_EM","pIC_LNPs_8h_LM","CpG_LNPs_2h_EM","CpG_LNPs_8h_LM"))
  Test[[j]][["V1"]]<-as.numeric(Test[[j]][["V1"]])
}


dir.create(paste0(output.dir,"Genelist_expression_plots/Paper_2024/"))

for (j in 1:length(Test)){
  p1<- ggplot(Test[[j]], aes(x=as.factor(Groups), y = V1, fill=as.factor(Groups) )) + 
    geom_bar(stat = "summary", fun = mean) +
    scale_fill_manual(values =c("cadetblue1",brewer.pal(10, "Paired")[1:4],brewer.pal(10, "Paired")[9:10],brewer.pal(10, "Paired")[5:8])) +
    theme(legend.position="none") +
    # theme(axis.text.x=element_blank(),
    #       axis.ticks.x=element_blank()) +
    scale_y_continuous(name = "Relative mRNA expression",
                       limits=c(-0.6,1.8),
                       breaks=c(-0.5,0,0.5,1,1.5)) +
    # scale_shape_manual(values=c(seq(0,9))) +
    labs(title=paste0(names(Test)[j]," genelist expression across various conditions"),
         x ="Samples", y = "Relative mRNA expression")
  
  p2 <- ggplot(Test[[j]], aes(x = as.factor(Groups_v2), y =V1)) +
    geom_point(aes(colour = factor(Groups), shape = factor(Groups)), size = 5) +
    # facet_wrap(~ Groups) +
    theme_light()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_y_continuous(name = "Relative mRNA expression",
                       limits=c(-0.6,1.8), #Had to increase limit otherwise point left off of graph! Always check warnings!
                       breaks=c(-0.5,0,0.5,1,1.5)) +
    scale_shape_manual(values=c(1,10,13,0,12,5,9,4,8,2,6)) + #seq(1,11)
    # scale_size_manual(values = rep(10,43)) +
    scale_color_manual(values = c("cadetblue1",brewer.pal(10, "Paired")[1:4],brewer.pal(10, "Paired")[9:10],brewer.pal(10, "Paired")[5:8])) +
    labs(title=paste0(names(Test)[j]," genelist expression various conditions"),
         x ="Samples", y = "Relative mRNA expression",
         colour = "Cell populations", shape = "Cell populations")
  
  pdf(file=paste0(output.dir,"Genelist_expression_plots/Paper_2024/Expression_plot_v2_",names(Test)[j],"_average.pdf"), height = 8, width = 10)
  print(p1)
  dev.off()
  
  pdf(file=paste0(output.dir,"Genelist_expression_plots/Paper_2024/Expression_plot_v2_",names(Test)[j],"_individual.pdf"), height = 8, width = 10)
  print(p2)
  dev.off()
}

#############################

dir.create(paste0(output.dir,"Feature_plots/Paper_2024/"))

## Featureplots of interesting marker genes and proteins
features<-c("CD43","CD117","Flt3","CD62L","Sell","XCR1","Xcr1","ESAM","Esam","CD103","Itgae","CD207-mh","Cd207",
            "CD197","Ccr7","Fscn1","Ccl22","CD63","Cd63","CD274","CD278","CD80","Cd80","Cd83","CD86","Cd86",
            "Cd69","Mki67","Tmem176a","Tmem176b","Il12b","Tnfrsf4","Cxcl9","Cxcl11",
            "Ms4a4c","Egr3","Hamp","Ccl5","Serpina3f","Serpina3g","H2-M2","Ifna2","Zmynd15",
            "Icosl","Scube3","Slco5a1","Ifit2","Isg20","Oas3","Ifi47",
            "Dnase1l3","Ch25h","Nlrp3","Ido1","Ido2","Kmo","Kyat3","Ern1","Josd1",
            "Srebf1","Srebf2","Nr1h2","Nr1h3","Abca1","Abcg1","Apol7c","Apoe",
            "Btla","Havcr2","Bag6","Lgals9","Casp4","Casp8","Smpdl3a",
            "Cxcl16") #Added Cxcl16 (29/02/24)

pdf(file=paste0(output.dir,"Feature_plots/Paper_2024/Feature_plot_markers_paper_2024_split_RNA_harmony_UMAP_",samplename,"_blue_grey.pdf"), height = 4, width = 40)
for (feature in features) {
  F1<-FeaturePlot(object = seuratObj, features =feature, cols = c("grey", "blue"), 
                  reduction = "RNA_harmony_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T,
                  split.by = "Condition", raster = T) & theme(legend.position = "right") #+ #wnn.umap
  # scale_color_viridis(option = "C")
  print(F1)
}
dev.off()

SCT_markers_WT<-read.xlsx("../../../RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Marker_lists/SCTmarkersList_SCTclus_SAM2and3_WT_subset_final_2021.xlsx",
                          sheet = "Early mature cDC1s")

ISG_signature <- list(c(SCT_markers_WT$gene[c(1,4,5,6,7,8,9,19,20,21,22)],"Ifit1","Ifit2","Ifit3"))

seuratObj <- AddModuleScore(object = seuratObj, features = ISG_signature, name = "ISG_signature_score")

F2<-FeaturePlot(object = seuratObj, features = "ISG_signature_score1", label = F, repel =T, order = T, cols = c("Yellow","Red"),
                reduction = "RNA_harmony_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, #wnn.umap
                split.by = "Condition", raster = T, keep.scale = "all") & theme(legend.position = "right") 

pdf(file = paste0(output.dir,"Feature_plots/Paper_2024/Feature_plot_markers_paper_2024_ISG_modulescore_split_RNA_harmony_UMAP_v2_",samplename,".pdf"), width = 40, height = 4)
F2
dev.off()


#############################

## Save objects
saveRDS(seuratObj, file = paste0(output.dir,"Robjects/seuratObj_",experiment,".rds"))
saveRDS(diagnostics, file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

## Read objects
seuratObj <- readRDS(file = paste0(output.dir,"Robjects/seuratObj_",experiment,".rds"))
diagnostics <- readRDS(file=paste0(output.dir,"Robjects/diagnostics_",experiment,".rds"))

#############################

## Export data for Single-cell.be (Use Test_Jade conda env!!)
seuratToParquetSimple(
  seuratObj,
  out = paste0(output.dir,"Robjects/Parquet_file_",experiment),
  assays = c("RNA","ADT"),
  metadata.columns = c("annotated_clusters_Muscat_paper_2024","annotated_clusters_Muscat_v2_paper_2024","annotated_clusters_paper_2024","Condition","MULTI_ID_merge"),
  metadata.names = c("annotated_clusters_Muscat_paper_2024","annotated_clusters_Muscat_paper_2024_split","annotated_clusters_paper_2024","Condition","MULTI_ID_merge"),
  reduction = "RNA_harmony_umap",
  chunk.size = 10000
)

############################################################################################################

## Create diet object for online tool (03/06/24)
seuratObj_diet<-DietSeurat(seuratObj, counts = T, data = T, scale.data = F,
                                    assays = c("RNA","ADT"), dimreducs = "RNA_harmony_umap", graphs = NULL)

## New metadata names
## c("annotated_clusters_Muscat_paper_2024","annotated_clusters_Muscat_v2_paper_2024","annotated_clusters_paper_2024","Condition","MULTI_ID_merge")
seuratObj_diet$Annotation_paper_detailed<-seuratObj_diet$annotated_clusters_paper_2024
seuratObj_diet$Annotation_paper<-seuratObj_diet$annotated_clusters_Muscat_paper_2024
seuratObj_diet$Annotation_paper_split_by_condition<-seuratObj_diet$annotated_clusters_Muscat_v2_paper_2024
seuratObj_diet$Sample<-seuratObj_diet$orig.ident
seuratObj_diet$Replicate_info<-seuratObj_diet$MULTI_ID_merge

## Update detailed annotation (Immature cDC1s and 2 Other cDC1s clusters) -> same as in paper_annotation
levels(seuratObj_diet$Annotation_paper_detailed)[5]<-"Immature_cDC1s"
levels(seuratObj_diet$Annotation_paper_detailed)[12]<-"Other_cDC1s_1"
levels(seuratObj_diet$Annotation_paper_detailed)[16]<-"Other_cDC1s_2"

DimPlot(seuratObj_diet, reduction = "RNA_harmony_umap", label = T, repel = T, group.by = "Condition", label.size = 1,
        cols =brewer.pal(10, "Paired")[c(2:4,9:10,5:8)]) + 
DimPlot(seuratObj_diet, reduction = "RNA_harmony_umap", label = T, repel = T, group.by = "Sample", label.size = 1,
          cols =brewer.pal(10, "Paired")[c(2,4,6,8,10,3,5,7,9)])  + NoLegend()

Colset<-brewer.pal(10, "Paired")[c(2:4,9:10,5:8)]
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

All<-c(levels(as.factor(seuratObj_diet$Annotation_paper)),
       levels(as.factor(seuratObj_diet$Condition)),levels(as.factor(seuratObj_diet$Annotation_paper_split_by_condition)),
       levels(as.factor(seuratObj_diet$Annotation_paper_detailed)),levels(as.factor(seuratObj_diet$Sample)),
       levels(as.factor(seuratObj_diet$Replicate_info)))
Color_info<-c(gg_color_hue(6),brewer.pal(10, "Paired")[c(2:4,9:10,5:8)],gg_color_hue(46),gg_color_hue(18),brewer.pal(10, "Paired")[c(2,4,6,8,10,3,5,7,9)],gg_color_hue(35))
Metadata_column<-c(rep("Annotation_paper",6),rep("Condition",9),rep("Annotation_paper_split_by_condition",46),rep("Annotation_paper_detailed",18),
                   rep("Sample",9),rep("Replicate_info",35))
Info_Kevin<-as.data.frame(cbind(All,Color_info,Metadata_column))

write.xlsx(Info_Kevin, file =paste0(output.dir, "Annotation/Info_Kevin_",experiment,".xlsx"))

########################

##### Save object
saveRDS(seuratObj_diet, file=paste0(output.dir,"Robjects/seuratObj_paper_diet",experiment,"_2024.rds"))

##### Read object
seuratObj_diet<-readRDS(file=paste0(output.dir,"Robjects/seuratObj_paper_diet",experiment,"_2024.rds"))