library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library('dorothea')
library('tidyr')
library('pheatmap')
library('tibble')
library('openxlsx')

library("future")
plan("multiprocess", workers = 6)

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

## Dorothea TF list
dir.create(paste0(sampleFolder,"results/Dorothea/"))
mm<-dorothea_mm

## Dorothea analysis
# Holland et al. (2020) showed that clustering the cells based on their TF activity profiles can also be very interesting. 
# Indeed, clustering cells using TF activity computed with VIPER and DoRothEA performs better than using the expression level of the same TFs. 
# In addition, it brings complementary information to the clusters based on transcriptomics profiles.
# 
# Here, we first run VIPER on DoRothEA’s regulons to obtain TFs activity, by using the wrapper function run_viper(). 
# This function can deal with different input types such as matrix, dataframe, ExpressionSet or even Seurat objects. 
# In case of a seurat object the function returns the same seurat object with an additonal assay called dorothea 
# containing the TF activities in the slot data.

## We read Dorothea Regulons for Human:
dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_mouse %>%
  dplyr::filter(confidence %in% c("A","B","C"))

mm_Rbpj<-mm[which(mm$tf == "Rbpj"),]
mm_Ahr<-mm[which(mm$tf == "Ahr"),]

write.xlsx(mm_Rbpj,paste0(sampleFolder,"results/Dorothea/mm_Dorothea_RBPJ_",sampleName,".xlsx"))
write.xlsx(mm_Ahr,paste0(sampleFolder,"results/Dorothea/mm_Dorothea_AHR_",sampleName,".xlsx"))

# ## Downsample seuratObj (otherwise memory shortage)
Idents(seuratObj)<-seuratObj$Condition
seuratObj.small <- subset(seuratObj, downsample = 2000)
rm(seuratObj)
gc()

## We compute Viper Scores 
seuratObj.small <- run_viper(seuratObj.small, regulon, assay_key = "RNA",
                       options = list(method = "scale", minsize = 4, 
                                      eset.filter = FALSE, cores = 4, 
                                      verbose = FALSE))

# We then apply Seurat to cluster the cells following the same protocol than above but using TF activity scores.

## We compute the Nearest Neighbours to perform cluster
DefaultAssay(object = seuratObj.small) <- "dorothea"
seuratObj.small <- ScaleData(seuratObj.small)
seuratObj.small <- RunPCA(seuratObj.small, features = rownames(seuratObj.small), verbose = FALSE)
seuratObj.small <- FindNeighbors(seuratObj.small, dims = 1:20, verbose = FALSE)
seuratObj.small <- FindClusters(seuratObj.small, resolution = 0.8, verbose = FALSE)

seuratObj.small <- RunUMAP(seuratObj.small, dims = 1:20, umap.method = "uwot", metric = "cosine")

## Save clusters
seuratObj.small@meta.data$dorothea_clusters<-seuratObj.small@meta.data$dorothea_snn_res.0.8

## Get new annotation
seuratObj.small$annotated_clusters_Muscat<-seuratObj.small$annotated_clusters
levels(seuratObj.small$annotated_clusters_Muscat)<-c('Early mature cDC1s',"Late mature cDC1s","Early mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Late mature cDC1s",
                                               "Proliferating cDC1s","Proliferating cDC1s","Late mature cDC1s",
                                               "Late mature cDC1s","Immature cDC1s","Other cDC1s","Late mature cDC1s",
                                               "Early mature cDC1s","Late mature cDC1s","Doublets","Late mature cDC1s","Early mature cDC1s")

seuratObj.small$annotated_clusters_Muscat_v2<-paste0(seuratObj.small$orig.ident,"_",seuratObj.small$annotated_clusters_Muscat)

DimPlot(seuratObj.small, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat")
DimPlot(seuratObj.small, reduction = "umap", group.by = "annotated_clusters_Muscat")
DimPlot(seuratObj.small, reduction = "umap", group.by = "orig.ident")
DimPlot(seuratObj.small, reduction = "umap", group.by = "annotated_clusters_Muscat_v2")
DimPlot(seuratObj.small, reduction = "umap", group.by = "dorothea_clusters")

## Save object
saveRDS(seuratObj.small, paste0(sampleFolder,"results/Robjects/seuratObj_small_Dorothea_",sampleName,".rds"))

## Marker analysis 
## New dorothea clusters
seuratObj.markers <- FindAllMarkers(seuratObj.small, only.pos = TRUE, min.pct = 0.10, 
                                    logfc.threshold = 0.25, verbose = FALSE)
table(seuratObj.markers$cluster)
saveRDS(seuratObj.markers, file=paste0(sampleFolder,"results/Robjects/Dorothea_markersList_Dorotheaclus_",sampleName,".rds"))

### Create list with markers
totalNrDorotheaclusters_Dorotheaclus<-max(as.numeric(names(table(seuratObj.markers$cluster))))
totalNrDorotheaclusters_DorotheaclusPlusOne<-totalNrDorotheaclusters_Dorotheaclus+1
DorotheamarkersList_Dorotheaclus<-list()

for(i in 1:totalNrDorotheaclusters_DorotheaclusPlusOne){
  clusterNr<-i-1
  
  tmp<-seuratObj.markers[seuratObj.markers$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  DorotheamarkersList_Dorotheaclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(DorotheamarkersList_Dorotheaclus)<-paste0("Dorotheacluster",0:totalNrDorotheaclusters_Dorotheaclus)

### Write to Excel
write.xlsx(DorotheamarkersList_Dorotheaclus, file =paste0(sampleFolder, "results/Dorothea/DorotheamarkersList_Dorotheaclus_",sampleName,".xlsx"))

## Create plots
U1<-DimPlot(seuratObj.small, reduction = "umap", group.by = "dorothea_clusters", label = TRUE, label.size = 8, pt.size = 0.5)
U2<-DimPlot(seuratObj.small, reduction = "umap", group.by = "orig.ident", label = TRUE, label.size = 5, pt.size = 0.5)
U3<-DimPlot(seuratObj.small, reduction = "umap", group.by = "annotated_clusters_Muscat", label = TRUE, label.size = 5, pt.size = 0.5)
U4<-DimPlot(seuratObj.small, reduction = "umap", group.by = "annotated_clusters_Muscat_v2", label = TRUE, label.size = 5, pt.size = 0.5)
U5<-DimPlot(seuratObj.small, reduction = "RNA_harmony_umap", group.by = "dorothea_clusters", label = TRUE, label.size = 8, pt.size = 0.5)
U6<-DimPlot(seuratObj.small, reduction = "RNA_harmony_umap", group.by = "orig.ident", label = TRUE, label.size = 5, pt.size = 0.5)
U7<-DimPlot(seuratObj.small, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat", label = TRUE, label.size = 8, pt.size = 0.5)
U8<-DimPlot(seuratObj.small, reduction = "RNA_harmony_umap", group.by = "annotated_clusters_Muscat_v2", label = TRUE, label.size = 5, pt.size = 0.5)

DefaultAssay(object = seuratObj.small) <- "dorothea"
F1<-FeaturePlot(object = seuratObj.small, features = c("Srebf1", "Srebf2","Nr1h2","Nr1h3"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
F2<-FeaturePlot(object = seuratObj.small, features = c("Srebf1", "Srebf2","Nr1h2","Nr1h3"), cols = c("grey", "blue"), 
                reduction = "RNA_harmony_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
F3<-FeaturePlot(object = seuratObj.small, features = "Nr1h2", cols = c("yellow", "red"), 
                reduction = "RNA_harmony_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
F4<-FeaturePlot(object = seuratObj.small, features = c("Fos", "Jun","Junb"), cols = c("yellow", "red"), 
                reduction = "RNA_harmony_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T)
Idents(seuratObj.small)<-paste0(seuratObj.small$Condition, "_", seuratObj.small$annotated_clusters_Muscat)
F5<-FeaturePlot(object = seuratObj.small, features = c("Fos", "Jun","Junb"), cols = c("yellow", "red"), 
                reduction = "RNA_umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5, order=T, label = T, repel = T, label.size = 2)

## Save plots
pdf(file=paste0(sampleFolder,"results/Dorothea/Dim_and_feature_plot_Dorothea_",sampleName,".pdf"), height = 14, width = 20)
U1
U2
U3
U4
F1
U5
U6
U7
U8
F2
dev.off()

pdf(file=paste0(sampleFolder,"results/Dorothea/Nr1h2_Dorothea_TF_activity_",sampleName,".pdf"), height = 7, width = 10)
F3
dev.off()

pdf(file=paste0(sampleFolder,"results/Dorothea/ERK_Dorothea_TF_activity_",sampleName,".pdf"), height = 13, width = 20)
F4
F5
dev.off()

##################################################################

## Read object
seuratObj.small <- readRDS(paste0(sampleFolder,"results/Robjects/seuratObj_small_Dorothea_",sampleName,".rds"))

## Subset for orig clusters
Idents(seuratObj.small)<-seuratObj.small@meta.data$annotated_clusters_Muscat_v2

seuratObj.small_subset<-subset(seuratObj.small, idents = levels(Idents(seuratObj.small))[c(1,2,6,11,16,23,
                                                                                           28,32,36,40)])
## Final RNA annotation
seuratObj.markers.2 <- FindAllMarkers(seuratObj.small_subset, only.pos = TRUE, min.pct = 0.10, 
                                      logfc.threshold = 0.25, verbose = FALSE, assay = "dorothea")
table(seuratObj.markers.2$cluster)
saveRDS(seuratObj.markers.2, file=paste0(sampleFolder,"results/Robjects/Dorothea_markersList_RNAclus_",sampleName,".rds"))

### Create list with markers
totalNrDorotheaclusters_RNAclus<-names(table(seuratObj.markers.2$cluster))
DorotheamarkersList_RNAclus<-list()

for(i in totalNrDorotheaclusters_RNAclus){
  tmp<-seuratObj.markers.2[seuratObj.markers.2$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_log2FC
  
  DorotheamarkersList_RNAclus[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
write.xlsx(DorotheamarkersList_RNAclus, file =paste0(sampleFolder, "results/Dorothea/DorotheamarkersList_RNAclus_",sampleName,".xlsx"))

##############################

## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df_subset <- GetAssayData(seuratObj.small_subset, slot = "scale.data", 
                                       assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()

viper_scores_df <- GetAssayData(seuratObj.small, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()

## We create a data frame containing the cells and their clusters
CellsOrigClusters <- data.frame(cell = names(Idents(seuratObj.small_subset)), 
                                cell_type = as.character(Idents(seuratObj.small_subset)),
                                check.names = F)

CellsDorotheaClusters <- data.frame(cell = rownames(seuratObj.small@meta.data), 
                                    cell_type = as.character(seuratObj.small@meta.data$dorothea_clusters),
                                    check.names = F)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_Origclusters <- viper_scores_df_subset  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsOrigClusters)

viper_scores_Dorotheaclusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsDorotheaClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_Origclusters_scores <- viper_scores_Origclusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

summarized_viper_Dorotheaclusters_scores <- viper_scores_Dorotheaclusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

# For visualization purposes, we select the 20 most variable TFs across clusters according to our scores.
## We select the 20 most variable TFs. (20*6 populations = 120)
highly_variable_tfs_Origclusters <- summarized_viper_Origclusters_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(500, -var) %>% #20000 for all TFs and -var for least variable
  distinct(tf)

highly_variable_tfs_Dorotheaclusters <- summarized_viper_Dorotheaclusters_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(500, var) %>%
  distinct(tf)

## We prepare the data for the plot
summarized_viper_Origclusters_scores_df <- summarized_viper_Origclusters_scores %>%
  semi_join(highly_variable_tfs_Origclusters, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 

palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

## Order on expr Late mature cDC1s
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df))
order_tfs<-order_tfs[order(order_tfs$`VBO004_Late mature cDC1s`, decreasing = T),]

my_breaks_Origclusters <- c(seq(min(summarized_viper_Origclusters_scores_df), 0, 
                                length.out=ceiling(palette_length/2) + 1),
                            seq(max(summarized_viper_Origclusters_scores_df)/palette_length, 
                                max(summarized_viper_Origclusters_scores_df), 
                                length.out=floor(palette_length/2)))

viper_hmap_Origclusters <- pheatmap(order_tfs,fontsize=14, #t(summarized_viper_Origclusters_scores_df)
                                    fontsize_row = 10, cluster_cols = F, cluster_rows = F,
                                    color=my_color, breaks = my_breaks_Origclusters, 
                                    main = "DoRothEA (Orig clusters)", angle_col = 45, 
                                    treeheight_col = 0,  border_color = NA) 

summarized_viper_Dorotheaclusters_scores_df <- summarized_viper_Dorotheaclusters_scores %>%
  semi_join(highly_variable_tfs_Dorotheaclusters, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks_Dorothea_clusters <- c(seq(min(summarized_viper_Dorotheaclusters_scores_df), 0, 
                                     length.out=ceiling(palette_length/2) + 1),
                                 seq(max(summarized_viper_Dorotheaclusters_scores_df)/palette_length, 
                                     max(summarized_viper_Dorotheaclusters_scores_df), 
                                     length.out=floor(palette_length/2)))

viper_hmap_Dorotheaclusters <- pheatmap(t(summarized_viper_Dorotheaclusters_scores_df),fontsize=14, 
                                        fontsize_row = 10,
                                        color=my_color, breaks = my_breaks_Dorothea_clusters, 
                                        main = "DoRothEA (Dorothea clusters)", angle_col = 45,
                                        treeheight_col = 0,  border_color = NA) 

pdf(file=paste0(sampleFolder,"results/Dorothea/Viper_heatmap_Orig_clusters_",sampleName,".pdf"), height = 15, width = 20)
viper_hmap_Origclusters
dev.off()

pdf(file=paste0(sampleFolder,"results/Dorothea/Viper_heatmap_Dorothea_clusters_",sampleName,".pdf"), height = 7, width = 10)
viper_hmap_Dorotheaclusters
dev.off()

##################################################################

## Extra analysis 2023: check all TFs
## Read object
seuratObj.small <- readRDS(paste0(sampleFolder,"results/Robjects/seuratObj_small_Dorothea_",sampleName,".rds"))

## Subset for orig clusters
Idents(seuratObj.small)<-seuratObj.small@meta.data$annotated_clusters_Muscat_v2

seuratObj.small_subset<-subset(seuratObj.small, idents = levels(Idents(seuratObj.small))[c(1,2,4,6,11,16,23,
                                                                                           28,32,36,40)]) #Also include SS immature cDC1s: 4

## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df_subset <- GetAssayData(seuratObj.small_subset, slot = "scale.data", #data for non-scaled
                                       assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()

## We create a data frame containing the cells and their clusters
CellsOrigClusters <- data.frame(cell = names(Idents(seuratObj.small_subset)), 
                                cell_type = as.character(Idents(seuratObj.small_subset)),
                                check.names = F)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_Origclusters <- viper_scores_df_subset  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsOrigClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_Origclusters_scores <- viper_scores_Origclusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

# For visualization purposes, we select the 20 most variable TFs across clusters according to our scores.
# Extra 2023: all TFs (check common ones)
highly_variable_tfs_Origclusters <- summarized_viper_Origclusters_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(20000, var) %>% #20000 for all TFs (and -var for least variable order)
  distinct(tf)

## We prepare the data for the plot
summarized_viper_Origclusters_scores_df <- summarized_viper_Origclusters_scores %>%
  semi_join(highly_variable_tfs_Origclusters, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 

# Change names for full list (2023)
rownames(summarized_viper_Origclusters_scores_df)<-c("SS_EM_cDC1s","SS_Immature_cDC1s","SS_LM_cDC1s","eLNPs_8h_LM_cDC1s","pIC_LNPs_8h_LM_cDC1s","CpG_LNPs_8h_LM_cDC1s","pIC_8h_LM_cDC1s",
                                                     "eLNPs_2h_EM_cDC1s","pIC_LNPs_2h_EM_cDC1s","CpG_LNPs_2h_EM_cDC1s","pIC_2h_EM_cDC1s")
summarized_viper_Origclusters_scores_df<-summarized_viper_Origclusters_scores_df[c(2,1,3,8,4,11,7,9,5,10,6),]

palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

## Order on expr Late mature cDC1s
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df))

my_breaks_Origclusters <- c(seq(min(summarized_viper_Origclusters_scores_df), 0, 
                                length.out=ceiling(palette_length/2) + 1),
                            seq(max(summarized_viper_Origclusters_scores_df)/palette_length, 
                                max(summarized_viper_Origclusters_scores_df), 
                                length.out=floor(palette_length/2)))

viper_hmap_Origclusters <- pheatmap(order_tfs,fontsize=14, #t(summarized_viper_Origclusters_scores_df)
                                    fontsize_row = 10, cluster_cols = F, cluster_rows = F,
                                    color=my_color, breaks = my_breaks_Origclusters, 
                                    main = "DoRothEA (Orig clusters)", angle_col = 90, #45
                                    treeheight_col = 0,  border_color = NA) 

pdf(file=paste0(sampleFolder,"results/Dorothea/Viper_heatmap_Orig_clusters_with_SS_imm_all_TFs_",sampleName,".pdf"), height = 45, width = 20) #273 TFs
viper_hmap_Origclusters
dev.off()

pdf(file=paste0(sampleFolder,"results/Dorothea/Viper_heatmap_Orig_clusters_with_SS_imm_all_TFs_not_scaled_",sampleName,".pdf"), height = 45, width = 20) #273 TFs
viper_hmap_Origclusters
dev.off()

##################################################################

## Extra analysis 2023: check patterned TFs
## Read object
seuratObj.small <- readRDS(paste0(sampleFolder,"results/Robjects/seuratObj_small_Dorothea_",sampleName,".rds"))

## Subset for orig clusters
Idents(seuratObj.small)<-seuratObj.small@meta.data$annotated_clusters_Muscat_v2

seuratObj.small_subset<-subset(seuratObj.small, idents = levels(Idents(seuratObj.small))[c(1,2,4,6,11,16,23,
                                                                                           28,32,36,40)]) #Also include SS immature cDC1s: 4

## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df_subset <- GetAssayData(seuratObj.small_subset, slot = "scale.data", #data for non-scaled
                                       assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()

## We create a data frame containing the cells and their clusters
CellsOrigClusters <- data.frame(cell = names(Idents(seuratObj.small_subset)), 
                                cell_type = as.character(Idents(seuratObj.small_subset)),
                                check.names = F)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_Origclusters <- viper_scores_df_subset  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsOrigClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_Origclusters_scores <- viper_scores_Origclusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

# For visualization purposes, we select the 20 most variable TFs across clusters according to our scores.
highly_variable_tfs_Origclusters <- summarized_viper_Origclusters_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(20000, var) %>% #20000 for all TFs (and -var for least variable order)
  distinct(tf)

## We prepare the data for the plot
summarized_viper_Origclusters_scores_df <- summarized_viper_Origclusters_scores %>%
  semi_join(highly_variable_tfs_Origclusters, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 

# Change names for full list (2023)
rownames(summarized_viper_Origclusters_scores_df)<-c("SS_EM_cDC1s","SS_Immature_cDC1s","SS_LM_cDC1s","eLNPs_8h_LM_cDC1s","pIC_LNPs_8h_LM_cDC1s","CpG_LNPs_8h_LM_cDC1s","pIC_8h_LM_cDC1s",
                                                     "eLNPs_2h_EM_cDC1s","pIC_LNPs_2h_EM_cDC1s","CpG_LNPs_2h_EM_cDC1s","pIC_2h_EM_cDC1s")
summarized_viper_Origclusters_scores_df<-summarized_viper_Origclusters_scores_df[c(2,1,3,8,4,11,7,9,5,10,6),]

# Extra 2023: patterned TFs
summarized_viper_Origclusters_scores_df_trans<-as.data.frame(t(summarized_viper_Origclusters_scores_df))
summarized_viper_Origclusters_scores_df_trans$TF <- rownames(summarized_viper_Origclusters_scores_df_trans)
Immuno_TFs_ON_LM = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s > SS_LM_cDC1s + 0.5 & pIC_LNPs_8h_LM_cDC1s > eLNPs_8h_LM_cDC1s + 0.5 & 
                                                                            CpG_LNPs_8h_LM_cDC1s > SS_LM_cDC1s + 0.5 & CpG_LNPs_8h_LM_cDC1s > eLNPs_8h_LM_cDC1s + 0.5) %>% pull(TF)
Homeo_TFs_ON_LM = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s < SS_LM_cDC1s - 0.5 & pIC_LNPs_8h_LM_cDC1s < eLNPs_8h_LM_cDC1s - 0.5 & 
                                                                           CpG_LNPs_8h_LM_cDC1s < SS_LM_cDC1s - 0.5 & CpG_LNPs_8h_LM_cDC1s < eLNPs_8h_LM_cDC1s - 0.5) %>% pull(TF)
Common_TFs_ON_LM_vs_SS_Imm = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 0.5 & SS_LM_cDC1s > SS_Immature_cDC1s + 0.5 & eLNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 0.5 & 
                                                                                   CpG_LNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 0.5) %>% pull(TF)
Common_TFs_OFF_LM_vs_SS_Imm = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s < SS_Immature_cDC1s - 0.5 & SS_LM_cDC1s < SS_Immature_cDC1s - 0.5 & eLNPs_8h_LM_cDC1s < SS_Immature_cDC1s - 0.5 & 
                                                                                     CpG_LNPs_8h_LM_cDC1s < SS_Immature_cDC1s - 0.5) %>% pull(TF)
Common_TFs_ON_EM_vs_SS_Imm = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_2h_EM_cDC1s > SS_Immature_cDC1s + 0.5 & SS_EM_cDC1s > SS_Immature_cDC1s + 0.5 & eLNPs_2h_EM_cDC1s > SS_Immature_cDC1s + 0.5 & 
                                                                                        CpG_LNPs_2h_EM_cDC1s > SS_Immature_cDC1s + 0.5) %>% pull(TF)
Common_TFs_OFF_EM_vs_SS_Imm = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_2h_EM_cDC1s < SS_Immature_cDC1s - 0.5 & SS_EM_cDC1s < SS_Immature_cDC1s - 0.5 & eLNPs_2h_EM_cDC1s < SS_Immature_cDC1s - 0.5 & 
                                                                                         CpG_LNPs_2h_EM_cDC1s < SS_Immature_cDC1s - 0.5) %>% pull(TF)
Common_TFs_ON_LM_vs_own_EM = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s > pIC_LNPs_2h_EM_cDC1s + 0.5 & SS_LM_cDC1s > SS_EM_cDC1s + 0.5 & eLNPs_8h_LM_cDC1s > eLNPs_2h_EM_cDC1s + 0.5 & 
                                                                                     CpG_LNPs_8h_LM_cDC1s > CpG_LNPs_2h_EM_cDC1s + 0.5) %>% pull(TF)
Common_TFs_OFF_LM_vs_own_EM = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s < pIC_LNPs_2h_EM_cDC1s - 0.5 & SS_LM_cDC1s < SS_EM_cDC1s - 0.5 & eLNPs_8h_LM_cDC1s < eLNPs_2h_EM_cDC1s - 0.5 & 
                                                                                        CpG_LNPs_8h_LM_cDC1s < CpG_LNPs_2h_EM_cDC1s - 0.5) %>% pull(TF)
# Extra August 2023 for Victor: stricter cutoffs
Immuno_TFs_ON_LM_strict = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s > SS_LM_cDC1s + 1 & pIC_LNPs_8h_LM_cDC1s > eLNPs_8h_LM_cDC1s + 1 & 
                                                                              CpG_LNPs_8h_LM_cDC1s > SS_LM_cDC1s + 1 & CpG_LNPs_8h_LM_cDC1s > eLNPs_8h_LM_cDC1s + 1) %>% pull(TF)
Homeo_TFs_ON_LM_strict = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s < SS_LM_cDC1s - 1 & pIC_LNPs_8h_LM_cDC1s < eLNPs_8h_LM_cDC1s - 1 & 
                                                                             CpG_LNPs_8h_LM_cDC1s < SS_LM_cDC1s - 1 & CpG_LNPs_8h_LM_cDC1s < eLNPs_8h_LM_cDC1s - 1) %>% pull(TF)
Common_TFs_ON_LM_vs_SS_Imm_strict = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 1 & SS_LM_cDC1s > SS_Immature_cDC1s + 1 & eLNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 1 & 
                                                                                        CpG_LNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 1) %>% pull(TF)
Common_TFs_OFF_LM_vs_SS_Imm_strict = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s < SS_Immature_cDC1s - 1 & SS_LM_cDC1s < SS_Immature_cDC1s - 1 & eLNPs_8h_LM_cDC1s < SS_Immature_cDC1s - 1 & 
                                                                                         CpG_LNPs_8h_LM_cDC1s < SS_Immature_cDC1s - 1) %>% pull(TF)

# Extra August 2023 for Victor: common TFs EM and LM without SS LM
Common_TFs_ON_LM_vs_SS_Imm_withoutSS = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 0.5 & eLNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 0.5 & 
                                                                                        CpG_LNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 0.5) %>% pull(TF)
Common_TFs_OFF_LM_vs_SS_Imm_withoutSS = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s < SS_Immature_cDC1s - 0.5 & eLNPs_8h_LM_cDC1s < SS_Immature_cDC1s - 0.5 & 
                                                                                         CpG_LNPs_8h_LM_cDC1s < SS_Immature_cDC1s - 0.5) %>% pull(TF)
Common_TFs_ON_EM_vs_SS_Imm_withoutSS = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_2h_EM_cDC1s > SS_Immature_cDC1s + 0.5 & eLNPs_2h_EM_cDC1s > SS_Immature_cDC1s + 0.5 & 
                                                                                        CpG_LNPs_2h_EM_cDC1s > SS_Immature_cDC1s + 0.5) %>% pull(TF)
Common_TFs_OFF_EM_vs_SS_Imm_withoutSS = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_2h_EM_cDC1s < SS_Immature_cDC1s - 0.5 & eLNPs_2h_EM_cDC1s < SS_Immature_cDC1s - 0.5 & 
                                                                                         CpG_LNPs_2h_EM_cDC1s < SS_Immature_cDC1s - 0.5) %>% pull(TF)

# Extra Feb 2024 paper: extra rule for common TFs
Immuno_TFs_ON_LM_paper = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s > SS_LM_cDC1s + 0.5 & pIC_LNPs_8h_LM_cDC1s > eLNPs_8h_LM_cDC1s + 0.5 & 
                                                                                  CpG_LNPs_8h_LM_cDC1s > SS_LM_cDC1s + 0.5 & CpG_LNPs_8h_LM_cDC1s > eLNPs_8h_LM_cDC1s + 0.5 &
                                                                                  pIC_LNPs_8h_LM_cDC1s  > 0 & CpG_LNPs_8h_LM_cDC1s > 0) %>% pull(TF)
Homeo_TFs_ON_LM_paper = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s < SS_LM_cDC1s - 0.5 & pIC_LNPs_8h_LM_cDC1s < eLNPs_8h_LM_cDC1s - 0.5 & 
                                                                                 CpG_LNPs_8h_LM_cDC1s < SS_LM_cDC1s - 0.5 & CpG_LNPs_8h_LM_cDC1s < eLNPs_8h_LM_cDC1s - 0.5 &
                                                                                 SS_LM_cDC1s  > 0 & eLNPs_8h_LM_cDC1s > 0) %>% pull(TF)
Common_TFs_ON_LM_vs_SS_Imm_paper = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 0.5 & SS_LM_cDC1s > SS_Immature_cDC1s + 0.5 & eLNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 0.5 & 
                                                                                            CpG_LNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 0.5 & SS_LM_cDC1s  > 0 & eLNPs_8h_LM_cDC1s > 0 &
                                                                                            pIC_LNPs_8h_LM_cDC1s  > 0 & CpG_LNPs_8h_LM_cDC1s > 0) %>% pull(TF)
Common_TFs_OFF_LM_vs_SS_Imm_paper = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s < SS_Immature_cDC1s - 0.5 & SS_LM_cDC1s < SS_Immature_cDC1s - 0.5 & eLNPs_8h_LM_cDC1s < SS_Immature_cDC1s - 0.5 & 
                                                                                                CpG_LNPs_8h_LM_cDC1s < SS_Immature_cDC1s - 0.5 & SS_LM_cDC1s  < 0 & eLNPs_8h_LM_cDC1s < 0 &
                                                                                               pIC_LNPs_8h_LM_cDC1s < 0 & CpG_LNPs_8h_LM_cDC1s < 0) %>% pull(TF)
Homeo_TFs_ON_LM_vs_SS_Imm_paper = summarized_viper_Origclusters_scores_df_trans %>% filter(SS_LM_cDC1s > SS_Immature_cDC1s + 0.5 & eLNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 0.5 & #pIC_LNPs_8h_LM_cDC1s < SS_Immature_cDC1s + 0.5 & 
                                                                                           pIC_LNPs_8h_LM_cDC1s  < 0 & CpG_LNPs_8h_LM_cDC1s < 0 & #CpG_LNPs_8h_LM_cDC1s < SS_Immature_cDC1s + 0.5 &
                                                                                           SS_LM_cDC1s  > 0 & eLNPs_8h_LM_cDC1s > 0) %>% pull(TF)
Immuno_TFs_ON_LM_vs_SS_Imm_paper = summarized_viper_Origclusters_scores_df_trans %>% filter(pIC_LNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 0.5 & #SS_LM_cDC1s < SS_Immature_cDC1s + 0.5 & eLNPs_8h_LM_cDC1s < SS_Immature_cDC1s + 0.5 & 
                                                                                             CpG_LNPs_8h_LM_cDC1s > SS_Immature_cDC1s + 0.5 & SS_LM_cDC1s  < 0 & eLNPs_8h_LM_cDC1s < 0 &
                                                                                             pIC_LNPs_8h_LM_cDC1s  > 0 & CpG_LNPs_8h_LM_cDC1s > 0) %>% pull(TF)


## Filter on TF lists
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Immuno_TFs_ON_LM]))
analysis<-"Immuno_TFs_ON_LM"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Homeo_TFs_ON_LM]))
analysis<-"Homeo_TFs_ON_LM"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_ON_LM_vs_SS_Imm]))
analysis<-"Common_TFs_ON_LM_vs_SS_Imm"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_OFF_LM_vs_SS_Imm]))
analysis<-"Common_TFs_OFF_LM_vs_SS_Imm"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_ON_EM_vs_SS_Imm]))
analysis<-"Common_TFs_ON_EM_vs_SS_Imm"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_OFF_EM_vs_SS_Imm, drop = F]))
analysis<-"Common_TFs_OFF_EM_vs_SS_Imm"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_ON_LM_vs_own_EM]))
analysis<-"Common_TFs_ON_LM_vs_own_EM"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_OFF_LM_vs_own_EM]))
analysis<-"Common_TFs_OFF_LM_vs_own_EM"

order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Immuno_TFs_ON_LM_strict]))
analysis<-"Immuno_TFs_ON_LM_strict"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Homeo_TFs_ON_LM_strict]))
analysis<-"Homeo_TFs_ON_LM_strict"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_ON_LM_vs_SS_Imm_strict]))
analysis<-"Common_TFs_ON_LM_vs_SS_Imm_strict"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_OFF_LM_vs_SS_Imm_strict]))
analysis<-"Common_TFs_OFF_LM_vs_SS_Imm_strict"

order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_ON_LM_vs_SS_Imm_withoutSS]))
analysis<-"Common_TFs_ON_LM_without_SS_vs_SS_Imm"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_OFF_LM_vs_SS_Imm_withoutSS]))
analysis<-"Common_TFs_OFF_LM_without_SS_vs_SS_Imm"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_ON_EM_vs_SS_Imm_withoutSS]))
analysis<-"Common_TFs_ON_EM_without_SS_vs_SS_Imm"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_OFF_EM_vs_SS_Imm_withoutSS]))
analysis<-"Common_TFs_OFF_EM_without_SS_vs_SS_Imm"

order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Immuno_TFs_ON_LM_paper]))
analysis<-"Immuno_TFs_ON_LM_paper"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Homeo_TFs_ON_LM_paper]))
analysis<-"Homeo_TFs_ON_LM_paper"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_ON_LM_vs_SS_Imm_paper]))
analysis<-"Common_TFs_ON_LM_vs_SS_Imm_paper"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Common_TFs_OFF_LM_vs_SS_Imm_paper]))
analysis<-"Common_TFs_OFF_LM_vs_SS_Imm_paper"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Homeo_TFs_ON_LM_vs_SS_Imm_paper]))
analysis<-"Homeo_TFs_ON_LM_vs_SS_Imm_paper"
order_tfs<-as.data.frame(t(summarized_viper_Origclusters_scores_df[,Immuno_TFs_ON_LM_vs_SS_Imm_paper]))
analysis<-"Immuno_TFs_ON_LM_vs_SS_Imm_paper"

# Heatmap
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks_Origclusters <- c(seq(min(summarized_viper_Origclusters_scores_df), 0, 
                                length.out=ceiling(palette_length/2) + 1),
                            seq(max(summarized_viper_Origclusters_scores_df)/palette_length, 
                                max(summarized_viper_Origclusters_scores_df), 
                                length.out=floor(palette_length/2)))

viper_hmap_Origclusters <- pheatmap(order_tfs,fontsize=14, #t(summarized_viper_Origclusters_scores_df)
                                    fontsize_row = 10, cluster_cols = F, cluster_rows = F,
                                    color=my_color, breaks = my_breaks_Origclusters, 
                                    main = "DoRothEA (Orig clusters)", angle_col = 90, #45
                                    treeheight_col = 0,  border_color = NA) 

pdf(file=paste0(sampleFolder,"results/Dorothea/Viper_heatmap_Orig_clusters_",analysis,"_",sampleName,".pdf"), height = 10, width = 20) #273 TFs
viper_hmap_Origclusters
dev.off()
