## Meta-analysis script
## Comparison of our LNP datasets vs other gene lists

setwd("/home/clintdn/VIB/DATA/Sophie/Meta_analysis/")

library(openxlsx)
library(stringr)
library(UpSetR)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggrepel)

###First letter upper case
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

##########################################################################################
##########################################################################################

## Load in data

## Baratin_data
## Isolate positive marker genes for mature DCs
Baratin_SS<-read.table("GSE49358.top.table_Baratin_SS_res_vs_mig.tsv",sep="\t", quote = "", fill = F, 
                       header=TRUE, stringsAsFactors=FALSE, check.names = FALSE) 
Baratin_SS_genelist<-Baratin_SS[which(Baratin_SS$adj.P.Val<0.05 & Baratin_SS$logFC < -1),"Gene.symbol"] #logFC< -1 because Res vs Mig comparison done on GEO2R
Baratin_SS_genelist<-unique(Baratin_SS_genelist[Baratin_SS_genelist != ""]) #Remove empty strings
for (i in 1:length(Baratin_SS_genelist)){ #Select first genesymbol of those separated by ///
  Baratin_SS_genelist[i]<-str_split(Baratin_SS_genelist[i],pattern = "///")[[1]][1]
}

####

## Mellman data
Mellman_12_18_36hr_vs_Control0<-read.table("GSE9241.top.table_Mellman_12hr+_CD_vs_control_0hr.tsv",sep="\t", quote = "", fill = F, 
                                           header=TRUE, stringsAsFactors=FALSE, check.names = FALSE) 
#Work with only logFC 3 cutoff like paper. Or use regular p-value and logFC 1??
Mellman_12_18_36hr_vs_Control0_genelist<-Mellman_12_18_36hr_vs_Control0[which(Mellman_12_18_36hr_vs_Control0$P.Value<0.05 & Mellman_12_18_36hr_vs_Control0$logFC < - 1),"Gene.symbol"] #logFC< -1 because Res vs Mig comparison done on GEO2R
Mellman_12_18_36hr_vs_Control0_genelist<-unique(Mellman_12_18_36hr_vs_Control0_genelist[Mellman_12_18_36hr_vs_Control0_genelist != ""]) #Remove empty strings
for (i in 1:length(Mellman_12_18_36hr_vs_Control0_genelist)){ #Select first genesymbol of those separated by ///
  Mellman_12_18_36hr_vs_Control0_genelist[i]<-str_split(Mellman_12_18_36hr_vs_Control0_genelist[i],pattern = "///")[[1]][1]
}

Mellman_12_18_36hr_vs_Control36<-read.table("GSE9241.top.table_Mellman_12hr+_CD_vs_control_36hr.tsv",sep="\t", quote = "", fill = F, 
                                            header=TRUE, stringsAsFactors=FALSE, check.names = FALSE) 
#Work with only logFC 3 cutoff like paper. Or use regular p-value and logFC 1??
Mellman_12_18_36hr_vs_Control36_genelist<-Mellman_12_18_36hr_vs_Control36[which(Mellman_12_18_36hr_vs_Control36$P.Value<0.05 & Mellman_12_18_36hr_vs_Control36$logFC > 1),"Gene.symbol"] #Correct direction on GEO2R
Mellman_12_18_36hr_vs_Control36_genelist<-unique(Mellman_12_18_36hr_vs_Control36_genelist[Mellman_12_18_36hr_vs_Control36_genelist != ""]) #Remove empty strings
for (i in 1:length(Mellman_12_18_36hr_vs_Control36_genelist)){ #Select first genesymbol of those separated by ///
  Mellman_12_18_36hr_vs_Control36_genelist[i]<-str_split(Mellman_12_18_36hr_vs_Control36_genelist[i],pattern = "///")[[1]][1]
}

####

## Cummings data
Cummings_data<-read.xlsx("Cummings_GSE85682_GFPminus_vs_GFPplus_CD103_annotated.xlsx")
Cummings_data_filtered<-Cummings_data[!is.na(Cummings_data$gene),]
#Work with lower logFC cutoff due to only 4hr post AC engulfment. Or use regular adj p-value and logFC 1?? Only 
Cummings_data_filtered[which(Cummings_data_filtered$adj.P.Val<0.05 & Cummings_data_filtered$logFC < -1),"gene"] #logFC< -1 because GFPminus vs GFPplus comparison done on GEO2R
Cummings_data_filtered_genelist<-Cummings_data_filtered[which(Cummings_data_filtered$P.Value<0.01 & Cummings_data_filtered$logFC < -1),"gene"] #logFC< -1 because GFPminus vs GFPplus comparison done on GEO2R

####

## Ardouin lists
## New way with DE analysis
Ardouin_homeo_spleen<-read.xlsx("results/Ardouin/Ardouin_summary_DEgenes_clint.xlsx", sheet = "Homeo_mat_vs_SS_imm_spleen")
Ardouin_homeo_spleen_genelist<-firstup(tolower(Ardouin_homeo_spleen[which(Ardouin_homeo_spleen$adj.P.Val < 0.05 & Ardouin_homeo_spleen$logFC > 1),"gene"])) 
Ardouin_homeo_lung<-read.xlsx("results/Ardouin/Ardouin_summary_DEgenes_clint.xlsx", sheet = "Homeo_mat_vs_SS_int_lung")
Ardouin_homeo_lung_genelist<-firstup(tolower(Ardouin_homeo_lung[which(Ardouin_homeo_lung$adj.P.Val < 0.05 & Ardouin_homeo_lung$logFC > 1),"gene"])) 
Ardouin_homeo_thymus<-read.xlsx("results/Ardouin/Ardouin_summary_DEgenes_clint.xlsx", sheet = "Homeo_mat_vs_SS_imm_thymus")
Ardouin_homeo_thymus_genelist<-firstup(tolower(Ardouin_homeo_thymus[which(Ardouin_homeo_thymus$adj.P.Val < 0.05 & Ardouin_homeo_thymus$logFC > 1),"gene"])) 
Ardouin_immuno_STAg_spleen<-read.xlsx("results/Ardouin/Ardouin_summary_DEgenes_clint.xlsx", sheet = "STAg_mat_vs_SS_imm_spleen") #Issue only 2 samples!!!! Very large list!!
Ardouin_immuno_STAg_spleen_genelist<-firstup(tolower(Ardouin_immuno_STAg_spleen[which(Ardouin_immuno_STAg_spleen$adj.P.Val < 0.05 & Ardouin_immuno_STAg_spleen$logFC > 1),"gene"])) 
Ardouin_immuno_MCMV_spleen<-read.xlsx("results/Ardouin/Ardouin_summary_DEgenes_clint.xlsx", sheet = "MCMV_mat_vs_SS_imm_spleen")
Ardouin_immuno_MCMV_spleen_genelist<-firstup(tolower(Ardouin_immuno_MCMV_spleen[which(Ardouin_immuno_MCMV_spleen$adj.P.Val < 0.05 & Ardouin_immuno_MCMV_spleen$logFC > 1),"gene"])) 
Ardouin_immuno_pIC_18hr_spleen<-read.xlsx("results/Ardouin/Ardouin_summary_DEgenes_clint.xlsx", sheet = "pIC_18hr_mat_vs_SS_imm_spleen")
Ardouin_immuno_pIC_18hr_spleen_genelist<-firstup(tolower(Ardouin_immuno_pIC_18hr_spleen[which(Ardouin_immuno_pIC_18hr_spleen$adj.P.Val < 0.05 & Ardouin_immuno_pIC_18hr_spleen$logFC > 1),"gene"])) 

####

# Maier Merad mouse list
Maier_Merad_CITEseq_lists<-read.xlsx("Maier_Merad_CITEseq_genelists.xlsx")
Maier_Merad_CITEseq_genelist<-str_split(Maier_Merad_CITEseq_lists[13,"X4"],pattern = ",")[[1]]  #Split on ,

# New 2024: Maier human list too!
Maier_Merad_CITEseq_human_genelist<-str_split(Maier_Merad_CITEseq_lists[16,"X4"],pattern = ",")[[1]]  #Split on ,

####

# Reanalyzed ImmGen data from the count table!! cDC1 and cDC2 list for lung specifically!
ImmGen_RNAseq_Lung_mig_vs_res_cDC1<-read.xlsx("results/ImmGen_DCs/summary_allgenes_clint_ImmGen_Bulk_RNA_Seq.xlsx", sheet = "Mig_cDC1_vs_Res_cDC1")
ImmGen_RNAseq_Lung_mig_vs_res_cDC1_genelist<-ImmGen_RNAseq_Lung_mig_vs_res_cDC1[which(ImmGen_RNAseq_Lung_mig_vs_res_cDC1$adj.P.Val<0.05 & ImmGen_RNAseq_Lung_mig_vs_res_cDC1$logFC>1),"gene"]

ImmGen_RNAseq_Lung_mig_vs_res_cDC2<-read.xlsx("results/ImmGen_DCs/summary_allgenes_clint_ImmGen_Bulk_RNA_Seq.xlsx", sheet = "Mig_cDC2_vs_Res_cDC2")
ImmGen_RNAseq_Lung_mig_vs_res_cDC2_genelist<-ImmGen_RNAseq_Lung_mig_vs_res_cDC2[which(ImmGen_RNAseq_Lung_mig_vs_res_cDC2$adj.P.Val<0.05 & ImmGen_RNAseq_Lung_mig_vs_res_cDC2$logFC>1),"gene"]

####

## Weckel scRNAseq lists
Weckel_Mouse_scRNAseq_cDC2_zsgreen_pos_vs_neg<-read.xlsx("results/Weckel/Markerlist_zsgreenpos_vs_neg_CD301b_cells_Mendeley.xlsx")
Weckel_Mouse_scRNAseq_cDC2_zsgreen_pos_vs_neg_genelist<-Weckel_Mouse_scRNAseq_cDC2_zsgreen_pos_vs_neg[which(Weckel_Mouse_scRNAseq_cDC2_zsgreen_pos_vs_neg$p_val_adj<0.01 & Weckel_Mouse_scRNAseq_cDC2_zsgreen_pos_vs_neg$avg_log2FC > 0.25),"gene"] #Wilcox

Weckel_Human_scRNAseq_cDC2_zsgreen_pos_vs_neg<-read.xlsx("results/Weckel/Markerlist_Foreskin_DCsubset_DC2_Zsgreen_pos_vs_neg_Mendeley.xlsx")
Weckel_Human_scRNAseq_cDC2_zsgreen_pos_vs_neg_genelist<-Weckel_Human_scRNAseq_cDC2_zsgreen_pos_vs_neg[which(Weckel_Human_scRNAseq_cDC2_zsgreen_pos_vs_neg$p_val_adj<0.01 & Weckel_Human_scRNAseq_cDC2_zsgreen_pos_vs_neg$avg_log2FC > 0.25),"gene"] #Wilcox

Weckel_Human_scRNAseq_cDC2_zsgreen_pos_maturation<-read.xlsx("results/Weckel/Markerlist_Foreskin_DCsubset_DC2_Zsgreenpos_maturation_Mendeley.xlsx")
Weckel_Human_scRNAseq_cDC2_zsgreen_pos_maturation_genelist<-Weckel_Human_scRNAseq_cDC2_zsgreen_pos_maturation[which(Weckel_Human_scRNAseq_cDC2_zsgreen_pos_maturation$p_val_adj<0.01 & Weckel_Human_scRNAseq_cDC2_zsgreen_pos_maturation$avg_log2FC > 0.25),"gene"] #Wilcox

Weckel_Human_scRNAseq_cDC2_zsgreen_neg_maturation<-read.xlsx("results/Weckel/Markerlist_Foreskin_DCsubset_DC2_Zsgreenneg_maturation_Mendeley.xlsx")
Weckel_Human_scRNAseq_cDC2_zsgreen_neg_maturation_genelist<-Weckel_Human_scRNAseq_cDC2_zsgreen_neg_maturation[which(Weckel_Human_scRNAseq_cDC2_zsgreen_neg_maturation$p_val_adj<0.01 & Weckel_Human_scRNAseq_cDC2_zsgreen_neg_maturation$avg_log2FC > 0.25),"gene"] #Wilcox

####

## Torow Peyer's Patches scRNA-seq data (only adult!)
## Switched to RNA assay for the DEA (instead of integrated!!)
# 1/ qDC1 (quiscent = immature) vs actDC1 (activated = mature) in PBS conditie = Homeostatic cDC1s
# 2/ qDC2 (quiscent = immature) vs actDC2 (activated = mature) in PBS conditie = Homeostatic cDC2s
# 3/ qDC1 van de PBS conditie vs actDC1 (activated = mature) van de R848 conditie = Immunogenic cDC1s
# 4/ qDC2 van de PBS conditie vs actDC2 (activated = mature) van de R848 conditie = Immunogenic cDC2s
Torow_Mouse_scRNAseq_actDC1_PBS_vs_qDC1_PBS<-read.xlsx("results/Torow_Peyer_patch/Maturation_markers_RNA_Torow_Peyer_patch.xlsx", sheet = "Homeostatic_cDC1s")
Torow_Mouse_scRNAseq_actDC1_PBS_vs_qDC1_PBS_genelist<-Torow_Mouse_scRNAseq_actDC1_PBS_vs_qDC1_PBS[which(Torow_Mouse_scRNAseq_actDC1_PBS_vs_qDC1_PBS$p_val_adj<0.01 & Torow_Mouse_scRNAseq_actDC1_PBS_vs_qDC1_PBS$avg_log2FC>0.25),"gene"]

Torow_Mouse_scRNAseq_actDC2_PBS_vs_qDC2_PBS<-read.xlsx("results/Torow_Peyer_patch/Maturation_markers_RNA_Torow_Peyer_patch.xlsx", sheet = "Homeostatic_cDC2s")
Torow_Mouse_scRNAseq_actDC2_PBS_vs_qDC2_PBS_genelist<-Torow_Mouse_scRNAseq_actDC2_PBS_vs_qDC2_PBS[which(Torow_Mouse_scRNAseq_actDC2_PBS_vs_qDC2_PBS$p_val_adj<0.01 & Torow_Mouse_scRNAseq_actDC2_PBS_vs_qDC2_PBS$avg_log2FC>0.25),"gene"]

Torow_Mouse_scRNAseq_actDC1_R848_vs_qDC1_PBS<-read.xlsx("results/Torow_Peyer_patch/Maturation_markers_RNA_Torow_Peyer_patch.xlsx", sheet = "Immunogenic_cDC1s")
Torow_Mouse_scRNAseq_actDC1_R848_vs_qDC1_PBS_genelist<-Torow_Mouse_scRNAseq_actDC1_R848_vs_qDC1_PBS[which(Torow_Mouse_scRNAseq_actDC1_R848_vs_qDC1_PBS$p_val_adj<0.01 & Torow_Mouse_scRNAseq_actDC1_R848_vs_qDC1_PBS$avg_log2FC>0.25),"gene"]

Torow_Mouse_scRNAseq_actDC2_R848_vs_qDC2_PBS<-read.xlsx("results/Torow_Peyer_patch/Maturation_markers_RNA_Torow_Peyer_patch.xlsx", sheet = "Immunogenic_cDC2s")
Torow_Mouse_scRNAseq_actDC2_R848_vs_qDC2_PBS_genelist<-Torow_Mouse_scRNAseq_actDC2_R848_vs_qDC2_PBS[which(Torow_Mouse_scRNAseq_actDC2_R848_vs_qDC2_PBS$p_val_adj<0.01 & Torow_Mouse_scRNAseq_actDC2_R848_vs_qDC2_PBS$avg_log2FC>0.25),"gene"]

####

## Silva-Sanchez lung data
# Bulk data
Silva_Sanchez_RNAseq_INT_vs_CD103<-read.xlsx("Raw_data/Silva-Sanchez_lung_neonatal_DCs/GSE112270_SilvaSanchez_differential_expression_180322.xlsx", sheet = "INT_vs_CD103")
Silva_Sanchez_RNAseq_INT_vs_CD103_genelist<-Silva_Sanchez_RNAseq_INT_vs_CD103[which(Silva_Sanchez_RNAseq_INT_vs_CD103$FDR<0.05 & Silva_Sanchez_RNAseq_INT_vs_CD103$logFC>1),"Gene"]

# scRNAseq
# 1/ Adult mature cDC1s vs rest adult cDC1s
# 2/ Adult mature cDC2s vs rest adult cDC2s
Silva_Sanchez_Mouse_scRNAseq_adult_mature_vs_other_cDC1s<-read.xlsx("results/Silva-Sanchez_lung_neonatal_DCs/Maturation_markers_Silva-Sanchez_neonatal_lung.xlsx", sheet = "Homeostatic_cDC1s")
Silva_Sanchez_Mouse_scRNAseq_adult_mature_vs_other_cDC1s_genelist<-Silva_Sanchez_Mouse_scRNAseq_adult_mature_vs_other_cDC1s[which(Silva_Sanchez_Mouse_scRNAseq_adult_mature_vs_other_cDC1s$p_val_adj<0.01 & Silva_Sanchez_Mouse_scRNAseq_adult_mature_vs_other_cDC1s$avg_log2FC>0.25),"gene"]

Silva_Sanchez_Mouse_scRNAseq_adult_mature_vs_other_cDC2s<-read.xlsx("results/Silva-Sanchez_lung_neonatal_DCs/Maturation_markers_Silva-Sanchez_neonatal_lung.xlsx", sheet = "Homeostatic_cDC2s")
Silva_Sanchez_Mouse_scRNAseq_adult_mature_vs_other_cDC2s_genelist<-Silva_Sanchez_Mouse_scRNAseq_adult_mature_vs_other_cDC2s[which(Silva_Sanchez_Mouse_scRNAseq_adult_mature_vs_other_cDC2s$p_val_adj<0.01 & Silva_Sanchez_Mouse_scRNAseq_adult_mature_vs_other_cDC2s$avg_log2FC>0.25),"gene"]

#####

## Brown mouse and human spleen
## Mature DCs contain a mix of cDC1s and cDC2s...
## Split up based on Irf8 expression (some overlap with S100a4 though...)
# 1/ Mouse spleen SS CCR7hi cDC1s vs cDC1
# 2/ Mouse spleen SS CCR7hi cDC2s vs c("cDC2 Tbet-","cDC2 Tbet+")
# 3/ CCR7+ cDC2 vs c("CLEC10A- cDC2","CLEC10A+ cDC2")
## Decision to leave out mouse spleen genelists: too messy to split up based on marker expression!!
# Brown_scRNAseq_Mouse_spleen_cDC1s<-read.xlsx("results/Brown_human_and_mouse_spleen_SS_and_melanoma/Maturation_markers_update_Brown_human_and_mouse_spleen_SS_and_melanoma.xlsx", sheet = "Mouse_spleen_cDC1s")
# Brown_scRNAseq_Mouse_spleen_cDC1s_genelist<-Brown_scRNAseq_Mouse_spleen_cDC1s[which(Brown_scRNAseq_Mouse_spleen_cDC1s$p_val_adj<0.01 & Brown_scRNAseq_Mouse_spleen_cDC1s$avg_log2FC>0.25),"gene"]
# 
# Brown_scRNAseq_Mouse_spleen_cDC2s<-read.xlsx("results/Brown_human_and_mouse_spleen_SS_and_melanoma/Maturation_markers_update_Brown_human_and_mouse_spleen_SS_and_melanoma.xlsx", sheet = "Mouse_spleen_cDC2s")
# Brown_scRNAseq_Mouse_spleen_cDC2s_genelist<-Brown_scRNAseq_Mouse_spleen_cDC2s[which(Brown_scRNAseq_Mouse_spleen_cDC2s$p_val_adj<0.01 & Brown_scRNAseq_Mouse_spleen_cDC2s$avg_log2FC>0.25),"gene"]

Brown_scRNAseq_Human_spleen_cDC2s<-read.xlsx("results/Brown_human_and_mouse_spleen_SS_and_melanoma/Maturation_markers_update_Brown_human_and_mouse_spleen_SS_and_melanoma.xlsx", sheet = "Human_spleen_cDC2s")
Brown_scRNAseq_Human_spleen_cDC2s_genelist<-Brown_scRNAseq_Human_spleen_cDC2s[which(Brown_scRNAseq_Human_spleen_cDC2s$p_val_adj<0.01 & Brown_scRNAseq_Human_spleen_cDC2s$avg_log2FC>0.25),"gene"]

####

## Cui cytokine paper (Immune dictionary) DC1s and mig DC1s
# 1/ seuratObj MigDC_2_PBS vs cDC1_1_PBS
# 2/ Polarization analysis paper -> cell state D for DC1s (primarily linked to TNFa cytokine)
Cui_scRNAseq_PBS_MigcDC1s_vs_RescDC1s<-read.xlsx("results/Cytokine_Cui/Maturation_markers_Cui_cytokine_mouse_LN.xlsx", sheet = "Homeostatic_cDC1s")
Cui_scRNAseq_PBS_MigcDC1s_vs_RescDC1s_genelist<-Cui_scRNAseq_PBS_MigcDC1s_vs_RescDC1s[which(Cui_scRNAseq_PBS_MigcDC1s_vs_RescDC1s$p_val_adj<0.01 & Cui_scRNAseq_PBS_MigcDC1s_vs_RescDC1s$avg_log2FC>0.25),"gene"]

# Cui_scRNAseq_TNF_cDC1_cell_state_D<-read.xlsx("Raw_data/New_meta_analysis_data_2024/Cui_cytokine_LN/Cytokine_paper_DEGs.xlsx", sheet = "cDC1")
# Cui_scRNAseq_TNF_cDC1_cell_state_D_genelist<-Cui_scRNAseq_TNF_cDC1_cell_state_D[which(Cui_scRNAseq_TNF_cDC1_cell_state_D$Polarization == "cDC1-d" & Cui_scRNAseq_TNF_cDC1_cell_state_D$P_adj<0.01 & Cui_scRNAseq_TNF_cDC1_cell_state_D$Avg_log2FC>0.25),"Gene"]

####

## Zillionis mouse and human lung cancer (seurat analysis!)
# 1a/ Human tDC3 (DC1s) vs tDC1
# 1b/ Human tDC3 (DC2s) vs tDC2
# 2a/ Mouse DC3 (DC1s) vs DC1
# 2b/ Mouse DC3 (DC2s) vs DC2

Zillionis_scRNAseq_human_cDC1_maturation<-read.xlsx("results/Zillionis_Human_and_mouse_Lung_cancer/Maturation_markers_updated_Seurat_processed_from_raw_Zillionis_lung_cancer_Seurat.xlsx", sheet = "Human_tumor_cDC1s")
Zillionis_scRNAseq_human_cDC1_maturation_genelist<-Zillionis_scRNAseq_human_cDC1_maturation[which(Zillionis_scRNAseq_human_cDC1_maturation$p_val_adj<0.01 & Zillionis_scRNAseq_human_cDC1_maturation$avg_log2FC>0.25),"gene"]

Zillionis_scRNAseq_human_cDC2_maturation<-read.xlsx("results/Zillionis_Human_and_mouse_Lung_cancer/Maturation_markers_updated_Seurat_processed_from_raw_Zillionis_lung_cancer_Seurat.xlsx", sheet = "Human_tumor_cDC2s")
Zillionis_scRNAseq_human_cDC2_maturation_genelist<-Zillionis_scRNAseq_human_cDC2_maturation[which(Zillionis_scRNAseq_human_cDC2_maturation$p_val_adj<0.01 & Zillionis_scRNAseq_human_cDC2_maturation$avg_log2FC>0.25),"gene"]

## Corrected to only tumor DCs on 30/05/24
Zillionis_scRNAseq_mouse_cDC1_maturation<-read.xlsx("results/Zillionis_Human_and_mouse_Lung_cancer/Maturation_markers_mouse_tumor_updated_Seurat_processed_from_raw_Zillionis_lung_cancer_Seurat.xlsx", sheet = "Mouse_tumor_cDC1s")
Zillionis_scRNAseq_mouse_cDC1_maturation_genelist<-Zillionis_scRNAseq_mouse_cDC1_maturation[which(Zillionis_scRNAseq_mouse_cDC1_maturation$p_val_adj<0.01 & Zillionis_scRNAseq_mouse_cDC1_maturation$avg_log2FC>0.25),"gene"]

Zillionis_scRNAseq_mouse_cDC2_maturation<-read.xlsx("results/Zillionis_Human_and_mouse_Lung_cancer/Maturation_markers_mouse_tumor_updated_Seurat_processed_from_raw_Zillionis_lung_cancer_Seurat.xlsx", sheet = "Mouse_tumor_cDC2s")
Zillionis_scRNAseq_mouse_cDC2_maturation_genelist<-Zillionis_scRNAseq_mouse_cDC2_maturation[which(Zillionis_scRNAseq_mouse_cDC2_maturation$p_val_adj<0.01 & Zillionis_scRNAseq_mouse_cDC2_maturation$avg_log2FC>0.25),"gene"]

####

## Massoni-Badosa Human tonsil paper 
# 1/ aDC1s vs DC1 mature
Massoni_Badosa_scRNAseq_human_tonsil_cDC1_maturation<-read.xlsx("results/Massoni-Badosa/Maturation_markers_Massoni-Badosa_human_tonsil.xlsx", sheet = "Homeostatic_cDC1s")
Massoni_Badosa_scRNAseq_human_tonsil_cDC1_maturation_genelist<-Massoni_Badosa_scRNAseq_human_tonsil_cDC1_maturation[which(Massoni_Badosa_scRNAseq_human_tonsil_cDC1_maturation$p_val_adj<0.01 & Massoni_Badosa_scRNAseq_human_tonsil_cDC1_maturation$avg_log2FC>0.25),"gene"]

#####

## Own lists: update LNP lists in 2024 (local adj pval -> global adj pval)!!!
Victor_RNAseq_WT_mig_vs_res<-read.xlsx("../RNA-seq_Victor/results/summary_allgenes_clint_june2020.xlsx", sheet = "WT_mig_vs_res")
Victor_RNAseq_WT_mig_vs_res_genelist<-Victor_RNAseq_WT_mig_vs_res[which(Victor_RNAseq_WT_mig_vs_res$adj.P.Val<0.05 & Victor_RNAseq_WT_mig_vs_res$logFC>1),"gene"]

Victor_CITEseq_cDC1s_LM_vs_Imm<-read.xlsx("../RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset/Marker_lists/Late_mature_vs_Immature_cDC1_SAM2and3_WT_subset_genelist.xlsx") #Wilcoxon
Victor_CITEseq_cDC1s_LM_vs_Imm_genelist<-Victor_CITEseq_cDC1s_LM_vs_Imm[which(Victor_CITEseq_cDC1s_LM_vs_Imm$p_val_adj<0.01 & Victor_CITEseq_cDC1s_LM_vs_Imm$avg_logFC>0.25),"gene"]

Victor_CITEseq_cDC2s_LM_vs_Imm<-read.xlsx("../RNA-seq_Sandra/CITEseq_Test/PROCESSED_DATA/SAM2and3_WT/results_subset2_v2/Marker_lists/Late_mature_vs_Immature_cDC2_SAM2and3_WT_subset2_v2_genelist.xlsx") #Wilcoxon
Victor_CITEseq_cDC2s_LM_vs_Imm_genelist<-Victor_CITEseq_cDC2s_LM_vs_Imm[which(Victor_CITEseq_cDC2s_LM_vs_Imm$p_val_adj<0.01 & Victor_CITEseq_cDC2s_LM_vs_Imm$avg_logFC>0.25),"gene"]

Victor_LNP_CITEseq_SS_mig_vs_res<-read.xlsx("../scRNA-seq_Victor/LNP_CITEseq_experiment/VBO_LNP_merge/results/Muscat_Immature_cross_condition/tbl_full_muscat_DESeq2_new_Steady_state_Late_mature_cDC1s-Steady_state_Immature_cDC1s_VBO_4-12.xlsx")
Victor_LNP_CITEseq_SS_mig_vs_res_genelist<-Victor_LNP_CITEseq_SS_mig_vs_res[which(Victor_LNP_CITEseq_SS_mig_vs_res$p_adj.glb<0.01 & Victor_LNP_CITEseq_SS_mig_vs_res$logFC > 1 & Victor_LNP_CITEseq_SS_mig_vs_res$baseMean > 50),"gene"]

Victor_LNP_CITEseq_eLNPs_mig_vs_SS_res<-read.xlsx("../scRNA-seq_Victor/LNP_CITEseq_experiment/VBO_LNP_merge/results/Muscat_Immature_cross_condition/tbl_full_muscat_DESeq2_new_eLNPs_8h_Late_mature_cDC1s-Steady_state_Immature_cDC1s_VBO_4-12.xlsx")
Victor_LNP_CITEseq_eLNPs_mig_vs_SS_res_genelist<-Victor_LNP_CITEseq_eLNPs_mig_vs_SS_res[which(Victor_LNP_CITEseq_eLNPs_mig_vs_SS_res$p_adj.glb<0.01 & Victor_LNP_CITEseq_eLNPs_mig_vs_SS_res$logFC > 1 & Victor_LNP_CITEseq_eLNPs_mig_vs_SS_res$baseMean > 50),"gene"]

Victor_LNP_CITEseq_pIC_LNPs_mig_vs_SS_res<-read.xlsx("../scRNA-seq_Victor/LNP_CITEseq_experiment/VBO_LNP_merge/results/Muscat_Immature_cross_condition/tbl_full_muscat_DESeq2_new_pIC_LNPs_8h_Late_mature_cDC1s-Steady_state_Immature_cDC1s_VBO_4-12.xlsx")
Victor_LNP_CITEseq_pIC_LNPs_mig_vs_SS_res_genelist<-Victor_LNP_CITEseq_pIC_LNPs_mig_vs_SS_res[which(Victor_LNP_CITEseq_pIC_LNPs_mig_vs_SS_res$p_adj.glb<0.01 & Victor_LNP_CITEseq_pIC_LNPs_mig_vs_SS_res$logFC > 1 & Victor_LNP_CITEseq_pIC_LNPs_mig_vs_SS_res$baseMean > 50),"gene"]

Victor_LNP_CITEseq_CpG_LNPs_mig_vs_SS_res<-read.xlsx("../scRNA-seq_Victor/LNP_CITEseq_experiment/VBO_LNP_merge/results/Muscat_Immature_cross_condition/tbl_full_muscat_DESeq2_new_CpG_LNPs_8h_Late_mature_cDC1s-Steady_state_Immature_cDC1s_VBO_4-12.xlsx")
Victor_LNP_CITEseq_CpG_LNPs_mig_vs_SS_res_genelist<-Victor_LNP_CITEseq_CpG_LNPs_mig_vs_SS_res[which(Victor_LNP_CITEseq_CpG_LNPs_mig_vs_SS_res$p_adj.glb<0.01 & Victor_LNP_CITEseq_CpG_LNPs_mig_vs_SS_res$logFC > 1 & Victor_LNP_CITEseq_CpG_LNPs_mig_vs_SS_res$baseMean > 50),"gene"]

Victor_LNP_CITEseq_pIC_alone_mig_vs_SS_res<-read.xlsx("../scRNA-seq_Victor/LNP_CITEseq_experiment/VBO_LNP_merge/results/Muscat_Immature_cross_condition/tbl_full_muscat_DESeq2_new_pIC_alone_8h_Late_mature_cDC1s-Steady_state_Immature_cDC1s_VBO_4-12.xlsx")
Victor_LNP_CITEseq_pIC_alone_mig_vs_SS_res_genelist<-Victor_LNP_CITEseq_pIC_alone_mig_vs_SS_res[which(Victor_LNP_CITEseq_pIC_alone_mig_vs_SS_res$p_adj.glb<0.01 & Victor_LNP_CITEseq_pIC_alone_mig_vs_SS_res$logFC > 1 & Victor_LNP_CITEseq_pIC_alone_mig_vs_SS_res$baseMean > 50),"gene"]

# Only include these in triwise and not in correlation heatmap!!
Victor_LNP_CITEseq_CpG_and_pIC_LNPs_mig_vs_SS_res_genelist<-unique(c(Victor_LNP_CITEseq_pIC_LNPs_mig_vs_SS_res_genelist,Victor_LNP_CITEseq_CpG_LNPs_mig_vs_SS_res_genelist)) #Take combination of all genes instead of intersect 1800 -> 2500!!
#Overlap with SS list
Victor_LNP_CITEseq_homeo_old<-setdiff(Victor_LNP_CITEseq_SS_mig_vs_res_genelist,Victor_LNP_CITEseq_CpG_and_pIC_LNPs_mig_vs_SS_res_genelist) #Only in SS
Victor_LNP_CITEseq_common_old<-intersect(Victor_LNP_CITEseq_SS_mig_vs_res_genelist,Victor_LNP_CITEseq_CpG_and_pIC_LNPs_mig_vs_SS_res_genelist) #Common between SS and CpG/pIC
Victor_LNP_CITEseq_immuno_old<-setdiff(Victor_LNP_CITEseq_CpG_and_pIC_LNPs_mig_vs_SS_res_genelist,Victor_LNP_CITEseq_SS_mig_vs_res_genelist) #Only in CpG/pIC

# Calculate top 200 genes for the three lists for axes triwise
rownames(Victor_LNP_CITEseq_SS_mig_vs_res)<-Victor_LNP_CITEseq_SS_mig_vs_res$gene
rownames(Victor_LNP_CITEseq_pIC_LNPs_mig_vs_SS_res)<-Victor_LNP_CITEseq_pIC_LNPs_mig_vs_SS_res$gene
rownames(Victor_LNP_CITEseq_CpG_LNPs_mig_vs_SS_res)<-Victor_LNP_CITEseq_CpG_LNPs_mig_vs_SS_res$gene

Homeo_filtered<-Victor_LNP_CITEseq_SS_mig_vs_res[Victor_LNP_CITEseq_homeo_old,]
Homeo_filtered<-Homeo_filtered[order(Homeo_filtered$p_adj.loc),]
Victor_LNP_CITEseq_homeo<-head(rownames(Homeo_filtered),200)

Common_filtered<-cbind(Victor_LNP_CITEseq_SS_mig_vs_res[Victor_LNP_CITEseq_common_old,],
                       Victor_LNP_CITEseq_pIC_LNPs_mig_vs_SS_res[Victor_LNP_CITEseq_common_old,],
                       Victor_LNP_CITEseq_CpG_LNPs_mig_vs_SS_res[Victor_LNP_CITEseq_common_old,])
colnames(Common_filtered)[1:10]<-paste0(colnames(Common_filtered)[1:10],"_SS")
colnames(Common_filtered)[11:20]<-paste0(colnames(Common_filtered)[11:20],"_pIC_LNP")
colnames(Common_filtered)[21:30]<-paste0(colnames(Common_filtered)[21:30],"_CpG_LNP")
Common_filtered$p_adj.loc<-(Common_filtered$p_adj.loc_SS+Common_filtered$p_adj.loc_pIC_LNP+Common_filtered$p_adj.loc_CpG_LNP)/3
Common_filtered<-Common_filtered[order(Common_filtered$p_adj.loc),]
Victor_LNP_CITEseq_common<-head(rownames(Common_filtered),200)

Immuno_filtered<-cbind(Victor_LNP_CITEseq_pIC_LNPs_mig_vs_SS_res[Victor_LNP_CITEseq_immuno_old,],
                       Victor_LNP_CITEseq_CpG_LNPs_mig_vs_SS_res[Victor_LNP_CITEseq_immuno_old,])
colnames(Immuno_filtered)[1:10]<-paste0(colnames(Immuno_filtered)[1:10],"_pIC_LNP")
colnames(Immuno_filtered)[11:20]<-paste0(colnames(Immuno_filtered)[11:20],"_CpG_LNP")
Immuno_filtered$p_adj.loc<-(Immuno_filtered$p_adj.loc_pIC_LNP+Immuno_filtered$p_adj.loc_CpG_LNP)/2
Immuno_filtered<-Immuno_filtered[order(Immuno_filtered$p_adj.loc),]
Victor_LNP_CITEseq_immuno<-head(rownames(Immuno_filtered),200)

#################################################################################################################

## Convert human to mouse symbols (BioMaRt) for human data

# Convert with BioMart
# Basic function to convert human to mouse gene names 
## Update: changed back to useMart and to archive (2023)

convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

Mellman_12_18_36hr_vs_Control0_mouse_genelist <- convertHumanGeneList(Mellman_12_18_36hr_vs_Control0_genelist)
Mellman_12_18_36hr_vs_Control36_mouse_genelist <- convertHumanGeneList(Mellman_12_18_36hr_vs_Control36_genelist)

Weckel_Human_scRNAseq_cDC2_zsgreen_pos_vs_neg_mouse_genelist<-convertHumanGeneList(Weckel_Human_scRNAseq_cDC2_zsgreen_pos_vs_neg_genelist)
Weckel_Human_scRNAseq_cDC2_zsgreen_neg_maturation_mouse_genelist<-convertHumanGeneList(Weckel_Human_scRNAseq_cDC2_zsgreen_neg_maturation_genelist)
Weckel_Human_scRNAseq_cDC2_zsgreen_pos_maturation_mouse_genelist<-convertHumanGeneList(Weckel_Human_scRNAseq_cDC2_zsgreen_pos_maturation_genelist)

# New 2024 human genelists!
Maier_Merad_CITEseq_human_converted_genelist<-convertHumanGeneList(Maier_Merad_CITEseq_human_genelist)
Brown_scRNAseq_Human_spleen_cDC2s_converted_genelist<-convertHumanGeneList(Brown_scRNAseq_Human_spleen_cDC2s_genelist)
Zillionis_scRNAseq_human_cDC1_maturation_converted_genelist<-convertHumanGeneList(Zillionis_scRNAseq_human_cDC1_maturation_genelist)
Zillionis_scRNAseq_human_cDC2_maturation_converted_genelist<-convertHumanGeneList(Zillionis_scRNAseq_human_cDC2_maturation_genelist)
Massoni_Badosa_scRNAseq_human_tonsil_cDC1_maturation_converted_genelist<-convertHumanGeneList(Massoni_Badosa_scRNAseq_human_tonsil_cDC1_maturation_genelist)

#################################################################################################################

## Calculate overlap between the lists

## Upset plot
###Modify fromList function!! https://github.com/hms-dbmi/UpSetR/issues/85
fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

###Intersection of UpsetR lists
get_intersect_members <- function (x, ...){
  require(dplyr)
  require(tibble)
  x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
  n <- names(x)
  x %>% rownames_to_column() -> x
  l <- c(...)
  a <- intersect(names(x), l)
  ar <- vector('list',length(n)+1)
  ar[[1]] <- x
  i=2
  for (item in n) {
    if (item %in% a){
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '>= 1')
        i <- i + 1
      }
    } else {
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '== 0')
        i <- i + 1
      }
    }
  }
  do.call(filter_, ar) %>% column_to_rownames() -> x
  return(x)
}

###############################
##### LNP as reference ####
###############################

## Final annotation Victor triwise (26/07/23)
Conversion_file_triwise<-read.xlsx("Paper/Info_genelists_paper.xlsx", sheet = "Triwise")
Conversion_file_triwise<-Conversion_file_triwise[,1:3]

de_gs_by_k_split_triwise <- lapply(Conversion_file_triwise$Object.annotation.in.R, function(x) {
  x <- eval(as.name(x))
})

names(de_gs_by_k_split_triwise)<-Conversion_file_triwise$Paper.annotation

## Final annotation Victor heatmap (26/07/23)
Conversion_file_heatmap<-read.xlsx("Paper/Info_genelists_paper.xlsx", sheet = "Heatmap")
Conversion_file_heatmap<-Conversion_file_heatmap[,1:3]

de_gs_by_k_split_heatmap <- lapply(Conversion_file_heatmap$Object.annotation.in.R, function(x) {
  x <- eval(as.name(x))
})

names(de_gs_by_k_split_heatmap)<-Conversion_file_heatmap$Paper.annotation

#################################################################################################################

## Save lists
saveRDS(de_gs_by_k_split_triwise,"Paper/Tibble_all_genelists_triwise_v9_with_top_200_LNP_lists.rds")
saveRDS(de_gs_by_k_split_heatmap,"Paper/Tibble_all_genelists_heatmap_v9_with_top_200_LNP_lists.rds")

## Read lists
de_gs_by_k_split_triwise<-readRDS("Paper/Tibble_all_genelists_triwise_v9_with_top_200_LNP_lists.rds")
de_gs_by_k_split_heatmap<-readRDS("Paper/Tibble_all_genelists_heatmap_v9_with_top_200_LNP_lists.rds")

## Write lists
de_gs_by_k_split_triwise_excel<-de_gs_by_k_split_triwise
names(de_gs_by_k_split_triwise_excel)<-c(1:41)
write.xlsx(de_gs_by_k_split_triwise_excel,"Paper/Tibble_all_genelists_triwise_v9_with_top_200_LNP_lists.xlsx")

## Extra paper: (full lists CpG-LNP, pIC--LNP and SS too)
de_gs_by_k_split_heatmap_excel<-de_gs_by_k_split_heatmap
names(de_gs_by_k_split_heatmap_excel)<-c(1:41)
write.xlsx(de_gs_by_k_split_heatmap_excel,"Paper/Tibble_all_genelists_heatmap_v9_with_all_full_lists.xlsx")

#################################################################################################################

## Check overlap
upset(fromList(de_gs_by_k_split_triwise), nsets = 100)
upset(fromList(de_gs_by_k_split_heatmap), nsets = 100)

## Better analysis
dir.create("results/")
dir.create("results/Upset_plots/")

analysis<-"Extra_v9_meta_analysis"

listInput_split_triwise <- de_gs_by_k_split_triwise
listInput_split_heatmap <- de_gs_by_k_split_heatmap

cairo_pdf(file=paste0("results/Upset_plots/Upset_plot_triwise_lists_with_top_200_LNP_lists_",analysis,".pdf"), width=15, height = 12)
upset(fromList(listInput_split_triwise), sets = names(de_gs_by_k_split_triwise), order.by = "freq")
dev.off()

cairo_pdf(file=paste0("results/Upset_plots/Upset_plot_heatmap_lists_with_top_200_LNP_lists_",analysis,".pdf"), width=15, height = 12)
upset(fromList(listInput_split_heatmap), sets = names(de_gs_by_k_split_heatmap), order.by = "freq")
dev.off()

## Create intersection dataframe
Intersect_df_split_triwise<-fromList(listInput_split_triwise)
Intersect_df_split_heatmap<-fromList(listInput_split_heatmap)

## Modify intersection dataframe for excel
Intersect_df_split_triwise_excel<-Intersect_df_split_triwise
Intersect_df_split_triwise_excel$Gene.symbol<-rownames(Intersect_df_split_triwise_excel)
Intersect_df_split_triwise_excel<-Intersect_df_split_triwise_excel[,c(ncol(Intersect_df_split_triwise_excel),1:(ncol(Intersect_df_split_triwise_excel)-1))] #Put genes first
write.xlsx(Intersect_df_split_triwise_excel, file = paste0("results/Upset_plots/Intersection_dataframe_triwise_lists_with_top_200_LNP_lists_",analysis,".xlsx"))

Intersect_df_split_heatmap_excel<-Intersect_df_split_heatmap
Intersect_df_split_heatmap_excel$Gene.symbol<-rownames(Intersect_df_split_heatmap_excel)
Intersect_df_split_heatmap_excel<-Intersect_df_split_heatmap_excel[,c(ncol(Intersect_df_split_heatmap_excel),1:(ncol(Intersect_df_split_heatmap_excel)-1))] #Put genes first
write.xlsx(Intersect_df_split_heatmap_excel, file = paste0("results/Upset_plots/Intersection_dataframe_heatmap_lists_with_top_200_LNP_lists_",analysis,".xlsx"))

#################################################################################################################

## Late addition to meta analysis

## Lennon bulk data
C_WT_vs_NC_WT_results<-read.xlsx("results/Lennon/summary_allgenes_clint_Lennon.xlsx", sheet = "C_WT_vs_NC_WT")
C_WT_vs_NC_WT_results_genelist<-C_WT_vs_NC_WT_results[which(C_WT_vs_NC_WT_results$adj.P.Val<0.05 & C_WT_vs_NC_WT_results$logFC>1),"gene"]
LPS_WT_vs_NC_WT_results<-read.xlsx("results/Lennon/summary_allgenes_clint_Lennon.xlsx", sheet = "LPS_WT_vs_NC_WT")
LPS_WT_vs_NC_WT_results_genelist<-LPS_WT_vs_NC_WT_results[which(LPS_WT_vs_NC_WT_results$adj.P.Val<0.05 & LPS_WT_vs_NC_WT_results$logFC>1),"gene"]

de_gs_by_k_split_triwise[[42]]<-C_WT_vs_NC_WT_results_genelist
de_gs_by_k_split_triwise[[43]]<-LPS_WT_vs_NC_WT_results_genelist

names(de_gs_by_k_split_triwise)[42]<-"Confinement-induced DC maturation (Alraies et al., 2024, Bulk RNA-seq, In-vitro BMDCs)"
names(de_gs_by_k_split_triwise)[43]<-"LPS-induced DC maturation (Alraies et al., 2024, Bulk RNA-seq, In-vitro BMDCs)"

## Also update the Zillionis genelists!!
## Changed to DC3 vs DC1/DC2 comparison. Don't split up DC3s!
Zillionis_scRNAseq_mouse_cDC1_maturation<-read.xlsx("results/Zillionis_Human_and_mouse_Lung_cancer/Maturation_markers_mouse_tumor_updated_DC3_Seurat_processed_from_raw_Zillionis_lung_cancer_Seurat.xlsx", sheet = "Mouse_tumor_DC3_vs_cDC1s")
Zillionis_scRNAseq_mouse_cDC1_maturation_genelist<-Zillionis_scRNAseq_mouse_cDC1_maturation[which(Zillionis_scRNAseq_mouse_cDC1_maturation$p_val_adj<0.01 & Zillionis_scRNAseq_mouse_cDC1_maturation$avg_log2FC>0.25),"gene"]

Zillionis_scRNAseq_mouse_cDC2_maturation<-read.xlsx("results/Zillionis_Human_and_mouse_Lung_cancer/Maturation_markers_mouse_tumor_updated_DC3_Seurat_processed_from_raw_Zillionis_lung_cancer_Seurat.xlsx", sheet = "Mouse_tumor_DC3_vs_cDC2s")
Zillionis_scRNAseq_mouse_cDC2_maturation_genelist<-Zillionis_scRNAseq_mouse_cDC2_maturation[which(Zillionis_scRNAseq_mouse_cDC2_maturation$p_val_adj<0.01 & Zillionis_scRNAseq_mouse_cDC2_maturation$avg_log2FC>0.25),"gene"]

de_gs_by_k_split_triwise[[38]]<-Zillionis_scRNAseq_mouse_cDC1_maturation_genelist
de_gs_by_k_split_triwise[[39]]<-Zillionis_scRNAseq_mouse_cDC2_maturation_genelist

## Add new nature comms paper Kim with LNP mRNA vaccine (03/10/24)
DC_mig_LNP_vs_DC_cDC1_LNP_results<-read.xlsx("/home/clintdn/VIB/DATA/Sophie/Meta_analysis/Raw_data/Kim_et_al_2024/Maturation_markers_muscle_dataset_nature_comms_LNP_mRNA_vaccine.xlsx", sheet = "DC_mig_LNP_vs_DC_cDC1_LNP")
DC_mig_LNP_vs_DC_cDC1_LNP_results_genelist<-DC_mig_LNP_vs_DC_cDC1_LNP_results[which(DC_mig_LNP_vs_DC_cDC1_LNP_results$p_val_adj<0.01 & DC_mig_LNP_vs_DC_cDC1_LNP_results$avg_log2FC>0.25),"gene"]
DC_mig_LNPmRNA_vs_DC_cDC1_LNP_results<-read.xlsx("/home/clintdn/VIB/DATA/Sophie/Meta_analysis/Raw_data/Kim_et_al_2024/Maturation_markers_muscle_dataset_nature_comms_LNP_mRNA_vaccine.xlsx", sheet = "DC_mig_LNPmRNA_vs_DC_cDC1_LNP")
DC_mig_LNPmRNA_vs_DC_cDC1_LNP_results_genelist<-DC_mig_LNPmRNA_vs_DC_cDC1_LNP_results[which(DC_mig_LNPmRNA_vs_DC_cDC1_LNP_results$p_val_adj<0.01 & DC_mig_LNPmRNA_vs_DC_cDC1_LNP_results$avg_log2FC>0.25),"gene"]
DC_mig_LNP_vs_DC_cDC2_LNP_results<-read.xlsx("/home/clintdn/VIB/DATA/Sophie/Meta_analysis/Raw_data/Kim_et_al_2024/Maturation_markers_muscle_dataset_nature_comms_LNP_mRNA_vaccine.xlsx", sheet = "DC_mig_LNP_vs_DC_cDC2_LNP")
DC_mig_LNP_vs_DC_cDC2_LNP_results_genelist<-DC_mig_LNP_vs_DC_cDC2_LNP_results[which(DC_mig_LNP_vs_DC_cDC2_LNP_results$p_val_adj<0.01 & DC_mig_LNP_vs_DC_cDC2_LNP_results$avg_log2FC>0.25),"gene"]
DC_mig_LNPmRNA_vs_DC_cDC2_LNP_results<-read.xlsx("/home/clintdn/VIB/DATA/Sophie/Meta_analysis/Raw_data/Kim_et_al_2024/Maturation_markers_muscle_dataset_nature_comms_LNP_mRNA_vaccine.xlsx", sheet = "DC_mig_LNPmRNA_vs_DC_cDC2_LNP")
DC_mig_LNPmRNA_vs_DC_cDC2_LNP_results_genelist<-DC_mig_LNPmRNA_vs_DC_cDC2_LNP_results[which(DC_mig_LNPmRNA_vs_DC_cDC2_LNP_results$p_val_adj<0.01 & DC_mig_LNPmRNA_vs_DC_cDC2_LNP_results$avg_log2FC>0.25),"gene"]

de_gs_by_k_split_triwise[[44]]<-DC_mig_LNP_vs_DC_cDC1_LNP_results_genelist
de_gs_by_k_split_triwise[[45]]<-DC_mig_LNPmRNA_vs_DC_cDC1_LNP_results_genelist
de_gs_by_k_split_triwise[[46]]<-DC_mig_LNP_vs_DC_cDC2_LNP_results_genelist
de_gs_by_k_split_triwise[[47]]<-DC_mig_LNPmRNA_vs_DC_cDC2_LNP_results_genelist

names(de_gs_by_k_split_triwise)[44]<-"LNP-induced mDCs vs cDC1 (Kim et al., 2024, scRNA-seq, Muscle)"
names(de_gs_by_k_split_triwise)[45]<-"mRNA_LNP-induced mDCs vs cDC1 (Kim et al., 2024, scRNA-seq, Muscle)"
names(de_gs_by_k_split_triwise)[46]<-"LNP-induced mDCs vs cDC2 (Kim et al., 2024, scRNA-seq, Muscle)"
names(de_gs_by_k_split_triwise)[47]<-"mRNA_LNP-induced mDCs vs cDC2 (Kim et al., 2024, scRNA-seq, Muscle)"

## Add other in-house signatures
de_gs_by_k_split_triwise[[48]]<-Victor_LNP_CITEseq_CpG_LNPs_mig_vs_SS_res_genelist
de_gs_by_k_split_triwise[[49]]<-Victor_LNP_CITEseq_pIC_LNPs_mig_vs_SS_res_genelist
de_gs_by_k_split_triwise[[50]]<-Victor_LNP_CITEseq_SS_mig_vs_res_genelist

names(de_gs_by_k_split_triwise)[48]<-"CpG-LNP-induced cDC1 maturation 8h (Rennen et al., This study, CITE-seq, Spleen)"
names(de_gs_by_k_split_triwise)[49]<-"pIC-LNP-induced cDC1 maturation 8h (Rennen et al., This study, CITE-seq, Spleen)"
names(de_gs_by_k_split_triwise)[50]<-"Homeostatic cDC1 maturation (Rennen et al., This study, CITE-seq, Spleen)"

## Also update the names to latest list of curated names
New_genelist_names<-read.xlsx("Paper/Suppl_Table_Meta_analysis_genelists_Lennon_check_Victor_170924_updated.xlsx", sheet = 1)
de_gs_by_k_split_triwise_updated<-de_gs_by_k_split_triwise
for (name in 1:length(de_gs_by_k_split_triwise_updated)){
  names(de_gs_by_k_split_triwise_updated)[1]
  for (name2 in 1:nrow(New_genelist_names)){
    if(New_genelist_names[name2,"Paper.annotation"] == names(de_gs_by_k_split_triwise_updated)[name]){
      names(de_gs_by_k_split_triwise_updated)[name]<-New_genelist_names[name2,"Nieuwe.naam"]
    }
  }
}

## Update other names with time info
names(de_gs_by_k_split_triwise_updated)[5]<-"pIC-induced cDC1 maturation 18h (Ardouin et al., 2016, Microarray, Spleen)"
names(de_gs_by_k_split_triwise_updated)[6]<-"STAg-induced cDC1 maturation 18h (Ardouin et al., 2016, Microarray, Spleen)"
names(de_gs_by_k_split_triwise_updated)[10]<-"eLNP-induced cDC1 maturation 8h (Rennen et al., This study, CITE-seq, Spleen)"  
names(de_gs_by_k_split_triwise_updated)[11]<-"pIC-induced cDC1 maturation 8h (Rennen et al., This study, CITE-seq, Spleen)"
names(de_gs_by_k_split_triwise_updated)[18]<-"Apoptotic cell-induced cDC1 signature (Cummings et al., 2016, Microarray, Small Intestine)"
names(de_gs_by_k_split_triwise_updated)[31]<-"R848-actDC1 8h (Torow et al., 2023, scRNA-seq, SI PP)"
names(de_gs_by_k_split_triwise_updated)[33]<-"R848-actDC2 8h (Torow et al., 2023, scRNA-seq, SI PP)"

## Double check
names(de_gs_by_k_split_triwise_updated)

##Update
de_gs_by_k_split_triwise<-de_gs_by_k_split_triwise_updated

## Reorder based on name author
de_gs_by_k_split_triwise_paper<-de_gs_by_k_split_triwise[c(42,43,1:6,7,8,14,16,18,44:47,21,24:33,37,38,39,10,11,48:50,9,12,13)]
Final_full_names<-names(de_gs_by_k_split_triwise_paper)
Final_abrev_names<-gsub(",.*,", ",", Final_full_names)


###################################################################################
##### 2D plot: Top 200 LNP overlapped lists as reference!!!
###################################################################################

# Calculate overlap lists (human lists and certain mouse lists removed for final plot paper)
Overlap_2D<-numeric()
Common_score<-numeric()
Gene_list_length<-numeric()
# de_gs_by_k_split_overlap_new<-de_gs_by_k_split_triwise[c(1:8,10,11,14:16,18:43)]
de_gs_by_k_split_overlap_new<-de_gs_by_k_split_triwise_paper[c(1:33)]
for (k in 1:length(de_gs_by_k_split_overlap_new)){
  counter<-k
  Overlap_2D[counter]<-(length(intersect(de_gs_by_k_split_triwise_paper[[38]],de_gs_by_k_split_overlap_new[[k]]))-length(intersect(de_gs_by_k_split_triwise_paper[[39]],de_gs_by_k_split_overlap_new[[k]])))
  Common_score[counter]<-length(intersect(de_gs_by_k_split_triwise_paper[[37]],de_gs_by_k_split_overlap_new[[k]]))
  Gene_list_length[counter]<-length(de_gs_by_k_split_overlap_new[[k]])
}
Overlap_2D
Gene_list_length
Common_score

Overlap_2D_df = as.matrix(cbind(Overlap_2D,Common_score,Gene_list_length))
rownames(Overlap_2D_df) <- names(de_gs_by_k_split_overlap_new)

Ordering1<-order(Overlap_2D_df[,1], decreasing = T)
Overlap_2D_df<-Overlap_2D_df[Ordering1,]
Colorset1<-Colorset[Ordering1]

Overlap_2D_df<-as.data.frame(Overlap_2D_df)
Overlap_2D_df$List_name<-rownames(Overlap_2D_df)

## Create plot
library(ggplot2)
library(ggrepel)
p <- ggplot(data=Overlap_2D_df, aes(x=Overlap_2D, y=log10(Gene_list_length))) + geom_point() + theme_minimal()

# add a column of NAs
Overlap_2D_df$diffexpressed <- "Undefined"
Overlap_2D_df$diffexpressed[Overlap_2D_df$Overlap_2D > 5] <- "Homeo"
Overlap_2D_df$diffexpressed[Overlap_2D_df$Overlap_2D < -5] <- "Immuno"

Overlap_2D_df$delabel <- NA
Overlap_2D_df$delabel[Overlap_2D_df$diffexpressed != "NO"] <- Overlap_2D_df$List_name[Overlap_2D_df$diffexpressed != "NO"]

# cDC1s vs cDC2s
Overlap_2D_df$CellType<- "cDCs"
Overlap_2D_df$CellType[grep("DC2",rownames(Overlap_2D_df))]<-"cDC2s"
Overlap_2D_df$CellType[grep("BMDCs",rownames(Overlap_2D_df))]<-"BMDCs"
Overlap_2D_df$CellType[grep("DC1",rownames(Overlap_2D_df))]<-"cDC1s"
Overlap_2D_df$CellType[grep("CD11b",rownames(Overlap_2D_df))]<-"cDC2s"
Overlap_2D_df$CellType[grep("CD103",rownames(Overlap_2D_df))]<-"cDC1s"

# Human vs Mouse
Overlap_2D_df$Species <- "Mouse"
Overlap_2D_df$Species[grep("Human",rownames(Overlap_2D_df))]<-"Human"

# Combo
Overlap_2D_df$Group<-paste0(Overlap_2D_df$Species,"_",Overlap_2D_df$CellType)

# colors
mycolors <- c("forestgreen", "firebrick1", "gray50")
names(mycolors) <- c("Homeo", "Immuno", "Undefined")

# New label parameter based on y-axis -> change to all!!
Overlap_2D_df$label_v2<-NA
Overlap_2D_df$label_v2[Overlap_2D_df$Gene_list_length >0] <- Overlap_2D_df$List_name[Overlap_2D_df$Gene_list_length > 0] #>100
# Replace names with shorter version
Overlap_2D_df$label_v2[Overlap_2D_df$Gene_list_length >0]<-Final_abrev_names[c(1:33)][Ordering1]

# Create new version without human datasets
Overlap_2D_df_without_human<-Overlap_2D_df[Overlap_2D_df$Species != "Human",]
colnames(Overlap_2D_df_without_human)[1]<-"Maturation Type"
colnames(Overlap_2D_df_without_human)[3]<-"Gene list length"
colnames(Overlap_2D_df_without_human)[5]<-"Score significance"

# Updated version 4 after feedback: Without species lists!!, color fixed based on score, label based on y-axis, bigger label, adapt limits x axis
cairo_pdf(file=paste0("results/2D_plot/Test_plot_top200_",analysis,"_1st_submission_v2.pdf"), width=15, height = 12)
ggplot(data=Overlap_2D_df_without_human, aes(x=`Maturation Type`, y=log10(`Gene list length`), col=`Score significance`, label=label_v2,shape=CellType)) + 
  geom_point(size = 2.5) + 
  theme_minimal() +
  geom_text_repel(size = 5) + #4
  scale_shape_manual(values=c(15,16,17,18))+
  scale_colour_manual(values = mycolors) +
  geom_vline(xintercept = 0) +
  xlim(-100, 100)
dev.off()
