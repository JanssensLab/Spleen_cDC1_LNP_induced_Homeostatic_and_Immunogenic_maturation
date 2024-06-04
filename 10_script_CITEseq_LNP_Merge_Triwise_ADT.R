## Script for analyzing LNP cDC1 CITE-seq project data
## Triwise analysis ADT

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

dir.create("VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor")

################################################################################
## TRIWISE COMPARING BETWEEN SAMPLES
################################################################################

# ################
# ### PREPARE DATA
# ################

## Read in muscat data
res<-readRDS("VBO_LNP_merge/results/Robjects_Thesis_Victor/res_Muscat_DeSeq2_new_DS_ADT_thesis_Victor_VBO_4-12.rds")
expression_matrix <- assays(varianceStabilizingTransformation(res[["data"]][["cDC1s"]]))[[1]] #all pops

saveRDS(expression_matrix,"VBO_LNP_merge/results/Robjects_Thesis_Victor/Expr_matrix_cDC1s_Muscat_DeSeq2_new_DS_ADT_thesis_Victor_VBO_4-12.rds")

##### Get columns
cols_eLNP_LM<-grep("eLNPs_8h_Late_mature_cDC1s",colnames(expression_matrix))
cols_SS_LM<-grep("Steady_state_Late_mature_cDC1s",colnames(expression_matrix))

cols_SS_Imm<-grep("Steady_state_Immature_cDC1s",colnames(expression_matrix))

cols_pIC_LNP_LM<-grep("pIC_LNPs_8h_Late_mature_cDC1s",colnames(expression_matrix))
cols_CpG_LNP_LM<-grep("CpG_LNPs_8h_Late_mature_cDC1s",colnames(expression_matrix))

##### Create expTableMean 
eLNP_LM_mean<-apply(expression_matrix[,cols_eLNP_LM],1,mean)
SS_LM_mean<-apply(expression_matrix[,cols_SS_LM],1,mean)

SS_Imm_mean<-apply(expression_matrix[,cols_SS_Imm],1,mean)

pIC_LNP_LM_mean<-apply(expression_matrix[,cols_pIC_LNP_LM],1,mean)
CpG_LNP_LM_mean<-apply(expression_matrix[,cols_CpG_LNP_LM],1,mean)

expTable_mean_v1<-cbind(SS_Imm_mean,eLNP_LM_mean,pIC_LNP_LM_mean)
colnames(expTable_mean_v1)<-c('SS_Imm','eLNP_LM','pIC_LNP_LM')

expTable_mean_v2<-cbind(SS_Imm_mean,eLNP_LM_mean,CpG_LNP_LM_mean)
colnames(expTable_mean_v2)<-c('SS_Imm','eLNP_LM','CpG_LNP_LM')

expTable_mean_v3<-cbind(SS_Imm_mean,SS_LM_mean,pIC_LNP_LM_mean)
colnames(expTable_mean_v3)<-c('SS_Imm','SS_LM','pIC_LNP_LM')

expTable_mean_v4<-cbind(SS_Imm_mean,SS_LM_mean,CpG_LNP_LM_mean)
colnames(expTable_mean_v4)<-c('SS_Imm','SS_LM','CpG_LNP_LM')

## Save exptables
write.xlsx(as.data.frame(expTable_mean_v1), "VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/ExpTable_mean_v1_ADT.xlsx", row.names = T)
write.xlsx(as.data.frame(expTable_mean_v2), "VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/ExpTable_mean_v2_ADT.xlsx", row.names = T)
write.xlsx(as.data.frame(expTable_mean_v3), "VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/ExpTable_mean_v3_ADT.xlsx", row.names = T)
write.xlsx(as.data.frame(expTable_mean_v4), "VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/ExpTable_mean_v4_ADT.xlsx", row.names = T)

######## CREATE BARYCOORDINATES
barycoords_v1 <- transformBarycentric(expTable_mean_v1)
barycoords_v2 <- transformBarycentric(expTable_mean_v2)
barycoords_v3 <- transformBarycentric(expTable_mean_v3)
barycoords_v4 <- transformBarycentric(expTable_mean_v4)

################################################################################
## TRIWISE CONTINUED FOR EITHER (first batch: only v2 and v3)
################################################################################

################
### DEG
################

## Work with lists with stricter cutoffs!! Otherwise too many DE genes
## clint suffix: p_adj.glb < 0.01, abs(logFC) > 1, baseMean > 50

eLNP_LM_vs_SS_Imm<-readRDS("VBO_LNP_merge/results/Robjects_Thesis_Victor/tbl_fil_muscat_DESeq2_new_DS_ADT_thesis_Victor_eLNPs_8h_Late_mature_cDC1s-Steady_state_Immature_cDC1s_VBO_4-12_clint.rds")

pIC_LNP_LM_vs_SS_Imm<-readRDS("VBO_LNP_merge/results/Robjects_Thesis_Victor/tbl_fil_muscat_DESeq2_new_DS_ADT_thesis_Victor_pIC_LNPs_8h_Late_mature_cDC1s-Steady_state_Immature_cDC1s_VBO_4-12_clint.rds")
pIC_LNP_LM_vs_eLNP_LM<-readRDS("VBO_LNP_merge/results/Robjects_Thesis_Victor/tbl_fil_muscat_DESeq2_new_DS_ADT_thesis_Victor_pIC_LNPs_8h_Late_mature_cDC1s-eLNPs_8h_Late_mature_cDC1s_VBO_4-12_clint.rds")

CpG_LNP_LM_vs_SS_Imm<-readRDS("VBO_LNP_merge/results/Robjects_Thesis_Victor/tbl_fil_muscat_DESeq2_new_DS_ADT_thesis_Victor_CpG_LNPs_8h_Late_mature_cDC1s-Steady_state_Immature_cDC1s_VBO_4-12_clint.rds")
CpG_LNP_LM_vs_eLNP_LM<-readRDS("VBO_LNP_merge/results/Robjects_Thesis_Victor/tbl_fil_muscat_DESeq2_new_DS_ADT_thesis_Victor_CpG_LNPs_8h_Late_mature_cDC1s-eLNPs_8h_Late_mature_cDC1s_VBO_4-12_clint.rds")

SS_LM_vs_SS_Imm<-readRDS("VBO_LNP_merge/results/Robjects_Thesis_Victor/tbl_fil_muscat_DESeq2_new_DS_ADT_thesis_Victor_Steady_state_Late_mature_cDC1s-Steady_state_Immature_cDC1s_VBO_4-12_clint.rds")

pIC_LNP_LM_vs_SS_LM<-readRDS("VBO_LNP_merge/results/Robjects_Thesis_Victor/tbl_fil_muscat_DESeq2_new_DS_ADT_thesis_Victor_pIC_LNPs_8h_Late_mature_cDC1s-Steady_state_Late_mature_cDC1s_VBO_4-12_clint.rds")

CpG_LNP_LM_vs_SS_LM<-readRDS("VBO_LNP_merge/results/Robjects_Thesis_Victor/tbl_fil_muscat_DESeq2_new_DS_ADT_thesis_Victor_CpG_LNPs_8h_Late_mature_cDC1s-Steady_state_Late_mature_cDC1s_VBO_4-12_clint.rds")

Gdiffexp_eLNP_LM_vs_SS_Imm <- c(eLNP_LM_vs_SS_Imm$`cDC1s`$gene)

Gdiffexp_pIC_LNP_LM_vs_SS_Imm <- c(pIC_LNP_LM_vs_SS_Imm$`cDC1s`$gene)
Gdiffexp_pIC_LNP_LM_vs_eLNP_LM <- c(pIC_LNP_LM_vs_eLNP_LM$`cDC1s`$gene)

Gdiffexp_CpG_LNP_LM_vs_SS_Imm <- c(CpG_LNP_LM_vs_SS_Imm$`cDC1s`$gene)
Gdiffexp_CpG_LNP_LM_vs_eLNP_LM <- c(CpG_LNP_LM_vs_eLNP_LM$`cDC1s`$gene)

Gdiffexp_SS_LM_vs_SS_Imm <- c(SS_LM_vs_SS_Imm$`cDC1s`$gene)

Gdiffexp_pIC_LNP_LM_vs_SS_LM <- c(pIC_LNP_LM_vs_SS_LM$`cDC1s`$gene)

Gdiffexp_CpG_LNP_LM_vs_SS_LM <- c(CpG_LNP_LM_vs_SS_LM$`cDC1s`$gene)


## Combine
All_DEG_v1<-unique(c(Gdiffexp_pIC_LNP_LM_vs_eLNP_LM,Gdiffexp_pIC_LNP_LM_vs_SS_Imm,Gdiffexp_eLNP_LM_vs_SS_Imm))
All_DEG_v2<-unique(c(Gdiffexp_CpG_LNP_LM_vs_eLNP_LM,Gdiffexp_CpG_LNP_LM_vs_SS_Imm,Gdiffexp_eLNP_LM_vs_SS_Imm))
All_DEG_v3<-unique(c(Gdiffexp_pIC_LNP_LM_vs_SS_LM,Gdiffexp_pIC_LNP_LM_vs_SS_Imm,Gdiffexp_SS_LM_vs_SS_Imm))
All_DEG_v4<-unique(c(Gdiffexp_CpG_LNP_LM_vs_SS_LM,Gdiffexp_CpG_LNP_LM_vs_SS_Imm,Gdiffexp_SS_LM_vs_SS_Imm))

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

Ardouin_C72_common<-read.xlsx("../../Meta_analysis/Ardouin_lists.xlsx", sheet = "C72")
Ardouin_C72_common_genelist<-Ardouin_C72_common$Gene.Symbol
Ardouin_C72_common_genelist<-firstup(tolower(Ardouin_C72_common_genelist)) # Mouse data, but all letters capitalized!!
Ardouin_C73_common<-read.xlsx("../../Meta_analysis/Ardouin_lists.xlsx", sheet = "C73")
Ardouin_C73_common_genelist<-Ardouin_C73_common$Gene.Symbol
Ardouin_C73_common_genelist<-firstup(tolower(Ardouin_C73_common_genelist)) # Mouse data, but all letters capitalized!!
Ardouin_C115_common<-read.xlsx("../../Meta_analysis/Ardouin_lists.xlsx", sheet = "C115")
Ardouin_C115_common_genelist<-Ardouin_C115_common$Gene.Symbol
Ardouin_C115_common_genelist<-firstup(tolower(Ardouin_C115_common_genelist)) # Mouse data, but all letters capitalized!!
ArdouinC72_73_115_common_genelist<-c(Ardouin_C72_common_genelist,Ardouin_C73_common_genelist,Ardouin_C115_common_genelist) #Combine as ref!

##########################
### GENES TO MARK ON PLOT
##########################

# GOI (GENES OF INTEREST) ARE TOP EXPRESSED OR OWN LIST

##########################
### GENERATE PLOT
##########################

p1 <- plotDotplot(barycoords_v1, Gdiffexp = Gdiffexp_pIC_LNP_LM_vs_eLNP_LM, Goi= genelist, rmax = 4)
pdf("VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/Plot_v1_ADT_Late_mature_cDC1s_pIC_LNP_LM_vs_eLNP_LM.pdf", width = 15, height = 15)
print(p1)
dev.off()

p2 <- plotDotplot(barycoords_v1, Gdiffexp = Gdiffexp_pIC_LNP_LM_vs_SS_Imm, Goi= genelist, rmax = 4)
pdf("VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/Plot_v1_ADT_Late_mature_cDC1s_pIC_LNP_LM_vs_SS_Imm.pdf", width = 15, height = 15)
print(p2)
dev.off()

p3 <- plotDotplot(barycoords_v1, Gdiffexp = Gdiffexp_eLNP_LM_vs_SS_Imm, Goi= genelist, rmax = 4)
pdf("VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/Plot_v1_ADT_Late_mature_cDC1s_eLNP_LM_vs_SS_Imm.pdf", width = 15, height = 15)
print(p3)
dev.off()

p4 <- plotDotplot(barycoords_v1, Gdiffexp = All_DEG_v1, Goi= genelist, rmax = 4)
pdf("VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/Plot_v1_ADT_All_DEG_three_way_comparison.pdf", width = 15, height = 15)
print(p4)
dev.off()

p5 <- plotDotplot(barycoords_v2, Gdiffexp = All_DEG_v2, Goi= genelist, rmax = 4)
pdf("VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/Plot_v2_ADT_All_DEG_three_way_comparison.pdf", width = 15, height = 15)
print(p5)
dev.off()

p6 <- plotDotplot(barycoords_v3, Gdiffexp = All_DEG_v3, Goi= genelist, rmax = 4)
pdf("VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/Plot_v3_ADT_All_DEG_three_way_comparison.pdf", width = 15, height = 15)
print(p6)
dev.off()

p7 <- plotDotplot(barycoords_v4, Gdiffexp = All_DEG_v4, Goi= genelist, rmax = 4)
pdf("VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/Plot_v4_ADT_All_DEG_three_way_comparison.pdf", width = 15, height = 15)
print(p7)
dev.off()

###Rose plot
p1r <- plotRoseplot(barycoords_v1, Gdiffexp = All_DEG_v1, showlabels = F)
pdf("VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/Roseplot_v1_ADT_All_DEG_three_way_comparison.pdf", width = 15, height = 15)
print(p1r)
dev.off()

# Interactive
## v1
p1i<-interactiveDotplot(expTable_mean_v1, Gdiffexp=c(All_DEG_v1), Goi= genelist, plotLocalEnrichment=FALSE, rmax=4)
print(p1i)
fileName<-paste0(getwd(),"/VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/interactivePlot_v1_ADT_All_DEGs_three_way_comparison.html")
saveWidget(p1i,file=fileName) ##needs full path!

####################

## Update: only label the genes which are outside inner circle (stronger logFC)
## Paper February 2024

# CHANGE SCALING SO THAT LABELS ARE VISIBLE
genelist<-All_DEG_v1

p <- plotDotplot(barycoords_v1, Gdiffexp = All_DEG_v1, Goi= genelist, rmax = 5)

p <- p + scale_y_continuous(expand = expansion(mult = 0.2)) + scale_x_continuous(expand = expansion(mult = 0.2))   #change scaling

#MARK GENES ON PLOT

# FROM FIND MAKERS -> based on r value (first version R > 1) -> 04/24 r > 0.99 to include Cd1d
NegXBarycoords <- subset(barycoords_v1[intersect(unlist(genelist),All_DEG_v1),],
                         barycoords_v1[intersect(unlist(genelist),All_DEG_v1),1]< 0 & barycoords_v1[intersect(unlist(genelist),All_DEG_v1),4]> 0.99)
PosXBarycoords <- subset(barycoords_v1[intersect(unlist(genelist),All_DEG_v1),],
                         barycoords_v1[intersect(unlist(genelist),All_DEG_v1),1]> 0 & barycoords_v1[intersect(unlist(genelist),All_DEG_v1),4]> 0.99)
p <- p +
  geom_text_repel(NegXBarycoords[],
                  mapping=aes(label=rownames(NegXBarycoords), x=NegXBarycoords[,1], y=NegXBarycoords[,2]),
                  force = 1, size=4, segment.alpha=0.5, segment.colour = '#999999',
                  nudge_x = (0-NegXBarycoords[,1]-4),
                  box.padding =0.25) +
  geom_text_repel(PosXBarycoords[,],
                  mapping=aes(label=rownames(PosXBarycoords), x=PosXBarycoords[,1], y=PosXBarycoords[,2]),
                  force = 1, size=4, segment.alpha=0.5, segment.colour = '#999999',
                  nudge_x = (0-PosXBarycoords[,1]+2),
                  box.padding =0.25)
print(p)

pdf("VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/Plot_v1_ADT_labeled_All_DEGs_three_way_comparison_Sofie_paper_new_filtering.pdf", width = 10, height = 10)
print(p)
dev.off()

###########################################################################################################
###########################################################################################################

################################################################################
######### TRIWISE PLOTS (colored lists)
################################################################################

########################################
##### Preparation
########################################

wantedColors7<-c(nodiffall="gray55",diffall="black",nodiffSet1="indianred1",diffSet1="red", 
                 nodiffSet2="limegreen",diffSet2="darkgreen", nodiffSet3="lightblue", diffSet3="darkblue", 
                 nodiffSet4="darkorchid1", diffSet4="darkmagenta", nodiffSet5="lightsalmon", diffSet5="darkorange", 
                 nodiffSet6="lightskyblue", diffSet6="turquoise4", nodiffSet7="bisque", diffSet7="burlywood4")

###Extra: create colored plots with gene lists! (only eBayes)
###Seperate the DE genes of the different comparisons
#set1 (red)=DE genes only present in pIC_LNP_LM_vs_SS_Imm
#set2 (green)=DE genes only present in eLNP_LM_vs_SS_Imm
#set3 (blue)=DE genes only present in pIC_LNP_LM_vs_eLNP_LM
#set4 (purple)=DE genes present in all 3 comparisons
#set5 (orange)=DE genes present in pIC_LNP_LM_vs_SS_Imm and eLNP_LM_vs_SS_Imm
#set6 (turquoise)=DE genes present in eLNP_LM_vs_SS_Imm and pIC_LNP_LM_vs_eLNP_LM
#set7 (brown)=DE genes present in pIC_LNP_LM_vs_SS_Imm and pIC_LNP_LM_vs_eLNP_LM

DEgenesGroup1 <- pIC_LNP_LM_vs_SS_Imm$cDC1s
rownames(DEgenesGroup1) <- DEgenesGroup1$gene
DEgenesGroup2 <- eLNP_LM_vs_SS_Imm$cDC1s
rownames(DEgenesGroup2) <- DEgenesGroup2$gene
DEgenesGroup3 <- pIC_LNP_LM_vs_eLNP_LM$cDC1s
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

#Take p adj global (extra strict cutoff used!!)
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
write.xlsx(geneset_list, "VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/triwisePlot_withColors_v1_ADT.xlsx")

p_col<-plotDotplot(barycoords_v1, Gdiffexp=All_DEG_v1, showlabels = F, rmax = 4,
                   Goi=list(Set1=genesSet1,Set2=genesSet2,Set3=genesSet3,Set4=genesSet4,
                            Set5=genesSet5,Set6=genesSet6,Set7=genesSet7), colorvalues=wantedColors7) +
  theme(legend.position = "none")

pdf("VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/Plot_v1_ADT_All_DEGs_colored.pdf", width = 15, height = 15)
print(p_col)
dev.off()

myList<-list("set1"=genesSet1,"set2"=genesSet2,"set3"=genesSet3,"set4"=genesSet4,"set5"=genesSet5,"set6"=genesSet6,"set7"=genesSet7)
saveRDS(myList,file="VBO_LNP_merge/results/Triwise/ADT_Thesis_Victor/listGenes_coloredTriwise_v1_ADT.rds")
