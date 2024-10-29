# Spleen_cDC1_LNP_induced_Homeostatic_and_Immunogenic_maturation

## About

Dendritic cells (DCs) are short-lived immune cells that continuously scan our body for the presence of foreign or self-antigens. As soon as they acquire antigen, they mature and start migrating to the lymph node to present the antigen to naïve T cells. Depending on how the antigen is perceived, as innocuous or as dangerous, they will mature in a homeostatic or immunogenic manner. So far, the field is lacking proper tools to distinguish between the two maturation states. Most maturation markers are shared and therefore inappropriate to use. Still, defining the proper maturation type is crucial as it determines how the DCs will instruct the T cells towards antigen expressing cells. In this study, we used a lipid nanoparticle (LNP)-based approach to steer DC maturation pathways in vivo. CITE-Seq analysis allowed us to design a panel of flow cytometry markers that can be used to reliably annotate the two maturation states. The data also revealed that the uptake of empty LNPs induces homeostatic DC maturation, indicating that LNPs have no intrinsic adjuvants activity.


## Overview scripts

Here's an overview of the various R scripts used in processing the CITE-Seq data in the manuscript Rennen et al.:
- 1_script_CITEseq_BareBones_RNA-ADT_HPCscript.R: Standard pipeline used on the High Performance Cluster for analyzing the RNA and ADT assays of the CITE-seq objects
- 2_script_CITEseq_LNP_Merge.R: Script for basic and downstream analysis of the merged cDC1 CITE-seq data of all conditions
- 3_script_CITEseq_LNP_Merge_muscat_RNA_part1.R: Script for muscat Differential State analysis of the merged CITE-seq RNA data part 1
- 4_script_CITEseq_LNP_Merge_muscat_RNA_part2.R: Script for muscat Differential State analysis of the merged CITE-seq RNA data part 2
- 5_script_CITEseq_LNP_Merge_muscat_RNA_part3.R: Script for muscat Differential State analysis of the merged CITE-seq RNA data part 3
- 6_script_CITEseq_LNP_Merge_muscat_ADT_part1.R: Script for muscat Differential State analysis of the merged CITE-seq ADT data part 1
- 7_script_CITEseq_LNP_Merge_muscat_ADT_part2.R: Script for muscat Differential State analysis of the merged CITE-seq ADT data part 2
- 8_script_CITEseq_LNP_Merge_Triwise_RNA_part1.R: Script for triwise analysis of the merged CITE-seq RNA data part 1
- 9_script_CITEseq_LNP_Merge_Triwise_RNA_part2.R: Script for triwise analysis of the merged CITE-seq RNA data part 2
- 10_script_CITEseq_LNP_Merge_Triwise_ADT.R: Script for triwise analysis of the merged CITE-seq ADT data
- 11_script_CITEseq_LNP_Merge_DoRothEA.R: Script for DoRothEA TF analysis of the merged CITE-seq data
- 12_Meta_analysis_script.R: Script for meta analysis of literature gene lists to determine type of DC maturation

## Citation

Sofie Rennen et al., Lipid nanoparticles steer dendritic cell maturation pathways in vivo.
