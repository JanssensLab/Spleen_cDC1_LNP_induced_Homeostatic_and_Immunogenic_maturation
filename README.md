# Spleen_cDC1_LNP_induced_Homeostatic_and_Immunogenic_maturation

Novel tools to deconvolute homeostatic from immunogenic mature DCs based on lipid nanoparticle-induced DC maturation.

## About

Dendritic cells (DCs) are short-lived immune cells that continuously roam our body in search for foreign or self-antigens. Upon acquisition of antigen, they mature and start migrating to the lymph node to present the antigen to naïve T cells. Depending on the context wherein the antigen is acquired, DCs will mature in a homeostatic or immunogenic manner. So far, the field is lacking proper tools to distinguish between the two maturation states. Most maturation markers are shared between the two states and therefore inappropriate to use. Still, defining the proper maturation type is crucial as it determines how the DCs will instruct the T cells towards antigen expressing cells. In this study, we used a lipid nanoparticle (LNP)-based approach to steer DC maturation pathways in vivo. CITE-seq analysis allowed us to design a panel of flow cytometry markers that reliably annotates the two DC maturation states, as validated in an infection and in a tumor model. Furthermore, the data corroborated that uptake of empty LNPs in DCs induces their homeostatic maturation, in contrast to uptake of mRNA-LNPs or TLR ligand-adjuvanted LNPs, leading to distinct effector T cell outputs. This reveals that LNPs themselves are not being decoded as “danger” by cDC1s, and that the cargo is essential to provide adjuvants activity, which is highly relevant for targeted design of LNP-based therapies.


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

Sofie Rennen et al., Lipid nanoparticles as a tool to dissect dendritic cell maturation pathways.
