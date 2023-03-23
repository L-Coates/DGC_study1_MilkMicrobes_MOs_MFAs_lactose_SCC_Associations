# Dairy Grand Challenge Study 1, 2017-2018
## This repository contains datasets resulting from analysis of raw bovine milk from cows in a randomized, cross-over feeding experiment consuming a high fiber/low starch diet and a low fiber/high starch diet.

### Funding and scientists
This study was funded by the USDA Agriculture Research Service. It is a collaboration among scientists at the ARS sites in Beltsville, MD, Madison, WI, Davis, CA, Wyndmoor, PA, University Park, PA, Grand Forks, ND, and Fort Collins, CO, and scientists at UC Davis. 

### Study aim
The aim of this study was to determine the relationships among microbes, fatty acids, lactose, 
oligosaccharides, and somatic cell counts in raw milk samples collected from Holstein dairy cows. 

### Study design
This study involves the use of our previously published datasets (please see https://doi.org/10.1093/cdn/nzac086; 10.1093/cdn/nzac033; https://doi.org/10.1093/cdn/nzab154; 10.1016/j.animal.2022.100599) to fulfill these secondary objectives: 

1) Determine relationships between milk microbes and milk oligosaccharides. 
2) Determine relationships between milk microbes and milk lactose. 
3) Determine relationships between milk microbes and milk fatty acids. 
4) Determine relationships between milk microbes and milk somatic cell count.
5) Determine relationships between milk somatic cell count and milk oligosaccharides. 
6) Determine relationships between milk somatic cell count and milk lactose. 
7) Determine relationships between milk somatic cell count and milk fatty acids.

### Folder and file descriptions
* picrust2_output_verbose -- Predicted milk microbial gene abundance. Output files from picrust2 plugin in Qiime2, version 2021.11.  
* BMO_Microbiota_Correlations_glmer_poisson.v2.R -- R code used to integrate and analyze bovine milk oligosaccharide data and bovine milk microbiota count data. 
* BMOs_SCC.R -- R code used to integrate and analyze bovine milk oligosaccharide data and bovine milk somatic cell count data. 
* Correlations_between_microbial_genera.R -- R code used to calculate correlations among milk microbial taxa and to generate a heatmap of these correlations.
* Dec 2020 DAPP 1 Fatty acid data.xlsx -- Milk fatty acid dataset. 
* KFK08 010218.xlsx -- Milk somatic cell count dataset from milk collected on January 2, 2018. 
* KFK08 010518.xlsx -- Milk somatic cell count dataset from milk collected on January 5, 2018. 
* KFK08 031218 031318 - March 16.xlsx -- Milk somatic cell count dataset from milk collected on March 12 and March 13, 2018. 
* KFK08 031918 032018 - March 23.xlsx -- Milk somatic cell count dataset from milk collected on March 19 and March 20, 2018. 
* KFK08 112017 112117 - Nov 24.xlsx -- Milk somatic cell count dataset from milk collected on November 20 and 21, 2017. 
* KFK08 112717 112817 - Dec 1.xlsx -- Milk somatic cell count dataset from milk collected on November 27 and 28, 2017. 
* Lactose_Microbiota_AlphaDiversity_Correlations.R -- R code used to integrate and analyze milk lactose data and milk microbiota count data. 
* Lactose_Microbiota_picrust2.R -- R code used to integrate and analyze milk lactose data and predicted milk microbial genes involved in lactose degradation. 
* Lactose_SCC.R -- R code used to integrate and analyze milk lactose data and milk somatic cell count data. 
* MFAs_microbiota_correlations_glmer_poisson.R -- R code used to integrate and analyze milk fatty acid data and milk microbial count data. 
* MilkFattyAcid_SCC_relationships.R -- R code used to integrate and analyze milk fatty acid data and milk somatic cell count data.
* PMT DVH_Dairy Grand Challenge coded FOSS composition data_101818.xlsx -- Milk lactose dataset. 
* USDA DGC Project 1 - BMO data (Barile Lab) - 11.17.2020 - SDD.xlsx -- Milk oligosaccharide dataset.
* feature-table-rarefied-at-507-samples-for-analysis.qza -- Qiime2 artifact for milk microbial counts rarefied at 507 reads, with a subset of samples used for subsequent analyses.
* feature-table-rarefied-at-507.qza --Qiime2 artifact for milk microbial counts rarefied at 507 reads.
* metadata.csv -- Metadata on the dairy cows in the study. 
* microbiota.SCC.LME.v2.R -- R code used to integrate and analyze milk microbial count data and milk somatic cell count data. 
* microbiota_sample_selection_rarefying_depth.R -- R code used to select milk microbial samples and calculate rarefying depth. 
* rep-seqs-rarefied-at-507-samples-for-analyses.qza -- Qiime2 artifact for milk microbial amplicon sequence variant sequences. 
* rooted-tree-rarefied507-samplesforanalyses.qza -- Qiime2 artifact for the rooted phylogenetic tree of the milk microbial amplicon sequence variants. 
* samples.for.microbiota.analyses.second.manuscript.txt -- List of milk microbial samples selected for analyses.
* taxonomy-PE-no-singletons-blanks-crmc-decontaminated-no-Bostaurus-mito-chloro-eukary.qza -- Qiime2 artifact for the taxonomic classifications of milk microbial amplicon sequence variants. 
