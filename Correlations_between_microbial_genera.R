#This script is the code used to generate a heat map of the correlations among microbial taxa. 

#Load needed libraries
library(devtools)
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
library(plotly)
library(vegan)
library(picante)
library(ggdendro)
library(grid)

#STEP 1:Reading in the metadata table
metadata  <- read.csv("../../../DairyGrandChallenge.Study1.2017-2018/ResearchMaterials.DairyGrandChallenge/Cleaned_For_MetadataFile/metadata.csv")

#keeping only the data that is needed. 
metadata<- metadata[,c(1:15)]
metadata <- metadata[-1,]

#STEP 2: Read in the table with number of reads per sample, before rarefaction. 
reads.per.sample <- read.csv("sample-frequency.csv", head=FALSE)
colnames(reads.per.sample) <- c("sample.id", "reads")
metadata.v2 <- merge(metadata, reads.per.sample, by="sample.id", all=FALSE)
dim(metadata.v2) #508 samples
length(unique(metadata.v2$cow))#111 cows

#STEP 3: Rename sample 5849-D1 to 5849-D2 because cow 5849 actually missed the milk collection
#on November 25, 2017, but milk was collected from this cow on November 26, 2017.
#The same was true for cow 6233. 
metadata.v3 <- metadata.v2
metadata.v3$sample.id <- gsub(pattern="5849-D1", replacement="5849-D2", x=metadata.v3$sample.id)
metadata.v3$sample.id <- gsub(pattern="6233-D1", replacement="6233-D2", x=metadata.v3$sample.id)
metadata.v3[metadata.v3$sample.id=="5849-D2","day"] <- 2
metadata.v3[metadata.v3$sample.id=="6233-D2","day"] <- 2

#STEP 4: Loading in the feature table that's been rarefied and filtered
library(qiime2R)
feature.table <- read_qza("feature-table-rarefied-at-507-samples-for-analysis.qza")
feature.table <- as.data.frame(feature.table$data)

#Rename 5849-D1 to 5849-D2
feature.table.v2 <- feature.table
library(dplyr)
feature.table.v2 <- rename(feature.table.v2, "5849-D2"="5849-D1")

#STEP 5: Filter the metadata samples down to what is the feature table
metadata.v4 <- metadata.v3[which(metadata.v3$sample.id %in% colnames(feature.table.v2)),]
dim(metadata.v4) #334 samples

#STEP 6: If a cow has more than one sample in a period then, select the sample with
#higher reads (before rarefying) for a cow in a period. 

ChooseSampleFromPeriod <-function(cow) {
    NoDuplicates <- data.frame()
    for(i in cow) {
        select.cow <- metadata.v4[metadata.v4$cow==i,]
        period0 <- select.cow[select.cow$period=="p_0",]
        if(length(period0$sample.id)>1) {
            period0<-period0[order(period0$reads, decreasing=TRUE),]
            Higher.Read.Sample <- period0[match(unique(period0$period), period0$period),]
            NoDuplicates <- rbind(NoDuplicates, Higher.Read.Sample)
        } else {
            NoDuplicates <- rbind(NoDuplicates, period0)
        }
        period1 <- select.cow[select.cow$period=="p_1",]
        if(length(period1$sample.id)>1) {
            period1<-period1[order(period1$reads, decreasing=TRUE),]
            Higher.Read.Sample <- period1[match(unique(period1$period), period1$period),]
            NoDuplicates <- rbind(NoDuplicates, Higher.Read.Sample)
        } else {
            NoDuplicates <- rbind(NoDuplicates, period1)
        }
        period2 <- select.cow[select.cow$period=="p_2",]
        if(length(period2$sample.id)>1) {
            period2<-period2[order(period2$reads, decreasing=TRUE),]
            Higher.Read.Sample <- period2[match(unique(period2$period), period2$period),]
            NoDuplicates <- rbind(NoDuplicates, Higher.Read.Sample)
        } else {
            NoDuplicates <- rbind(NoDuplicates, period2)
        }
    }
    return(NoDuplicates)
}

#Run the function "ChooseSampleFromPeriod"
cows <-unique(metadata.v4$cow)
metadata.v5 <- ChooseSampleFromPeriod(cows) 
table(metadata.v5$sample.id)
length(metadata.v5$sample.id) #191 samples
length(unique(metadata.v5$cow)) #66 cows

#STEP 7: Keep only the cows with a sample in each period. 
cows.three.samples <- data.frame(table(metadata.v5$cow))
cows.three.samples <- cows.three.samples[cows.three.samples$Freq==3,]
cows.three.samples <- cows.three.samples$Var1
length(cows.three.samples) #60 cows
metadata.v6 <- metadata.v5[c(which(metadata.v5$cow %in% cows.three.samples)),]
dim(metadata.v6) #180 samples

#STEP 8: Read in taxonomy
taxonomy <- read_qza("taxonomy-PE-no-singletons-blanks-crmc-decontaminated-no-Bostaurus-mito-chloro-eukary.qza")

#change it to the way phyloseq wants it 
taxonomy2 <- as.data.frame(taxonomy$data) %>% column_to_rownames("Feature.ID") %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% as.matrix()
#remove "Confidence" column
taxonomy2 <- taxonomy2[,-8]

#trim taxonomy to include only ASVs in the ASV table
taxonomy3 <- taxonomy2[c(rownames(taxonomy2)%in%rownames(feature.table.v2)),]

#Replace all empty levels with NA for the sake of consistency
#This prevents inaccurate glomming later
taxonomy3[,grep("Species", colnames(taxonomy3))][grep("s__$", taxonomy3[,grep("Species", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Genus", colnames(taxonomy3))][grep("g__$", taxonomy3[,grep("Genus", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Family", colnames(taxonomy3))][grep("f__$", taxonomy3[,grep("Family", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Order", colnames(taxonomy3))][grep("o__$", taxonomy3[,grep("Order", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Class", colnames(taxonomy3))][grep("c__$", taxonomy3[,grep("Class", colnames(taxonomy3))])] <- NA
taxonomy3[,grep("Phylum", colnames(taxonomy3))][grep("p__$", taxonomy3[,grep("Phylum", colnames(taxonomy3))])] <- NA

#STEP 9: Read in rooted tree. 
tree <- read_qza("rooted-tree-rarefied507-samplesforanalyses.qza")

#STEP 10: build a phyloseq object
rownames(metadata.v6)=NULL
phylobj <- phyloseq(otu_table(feature.table.v2, taxa_are_rows=T),tax_table(taxonomy3), phy_tree(tree$data), sample_data(metadata.v6 %>% as.data.frame() %>% column_to_rownames("sample.id")))

#keep only the taxa that have at least one count
phylobj.filtered <- filter_taxa(phylobj, function(x) sum(x)>0, TRUE)

phylobj.genus <- tax_glom(phylobj.filtered, 'Genus', NArm=F)
dim(otu_table(phylobj.genus)) #There are 1166 genera and  samples

#Create a list of unique taxa names that can be assigned to the features within the phyloseq object
taxa.names <- tax_table(phylobj.genus)[,6]

#replace empty names with family level identities
taxa.names2 <- tax_table(phylobj.genus)[,5] #make a list of family level designations
taxa.names[c(which(is.na(taxa.names)), grep("g__$", taxa.names))] <-taxa.names2[c(which(is.na(taxa.names)),grep("g__$", taxa.names))]

#replace empty names with order level identities
taxa.names3 <-tax_table(phylobj.genus)[,4] #make a list of order level designations
taxa.names[c(which(is.na(taxa.names2)), grep( "f__$", taxa.names))] <- taxa.names3[c(which(is.na(taxa.names2)), grep("f__$", taxa.names))]

#replace empty names with class level identities
taxa.names4 <- tax_table(phylobj.genus)[,3] #make a list of class level designations
taxa.names[c(which(is.na(taxa.names3)), grep("o__$", taxa.names))] <-taxa.names4[c(which(is.na(taxa.names3)),grep("o__$", taxa.names))]

#replace empty names with phylum level identities
taxa.names5 <- tax_table(phylobj.genus)[,2] #make a list of phylum level designations
taxa.names[c(which(is.na(taxa.names4)),grep("c__$", taxa.names))] <-taxa.names5[c(which(is.na(taxa.names4)),grep("c__$", taxa.names))]

#replace empty names with kingdom level identities
taxa.names6 <- tax_table(phylobj.genus)[,1] #make a list of phylum level designations
taxa.names[c(which(is.na(taxa.names5)),grep("p__$", taxa.names))] <-taxa.names6[c(which(is.na(taxa.names5)),grep("p__$", taxa.names))]

#remove leading white space
taxa.names <- gsub(" ","", taxa.names)

#change the taxa designations to names compatible with R
taxa.names <- make.names(taxa.names, unique = TRUE)

#Make all of the feature names within the phyloseq object correspond to the genus level designation
taxa_names(phylobj.genus) <- taxa.names

#STEP 11: Preparing the phyloseq genus level object for ggplot and statistics
melted.phylgenus <- psmelt(phylobj.genus)

#confirm that the set of 60 cows and 180 samples were retained
length(unique(melted.phylgenus$cow)) #60 cows
length(unique(melted.phylgenus$Sample)) #180 samples

#STEP 12: Filtering down to the list of 44 genera that were included in any one of the set
#of analyses among milk components, or so, and microbiota. 

#read in set of genera that were analyzed across the different models
set.of.44genera <- read.csv("Genera_tested_for_relationships_with_bioactives.csv")
set.of.44genera <- set.of.44genera$Genus

melted.phylgenus.v2 <- melted.phylgenus[c(which(melted.phylgenus$OTU %in% set.of.44genera)),]
length(unique(melted.phylgenus.v2$OTU)) #44 genera

#STEP 13: Run a function to calculate Spearman rank correlation coefficient on the genus counts

#turn into wide format
melted.phylgenus.v2.wide <- pivot_wider(data=melted.phylgenus.v2, id_cols=Sample, names_from=OTU, values_from=Abundance)

Spearman.function <-function(genus1, genus2) {
    correlations <- data.frame()
    for(i in genus1) {
        x.vector <- melted.phylgenus.v2.wide[[i]]
        for(j in genus2){
            y.vector <- melted.phylgenus.v2.wide[[j]]
            corr.output <-cor.test(x=x.vector,y=y.vector, method="spearman")
            corr.output.2 <- data.frame(genus_1=i, genus_2=j, estimate=corr.output[["estimate"]])
            correlations <- rbind(correlations, corr.output.2)
        } 
    }
    return(correlations)
}

genus1<-colnames(melted.phylgenus.v2.wide[,2:45])
genus2 <- colnames(melted.phylgenus.v2.wide[,2:45])
spearman.correlations <- Spearman.function(genus1, genus2)

#STEP 14: Create heat map of correlations among milk microbial genera

#remove the "f__" and "g__" designations in front of the taxon name
spearman.correlations$genus_1 <- gsub(pattern="g__", replacement="", x=spearman.correlations$genus_1)
spearman.correlations$genus_1 <- gsub(pattern="f__", replacement="", x=spearman.correlations$genus_1)
spearman.correlations$genus_2 <- gsub(pattern="g__", replacement="", x=spearman.correlations$genus_2)
spearman.correlations$genus_2 <- gsub(pattern="f__", replacement="", x=spearman.correlations$genus_2)

#change to wide format
spearman.correlations.wide <- pivot_wider(data=spearman.correlations, id_cols=genus_1, names_from=genus_2, values_from=estimate)
spearman.correlations.wide <- as.data.frame(spearman.correlations.wide)
rownames(spearman.correlations.wide)<-spearman.correlations.wide$genus_1
spearman.correlations.wide.v2 <- spearman.correlations.wide[,-1]

#get phylogenetic tree for cladogram of the genera for the heatmap
library(phylogram)
phylobj.44genera <- prune_taxa(x=phylobj.genus, taxa=genus1)#cut down to only the 44 genera t 
taxa_names(phylobj.44genera)

phylogenetic.tree <- phy_tree(phylobj.44genera)
dendrogram <- as.dendrogram.phylo(phylogenetic.tree)
cladogram <- as.cladogram(as.dendrogram.phylo(phylogenetic.tree))
plot(cladogram)
#re-order the genera in the heatmap by the order of genera in the phylogenetic tree
genera.in.phy.tree <- taxa_names(phylogenetic.tree)
genera.in.phy.tree <- gsub(pattern="g__", replacement="", genera.in.phy.tree)
genera.in.phy.tree <- gsub(pattern="f__", replacement="", genera.in.phy.tree)
spearman.correlations.wide.v3 <- arrange(spearman.correlations.wide,factor(genus_1, levels = genera.in.phy.tree))
spearman.correlations.wide.v3 <- spearman.correlations.wide.v3[,-1]
spearman.correlations.wide.v4 <- spearman.correlations.wide.v3[,c(genera.in.phy.tree)]
heatmap(x=as.matrix(spearman.correlations.wide.v4), scale="none", Rowv = as.cladogram(as.dendrogram.phylo(phylogenetic.tree)), Colv = "Rowv")

library(ComplexHeatmap)
Heatmap(matrix=as.matrix(spearman.correlations.wide.v4), cluster_rows= as.cladogram(as.dendrogram.phylo(phylogenetic.tree)), cluster_columns = as.cladogram(as.dendrogram.phylo(phylogenetic.tree)), show_heatmap_legend = TRUE, heatmap_legend_param = list(title="Spearman's rho"),row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize=10), row_names_max_width=unit(12,"cm"), column_names_max_height = unit(12,"cm"), heatmap_width = unit(1,"npc"), heatmap_height = unit(1,"npc"))
#saved the heatmap as "microbiota.SpearmanCorrelation.heatmap.jpg"

#the phylogenetic tree doesn't agree with the taxonomy, so now generating 
#a cladogram based on the taxonomic assignment, instead of the phylogenetic tree

taxa.table <- as.data.frame(tax_table(phylobj.44genera))
#get rid of species column since it's empty anyway
taxa.table <- taxa.table[,-7]
taxa.table$Genus <- rownames(taxa.table)
#for the unknown genera, fill in the blank with the family name
taxa.table[is.na(taxa.table$Genus),6] <- taxa.table[is.na(taxa.table$Genus),][["Family"]]

#remove "g_" and "f_" prefix from genus names 
taxa.table$Genus <- gsub(taxa.table$Genus, pattern="f__|g__", replacement="")

#change all taxa name columns into factors
taxa.table <- data.frame(lapply(taxa.table, FUN=as.factor))

library(ape)
#set the taxa order 
taxa.order <- ~Kingdom/Phylum/Class/Order/Family/Genus
taxa.tree <- as.phylo.formula(taxa.order, taxa.table, collapse=F)
taxa.tree$edge.length <- rep(1, nrow(taxa.tree$edge))
plot(tree.1, node.label=TRUE)

#re-order the genera in the heatmap by the order of genera in the taxonomic cladogram tree
genera.taxa.tree <- taxa.tree$tip.label
spearman.correlations.wide.v5 <- arrange(spearman.correlations.wide,factor(genus_1, levels = genera.taxa.tree))
spearman.correlations.wide.v5 <- spearman.correlations.wide.v5[,-1]
spearman.correlations.wide.v6 <- spearman.correlations.wide.v5[,c(genera.taxa.tree)]
heatmap(x=as.matrix(spearman.correlations.wide.v6), scale="none", Rowv = taxa.tree, Colv = "Rowv")

Heatmap(matrix=as.matrix(spearman.correlations.wide.v6),cluster_rows= taxa.tree, cluster_columns = taxa.tree, show_heatmap_legend = TRUE, heatmap_legend_param = list(title="Spearman's rho"),row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize=10), row_names_max_width=unit(12,"cm"), column_names_max_height = unit(12,"cm"), heatmap_width = unit(1,"npc"), heatmap_height = unit(1,"npc"))
#saved as "microbiota.SpearmanCorrelation.heatmap.TaxonomyCladogram.jpeg"
