#This file contains the scripts/steps taken to analyze picrust2 output

#loading libraries
library(qiime2R)
library(phyloseq)
library(lme4)
library(lattice)
library(lmerTest)

#STEP 1: Read in predicted enzyme abundance (this dataset already has had the
#antibiotic-treated cows removed and has already had duplicate samples selected)
EC.metagenome <- read_qza("picrust2_output_verbose/ec_metagenome.qza")
EC.metagenome <- data.frame(EC.metagenome[["data"]])
EC.metagenome.v2 <- data.frame(t(EC.metagenome))
EC.metagenome.v2$sample_id <- rownames(EC.metagenome.v2)

#STEP 2: Calculate the abundance of beta galactosidase, EC: 3.2.1.23,
#relative to all other enzyme abundances, for each sample.

Function.Bgalactosidase.RelativeAbundance<-function(sample){
  galactosidase.dataframe<-matrix(nrow=1, ncol=1)
  for(i in sample){
    EC.metagenome.subset <- EC.metagenome.v2[EC.metagenome.v2$sample_id==i,]
    EC.metagenome.subset <- EC.metagenome.subset[,-c(2413)]
    galactosidase.relabund<-EC.metagenome.subset[["EC.3.2.1.23"]]/(rowSums(EC.metagenome.subset))
    galactosidase.dataframe<-rbind(galactosidase.dataframe, galactosidase.relabund)
  }
  print(galactosidase.dataframe)
}

sample <- EC.metagenome.v2$sample_id
beta_galactosidase_rel_abund <- Function.Bgalactosidase.RelativeAbundance(sample)
beta_galactosidase_rel_abund <- data.frame(beta_galactosidase_rel_abund)
beta_galactosidase_rel_abund <- beta_galactosidase_rel_abund[-1,]
beta_galactosidase_rel_abund <- cbind(beta_galactosidase_rel_abund, sample)
beta_galactosidase_rel_abund <- data.frame(beta_galactosidase_rel_abund)

#Make the relative abundance numeric
beta_galactosidase_rel_abund$beta_galactosidase_rel_abund <- as.numeric(beta_galactosidase_rel_abund$beta_galactosidase_rel_abund)

#STEP 3:Calculate the combined relative abundance of beta-galactosidase 
#AND 6-phospho-beta-galactosidase (EC 3.2.1.85).

Function.phosphoBgalactosidase.RelativeAbundance<-function(sample){
    galactosidase.dataframe<-matrix(nrow=1, ncol=1)
    for(i in sample){
        EC.metagenome.subset <- EC.metagenome.v2[EC.metagenome.v2$sample_id==i,]
        EC.metagenome.subset <- EC.metagenome.subset[,-c(2413)]
        galactosidase.relabund<-(EC.metagenome.subset[["EC.3.2.1.23"]]+EC.metagenome.subset[["EC.3.2.1.85"]])/(rowSums(EC.metagenome.subset))
        galactosidase.dataframe<-rbind(galactosidase.dataframe, galactosidase.relabund)
    }
    print(galactosidase.dataframe)
}

sample <- EC.metagenome.v2$sample_id
phospho_beta_galactosidase_rel_abund <- Function.phosphoBgalactosidase.RelativeAbundance(sample)
phospho_beta_galactosidase_rel_abund <- data.frame(phospho_beta_galactosidase_rel_abund)
phospho_beta_galactosidase_rel_abund <- phospho_beta_galactosidase_rel_abund[-1,]
phospho_beta_galactosidase_rel_abund <- cbind(phospho_beta_galactosidase_rel_abund, sample)
phospho_beta_galactosidase_rel_abund <- data.frame(phospho_beta_galactosidase_rel_abund)

#Make the relative abundance numeric
phospho_beta_galactosidase_rel_abund$phospho_beta_galactosidase_rel_abund <- as.numeric(phospho_beta_galactosidase_rel_abund$phospho_beta_galactosidase_rel_abund)

#STEP 4:Reading in the metadata table
metadata  <- read.csv("../../../DairyGrandChallenge.Study1.2017-2018/ResearchMaterials.DairyGrandChallenge/Cleaned_For_MetadataFile/metadata.csv")

#keeping only the data that is needed. 
metadata<- metadata[,c(1:15)]
metadata <- metadata[-1,]


#STEP 5: Merge the metadata and beta-galactosidase and phospho-betagalactosidase datasets
beta_galactosidase_rel_abund <- rename(beta_galactosidase_rel_abund, "sample.id"="sample")
beta_galactosidase_rel_abund$sample.id <- gsub(beta_galactosidase_rel_abund$sample.id, pattern="X", replacement="")
beta_galactosidase_rel_abund$sample.id <- gsub(beta_galactosidase_rel_abund$sample.id, pattern="[.]", replacement="-")

phospho_beta_galactosidase_rel_abund <- rename(phospho_beta_galactosidase_rel_abund, "sample.id"="sample")
phospho_beta_galactosidase_rel_abund$sample.id <- gsub(phospho_beta_galactosidase_rel_abund$sample.id, pattern="X", replacement="")
phospho_beta_galactosidase_rel_abund$sample.id <- gsub(phospho_beta_galactosidase_rel_abund$sample.id, pattern="[.]", replacement="-")

metadata.v2 <- merge(metadata, beta_galactosidase_rel_abund, by="sample.id", all=FALSE)
dim(metadata.v2) #all 334 samples retained
metadata.v2 <- merge(metadata.v2, phospho_beta_galactosidase_rel_abund, by="sample.id", all=FALSE)
dim(metadata.v2) #all 334 samples retained. 

#STEP 6: Read in the number of reads per sample, before the samples were rarefied. 
#And then merge with the metadata file
reads.per.sample <- read.csv("sample-frequency.csv", head=FALSE)
colnames(reads.per.sample) <- c("sample.id", "reads")

metadata.v3 <- merge(metadata.v2, reads.per.sample, by="sample.id", all=FALSE)
dim(metadata.v3) #all 334 samples retained

#remove "-B" from sample names which indicates the sample was a duplicate
#note: it's ok to remove "-B" because all samples were already deduplicated prior
#to rarefaction. 
metadata.v3$sample.id <- gsub(pattern="-B", x=metadata.v3$sample.id, replacement="")


#STEP 7: Read in Peggy Tomasula's macronutrient data (which contains the lactose measurements)
macronutrient <- readxl::read_xlsx("PMT DVH_Dairy Grand Challenge coded FOSS composition data_101818.xlsx")

#Looking at the first few lines of the file
head(macronutrient, 10) #The intended column headers are on the eighth row of the excel file

#In order to set the correct the column headers set, need to read in the dataset again but without the first row being used as the column headers
macronutrient <- readxl::read_xlsx("PMT DVH_Dairy Grand Challenge coded FOSS composition data_101818.xlsx", col_names=FALSE)

#Set the column headers to the values on the eighth row
colnames(macronutrient)<-macronutrient[8,]
#Remove the first eight rows
macronutrient <- macronutrient[-c(1:8),]

dim(macronutrient) #there are 76 observations and 46 variables in the table
colnames(macronutrient)
head(macronutrient) #The sample IDs are listed under the "Day" columns
#I can tell from the column names that the observations for days 2,3,4,5,6 are listed to the right of day 1 observations instead of below
#I can also see that the first two columns appear to be empty
sum(is.na(macronutrient[,1])) #all entries in the first column are empty/NA
sum(is.na(macronutrient[,2])) #all entries in the second column are empty/NA

#Remove the first two columns
macronutrient <- macronutrient[,-c(1:2)]
colnames(macronutrient) #no longer see those empty column headers

#Remove the "Period 1", "Period 2", "Cow ID", "Treatment 0", "treatment 1", "Treatment 2", "Stolen from", and "Location" columns because these variables will be retrieved from the Kalscheur lab data or at not needed in the compiled dataset
macronutrient <- macronutrient[,-c(1:4, 17, 30, 43, 44)]
head(macronutrient)
colnames(macronutrient)

#Need to subset the by day (i.e. take subset chunks of columns), then add back together under the same column headers so that there are only 
macronutrient.Day1 <- macronutrient[,1:6]
macronutrient.Day2 <- macronutrient[,7:12]
macronutrient.Day3 <- macronutrient[,13:18]
macronutrient.Day4 <- macronutrient[,19:24]
macronutrient.Day5 <- macronutrient[,25:30]
macronutrient.Day6 <- macronutrient[,31:36]

#Change the "Day" column headers in the data subsets to "sample_id"
macronutrient.Day1 <- rename(macronutrient.Day1, sample.id=`Day 1`)
macronutrient.Day2 <- rename(macronutrient.Day2, sample.id=`Day 2`)
macronutrient.Day3 <- rename(macronutrient.Day3, sample.id=`Day 3`)
macronutrient.Day4 <- rename(macronutrient.Day4, sample.id=`Day 4`)
macronutrient.Day5 <- rename(macronutrient.Day5, sample.id=`Day  5`)
macronutrient.Day6 <- rename(macronutrient.Day6, sample.id=`Day 6`)

#All of the column headers should be the same now for the subsets fo the macronutrient data
#Checking to see that they are the same..
colnames(macronutrient.Day1)==colnames(macronutrient.Day2) #All column headers are the same
colnames(macronutrient.Day1)==colnames(macronutrient.Day3) #All column headers are the same
colnames(macronutrient.Day1)==colnames(macronutrient.Day4) #All column headers are the same
colnames(macronutrient.Day1)==colnames(macronutrient.Day5) #All column headers are the same
colnames(macronutrient.Day1)==colnames(macronutrient.Day6) #All column headers are the same

#Now add the subsets together by rows under the same set of 6 column headers
macronutrient <- rbind(macronutrient.Day1, macronutrient.Day2, macronutrient.Day3, macronutrient.Day4, macronutrient.Day5, macronutrient.Day6)
dim(macronutrient) #456 observations and 6 variables (this is what's expected with 76 observations in each subset)

#Check for NA values and duplicate values in the sample_id column
sum(is.na(macronutrient$sample.id)) #the output is zero, so no NA values among the sample_id column
macronutrient.distinct.col1 <- unique(macronutrient$sample.id)
length(macronutrient.distinct.col1) #456 unique observations in sample_id so there are no duplicate sample IDs

#Remove the destination information ("_Wyn") from the sample IDs so that they can be matched and merged to other datasets
macronutrient$sample.id <-gsub(pattern="_Wyn", x=macronutrient$sample.id, replacement="") 
head(macronutrient)

#Add a column specifying the units of measure for the macronutrient data
macronutrient$nut_Units <-"gram/100 gram of milk "

#STEP 8: Merge lactose data with metadata. 
#Since the dates are in the sample IDs for macronutrients dataset, 
#I will change them to D1, D2, D3, D4, D5, and D6, so that they match up with 
#the metadata sampe.id's 
macronutrient$sample.id <- gsub(replacement="-D1", pattern="_171125", macronutrient$sample.id)
macronutrient$sample.id <- gsub(replacement="-D2", pattern="_171126", macronutrient$sample.id)
macronutrient$sample.id <- gsub(replacement="-D3", pattern="_180104", macronutrient$sample.id)
macronutrient$sample.id <- gsub(replacement="-D4", pattern="_180105", macronutrient$sample.id)
macronutrient$sample.id <- gsub(replacement="-D5", pattern="_180314", macronutrient$sample.id)
macronutrient$sample.id <- gsub(replacement="-D6", pattern="_180315", macronutrient$sample.id)

metadata.v4 <- merge(metadata.v3, macronutrient, by="sample.id", all=FALSE)
dim(metadata.v4) #all 334 samples retained. 


#STEP 9: Rename sample 5849-D1 to 5849-D2 because cow 5849 actually missed
#the milk collection on November 25, 2017, but milk was collected from this cow
#on November 26, 2017. The same was true for cow 6233. 
metadata.v4$sample.id <- gsub(pattern="5849-D1", replacement="5849-D2", x=metadata.v4$sample.id)
metadata.v4$sample.id <- gsub(pattern="6233-D1", replacement="6233-D2", x=metadata.v4$sample.id)
metadata.v4[metadata.v4$sample.id=="5849-D2","day"] <- 2
metadata.v4[metadata.v4$sample.id=="6233-D2","day"] <- 2

#STEP 10: Remove cows with less than one sample per period.
cows.pd0 <- unique(metadata.v4[metadata.v4$period=="p_0",][["cow"]])
cows.pd1 <- unique(metadata.v4[metadata.v4$period=="p_1",][["cow"]])
cows.pd2 <- unique(metadata.v4[metadata.v4$period=="p_2",][["cow"]])
cows.to.keep <- cows.pd0[which(cows.pd0 %in% cows.pd1)]
cows.to.keep.v2 <- cows.to.keep[which(cows.to.keep %in% cows.pd2)]
length(cows.to.keep.v2)#60 cows
metadata.v5 <- metadata.v4[which(metadata.v4$cow %in% cows.to.keep.v2),]
dim(metadata.v5) #313 samples
length(unique(metadata.v5$cow)) #60 cows

#STEP 11: Some cows have more than one sample in each period, 
#so selecting the sample from the period that has the higher # of reads 
#(before rarefying)

ChooseSampleFromPeriod <-function(cow) {
    NoDuplicates <- data.frame()
    for(i in cow) {
        select.cow <- metadata.v5[metadata.v5$cow==i,]
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
cows <-unique(metadata.v5$cow)
metadata.v6 <- ChooseSampleFromPeriod(cows) 
length(metadata.v6$sample.id) #180 samples
length(unique(metadata.v6$cow)) #60 cows

#STEP 12: Remove cows that the Barile lab identified as lactose outliers (cow #6229 & #5651)
metadata.v7 <- metadata.v6[metadata.v6$cow!=5651,]
metadata.v7 <- metadata.v7[metadata.v7$cow!=6229,]
dim(metadata.v7) #174 samples
length(unique(metadata.v7$cow)) #58 cows

#STEP 13: Linear mixed effect modeling
hist(metadata.v7$beta_galactosidase_rel_abund) #slightly right skewed
hist(metadata.v7$phospho_beta_galactosidase_rel_abund) #slightly right skewed
metadata.v7$Lactose <- as.numeric(metadata.v7$Lactose)
hist(metadata.v7$Lactose) #a little left skewed

summary(lmer(formula=beta_galactosidase_rel_abund~Lactose+(1|cow), data=metadata.v7)) #singular fit warning, not significant
allFit(lmer(formula=beta_galactosidase_rel_abund~Lactose+(1|cow), data=metadata.v7))#singular fit warning, not significant
summary(lmer(formula=phospho_beta_galactosidase_rel_abund~Lactose+(1|cow), data=metadata.v7))#singular fit warning, not significant
allFit(lmer(formula=phospho_beta_galactosidase_rel_abund~Lactose+(1|cow), data=metadata.v7))#singular fit warning, not significant

#STEP 14: Repeated measures correlation
metadata.v7$cow <- as.factor(metadata.v7$cow)
rmcorr(participant=cow, measure1=Lactose, measure2=beta_galactosidase_rel_abund, dataset=metadata.v7) #not signficant
rmcorr(participant=cow, measure1=Lactose, measure2=phospho_beta_galactosidase_rel_abund, dataset=metadata.v7) #not significant

#STEP 15: Performing a square root transformation on the b-galactosidase, 
#phospho-b-galactosidase, and lactose abundances to make data more normally distributed. 
metadata.v7$sqrt.beta_galactosidase_rel_abund <- sqrt(metadata.v7$beta_galactosidase_rel_abund)
hist(metadata.v7$sqrt.beta_galactosidase_rel_abund) #normal distribution now

metadata.v7$sqrt.phospho_beta_galactosidase_rel_abund <- sqrt(metadata.v7$phospho_beta_galactosidase_rel_abund)
hist(metadata.v7$sqrt.phospho_beta_galactosidase_rel_abund) #normal distribution now

metadata.v7$sqrt.Lactose <- sqrt(metadata.v7$Lactose)
hist(metadata.v7$sqrt.Lactose) #not as skewed now

#STEP 16: Redoing linear mixed effects modeling and Repeated measures correlation
#with the square root transformed values
summary(lmer(sqrt.beta_galactosidase_rel_abund~sqrt.Lactose+(1|cow), data=metadata.v7)) #singular fit warning,not significant
allFit(lmer(sqrt.beta_galactosidase_rel_abund~sqrt.Lactose+(1|cow), data=metadata.v7)) #singular fit warning,not significant
summary(lmer(sqrt.beta_galactosidase_rel_abund~sqrt.Lactose+(Lactose|cow), data=metadata.v7)) #singular fit warning, not significant
allFit(lmer(sqrt.beta_galactosidase_rel_abund~sqrt.Lactose+(Lactose|cow), data=metadata.v7)) #singular fit warning, not significant

rmcorr(participant=cow, measure1=sqrt.Lactose, measure2=sqrt.beta_galactosidase_rel_abund, dataset=metadata.v7) #not signficant
rmcorr(participant=cow, measure1=sqrt.Lactose, measure2=sqrt.phospho_beta_galactosidase_rel_abund, dataset=metadata.v7) #not significant
