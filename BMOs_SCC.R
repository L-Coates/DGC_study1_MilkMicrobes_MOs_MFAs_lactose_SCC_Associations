#This file describes the steps taken to analyze the relationships
#between BMOs and SCC. 

#Loading the necessary libraries 
library(dplyr)
library(lmerTest)
library(lme4)
library(ggplot2)
library(vegan)
library(rmcorr)

#STEP 1: Reading in datasets and merging datasets. 

##Reading in the AgSource data for the SCC measures from 11/21/17, 11/28/17, 01/02/18, 01/05/18, 03/13/18, and 03/19/18 (which are the same agsource SCC measures that were used for the analysis of SCC & microbiota data.)
Nov.21.17 <-readxl::read_xlsx("../../ResearchMaterials.DairyGrandChallenge/Materials.FromKableLab/DairyGrandChallenge_SCCdata/KFK08 112017 112117 - Nov 24.xlsx", col_names=FALSE, col_types=c("date", "text", "skip", "skip", "skip", "skip", "numeric", "skip", "text", "skip", "skip", "skip", "skip"))
Nov.28.17  <-readxl::read_xlsx("../../ResearchMaterials.DairyGrandChallenge/Materials.FromKableLab/DairyGrandChallenge_SCCdata/KFK08 112717 112817 - Dec 1.xlsx", col_names=FALSE, col_types=c("date", "text", "skip", "skip", "skip", "skip", "numeric", "skip", "text", "skip", "skip", "skip", "skip"))
Jan.02.18  <-readxl::read_xlsx("../../ResearchMaterials.DairyGrandChallenge/Materials.FromKableLab/DairyGrandChallenge_SCCdata/KFK08 010218.xlsx", col_names=FALSE, col_types=c("date", "text", "skip", "skip", "skip", "skip", "numeric", "skip", "text", "skip", "skip", "skip", "skip"))
Jan.05.18  <-readxl::read_xlsx("../../ResearchMaterials.DairyGrandChallenge/Materials.FromKableLab/DairyGrandChallenge_SCCdata/KFK08 010518.xlsx", col_names=FALSE, col_types=c("date", "text", "skip", "skip", "skip", "skip", "numeric", "skip", "text", "skip", "skip", "skip", "skip"))
March.13.18  <-readxl::read_xlsx("../../ResearchMaterials.DairyGrandChallenge/Materials.FromKableLab/DairyGrandChallenge_SCCdata/KFK08 031218 031318 - March 16.xlsx", col_names=FALSE, col_types=c("date", "text", "skip", "skip", "skip", "skip", "numeric", "skip", "text", "skip", "skip", "skip", "skip", "skip", "skip", "skip", "skip"))
March.19.18  <-readxl::read_xlsx("../../ResearchMaterials.DairyGrandChallenge/Materials.FromKableLab/DairyGrandChallenge_SCCdata/KFK08 031918 032018 - March 23.xlsx", col_names=FALSE, col_types=c("date", "text", "skip", "skip", "skip", "skip", "numeric", "skip", "text", "skip", "skip", "skip", "skip"))

#Renaming columns
colnames(Nov.21.17)<- c("date", "cow", "SCCx1000.milking1.beforesample", "milking.number")
colnames(Nov.28.17)<- c("date", "cow", "SCCx1000.milking1.aftersample", "milking.number")
colnames(Jan.02.18)<- c("date", "cow", "SCCx1000.milking1.beforesample", "milking.number")
colnames(Jan.05.18)<- c("date", "cow", "SCCx1000.milking1.aftersample", "milking.number")
colnames(March.13.18)<- c("date", "cow", "SCCx1000.milking1.beforesample", "milking.number")
colnames(March.19.18)<- c("date", "cow", "SCCx1000.milking1.aftersample", "milking.number")

#Removing the first nine rows which are empty. 
Nov.21.17 <- Nov.21.17[-c(1:9),]
Nov.28.17 <- Nov.28.17[-c(1:9),]
Jan.02.18 <- Jan.02.18[-c(1:9),]
Jan.05.18 <- Jan.05.18[-c(1:9),]
March.13.18 <- March.13.18[-c(1:9),]
March.19.18 <- March.19.18[-c(1:9),]

#Selecting only measurements for November 21, 2017
Nov.21.17$date <-as.Date(Nov.21.17$date)
Nov.21.17 <- Nov.21.17[Nov.21.17$date=="2017-11-21",]
length(unique(Nov.21.17$cow)) #76 cows

#Selecting only measurements for Nov. 28, 2017
Nov.28.17$date <-as.Date(Nov.28.17$date)
Nov.28.17 <- Nov.28.17[Nov.28.17$date=="2017-11-28",]
length(unique(Nov.28.17$cow)) #76 cows

#selecting only measurements for Jan 2, 2018
Jan.02.18$date <-as.Date(Jan.02.18$date)
Jan.02.18 <- Jan.02.18[Jan.02.18$date=="2018-01-02",]
length(unique(Jan.02.18$cow)) #76 cows

#selecting only measurements for Jan 5, 2018
Jan.05.18$date <-as.Date(Jan.05.18$date)
Jan.05.18 <- Jan.05.18[Jan.05.18$date=="2018-01-05",]
length(unique(Jan.05.18$cow)) #76 cows

#selecting only measurements for march 13, 2018
March.13.18$date <-as.Date(March.13.18$date)
March.13.18 <- March.13.18[March.13.18$date=="2018-03-13",]
length(unique(March.13.18$cow)) #74 cows

#selecting only measurements for march 19, 2018
March.19.18$date <-as.Date(March.19.18$date)
March.19.18 <- March.19.18[March.19.18$date=="2018-03-19",]
length(unique(March.19.18$cow)) #74 cows

##Selecting only the M1 measurements (these are from the first milking of the day) from Nov 21st (this is the most complete dataset taken closest to the date of collection of the milk samples used for microbiota analysis.)
Nov.21.17 <- Nov.21.17[Nov.21.17$milking.number=="M1",]
Nov.21.17 <- Nov.21.17[, c("cow", "SCCx1000.milking1.beforesample")]
length(unique(Nov.21.17$cow)) 

Nov.28.17 <- Nov.28.17[Nov.28.17$milking.number=="M1",]
Nov.28.17 <- Nov.28.17[, c("cow", "SCCx1000.milking1.aftersample")]
length(unique(Nov.28.17$cow)) 

Jan.02.18 <- Jan.02.18[Jan.02.18$milking.number=="M1",]
Jan.02.18 <- Jan.02.18[, c("cow", "SCCx1000.milking1.beforesample")]
length(unique(Jan.02.18$cow)) 

Jan.05.18 <- Jan.05.18[Jan.05.18$milking.number=="M1",]
Jan.05.18 <- Jan.05.18[, c("cow", "SCCx1000.milking1.aftersample")]
length(unique(Jan.05.18$cow)) 

March.13.18 <- March.13.18[March.13.18$milking.number=="M1",]
March.13.18 <- March.13.18[, c("cow", "SCCx1000.milking1.beforesample")]
length(unique(March.13.18$cow)) 

March.19.18 <- March.19.18[March.19.18$milking.number=="M1",]
March.19.18 <- March.19.18[, c("cow", "SCCx1000.milking1.aftersample")]
length(unique(March.19.18$cow))

#reading in the BMO file
BMO <- readxl::read_xlsx("USDA DGC Project 1 - BMO data (Barile Lab) - 11.17.2020 - SDD.xlsx", sheet=1)

#removing a few columns that won't be used
BMO <- BMO[,-c(2:7)]

#reformatting the sample IDs
BMO$`Sample ID 1` <- gsub(pattern="_Davis", x=BMO$`Sample ID 1`, replacement="")
BMO$`Sample ID 1` <- gsub(pattern="_171125", x=BMO$`Sample ID 1`, replacement="-D1")
BMO$`Sample ID 1` <- gsub(pattern="_171126", x=BMO$`Sample ID 1`, replacement="-D2")
BMO$`Sample ID 1` <- gsub(pattern="_180104", x=BMO$`Sample ID 1`, replacement="-D3")
BMO$`Sample ID 1` <- gsub(pattern="_180105", x=BMO$`Sample ID 1`, replacement="-D4")
BMO$`Sample ID 1` <- gsub(pattern="_180314", x=BMO$`Sample ID 1`, replacement="-D5")
BMO$`Sample ID 1` <- gsub(pattern="_180315", x=BMO$`Sample ID 1`, replacement="-D6")

BMO<-rename(BMO, sample.id = `Sample ID 1`)

#remove BMO 3_3_1_0_0 in column 11 because it has "NA" values
BMO <- BMO[,-11]
dim(BMO) #348 rows, 20 colums

#Put in a column "cow" and column "period"
BMO$cow <- gsub(BMO$sample.id, pattern="-D[1-6]", replacement="")
BMO$period <- BMO$sample.id
BMO$period <- gsub(BMO$period, pattern="....-D[1-2]", replacement="0")
BMO$period <- gsub(BMO$period, pattern="....-D[3-4]", replacement="1")
BMO$period <- gsub(BMO$period, pattern="....-D[5-6]", replacement="2")
length(unique(BMO$cow)) #62 cows
length(BMO$sample.id) #348 samples

#Merge the BMO and SCC data by period. 
BMO.pd0 <- BMO[BMO$period=="0",]
BMO.pd1 <- BMO[BMO$period=="1",]
BMO.pd2 <- BMO[BMO$period=="2",]

SCC.pd0 <- merge(Nov.21.17, Nov.28.17, by="cow", all=FALSE)
BMO.pd0 <- merge(BMO.pd0, SCC.pd0, by="cow", all=FALSE)
length(unique(BMO.pd0$cow)) #61 cows

SCC.pd1 <- merge(Jan.02.18, Jan.05.18, by="cow", all=FALSE)
BMO.pd1 <- merge(BMO.pd1, SCC.pd1, by="cow", all=FALSE)
length(unique(BMO.pd1$cow)) #58 cows

SCC.pd2 <- merge(March.13.18, March.19.18, by="cow", all=FALSE)
BMO.pd2 <- merge(BMO.pd2, SCC.pd2, by="cow", all=FALSE)
length(unique(BMO.pd2$cow)) #55 cows

BMO.v2 <- rbind(BMO.pd0, BMO.pd1, BMO.pd2)
length(unique(BMO.v2$cow)) #62 cows
length(BMO.v2$sample.id) #348 samples

#STEP 2: Filtering dataset and preparing for correlation analysis 

#Remove cows that don't have a sample from each period and select one sample per 
#period per cow since samples within the same period were combined for BMO 
#analysis and since SCC measures don't match up to the same days as the BMO milk
#sample days. 
samples.per.cow<-data.frame(table(BMO.v2$cow))
colnames(samples.per.cow) <- c("cow", "samples")
cows.to.keep <- samples.per.cow[samples.per.cow$samples==6,]$cow
length(cows.to.keep) #54 cows
BMO.v3 <- BMO.v2[c(BMO.v2$cow%in%cows.to.keep),]
BMO.v3 <- BMO.v3[grep(pattern="...-D1|....-D3|....-D5", BMO.v3$sample.id),]
length(unique(BMO.v3$cow)) #54 cows

#Multiply SCC by 1000 since it's currently in units of thousands, and then
#log transform to achieve more normal distribution. 
BMO.v3$LogSCC.milking1.beforesample <- log10(1000*(BMO.v3$SCCx1000.milking1.beforesample))
BMO.v3$LogSCC.milking1.aftersample <- log10(1000*(BMO.v3$SCCx1000.milking1.aftersample))

#confirming that none of the antibiotic-treated cows are in the BMO dataset
antibiotic.cows <- c("5020", "5212", "5297", "5405", "5697", "6076", "6215", "6232", "6233", "6241")
BMO.v3$cow %in% antibiotic.cows #None are in the BMO dataset

#Calculating BMO alpha diversity (shannon and inverse simpson)
BMO.v3[,3:21]<-lapply(BMO.v3[,3:21], as.numeric)
BMO.v4 <- BMO.v3
BMO.v4$BMO.diversity.inversesimpson<-diversity(x=BMO.v4[,3:21], index="invsimpson")
BMO.v4$BMO.diversity.shannon<-diversity(x=BMO.v4[,3:21], index="shannon")

#Change BMO column names to those that can be recognized in a lmer function
#in a for loop. 
colnames(BMO.v4)[3:21]<-gsub(pattern="'", replacement="", colnames(BMO.v4[3:21]))
colnames(BMO.v4)[3:21]<-gsub(pattern="-", replacement="_", colnames(BMO.v4[3:21]))
colnames(BMO.v4)[3:21]<-paste0("X", colnames(BMO.v4[3:21]))

#STEP 3: Running correlation analyses between BMOs and logSCC, and BMO
#alpha diversity and logSCC. 

#First, since some of the BMOs are skewed, the BMOs will be log-transformed
BMO.v4[,3:21] <- log10(BMO.v4[,3:21])

#Function to run a linear mixed effect model testing log(SCC) as a predictor 
#of each BMO with the intercept varying by cow.
Function.BMO.SCC.v1 <- function(BMO){
  coefficient.pvalue <- data.frame(coefficient=NA, pvalue=NA, BMO=NA)
  for(i in BMO){
    form<-as.formula(paste0(i,"~LogSCC.milking1.beforesample+(1|cow)"))
    lmer.output<-try(summary(lmer(form, data=BMO.v4)), silent=TRUE)
    lmer.output2<-try(data.frame(coefficient=lmer.output[["coefficients"]][2,1], pvalue=lmer.output[["coefficients"]][2,5], BMO=i), silent=TRUE)
    coefficient.pvalue <- rbind(coefficient.pvalue, lmer.output2)  
  }
  print(coefficient.pvalue)
}

BMO.names <- colnames(BMO.v4[3:21])
BMO.SCC.correlations.v1<- Function.BMO.SCC.v1(BMO.names)
BMO.SCC.correlations.v1<- BMO.SCC.correlations.v1[-1,]
BMO.SCC.correlations.v1$fdr <- p.adjust(BMO.SCC.correlations.v1$pvalue, method="fdr")
#After multiple test correction, only BMO 8_0_0_0_0 had a relationship with logSCC. 

#Function to run a linear mixed effect model testing each BMO as a predictor 
#and logSCC as the response variable, with the intercept varying by cow.
Function.BMO.SCC.v2 <- function(BMO){
coefficient.pvalue <- data.frame(coefficient=NA, pvalue=NA, BMO=NA)
for(i in BMO){
  form<-as.formula(paste0("LogSCC.milking1.aftersample~",i,"+(1|cow)"))
  lmer.output<-try(summary(lmer(form, data=BMO.v4)), silent=TRUE)
  lmer.output2<-try(data.frame(coefficient=lmer.output[["coefficients"]][2,1], pvalue=lmer.output[["coefficients"]][2,5], BMO=i), silent=TRUE)
  coefficient.pvalue <- rbind(coefficient.pvalue, lmer.output2)  
}
print(coefficient.pvalue)
}

BMO.SCC.correlations.v2<- Function.BMO.SCC.v2(BMO.names)
BMO.SCC.correlations.v2<- BMO.SCC.correlations.v2[-1,]
BMO.SCC.correlations.v2$fdr <- p.adjust(BMO.SCC.correlations.v2$pvalue, method="fdr") #positive relationship between logSCC and 6'-SL, after FDR correction.  


#Testing for BMO alpha diversity predictor of logSCC
summary(lmer(LogSCC.milking1.aftersample~BMO.diversity.shannon+(1|cow), data=BMO.v4)) #not significant predictor
summary(lmer(LogSCC.milking1.aftersample~BMO.diversity.inversesimpson+(1|cow), data=BMO.v4)) #not significant predictor
summary(lmer(LogSCC.milking1.aftersample~BMO.diversity.shannon+(BMO.diversity.shannon|cow), data=BMO.v4)) #singular fit warning
summary(lmer(LogSCC.milking1.aftersample~BMO.diversity.inversesimpson+(BMO.diversity.inversesimpson|cow), data=BMO.v4)) #warning -- failure to converge
#following up with repeat measures correlation
hist(BMO.v4$LogSCC.milking1.aftersample)
hist(BMO.v4$BMO.diversity.inversesimpson) #normal distribution
hist(BMO.v4$BMO.diversity.shannon) #a little left skewed

BMO.v4$cow <- as.factor(BMO.v4$cow)
rmcorr(measure1=LogSCC.milking1.beforesample, measure2=BMO.diversity.shannon, participant=cow, dataset=BMO.v4)
rmcorr(measure1=LogSCC.milking1.beforesample, measure2=BMO.diversity.inversesimpson, participant=cow, dataset=BMO.v4)

#testing for logSCC predictor of BMO alpha diversity
summary(lmer(BMO.diversity.shannon~LogSCC.milking1.beforesample+(1|cow), data=BMO.v4)) #not significant predictor
summary(lmer(BMO.diversity.inversesimpson~LogSCC.milking1.beforesample+(1|cow), data=BMO.v4)) #not significant predictor
summary(lmer(BMO.diversity.shannon~LogSCC.milking1.beforesample+(LogSCC.milking1.beforesample|cow), data=BMO.v4)) #warning -- failure to converge
summary(lmer(BMO.diversity.inversesimpson~LogSCC.milking1.beforesample+(LogSCC.milking1.beforesample|cow), data=BMO.v4)) #singular fit warning
#following up with repeat measures correlation
hist(BMO.v4$LogSCC.milking1.beforesample) #somewhat normal distribution
rmcorr(measure1=LogSCC.milking1.aftersample, measure2=BMO.diversity.shannon, participant=cow, dataset=BMO.v4) #no significant correlation
rmcorr(measure1=LogSCC.milking1.aftersample, measure2=BMO.diversity.inversesimpson, participant=cow, dataset=BMO.v4) #no significant correlation


#Now, running repeated measures correlation
#Trying repeated measures correlation (RMCORR)
Function.RMCORR.BMO.SCC.v1 <- function(BMOs){
  coefficient.pvalue <- data.frame(coefficient=NA, pvalue=NA, BMO=NA)
  for(i in BMOs){
    rmcorr.output<-try(rmcorr(measure1=i, measure2=LogSCC.milking1.beforesample, participant=cow, dataset=BMO.v4), silent=TRUE)
    rmcorr.output2<-try(data.frame(coefficient=rmcorr.output[["r"]], pvalue=rmcorr.output[["p"]], BMO=i), silent=TRUE)
    coefficient.pvalue <- rbind(coefficient.pvalue, rmcorr.output2)  
  }
  print(coefficient.pvalue)
}

RMCORR.BMO.SCC.correlations.BeforeSample <- Function.RMCORR.BMO.SCC.v1(BMO.names)
RMCORR.BMO.SCC.correlations.BeforeSample.v2 <- RMCORR.BMO.SCC.correlations.BeforeSample[-1,]
RMCORR.BMO.SCC.correlations.BeforeSample.v2$FDR <- p.adjust(p=RMCORR.BMO.SCC.correlations.BeforeSample.v2$pvalue, method="fdr") #no significant relationships

Function.RMCORR.BMO.SCC.v2 <- function(BMOs){
  coefficient.pvalue <- data.frame(coefficient=NA, pvalue=NA, BMO=NA)
  for(i in BMOs){
    rmcorr.output<-try(rmcorr(measure1=i, measure2=LogSCC.milking1.aftersample, participant=cow, dataset=BMO.v4), silent=TRUE)
    rmcorr.output2<-try(data.frame(coefficient=rmcorr.output[["r"]], pvalue=rmcorr.output[["p"]], BMO=i), silent=TRUE)
    coefficient.pvalue <- rbind(coefficient.pvalue, rmcorr.output2)  
  }
  print(coefficient.pvalue)
}

RMCORR.BMO.SCC.correlations.AfterSample <- Function.RMCORR.BMO.SCC.v2(BMO.names)
RMCORR.BMO.SCC.correlations.AfterSample.v2 <- RMCORR.BMO.SCC.correlations.AfterSample[-1,]
RMCORR.BMO.SCC.correlations.AfterSample.v2$FDR <- p.adjust(p=RMCORR.BMO.SCC.correlations.AfterSample.v2$pvalue, method="fdr") #no significant relationships

#STEP 4: Graph significant correlation(s) and determine the fit of the model. 

#graph with the regression line by cow
results.6SL <- lmer(LogSCC.milking1.aftersample~X6_SL+(1|cow), data=BMO.v4)
ggplot(aes(x=X6_SL, y=LogSCC.milking1.aftersample), data=BMO.v4)+geom_point()+ylab("log(SCC)")+xlab("log(6'-sialyllactose)")+geom_line(aes(y=predict(results.6SL), group=cow))
ggsave(filename="log6SL.vs.LogSCC.regressionline.jpeg", device="jpeg")

#graph with the regression line by cow
results.80000 <- lmer(X8_0_0_0_0~LogSCC.milking1.beforesample+(1|cow), data=BMO.v4)
ggplot(aes(y=X8_0_0_0_0, x=LogSCC.milking1.beforesample), data=BMO.v4)+geom_point()+ylab("log(8_0_0_0_0)")+xlab("log(SCC)")+geom_line(aes(y=predict(results.80000), group=cow))+annotate("text",x=6.25, y=-0.5, label="0.137")
ggsave(filename="LogSCC.vs.8_0_0_0_0.regressionline.jpeg", device="jpeg", width=10, height=5)

#determine if the models are good fits
#Do the models meet the assumptions?
#Plotting the fitted versus residuals
plot(fitted(results.6SL), resid(results.6SL))
plot(fitted(results.80000), residuals(results.80000))

#The points showed a consistent and fairly even spread around the residual=0 line
#so this suggests that the assumptions of linearity and equal variances has not been violated.

#Looking for normality of residuals with a histogram
hist(resid(results.6SL)) #residuals look like they have fairly normal distribution
hist(residuals(results.80000)) #close to normal distribution. 
#Looking at QQ plot to confirm that the observations come from a normal distribution
qqmath(results.6SL) #observations mostly fall on the line
qqmath(results.80000) #observations mostly fall on the line
