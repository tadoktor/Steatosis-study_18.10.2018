
###################################################################
####Author:  Tatyana Doktorova
####Version: 1.0 08-28-2018
####Description: Differential gene list expression script: 
####            - filtering according to p-value and fold change
####
#####################################################################


library("biomaRt")
library("dplyr")


Root_location <- "C:\\Users\\tatyana\\Documents\\Steatosis"


##Output directory
outDir <- paste(Root_location, "/IntermediateData/", sep='')

#Generation of Steatosis list of genes

TG_Gates_rat_log2fold_steatosis<- read.table("file:///C:/Users/tatyana/Documents/Advance project/Oana/TGG_data_log2fold_Steatotic_rat.txt", header=TRUE)
steatosis_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Advance project/Oana/TGG_in_Ste_samples.txt", header=TRUE)


steatosis_all<-merge(TG_Gates_rat_log2fold_steatosis,steatosis_sampleID, by.x="sampleId",by.y="sampleId")

steatosis_filtering<-steatosis_all[!is.na(steatosis_all$pvalue),]

steatosis_pval0.05<- steatosis_filtering[steatosis_filtering$pvalue<0.001,]
steatosis_pval0.05_FC_up_1<-steatosis_pval0.05[steatosis_pval0.05$ value > 1,]
steatosis_pval0.05_FC_down_1<-steatosis_pval0.05[steatosis_pval0.05$ value < -1,]
steatosis_significant<- rbind(steatosis_pval0.05_FC_down_1,steatosis_pval0.05_FC_up_1)
write.csv(steatosis_significant, file="steatosis_specific.csv")

#Generation of non-steatosis list of genes
TG_Gates_rat_log2fold_nonsteatosis<- read.table("file:///C:/Users/tatyana/Documents/Advance project/Oana/TGG_nonsteatosis_data_log2fold_rat.txt", header=TRUE)
nonsteatosis_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Advance project/Oana/TGG_in_nonsteatosis_samples.txt", header=TRUE)

nonsteatosis_all<-merge(TG_Gates_rat_log2fold_nonsteatosis,nonsteatosis_sampleID, by.x="sampleId",by.y="sampleId")

nonsteatosis_filtering<-nonsteatosis_all[!is.na(nonsteatosis_all$pvalue),]
                                             
nonsteatosis_pval0.05<- nonsteatosis_filtering[nonsteatosis_filtering$pvalue<0.001,]
nonsteatosis_pval0.05_FC_up_1<-nonsteatosis_pval0.05[nonsteatosis_pval0.05$ value > 1,]
nonsteatosis_pval0.05_FC_down_1<-nonsteatosis_pval0.05[nonsteatosis_pval0.05$ value < -1,]
nonsteatosis_significant<- rbind(nonsteatosis_pval0.05_FC_down_1,nonsteatosis_pval0.05_FC_up_1)
write.csv(nonsteatosis_significant, file="nonsteatosis_specific.csv")




