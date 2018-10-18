
library("biomaRt")
library("dplyr")
library(collapsibleTree)

#Generation of GTX list of genes

TG_Gates_human_log2fold_GTX<- read.table("file:///C:/Users/tatyana/Documents/Advance project/TGG_HCC_data_log2fold_human_GTX.txt", header=TRUE)
GTX_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Advance project/TGG_in_HCC_samples_GTX.txt", header=TRUE)
HCC_Toxcast_CTD<-read.csv("file:///C:/Users/tatyana/Documents/Advance project/HCC_Toxcast_CTD.csv", header=TRUE)
Unique_genes = HCC_Toxcast_CTD %>% distinct(lhs)

GTX_all<-merge(TG_Gates_human_log2fold_GTX,GTX_sampleID, by.x="sampleId",by.y="sampleId")

GTX_filtering<-GTX_all[!is.na(GTX_all$pvalue),]

GTX_pval0.05<- GTX_filtering[GTX_filtering$pvalue<0.05,]
GTX_pval0.05_FC_up_1<-GTX_pval0.05[GTX_pval0.05$ value > 0.58,]
GTX_pval0.05_FC_down_1<-GTX_pval0.05[GTX_pval0.05$ value < -0.58,]
GTX_significant<- rbind(GTX_pval0.05_FC_down_1,GTX_pval0.05_FC_up_1)
write.csv(GTX_significant, file="GTX_specific.csv")

#Generation of NGTX list of genes
TG_Gates_human_log2fold_NGTX<- read.table("file:///C:/Users/tatyana/Documents/Advance project/TGG_HCC_data_log2fold_human_NGTX.txt", header=TRUE)
NGTX_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Advance project/TGG_in_HCC_samples_NGTX.txt", header=TRUE)

NGTX_all<-merge(TG_Gates_human_log2fold_NGTX,NGTX_sampleID, by.x="sampleId",by.y="sampleId")

NGTX_filtering<-NGTX_all[!is.na(NGTX_all$pvalue),]
                                             
NGTX_pval0.05<- NGTX_filtering[NGTX_filtering$pvalue<0.05,]
NGTX_pval0.05_FC_up_1<-NGTX_pval0.05[NGTX_pval0.05$ value > 0.58,]
NGTX_pval0.05_FC_down_1<-NGTX_pval0.05[NGTX_pval0.05$ value < -0.58,]
NGTX_significant<- rbind(NGTX_pval0.05_FC_down_1,NGTX_pval0.05_FC_up_1)
write.csv(NGTX_significant, file="NGTX_specific.csv")

#Generation of NC list of genes
TG_Gates_human_log2fold_NC<- read.table("file:///C:/Users/tatyana/Documents/Advance project/TGG_HCC_data_log2fold_human_NC.txt", header=TRUE)
NC_sampleID<- read.table("file:///C:/Users/tatyana/Documents/Advance project/TGG_in_HCC_samples_NC.txt", header=TRUE)
#merge the two sources and remove NA values
NC_all<-merge(TG_Gates_human_log2fold_NC,NC_sampleID, by.x="sampleId",by.y="sampleId")

NC_filtering<-NC_all[!is.na(NC_all$pvalue),]

NC_pval0.05<- NC_filtering[NC_filtering$pvalue<0.05,]
NC_pval0.05_FC_up_1<-NC_pval0.05[NC_pval0.05$ value > 0.58,]
NC_pval0.05_FC_down_1<-NC_pval0.05[NC_pval0.05$ value < -0.58,]
NC_significant<- rbind(NC_pval0.05_FC_down_1,NC_pval0.05_FC_up_1)
write.csv(NC_significant, file="NC_specific.csv")


