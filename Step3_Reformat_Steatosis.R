####################################################################
####Author: Noffisat Oki, Tatyana Doktorova
####Version: 2.0 08-28-2018
####Description: 1) Compound Subset script:
####             Subsets chemical-gene data to filter for only genes active
####             in at least 3 chemicals 
####            2)re-annotation
####            3) union/difference analysis
####Notes: The path to the 'ADVANCE_Project' project folder should be modified to user specified 
####       locations for proper execution of the script. The 'Root_location' object on line 20
####       should be used for this purpose.
####Potential issues:
#####################################################################

library(plyr)
library(data.table)
library(reshape2)

##The root location folder should be specified here. 
#This is the main folder where any subfolders from which data will be read from or written to are located
Root_location <- "C:\\Users\\tatyana\\Documents\\Steatosis"


##Output directory
outDir <- paste(Root_location, "/IntermediateData/", sep='')

################################# Chemical X Gene data_in vivo_high dose ##################
steatosis_chemgene_data.long <- read.csv(paste(Root_location,"/InputFiles/steatosis_specific.csv", sep=''), header=TRUE, sep=',', strip.white=TRUE, fill=TRUE)#, comment.char='', quote='', fill=TRUE, strip.white=TRUE , colClasses="character")
steatosis_chemgene_data.long <- steatosis_chemgene_data.long[steatosis_chemgene_data.long$cellType == 'in vivo',]
steatosis_chemgene_data.long <- steatosis_chemgene_data.long[steatosis_chemgene_data.long$doseLevel == 'High',]
##reformatting the data
steatosis_chemgene_data.wide <- dcast(steatosis_chemgene_data.long, TGG_compoundName ~ assayId , value.var="geneSymbols")
recoding_steatosis_chemgene_data.wide <- (steatosis_chemgene_data.wide)
recoding_steatosis_chemgene_data.wide[,2:ncol(recoding_steatosis_chemgene_data.wide)][recoding_steatosis_chemgene_data.wide[,2:ncol(recoding_steatosis_chemgene_data.wide)]>1] <- 1
rownames(recoding_steatosis_chemgene_data.wide) <- recoding_steatosis_chemgene_data.wide[,1]
recoding_steatosis_chemgene_data.wide <- ((recoding_steatosis_chemgene_data.wide[,-1]))

Active_steatosis_dataset<- recoding_steatosis_chemgene_data.wide
#remove columns which are all '0' indicating no effect for an assay, gene or disease
Active_steatosis_dataset<- Active_steatosis_dataset[,!apply((Active_steatosis_dataset),2,function(x) sum(abs(x), na.rm=TRUE) < 3)]
chemicals <- as.data.frame(rownames(Active_steatosis_dataset))
colnames(chemicals) <- "chemical"
Probes <- as.data.frame(colnames(Active_steatosis_dataset))
Active_steatosis_dataset <- cbind(chemicals, Active_steatosis_dataset)
Active_steatosis_dataset_T<-t(Active_steatosis_dataset)
Active_steatosis_dataset<-cbind(Row.Names=rownames(Active_steatosis_dataset_T),Active_steatosis_dataset_T)

write.table(Active_steatosis_dataset,file=file.path(outDir, "Active_steatosis_dataset.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
write.table(Probes,file=file.path(outDir, "Active_steatosis_probes.csv"), sep=',', col.names=FALSE, row.names=FALSE,quote=FALSE)
   
# re-annotation to ensembl ID                            
library("biomaRt")
ensembl <- useMart("ensembl")
ensembl <- useDataset("rnorvegicus_gene_ensembl",mart=ensembl)

mapping_steatosis <- getBM(attributes = c("affy_rat230_2",
                                    "ensembl_gene_id", "description", "rgd_symbol"), filters = "affy_rat230_2",
                     values = Probes, mart = ensembl)


steatosis_GL_final_High <- mapping_steatosis[!duplicated(mapping_steatosis[,4]),] 



################################# Chemical X Gene data_in vivo_low dose ##################
steatosis_chemgene_data.long_x <- read.csv(paste(Root_location,"/InputFiles/steatosis_specific.csv", sep=''), header=TRUE, sep=',', strip.white=TRUE, fill=TRUE)#, comment.char='', quote='', fill=TRUE, strip.white=TRUE , colClasses="character")
steatosis_chemgene_data.long_x <- steatosis_chemgene_data.long_x[steatosis_chemgene_data.long_x$cellType == 'in vivo',]
steatosis_chemgene_data.long_x <- steatosis_chemgene_data.long_x[steatosis_chemgene_data.long_x$doseLevel == 'Low',]
##reformatting the data
steatosis_chemgene_data.wide <- dcast(steatosis_chemgene_data.long, TGG_compoundName ~ assayId , value.var="geneSymbols")
recoding_steatosis_chemgene_data.wide <- (steatosis_chemgene_data.wide)
recoding_steatosis_chemgene_data.wide[,2:ncol(recoding_steatosis_chemgene_data.wide)][recoding_steatosis_chemgene_data.wide[,2:ncol(recoding_steatosis_chemgene_data.wide)]>1] <- 1
rownames(recoding_steatosis_chemgene_data.wide) <- recoding_steatosis_chemgene_data.wide[,1]
recoding_steatosis_chemgene_data.wide <- ((recoding_steatosis_chemgene_data.wide[,-1]))
#write.table(recoding_GTX_chemgene_data.wide,file=file.path(outDir, "GTX_chem_gene_data_wide.csv"), sep=',', col.names=TRUE, row.names=TRUE,quote=FALSE)

Active_steatosis_dataset<- recoding_steatosis_chemgene_data.wide
#remove columns which are all '0' indicating no effect for an assay, gene or disease
Active_steatosis_dataset<- Active_steatosis_dataset[,!apply((Active_steatosis_dataset),2,function(x) sum(abs(x), na.rm=TRUE) < 3)]
chemicals <- as.data.frame(rownames(Active_steatosis_dataset))
colnames(chemicals) <- "chemical"
Probes <- as.data.frame(colnames(Active_steatosis_dataset))
Active_steatosis_dataset <- cbind(chemicals, Active_steatosis_dataset)
Active_steatosis_dataset_T<-t(Active_steatosis_dataset)
Active_steatosis_dataset<-cbind(Row.Names=rownames(Active_steatosis_dataset_T),Active_steatosis_dataset_T)

write.table(Active_steatosis_dataset,file=file.path(outDir, "Active_steatosis_dataset.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
write.table(Probes,file=file.path(outDir, "Active_steatosis_probes.csv"), sep=',', col.names=FALSE, row.names=FALSE,quote=FALSE)

# re-annotation to ensembl ID  
mapping_steatosis <- getBM(attributes = c("affy_rat230_2",
                                          "ensembl_gene_id", "description", "rgd_symbol"), filters = "affy_rat230_2",
                           values = Probes, mart = ensembl)


steatosis_GL_final_Low <- mapping_steatosis[!duplicated(mapping_steatosis[,4]),] 

################################# non-steatosis Chemical X Gene data_High_dose ##################


nonsteatosis_chemgene_data.long <- read.csv(paste(Root_location,"/InputFiles/nonsteatosis_specific.csv", sep=''), header=TRUE, sep=',', strip.white=TRUE, fill=TRUE)#, comment.char='', quote='', fill=TRUE, strip.white=TRUE , colClasses="character")
nonsteatosis_chemgene_data.long <- nonsteatosis_chemgene_data.long[nonsteatosis_chemgene_data.long$cellType == 'in vivo',]
nonsteatosis_chemgene_data.long <- nonsteatosis_chemgene_data.long[nonsteatosis_chemgene_data.long$doseLevel == 'High',]


##reformatting the data
nonsteatosis_chemgene_data.wide <- dcast(nonsteatosis_chemgene_data.long, TGG_compoundName ~ assayId , value.var="geneSymbols")
recoding_nonsteatosis_chemgene_data.wide <- (nonsteatosis_chemgene_data.wide)
recoding_nonsteatosis_chemgene_data.wide[,2:ncol(recoding_nonsteatosis_chemgene_data.wide)][recoding_nonsteatosis_chemgene_data.wide[,2:ncol(recoding_nonsteatosis_chemgene_data.wide)]>1] <- 1
rownames(recoding_nonsteatosis_chemgene_data.wide) <- recoding_nonsteatosis_chemgene_data.wide[,1]
recoding_nonsteatosis_chemgene_data.wide <- ((recoding_nonsteatosis_chemgene_data.wide[,-1]))

Active_nonsteatosis_dataset<- recoding_nonsteatosis_chemgene_data.wide
#remove columns which are all '0' indicating no effect for an assay, gene or disease
Active_nonsteatosis_dataset<- Active_nonsteatosis_dataset[,!apply((Active_nonsteatosis_dataset),2,function(x) sum(abs(x), na.rm=TRUE) < 3)]
chemicals <- as.data.frame(rownames(Active_nonsteatosis_dataset))
colnames(chemicals) <- "chemical"
nonsteatosis_Probes <- as.data.frame(colnames(Active_nonsteatosis_dataset))
Active_nonsteatosis_dataset <- cbind(chemicals, Active_nonsteatosis_dataset)
Active_nonsteatosis_dataset_T<-t(Active_nonsteatosis_dataset)
Active_nonsteatosis_dataset<-cbind(Row.Names=rownames(Active_nonsteatosis_dataset_T),Active_nonsteatosis_dataset_T)

write.table(Active_nonsteatosis_dataset,file=file.path(outDir, "Active_nonsteatosis_dataset.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
write.table(nonsteatosis_Probes,file=file.path(outDir, "Active_nonsteatosis_probes.csv"), sep=',', col.names=FALSE, row.names=FALSE,quote=FALSE)


# re-annotation to ensembl ID  

mapping_nonsteatosis <- getBM(attributes = c("affy_rat230_2",
                                             "ensembl_gene_id", "description", "rgd_symbol"), filters = "affy_rat230_2",
                              values = nonsteatosis_Probes, mart = ensembl)


nonsteatotic_GL_final_High <- mapping_nonsteatosis[!duplicated(mapping_nonsteatosis[,3]),] 



################################# NC Chemical X Gene data ##################

################################# non-steatosis Chemical X Gene data_Low_dose ##################


nonsteatosis_chemgene_data.long <- read.csv(paste(Root_location,"/InputFiles/nonsteatosis_specific.csv", sep=''), header=TRUE, sep=',', strip.white=TRUE, fill=TRUE)#, comment.char='', quote='', fill=TRUE, strip.white=TRUE , colClasses="character")
nonsteatosis_chemgene_data.long <- nonsteatosis_chemgene_data.long[nonsteatosis_chemgene_data.long$cellType == 'in vivo',]
nonsteatosis_chemgene_data.long <- nonsteatosis_chemgene_data.long[nonsteatosis_chemgene_data.long$doseLevel == 'Low',]


##reformatting the data
nonsteatosis_chemgene_data.wide <- dcast(nonsteatosis_chemgene_data.long, TGG_compoundName ~ assayId , value.var="geneSymbols")
recoding_nonsteatosis_chemgene_data.wide <- (nonsteatosis_chemgene_data.wide)
recoding_nonsteatosis_chemgene_data.wide[,2:ncol(recoding_nonsteatosis_chemgene_data.wide)][recoding_nonsteatosis_chemgene_data.wide[,2:ncol(recoding_nonsteatosis_chemgene_data.wide)]>1] <- 1
rownames(recoding_nonsteatosis_chemgene_data.wide) <- recoding_nonsteatosis_chemgene_data.wide[,1]
recoding_nonsteatosis_chemgene_data.wide <- ((recoding_nonsteatosis_chemgene_data.wide[,-1]))

Active_nonsteatosis_dataset<- recoding_nonsteatosis_chemgene_data.wide
#remove columns which are all '0' indicating no effect for an assay, gene or disease
Active_nonsteatosis_dataset<- Active_nonsteatosis_dataset[,!apply((Active_nonsteatosis_dataset),2,function(x) sum(abs(x), na.rm=TRUE) < 3)]
chemicals <- as.data.frame(rownames(Active_nonsteatosis_dataset))
colnames(chemicals) <- "chemical"
nonsteatosis_Probes <- as.data.frame(colnames(Active_nonsteatosis_dataset))
Active_nonsteatosis_dataset <- cbind(chemicals, Active_nonsteatosis_dataset)
Active_nonsteatosis_dataset_T<-t(Active_nonsteatosis_dataset)
Active_nonsteatosis_dataset<-cbind(Row.Names=rownames(Active_nonsteatosis_dataset_T),Active_nonsteatosis_dataset_T)

write.table(Active_nonsteatosis_dataset,file=file.path(outDir, "Active_nonsteatosis_dataset.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
write.table(nonsteatosis_Probes,file=file.path(outDir, "Active_nonsteatosis_probes.csv"), sep=',', col.names=FALSE, row.names=FALSE,quote=FALSE)


# re-annotation to ensembl ID  

mapping_nonsteatosis <- getBM(attributes = c("affy_rat230_2",
                                             "ensembl_gene_id", "description", "rgd_symbol"), filters = "affy_rat230_2",
                              values = nonsteatosis_Probes, mart = ensembl)


nonsteatosis_GL_final_Low <- mapping_nonsteatosis[!duplicated(mapping_nonsteatosis[,3]),] 


#Extracting genes unique to the steatosis group

steatosis_GL_High<-unique(subset(steatosis_GL_final_High, select=c("rgd_symbol")))
nonsteatosis_GL_High<-unique(subset(nonsteatotic_GL_final_High, select=c("rgd_symbol")))

steatosis_GL_Low<-unique(subset(steatosis_GL_final_Low, select=c("rgd_symbol")))
nonsteatosis_GL_Low<-unique(subset(nonsteatosis_GL_final_Low, select=c("rgd_symbol")))


steatotic_specific_High<-setdiff(steatosis_GL_High, nonsteatosis_GL_High)


steatotic_specific_Low<-setdiff(steatosis_GL_Low, nonsteatosis_GL_Low)

write.table(steatotic_specific_High,file=file.path(outDir, "Steatosis specific High Dose.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
write.table(steatotic_specific_Low,file=file.path(outDir, "Steatosis specific Low Dose.csv"), sep=',', col.names=TRUE, row.names=FALSE,quote=FALSE)
