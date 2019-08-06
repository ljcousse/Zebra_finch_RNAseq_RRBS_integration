###################################
###################################
##                               ##
## DE_analysis_statisticsPaper.R ##
##                               ##
###################################
###################################



##############
# Clustering #
##############


## Load libraries
################

# library("xlsx")
library("UniProt.ws")
library("stringr")
library("edgeR")
library("tidyverse")
library("stageR")
library("biomaRt")


## Load data
##########

fulldata_output <- read.table(file="FullData_output.txt",header=T,sep="\t",stringsAsFactors = F)

rnadata<-read.table("Plastoscine.txt",header=T,sep="\t",quote="")
rna_annot <- read.table("Plastoscine_annot.txt",header=F,sep="\t",quote="")
row.names(rnadata) <- unlist(rna_annot)

f1 <- "1dph,1dph,1dph,1dph,1dph,1dph,20dph,20dph,20dph,20dph,20dph,20dph,65dph,65dph,65dph,65dph,65dph,65dph,adult,adult,adult,adult,adult,adult"
f2 <- "F,F,F,M,M,M,F,F,F,M,M,M,F,F,F,M,M,M,F,F,F,M,M,M"
f1 <- unlist(strsplit(f1,","))
f2 <- unlist(strsplit(f2,","))

y <- DGEList(rnadata)
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep,  ,keep.lib.sizes=FALSE]
summary(keep)


## Clustering
##########

m_matrix <- cpm(y,log=T)
cn <- paste(f1,f2,sep=",")
colnames(m_matrix) <- cn
print(cn)

correl <- 1- abs(cov(m_matrix))
distance <- as.dist(correl)
hclust_complete <- hclust(distance, method = "complete")
dend_complete <- as.dendrogram(hclust_complete)

jpeg("Clustering.jpg")
par(xpd=TRUE)
par(mar=c(2.5,2.5,2.5,11))
Bclust<-plot(hclust_complete)
# legend(20,0,legend=c(paste(colnames(RRBS_repl),"is",cn)),ncol = 1)
dev.off()


##############
# Clustering #
##############


# Age
##########

print(sum((fulldata_output$Pvalue.age<0.05) & !is.na(fulldata_output$Pvalue.age<0.05)))


# Sex
############

# All timepoints
print(sum((!is.na(fulldata_output$Pvalue.1dph) & (fulldata_output$Pvalue.1dph<0.05)) | 
            (!is.na(fulldata_output$Pvalue.20dph) & (fulldata_output$Pvalue.20dph<0.05)) |
            (!is.na(fulldata_output$Pvalue.65dph) & (fulldata_output$Pvalue.65dph<0.05)) |
            (!is.na(fulldata_output$Pvalue.adult) & (fulldata_output$Pvalue.adult<0.05))))
genes_all <- genes[(genes$ensembl_gene_id%in%list_1dph) |
                     (genes$ensembl_gene_id%in%list_20dph) |
                     (genes$ensembl_gene_id%in%list_65dph) |
                     (genes$ensembl_gene_id%in%list_adult),]
dim(genes_all)
barplot(table(genes_all$chromosome_name),las=2)
sum(grepl("Z",genes_all$chromosome_name))/length(genes_all$chromosome_name)


fulldata_output_allsex <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_all$ensembl_gene_id,]


# Joint between all timepoints
print(sum((!is.na(fulldata_output$Pvalue.1dph) & (fulldata_output$Pvalue.1dph<0.05)) &
            (!is.na(fulldata_output$Pvalue.20dph) & (fulldata_output$Pvalue.20dph<0.05)) &
            (!is.na(fulldata_output$Pvalue.65dph) & (fulldata_output$Pvalue.65dph<0.05)) &
            (!is.na(fulldata_output$Pvalue.adult) & (fulldata_output$Pvalue.adult<0.05))))
genes_joint <- genes[(genes$ensembl_gene_id%in%list_1dph) &
                       (genes$ensembl_gene_id%in%list_20dph) &
                       (genes$ensembl_gene_id%in%list_65dph) &
                       (genes$ensembl_gene_id%in%list_adult),]
genes_joint <- genes_joint[!duplicated(genes_joint$ensembl_gene_id),]
dim(genes_joint)
barplot(table(genes_joint$chromosome_name),las=2)
sum(grepl("Z",genes_joint$chromosome_name))/length(genes_joint$chromosome_name)

fulldata_output_joint <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_joint$ensembl_gene_id,]
dim(fulldata_output_joint)
sum(rowSums(abs(sign(fulldata_output_joint[,grepl("logFC",colnames(fulldata_output_joint)) & !grepl("age",colnames(fulldata_output_joint))])))==4)

# Time point specific genes
genes_1dph <- genes[genes$ensembl_gene_id%in%list_1dph[!(list_1dph%in%c(list_adult,list_20dph,list_65dph))],]
genes_1dph <- genes_1dph[!duplicated(genes_1dph$ensembl_gene_id),]
dim(genes_1dph)
barplot(table(genes_1dph$chromosome_name),las=2)
sum(grepl("Z",genes_1dph$chromosome_name))/length(genes_1dph$chromosome_name)
fulldata_output_1dph <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_1dph$ensembl_gene_id,]
dim(fulldata_output_1dph)

genes_20dph <- genes[genes$ensembl_gene_id%in%list_20dph[!(list_20dph%in%c(list_1dph,list_adult,list_65dph))],]
genes_20dph <- genes_20dph[!duplicated(genes_20dph$ensembl_gene_id),]
dim(genes_20dph)
barplot(table(genes_20dph$chromosome_name),las=2)
sum(grepl("Z",genes_20dph$chromosome_name))/length(genes_20dph$chromosome_name)
fulldata_output_20dph <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_20dph$ensembl_gene_id,]
dim(fulldata_output_20dph)

genes_65dph <- genes[genes$ensembl_gene_id%in%list_65dph[!(list_65dph%in%c(list_1dph,list_20dph,list_adult))],]
genes_65dph <- genes_65dph[!duplicated(genes_65dph$ensembl_gene_id),]
dim(genes_65dph)
barplot(table(genes_65dph$chromosome_name),las=2)
sum(grepl("Z",genes_65dph$chromosome_name))/length(genes_65dph$chromosome_name)
fulldata_output_65dph <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_65dph$ensembl_gene_id,]
dim(fulldata_output_65dph)

genes_adult <- genes[genes$ensembl_gene_id%in%list_adult[!(list_adult%in%c(list_1dph,list_20dph,list_65dph))],]
genes_adult <- genes_adult[!duplicated(genes_adult$ensembl_gene_id),]
dim(genes_adult)
barplot(table(genes_adult$chromosome_name),las=2)
sum(grepl("Z",genes_adult$chromosome_name))/length(genes_adult$chromosome_name)
fulldata_output_adult <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_adult$ensembl_gene_id,]
dim(fulldata_output_adult)

head(fulldata_output[fulldata_output$ensembl.gene.id%in%genes_1dph$ensembl_gene_id,])
head(fulldata_output[fulldata_output$ensembl.gene.id%in%genes_20dph$ensembl_gene_id,])
head(fulldata_output[fulldata_output$ensembl.gene.id%in%genes_65dph$ensembl_gene_id,])
head(fulldata_output[fulldata_output$ensembl.gene.id%in%genes_adult$ensembl_gene_id,])


print(sum((fulldata_output$Pvalue.age<0.05) & !is.na(fulldata_output$Pvalue.age<0.05)))


## TF

annot <- read.table("MATRIX.txt",stringsAsFactors = F,fill=T)
print(head(annot))
species <- read.table("MATRIX_SPECIES.txt",stringsAsFactors = F,fill=T)
print(head(species))
species <- species[!duplicated(species$V1),]
annot2 <- read.table("MATRIX_ANNOTATION.txt",stringsAsFactors = F,sep="\t",quote='')
print(head(annot2))
print(dim(annot2))
print(str(annot2))
class <- annot2[grepl("class",annot2$V2),c(1,3)]
print(head(class))
class <- class[!duplicated(class$V1),]
family <- annot2[grepl("family",annot2$V2),c(1,3)]
print(head(family))
family <- family[!duplicated(family$V1),]
medline <- annot2[grepl("medline",annot2$V2),c(1,3)]
print(head(medline))
medline <- medline[!duplicated(medline$V1),]
print(sum(duplicated(annot$V1)))
annot_new <- as.data.frame(cbind(rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot))),stringsAsFactors=F)
colnames(annot_new) <- c("ID","M_ID","Gene_Symbol","Species","Family","Class","Medline")
for (i in 1:nrow(annot)){
  annot_new$ID[i] <- annot[i,1]
  annot_new$M_ID[i] <- annot[i,3]
  annot_new$Gene_Symbol[i] <- annot[i,5]
  if (!identical(species$V2[species$V1==annot$V1[i]],character(0))){
    annot_new$Species[i] <- species$V2[species$V1==annot$V1[i]]
  }
  if (!identical(family$V3[family$V1==annot$V1[i]],character(0))){
    annot_new$Family[i] <- family$V3[family$V1==annot$V1[i]]
  }
  if (!identical(class$V3[class$V1==annot$V1[i]],character(0))){
    annot_new$Class[i] <- class$V3[class$V1==annot$V1[i]]
  }
  if (!identical(medline$V3[medline$V1==annot$V1[i]],character(0))){
    annot_new$Medline[i] <-medline$V3[medline$V1==annot$V1[i]]
  }
}
annot <- annot_new

annot_human <- annot[annot$Species==9606,]

fulldata_output_jaspar <- fulldata_output[(fulldata_output$gene.symbol%in%annot_human$Gene_Symbol) | (fulldata_output$uniprot.annot%in%annot_human$Gene_Symbol),]
write.xlsx(fulldata_output_jaspar,file="Table_S9_RNAseq_output_TF.xlsx",row.names = F)


# check TF analysis

JASPAR_results_1dph <- read.table("JASPAR_TF_1dph_comp_h_filtered.txt",sep="\t",header=T,stringsAsFactors = F)
JASPAR_results_1dph_sign <- JASPAR_results_1dph[JASPAR_results_1dph$FDR<0.05,]
JASPAR_results_20dph <- read.table("JASPAR_TF_20dph_comp_h_filtered.txt",sep="\t",header=T,stringsAsFactors = F)
JASPAR_results_20dph_sign <- JASPAR_results_20dph[JASPAR_results_20dph$FDR<0.05,]
JASPAR_results_65dph <- read.table("JASPAR_TF_65dph_comp_h_filtered.txt",sep="\t",header=T,stringsAsFactors = F)
JASPAR_results_65dph_sign <- JASPAR_results_65dph[JASPAR_results_65dph$FDR<0.05,]
JASPAR_results_adult <- read.table("JASPAR_TF_adult_comp_h_filtered.txt",sep="\t",header=T,stringsAsFactors = F)
JASPAR_results_adult_sign <- JASPAR_results_adult[JASPAR_results_adult$FDR<0.05,]
JASPAR_results_age <- read.table("JASPAR_TF_age_comp_h_filtered.txt",sep="\t",header=T,stringsAsFactors = F)
JASPAR_results_age_sign <- JASPAR_results_age[JASPAR_results_age$FDR<0.05,]

TF_all_sex <- unique(c(JASPAR_results_1dph_sign$TF,JASPAR_results_20dph_sign$TF,JASPAR_results_65dph_sign$TF,JASPAR_results_adult_sign$TF))
TF_all_sex[!(TF_all_sex%in%JASPAR_results_age_sign$TF)]

# Only DE on 20dph and 65dph
selection <- (out2excel$Pvalue.20dph<0.05) & (out2excel$Pvalue.65dph<0.05) & 
  (out2excel$Pvalue.1dph>0.05) & (out2excel$Pvalue.adult>0.05) & 
  (!is.na(out2excel$Pvalue.20dph)) & (!is.na(out2excel$Pvalue.65dph)) & 
  (!is.na(out2excel$Pvalue.1dph)) & (!is.na(out2excel$Pvalue.adult)) 
dim(out2excel[selection,c(3,grep("Pvalue",colnames(out2excel)))])
out2excel[selection,c(3,grep("Pvalue",colnames(out2excel)))]

# Only DE on 20dph
selection <- (out2excel$Pvalue.20dph<0.05) & (out2excel$Pvalue.65dph>0.05) & 
  (out2excel$Pvalue.1dph>0.05) & (out2excel$Pvalue.adult>0.05) & 
  (!is.na(out2excel$Pvalue.20dph)) & (!is.na(out2excel$Pvalue.65dph)) & 
  (!is.na(out2excel$Pvalue.1dph)) & (!is.na(out2excel$Pvalue.adult)) 
dim(out2excel[selection,c(3,grep("Pvalue",colnames(out2excel)))])
out2excel[selection,c(3,grep("Pvalue",colnames(out2excel)))]

# Only DE on 65dph
selection <- (out2excel$Pvalue.20dph>0.05) & (out2excel$Pvalue.65dph<0.05) & 
  (out2excel$Pvalue.1dph>0.05) & (out2excel$Pvalue.adult>0.05) & 
  (!is.na(out2excel$Pvalue.20dph)) & (!is.na(out2excel$Pvalue.65dph)) & 
  (!is.na(out2excel$Pvalue.1dph)) & (!is.na(out2excel$Pvalue.adult)) 
dim(out2excel[selection,c(3,grep("Pvalue",colnames(out2excel)))])
out2excel[selection,c(3,grep("Pvalue",colnames(out2excel)))]

# Only DE on adult
selection <- (out2excel$Pvalue.20dph>0.05) & (out2excel$Pvalue.65dph<0.05) & 
  (out2excel$Pvalue.1dph>0.05) & (out2excel$Pvalue.adult<0.05) & 
  (!is.na(out2excel$Pvalue.20dph)) & (!is.na(out2excel$Pvalue.65dph)) & 
  (!is.na(out2excel$Pvalue.1dph)) & (!is.na(out2excel$Pvalue.adult)) 
dim(out2excel[selection,c(3,grep("Pvalue",colnames(out2excel)))])
out2excel[selection,c(3,grep("Pvalue",colnames(out2excel)))]


chisq.test(cbind(c(sum(genes_adult$chromosome_name=="Z"),length(genes_adult$chromosome_name)),
                 c(sum(genes_65dph$chromosome_name=="Z"),length(genes_65dph$chromosome_name))))
chisq.test(cbind(c(sum(genes_adult$chromosome_name=="Z"),length(genes_adult$chromosome_name)),
                 c(sum(genes_20dph$chromosome_name=="Z"),length(genes_20dph$chromosome_name))))
chisq.test(cbind(c(sum(genes_adult$chromosome_name=="Z"),length(genes_adult$chromosome_name)),
                 c(sum(genes_1dph$chromosome_name=="Z"),length(genes_1dph$chromosome_name))))


# # Heatmap
# ###################
# 
# rnadata<-read.table("Plastoscine.txt",header=T,sep="\t",quote="")
# rna_annot <- read.table("Plastoscine_annot.txt",header=F,sep="\t",quote="")
# row.names(rnadata) <- unlist(rna_annot)
# 
# f1 <- "1dph,1dph,1dph,1dph,1dph,1dph,20dph,20dph,20dph,20dph,20dph,20dph,65dph,65dph,65dph,65dph,65dph,65dph,adult,adult,adult,adult,adult,adult"
# f2 <- "F,F,F,M,M,M,F,F,F,M,M,M,F,F,F,M,M,M,F,F,F,M,M,M"
# f1 <- unlist(strsplit(f1,","))
# f2 <- unlist(strsplit(f2,","))
# 
# colnames(rnadata) <- paste(f1,f2,sep=",")
# 
# rnadata <- rnadata[rowMeans(rnadata)>5,]
# 
select = order(apply(y$counts,1,var),decreasing=TRUE)[1:10000]
# 
# y <- DGEList(rnadata)
# y <- calcNormFactors(y)
m_matrix <- cpm(y,log=T)
# 
# 
library("RColorBrewer")
library("gplots")

hmcol = colorRampPalette(brewer.pal(9, "RdBu"))(100)
cor_t <- 1- abs(cor(t(m_matrix[select,])))
cor <- 1- abs(cor(m_matrix[select,]))
distancet <- as.dist(cor_t)
distance <- as.dist(cor)
hclust_complete <- hclust(distancet, method = "complete")
dend_tcomplete <- as.dendrogram(hclust_complete)
hclus_complete <- hclust(distance, method = "complete")
dend_complete <- as.dendrogram(hclus_complete)
jpeg("Heatmap.jpg")
heatmap.2(m_matrix[select,], col = hmcol,Rowv=dend_tcomplete, Colv=dend_complete, trace="non", margin=c(5,1),labRow = "")
dev.off()


# Time point genes
genes_1dph_all <- genes[genes$ensembl_gene_id%in%list_1dph,]
genes_1dph_all <- genes_1dph_all[!duplicated(genes_1dph_all$ensembl_gene_id),]
dim(genes_1dph_all)
barplot(table(genes_1dph_all$chromosome_name),las=2)
sum(grepl("Z",genes_1dph_all$chromosome_name))/length(genes_1dph_all$chromosome_name)
fulldata_output_1dph_all <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_1dph_all$ensembl_gene_id,]
dim(fulldata_output_1dph_all)

genes_20dph_all <- genes[genes$ensembl_gene_id%in%list_20dph,]
genes_20dph_all <- genes_20dph_all[!duplicated(genes_20dph_all$ensembl_gene_id),]
dim(genes_20dph_all)
barplot(table(genes_20dph_all$chromosome_name),las=2)
sum(grepl("Z",genes_20dph_all$chromosome_name))/length(genes_20dph_all$chromosome_name)
fulldata_output_20dph_all <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_20dph_all$ensembl_gene_id,]
dim(fulldata_output_20dph_all)

genes_65dph_all <- genes[genes$ensembl_gene_id%in%list_65dph,]
genes_65dph_all <- genes_65dph_all[!duplicated(genes_65dph_all$ensembl_gene_id),]
dim(genes_65dph_all)
barplot(table(genes_65dph_all$chromosome_name),las=2)
sum(grepl("Z",genes_65dph_all$chromosome_name))/length(genes_65dph_all$chromosome_name)
fulldata_output_65dph_all <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_65dph_all$ensembl_gene_id,]
dim(fulldata_output_65dph_all)

genes_adult_all <- genes[genes$ensembl_gene_id%in%list_adult,]
genes_adult_all <- genes_adult_all[!duplicated(genes_adult_all$ensembl_gene_id),]
dim(genes_adult_all)
barplot(table(genes_adult_all$chromosome_name),las=2)
sum(grepl("Z",genes_adult_all$chromosome_name))/length(genes_adult_all$chromosome_name)
fulldata_output_adult_all <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_adult_all$ensembl_gene_id,]
dim(fulldata_output_adult_all)

pvals <- NULL
pvals <- c(pvals,chisq.test(cbind(c(sum(grepl("Z",genes_20dph_all$chromosome_name)),sum(!grepl("Z",genes_20dph_all$chromosome_name))),
                                  c(sum(grepl("Z",genes_1dph_all$chromosome_name)),sum(!grepl("Z",genes_1dph_all$chromosome_name)))))$p.value)
pvals <- c(pvals,chisq.test(cbind(c(sum(grepl("Z",genes_65dph_all$chromosome_name)),sum(!grepl("Z",genes_65dph_all$chromosome_name))),
                                  c(sum(grepl("Z",genes_1dph_all$chromosome_name)),sum(!grepl("Z",genes_1dph_all$chromosome_name)))))$p.value)
pvals <- c(pvals,chisq.test(cbind(c(sum(grepl("Z",genes_adult_all$chromosome_name)),sum(!grepl("Z",genes_adult_all$chromosome_name))),
                                  c(sum(grepl("Z",genes_1dph_all$chromosome_name)),sum(!grepl("Z",genes_1dph_all$chromosome_name)))))$p.value)
pvals <- c(pvals,chisq.test(cbind(c(sum(grepl("Z",genes_65dph_all$chromosome_name)),sum(!grepl("Z",genes_65dph_all$chromosome_name))),
                                  c(sum(grepl("Z",genes_20dph_all$chromosome_name)),sum(!grepl("Z",genes_20dph_all$chromosome_name)))))$p.value)
pvals <- c(pvals,chisq.test(cbind(c(sum(grepl("Z",genes_adult_all$chromosome_name)),sum(!grepl("Z",genes_adult_all$chromosome_name))),
                                  c(sum(grepl("Z",genes_20dph_all$chromosome_name)),sum(!grepl("Z",genes_20dph_all$chromosome_name)))))$p.value)
pvals <- c(pvals,chisq.test(cbind(c(sum(grepl("Z",genes_adult_all$chromosome_name)),sum(!grepl("Z",genes_adult_all$chromosome_name))),
                                  c(sum(grepl("Z",genes_65dph_all$chromosome_name)),sum(!grepl("Z",genes_65dph_all$chromosome_name)))))$p.value)


sum(fulldata_output_1dph_all$logFC.1dph>0)

sum(fulldata_output_allsex$logFC.adult[fulldata_output_allsex$ensembl.gene.id%in%genes_all$ensembl_gene_id[grepl("Z",genes_all$chromosome_name)]]>0)/length(fulldata_output_allsex$logFC.adult[fulldata_output_allsex$ensembl.gene.id%in%genes_all$ensembl_gene_id[grepl("Z",genes_all$chromosome_name)]])
sum(fulldata_output_joint$logFC.1dph[fulldata_output_joint$ensembl.gene.id%in%genes_joint$ensembl_gene_id[!grepl("Z",genes_joint$chromosome_name)]]>0)/nrow(fulldata_output_joint)


mat <- NULL
mat <- rbind(mat,c(sum(fulldata_output_1dph_all$logFC.1dph[fulldata_output_1dph_all$ensembl.gene.id%in%genes_1dph_all$ensembl_gene_id[grepl("Z",genes_1dph_all$chromosome_name)]]>0),nrow(fulldata_output_1dph_all[fulldata_output_1dph_all$ensembl.gene.id%in%genes_1dph_all$ensembl_gene_id[grepl("Z",genes_1dph_all$chromosome_name)],])))
mat <- rbind(mat,c(sum(fulldata_output_20dph_all$logFC.20dph[fulldata_output_20dph_all$ensembl.gene.id%in%genes_20dph_all$ensembl_gene_id[grepl("Z",genes_20dph_all$chromosome_name)]]>0),nrow(fulldata_output_20dph_all[fulldata_output_20dph_all$ensembl.gene.id%in%genes_20dph_all$ensembl_gene_id[grepl("Z",genes_20dph_all$chromosome_name)],])))
mat <- rbind(mat,c(sum(fulldata_output_65dph_all$logFC.65dph[fulldata_output_65dph_all$ensembl.gene.id%in%genes_65dph_all$ensembl_gene_id[grepl("Z",genes_65dph_all$chromosome_name)]]>0),nrow(fulldata_output_65dph_all[fulldata_output_65dph_all$ensembl.gene.id%in%genes_65dph_all$ensembl_gene_id[grepl("Z",genes_65dph_all$chromosome_name)],])))
mat <- rbind(mat,c(sum(fulldata_output_adult_all$logFC.adult[fulldata_output_adult_all$ensembl.gene.id%in%genes_adult_all$ensembl_gene_id[grepl("Z",genes_adult_all$chromosome_name)]]>0),nrow(fulldata_output_adult_all[fulldata_output_adult_all$ensembl.gene.id%in%genes_adult_all$ensembl_gene_id[grepl("Z",genes_adult_all$chromosome_name)],])))


mat <- NULL
mat <- rbind(mat,c(sum(fulldata_output_1dph$logFC.1dph[fulldata_output_1dph$ensembl.gene.id%in%genes_1dph$ensembl_gene_id[grepl("Z",genes_1dph$chromosome_name)]]>0),nrow(fulldata_output_1dph[fulldata_output_1dph$ensembl.gene.id%in%genes_1dph$ensembl_gene_id[grepl("Z",genes_1dph$chromosome_name)],])))
mat <- rbind(mat,c(sum(fulldata_output_20dph$logFC.20dph[fulldata_output_20dph$ensembl.gene.id%in%genes_20dph$ensembl_gene_id[grepl("Z",genes_20dph$chromosome_name)]]>0),nrow(fulldata_output_20dph[fulldata_output_20dph$ensembl.gene.id%in%genes_20dph$ensembl_gene_id[grepl("Z",genes_20dph$chromosome_name)],])))
mat <- rbind(mat,c(sum(fulldata_output_65dph$logFC.65dph[fulldata_output_65dph$ensembl.gene.id%in%genes_65dph$ensembl_gene_id[grepl("Z",genes_65dph$chromosome_name)]]>0),nrow(fulldata_output_65dph[fulldata_output_65dph$ensembl.gene.id%in%genes_65dph$ensembl_gene_id[grepl("Z",genes_65dph$chromosome_name)],])))
mat <- rbind(mat,c(sum(fulldata_output_adult$logFC.adult[fulldata_output_adult$ensembl.gene.id%in%genes_adult$ensembl_gene_id[grepl("Z",genes_adult$chromosome_name)]]>0),nrow(fulldata_output_adult[fulldata_output_adult$ensembl.gene.id%in%genes_adult$ensembl_gene_id[grepl("Z",genes_adult$chromosome_name)],])))

mat <- NULL
mat <- rbind(mat,c(sum(fulldata_output_1dph_all$logFC.1dph[fulldata_output_1dph_all$ensembl.gene.id%in%genes_1dph_all$ensembl_gene_id[!grepl("Z",genes_1dph_all$chromosome_name)]]>0),nrow(fulldata_output_1dph_all[fulldata_output_1dph_all$ensembl.gene.id%in%genes_1dph_all$ensembl_gene_id[!grepl("Z",genes_1dph_all$chromosome_name)],])))
mat <- rbind(mat,c(sum(fulldata_output_20dph_all$logFC.20dph[fulldata_output_20dph_all$ensembl.gene.id%in%genes_20dph_all$ensembl_gene_id[!grepl("Z",genes_20dph_all$chromosome_name)]]>0),nrow(fulldata_output_20dph_all[fulldata_output_20dph_all$ensembl.gene.id%in%genes_20dph_all$ensembl_gene_id[!grepl("Z",genes_20dph_all$chromosome_name)],])))
mat <- rbind(mat,c(sum(fulldata_output_65dph_all$logFC.65dph[fulldata_output_65dph_all$ensembl.gene.id%in%genes_65dph_all$ensembl_gene_id[!grepl("Z",genes_65dph_all$chromosome_name)]]>0),nrow(fulldata_output_65dph_all[fulldata_output_65dph_all$ensembl.gene.id%in%genes_65dph_all$ensembl_gene_id[!grepl("Z",genes_65dph_all$chromosome_name)],])))
mat <- rbind(mat,c(sum(fulldata_output_adult_all$logFC.adult[fulldata_output_adult_all$ensembl.gene.id%in%genes_adult_all$ensembl_gene_id[!grepl("Z",genes_adult_all$chromosome_name)]]>0),nrow(fulldata_output_adult_all[fulldata_output_adult_all$ensembl.gene.id%in%genes_adult_all$ensembl_gene_id[!grepl("Z",genes_adult_all$chromosome_name)],])))
