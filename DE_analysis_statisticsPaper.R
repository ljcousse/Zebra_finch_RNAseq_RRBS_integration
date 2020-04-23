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

list_1dph <- read.table("GeneLists/ListGenes_1dph.txt",header=F,stringsAsFactors = F)$V1
list_20dph <- read.table("GeneLists/ListGenes_20dph.txt",header=F,stringsAsFactors = F)$V1
list_65dph <- read.table("GeneLists/ListGenes_65dph.txt",header=F,stringsAsFactors = F)$V1
list_adult <- read.table("GeneLists/ListGenes_adult.txt",header=F,stringsAsFactors = F)$V1

if (!file.exists("Genes.Rda")){
  ensembl=useMart("ensembl",host="http://jan2019.archive.ensembl.org")
  ensembl = useDataset("tguttata_gene_ensembl",mart=ensembl)
  genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','transcription_start_site'),mart = ensembl)
  print(dim(genes))
  genes <- genes[!duplicated(genes$ensembl_gene_id),]
  print(dim(genes))
  save(genes,file="Genes.Rda")
} else {
  load(file="Genes.Rda")
  genes <- genes[!duplicated(genes$ensembl_gene_id),]
  print(dim(genes))
}
genes <- genes[,1:3]

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



tiff("Clustering.tiff",units="in", width=10, height=10,res=300,family = "sans")
par(xpd=TRUE)
par(mar=c(2.5,2.5,2.5,11))
Bclust<-plot(hclust_complete,cex=2.5,lwd=3,main="",yaxt="n")
# legend(20,0,legend=c(paste(colnames(RRBS_repl),"is",cn)),ncol = 1)
dev.off()


## Amount DE
##########

# Age #

print(sum((fulldata_output$Pvalue.age<0.05) & !is.na(fulldata_output$Pvalue.age)))

# Sex #

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


