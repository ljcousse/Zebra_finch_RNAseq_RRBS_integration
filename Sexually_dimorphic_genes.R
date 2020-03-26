################################
################################
##                            ##
## Sexually_dimorphic_genes.R ##
##                            ##
################################
################################


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

# Load in ensembl annotation
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

## Calculate percentage of Z-linked genes
############## 


# All sex genes #

print(sum((!is.na(fulldata_output$Pvalue.1dph) & (fulldata_output$Pvalue.1dph<0.05)) |
            (!is.na(fulldata_output$Pvalue.20dph) & (fulldata_output$Pvalue.20dph<0.05)) |
            (!is.na(fulldata_output$Pvalue.65dph) & (fulldata_output$Pvalue.65dph<0.05)) |
            (!is.na(fulldata_output$Pvalue.adult) & (fulldata_output$Pvalue.adult<0.05))))
genes_allsex <- genes[(genes$ensembl_gene_id%in%list_1dph) |
                     (genes$ensembl_gene_id%in%list_20dph) |
                     (genes$ensembl_gene_id%in%list_65dph) |
                     (genes$ensembl_gene_id%in%list_adult),]
genes_allsex <- genes_allsex[!duplicated(genes_allsex$ensembl_gene_id),]
dim(genes_allsex)
barplot(table(genes_allsex$chromosome_name),las=2)
sum(grepl("Z",genes_allsex$chromosome_name))/length(genes_allsex$chromosome_name)

fulldata_output_allsex <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_allsex$ensembl_gene_id,]
dim(fulldata_output_allsex)
sum(rowSums(abs(sign(fulldata_output_allsex[,grepl("logFC",colnames(fulldata_output_allsex)) & !grepl("age",colnames(fulldata_output_allsex))])))==4)


# Joint between all timepoints #

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


# Time point specific genes #

# 1dph
genes_1dph <- genes[genes$ensembl_gene_id%in%list_1dph[!(list_1dph%in%c(list_adult,list_20dph,list_65dph))],]
genes_1dph <- genes_1dph[!duplicated(genes_1dph$ensembl_gene_id),]
dim(genes_1dph)
barplot(table(genes_1dph$chromosome_name),las=2)
sum(grepl("Z",genes_1dph$chromosome_name))/length(genes_1dph$chromosome_name)
fulldata_output_1dph <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_1dph$ensembl_gene_id,]
dim(fulldata_output_1dph)

# 20dph
genes_20dph <- genes[genes$ensembl_gene_id%in%list_20dph[!(list_20dph%in%c(list_1dph,list_adult,list_65dph))],]
genes_20dph <- genes_20dph[!duplicated(genes_20dph$ensembl_gene_id),]
dim(genes_20dph)
barplot(table(genes_20dph$chromosome_name),las=2)
sum(grepl("Z",genes_20dph$chromosome_name))/length(genes_20dph$chromosome_name)
fulldata_output_20dph <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_20dph$ensembl_gene_id,]
dim(fulldata_output_20dph)

# 65dph
genes_65dph <- genes[genes$ensembl_gene_id%in%list_65dph[!(list_65dph%in%c(list_1dph,list_20dph,list_adult))],]
genes_65dph <- genes_65dph[!duplicated(genes_65dph$ensembl_gene_id),]
dim(genes_65dph)
barplot(table(genes_65dph$chromosome_name),las=2)
sum(grepl("Z",genes_65dph$chromosome_name))/length(genes_65dph$chromosome_name)
fulldata_output_65dph <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_65dph$ensembl_gene_id,]
dim(fulldata_output_65dph)

# adult
genes_adult <- genes[genes$ensembl_gene_id%in%list_adult[!(list_adult%in%c(list_1dph,list_20dph,list_65dph))],]
genes_adult <- genes_adult[!duplicated(genes_adult$ensembl_gene_id),]
dim(genes_adult)
barplot(table(genes_adult$chromosome_name),las=2)
sum(grepl("Z",genes_adult$chromosome_name))/length(genes_adult$chromosome_name)
fulldata_output_adult <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_adult$ensembl_gene_id,]
dim(fulldata_output_adult)


# Time point genes #

# 1dph
genes_1dph_all <- genes[genes$ensembl_gene_id%in%list_1dph,]
genes_1dph_all <- genes_1dph_all[!duplicated(genes_1dph_all$ensembl_gene_id),]
dim(genes_1dph_all)
barplot(table(genes_1dph_all$chromosome_name),las=2)
sum(grepl("Z",genes_1dph_all$chromosome_name))/length(genes_1dph_all$chromosome_name)
fulldata_output_1dph_all <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_1dph_all$ensembl_gene_id,]
dim(fulldata_output_1dph_all)

# 20dph
genes_20dph_all <- genes[genes$ensembl_gene_id%in%list_20dph,]
genes_20dph_all <- genes_20dph_all[!duplicated(genes_20dph_all$ensembl_gene_id),]
dim(genes_20dph_all)
barplot(table(genes_20dph_all$chromosome_name),las=2)
sum(grepl("Z",genes_20dph_all$chromosome_name))/length(genes_20dph_all$chromosome_name)
fulldata_output_20dph_all <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_20dph_all$ensembl_gene_id,]
dim(fulldata_output_20dph_all)

# 65dph
genes_65dph_all <- genes[genes$ensembl_gene_id%in%list_65dph,]
genes_65dph_all <- genes_65dph_all[!duplicated(genes_65dph_all$ensembl_gene_id),]
dim(genes_65dph_all)
barplot(table(genes_65dph_all$chromosome_name),las=2)
sum(grepl("Z",genes_65dph_all$chromosome_name))/length(genes_65dph_all$chromosome_name)
fulldata_output_65dph_all <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_65dph_all$ensembl_gene_id,]
dim(fulldata_output_65dph_all)

# adult
genes_adult_all <- genes[genes$ensembl_gene_id%in%list_adult,]
genes_adult_all <- genes_adult_all[!duplicated(genes_adult_all$ensembl_gene_id),]
dim(genes_adult_all)
barplot(table(genes_adult_all$chromosome_name),las=2)
sum(grepl("Z",genes_adult_all$chromosome_name))/length(genes_adult_all$chromosome_name)
fulldata_output_adult_all <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_adult_all$ensembl_gene_id,]
dim(fulldata_output_adult_all)


# Arising after 1dph #

n_Zafter1dph <- sum(grepl("Z",genes_allsex$chromosome_name[!(genes_allsex$ensembl_gene_id%in%genes_1dph_all$ensembl_gene_id)]))
n_Aafter1dph <- sum(!grepl("Z",genes_allsex$chromosome_name[!(genes_allsex$ensembl_gene_id%in%genes_1dph_all$ensembl_gene_id)]))

n_Aafter1dph/(n_Zafter1dph+n_Aafter1dph)


# Net result #

# 1dph 
sum(grepl("Z",genes_1dph_all$chromosome_name))/(sum(grepl("Z",genes_1dph_all$chromosome_name))+sum(!grepl("Z",genes_1dph_all$chromosome_name)))

# 20dph 
sum(grepl("Z",genes_20dph_all$chromosome_name))/(sum(grepl("Z",genes_20dph_all$chromosome_name))+sum(!grepl("Z",genes_20dph_all$chromosome_name)))

# 65dph 
sum(grepl("Z",genes_65dph_all$chromosome_name))/(sum(grepl("Z",genes_65dph_all$chromosome_name))+sum(!grepl("Z",genes_65dph_all$chromosome_name)))

# adult 
sum(grepl("Z",genes_adult_all$chromosome_name))/(sum(grepl("Z",genes_adult_all$chromosome_name))+sum(!grepl("Z",genes_adult_all$chromosome_name)))

# significance
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

adj.pvals <- pvals*6
adj.pvals[adj.pvals>1] <- 1
adj.pvals

# % Z-linked chi sq test #

mat_Z <- NULL
mat_Z <- rbind(mat_Z,c(sum(fulldata_output_1dph_all$logFC.1dph[fulldata_output_1dph_all$ensembl.gene.id%in%genes_1dph_all$ensembl_gene_id[grepl("Z",genes_1dph_all$chromosome_name)]]>0),
                       sum(fulldata_output_1dph_all$logFC.1dph[fulldata_output_1dph_all$ensembl.gene.id%in%genes_1dph_all$ensembl_gene_id[grepl("Z",genes_1dph_all$chromosome_name)]]<0)))
mat_Z <- rbind(mat_Z,c(sum(fulldata_output_20dph_all$logFC.20dph[fulldata_output_20dph_all$ensembl.gene.id%in%genes_20dph_all$ensembl_gene_id[grepl("Z",genes_20dph_all$chromosome_name)]]>0),
                       sum(fulldata_output_20dph_all$logFC.20dph[fulldata_output_20dph_all$ensembl.gene.id%in%genes_20dph_all$ensembl_gene_id[grepl("Z",genes_20dph_all$chromosome_name)]]<0)))
mat_Z <- rbind(mat_Z,c(sum(fulldata_output_65dph_all$logFC.65dph[fulldata_output_65dph_all$ensembl.gene.id%in%genes_65dph_all$ensembl_gene_id[grepl("Z",genes_65dph_all$chromosome_name)]]>0),
                       sum(fulldata_output_65dph_all$logFC.65dph[fulldata_output_65dph_all$ensembl.gene.id%in%genes_65dph_all$ensembl_gene_id[grepl("Z",genes_65dph_all$chromosome_name)]]<0)))
mat_Z <- rbind(mat_Z,c(sum(fulldata_output_adult_all$logFC.adult[fulldata_output_adult_all$ensembl.gene.id%in%genes_adult_all$ensembl_gene_id[grepl("Z",genes_adult_all$chromosome_name)]]>0),
                       sum(fulldata_output_adult_all$logFC.adult[fulldata_output_adult_all$ensembl.gene.id%in%genes_adult_all$ensembl_gene_id[grepl("Z",genes_adult_all$chromosome_name)]]<0)))

sum(mat_Z[,1])/sum(mat_Z)
chisq.test(mat_Z)


mat_A <- NULL
mat_A <- rbind(mat_A,c(sum(fulldata_output_1dph_all$logFC.1dph[fulldata_output_1dph_all$ensembl.gene.id%in%genes_1dph_all$ensembl_gene_id[!grepl("Z",genes_1dph_all$chromosome_name)]]>0),
                       sum(fulldata_output_1dph_all$logFC.1dph[fulldata_output_1dph_all$ensembl.gene.id%in%genes_1dph_all$ensembl_gene_id[!grepl("Z",genes_1dph_all$chromosome_name)]]<0)))
mat_A <- rbind(mat_A,c(sum(fulldata_output_20dph_all$logFC.20dph[fulldata_output_20dph_all$ensembl.gene.id%in%genes_20dph_all$ensembl_gene_id[!grepl("Z",genes_20dph_all$chromosome_name)]]>0),
                       sum(fulldata_output_20dph_all$logFC.20dph[fulldata_output_20dph_all$ensembl.gene.id%in%genes_20dph_all$ensembl_gene_id[!grepl("Z",genes_20dph_all$chromosome_name)]]<0)))
mat_A <- rbind(mat_A,c(sum(fulldata_output_65dph_all$logFC.65dph[fulldata_output_65dph_all$ensembl.gene.id%in%genes_65dph_all$ensembl_gene_id[!grepl("Z",genes_65dph_all$chromosome_name)]]>0),
                       sum(fulldata_output_65dph_all$logFC.65dph[fulldata_output_65dph_all$ensembl.gene.id%in%genes_65dph_all$ensembl_gene_id[!grepl("Z",genes_65dph_all$chromosome_name)]]<0)))
mat_A <- rbind(mat_A,c(sum(fulldata_output_adult_all$logFC.adult[fulldata_output_adult_all$ensembl.gene.id%in%genes_adult_all$ensembl_gene_id[!grepl("Z",genes_adult_all$chromosome_name)]]>0),
                       sum(fulldata_output_adult_all$logFC.adult[fulldata_output_adult_all$ensembl.gene.id%in%genes_adult_all$ensembl_gene_id[!grepl("Z",genes_adult_all$chromosome_name)]]<0)))

sum(mat_A[1:3,1])/sum(mat_A[1:3,])
mat_A[4,1]/sum(mat_A[4,])
pvals <- NULL
pvals <- c(pvals,chisq.test(cbind(mat_A[1,],mat_A[2,]))$p.value)
pvals <- c(pvals,chisq.test(cbind(mat_A[1,],mat_A[3,]))$p.value)
pvals <- c(pvals,chisq.test(cbind(mat_A[1,],mat_A[4,]))$p.value)
pvals <- c(pvals,chisq.test(cbind(mat_A[2,],mat_A[3,]))$p.value)
pvals <- c(pvals,chisq.test(cbind(mat_A[2,],mat_A[4,]))$p.value)
pvals <- c(pvals,chisq.test(cbind(mat_A[3,],mat_A[4,]))$p.value)

adj.pvals <- pvals*6
adj.pvals[adj.pvals>1] <- 1
adj.pvals

# Amount DE sex genes covered in RRBS #
rrbs_1dph_RNAseq <- read.table('RRBS_1dph_RNAseq.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_20dph_RNAseq <- read.table('RRBS_20dph_RNAseq.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_65dph_RNAseq <- read.table('RRBS_65dph_RNAseq.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_adult_RNAseq <- read.table('RRBS_adult_RNAseq.txt',header=T,sep="\t",stringsAsFactors = F)

RRBS_sex_covered <- c(unlist(rrbs_1dph_RNAseq[,grepl("gene",colnames(rrbs_1dph_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_1dph_RNAseq))]),
                      unlist(rrbs_20dph_RNAseq[,grepl("gene",colnames(rrbs_20dph_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_20dph_RNAseq))]),
                      unlist(rrbs_65dph_RNAseq[,grepl("gene",colnames(rrbs_65dph_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_65dph_RNAseq))]),
                      unlist(rrbs_adult_RNAseq[,grepl("gene",colnames(rrbs_adult_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_adult_RNAseq))]))
RRBS_sex_covered <- unique(RRBS_sex_covered[!is.na(RRBS_sex_covered)])
length(RRBS_sex_covered)


## Male and Female overexpressed genes
#######################################"

# Data
fulldata_output <- read.table("FullData_output.txt",header=T,sep="\t")

# help function
getZorA <- function(genelist){
  # print(head(genelist))
  ix <- rank(as.character(genelist))
  genelist_sorted <- sort(as.character(genelist))
  chr_sorted <- genes$chromosome_name[genes$ensembl_gene_id%in%genelist_sorted]
  # print(head(genelist_sorted[ix]))
  chr <- chr_sorted[ix]
  return(chr)
}

# Load in ensembl annotation
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

# Genelists
timepoints <- c("1dph","20dph","65dph","adult")
diff_expressed <- matrix(rep(0,8),ncol=4)
colnames(diff_expressed) <- timepoints
rownames(diff_expressed) <- c("over in M","under in M")
diff_expressed_Z <- diff_expressed
diff_expressed_A <- diff_expressed
counts_ZvsA <- diff_expressed
rownames(counts_ZvsA) <- c("Z","A")
genes_male <- NULL
genes_female <- NULL
genes_male_A <- NULL
genes_female_A <- NULL
genes_male_Z <- NULL
genes_female_Z <- NULL
for (i in 1:length(timepoints)){
  mask <- !is.na(fulldata_output[,grepl("Pvalue",colnames(fulldata_output)) & grepl(timepoints[i],colnames(fulldata_output))]) &
    (fulldata_output[,grepl("Pvalue",colnames(fulldata_output)) & grepl(timepoints[i],colnames(fulldata_output))]<0.05)
  fulldata_output_part <- fulldata_output[mask,]
  print(timepoints[i])
  counts_ZvsA[1,i] <- sum(grepl("Z",getZorA(fulldata_output_part$ensembl.gene.id)))
  counts_ZvsA[2,i] <- sum(!grepl("Z",getZorA(fulldata_output_part$ensembl.gene.id)))
  print(sum(grepl("Z",getZorA(fulldata_output_part$ensembl.gene.id)))/length(fulldata_output_part$ensembl.gene.id))
  
  print(dim(fulldata_output_part))
  logFC_part <- fulldata_output_part[,grepl("logFC",colnames(fulldata_output_part)) & grepl(timepoints[i],colnames(fulldata_output_part))]
  print(length(logFC_part))
  genes_male <- c(genes_male,as.character(fulldata_output_part$ensembl.gene.id[logFC_part>0]))
  genes_female <- c(genes_female,as.character(fulldata_output_part$ensembl.gene.id[logFC_part<0]))
  genes_male_A <- c(genes_male_A,as.character(fulldata_output_part$ensembl.gene.id[(logFC_part>0) & !(getZorA(fulldata_output_part$ensembl.gene.id)%in%c("Z","Z_random"))]))
  genes_female_A <- c(genes_female_A,as.character(fulldata_output_part$ensembl.gene.id[(logFC_part<0) & !(getZorA(fulldata_output_part$ensembl.gene.id)%in%c("Z","Z_random"))]))
  genes_male_Z <- c(genes_male_Z,as.character(fulldata_output_part$ensembl.gene.id[(logFC_part>0) & (getZorA(fulldata_output_part$ensembl.gene.id)%in%c("Z","Z_random"))]))
  genes_female_Z <- c(genes_female_Z,as.character(fulldata_output_part$ensembl.gene.id[(logFC_part<0) & (getZorA(fulldata_output_part$ensembl.gene.id)%in%c("Z","Z_random"))]))
  diff_expressed[1,i] <- sum(logFC_part>0)
  diff_expressed[2,i] <- sum(logFC_part<0)
  diff_expressed_A[1,i] <- sum(logFC_part[!(getZorA(fulldata_output_part$ensembl.gene.id)%in%c("Z","Z_random"))]>0)
  diff_expressed_A[2,i] <- sum(logFC_part[!(getZorA(fulldata_output_part$ensembl.gene.id)%in%c("Z","Z_random"))]<0)
  diff_expressed_Z[1,i] <- sum(logFC_part[getZorA(fulldata_output_part$ensembl.gene.id)%in%c("Z","Z_random")]>0)
  diff_expressed_Z[2,i] <- sum(logFC_part[getZorA(fulldata_output_part$ensembl.gene.id)%in%c("Z","Z_random")]<0)
}
counts_ZvsA
diff_expressed
diff_expressed_A
diff_expressed_Z

sum(diff_expressed_Z[1,])/sum(c(diff_expressed_Z[1,],diff_expressed_A[1,]))
sum(diff_expressed_Z[2,])/sum(c(diff_expressed_Z[2,],diff_expressed_A[2,]))

genes_female <- unique(genes_female)
genes_male <- unique(genes_male)
genes_female_A <- unique(genes_female_A)
genes_male_A <- unique(genes_male_A)
genes_female_Z <- unique(genes_female_Z)
genes_male_Z <- unique(genes_male_Z)

length(genes_male)
length(genes_female)
length(genes_male_A)
length(genes_female_A)
length(genes_male_Z)
length(genes_female_Z)

genes_male_Z <- genes_male_Z[!(genes_male_Z%in%genes_female_Z)]

write.table(genes_male_Z,"GeneLists/ListGenes_MaleOver.txt",col.names = F,row.names = F)
write.table(genes_female_Z,"GeneLists/ListGenes_FemaleOver.txt",col.names = F,row.names = F)



