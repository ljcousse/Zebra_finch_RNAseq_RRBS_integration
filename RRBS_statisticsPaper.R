#############################
#############################
##                         ##
## RRBS_statistics_Paper.R ##
##                         ##
#############################
#############################

library("edgeR")
library("BiSeq")
library("biomaRt")


##################
## Load in data ##
##################

## RNAseq
fulldata_output <- read.table(file="FullData_output.txt",header=T,sep="\t",stringsAsFactors = F)

## CpG
rrbs_1dph_full <- read.table('RRBS_1dph_full.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_1dph_RNAseq <- read.table('RRBS_1dph_RNAseq.txt',header=T,sep="\t",stringsAsFactors = F)

rrbs_20dph_full <- read.table('RRBS_20dph_full.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_20dph_RNAseq <- read.table('RRBS_20dph_RNAseq.txt',header=T,sep="\t",stringsAsFactors = F)

rrbs_65dph_full <- read.table('RRBS_65dph_full.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_65dph_RNAseq <- read.table('RRBS_65dph_RNAseq.txt',header=T,sep="\t",stringsAsFactors = F)

rrbs_adult_full <- read.table('RRBS_adult_full.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_adult_RNAseq <- read.table('RRBS_adult_RNAseq.txt',header=T,sep="\t",stringsAsFactors = F)

rrbs_age_full <- read.table('RRBS_age_full.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_age_RNAseq <- read.table('RRBS_age_RNAseq.txt',header=T,sep="\t",stringsAsFactors = F)


## CpH
rrbs_nonCpG_1dph_full <- read.table('RRBS_1dph_full_nonCpG.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_nonCpG_1dph_RNAseq <- read.table('RRBS_1dph_RNAseq_nonCpG.txt',header=T,sep="\t",stringsAsFactors = F)

rrbs_nonCpG_20dph_full <- read.table('RRBS_20dph_full_nonCpG.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_nonCpG_20dph_RNAseq <- read.table('RRBS_20dph_RNAseq_nonCpG.txt',header=T,sep="\t",stringsAsFactors = F)

rrbs_nonCpG_65dph_full <- read.table('RRBS_65dph_full_nonCpG.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_nonCpG_65dph_RNAseq <- read.table('RRBS_65dph_RNAseq_nonCpG.txt',header=T,sep="\t",stringsAsFactors = F)

rrbs_nonCpG_adult_full <- read.table('RRBS_adult_full_nonCpG.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_nonCpG_adult_RNAseq <- read.table('RRBS_adult_RNAseq_nonCpG.txt',header=T,sep="\t",stringsAsFactors = F)

rrbs_nonCpG_age_full <- read.table('RRBS_age_full_nonCpG.txt',header=T,sep="\t",stringsAsFactors = F)
rrbs_nonCpG_age_RNAseq <- read.table('RRBS_age_RNAseq_nonCpG.txt',header=T,sep="\t",stringsAsFactors = F)


########################
## General Statistics ##
########################

## CpG
print("CpG")
print("Full")
length(unique(paste(rrbs_1dph_full$chr[rrbs_1dph_full$adj.P.Val<0.05],rrbs_1dph_full$locations,sep="_")[rrbs_1dph_full$adj.P.Val<0.05]))
length(unique(paste(rrbs_20dph_full$chr[rrbs_20dph_full$adj.P.Val<0.05],rrbs_20dph_full$locations,sep="_")[rrbs_20dph_full$adj.P.Val<0.05]))
length(unique(paste(rrbs_65dph_full$chr[rrbs_65dph_full$adj.P.Val<0.05],rrbs_65dph_full$locations,sep="_")[rrbs_65dph_full$adj.P.Val<0.05]))
length(unique(paste(rrbs_adult_full$chr[rrbs_adult_full$adj.P.Val<0.05],rrbs_adult_full$locations,sep="_")[rrbs_adult_full$adj.P.Val<0.05]))
length(unique(paste(rrbs_age_full$chr[rrbs_age_full$adj.P.Val<0.05],rrbs_age_full$locations,sep="_")[rrbs_age_full$adj.P.Val<0.05]))

n_DM_full <- length(c(unique(paste(rrbs_1dph_full$chr[rrbs_1dph_full$adj.P.Val<0.05],rrbs_1dph_full$locations,sep="_")[rrbs_1dph_full$adj.P.Val<0.05]),
                      unique(paste(rrbs_20dph_full$chr[rrbs_20dph_full$adj.P.Val<0.05],rrbs_20dph_full$locations,sep="_")[rrbs_20dph_full$adj.P.Val<0.05]),
                      unique(paste(rrbs_65dph_full$chr[rrbs_65dph_full$adj.P.Val<0.05],rrbs_65dph_full$locations,sep="_")[rrbs_65dph_full$adj.P.Val<0.05]),
                      unique(paste(rrbs_adult_full$chr[rrbs_adult_full$adj.P.Val<0.05],rrbs_adult_full$locations,sep="_")[rrbs_adult_full$adj.P.Val<0.05]),
                      unique(paste(rrbs_age_full$chr[rrbs_age_full$adj.P.Val<0.05],rrbs_age_full$locations,sep="_")[rrbs_age_full$adj.P.Val<0.05])))

print("RNA_seq")
length(unique(paste(rrbs_1dph_RNAseq$chr[rrbs_1dph_RNAseq$adj.P.Val<0.05],rrbs_1dph_RNAseq$locations,sep="_")[rrbs_1dph_RNAseq$adj.P.Val<0.05]))
length(unique(paste(rrbs_20dph_RNAseq$chr[rrbs_20dph_RNAseq$adj.P.Val<0.05],rrbs_20dph_RNAseq$locations,sep="_")[rrbs_20dph_RNAseq$adj.P.Val<0.05]))
length(unique(paste(rrbs_65dph_RNAseq$chr[rrbs_65dph_RNAseq$adj.P.Val<0.05],rrbs_65dph_RNAseq$locations,sep="_")[rrbs_65dph_RNAseq$adj.P.Val<0.05]))
length(unique(paste(rrbs_adult_RNAseq$chr[rrbs_adult_RNAseq$adj.P.Val<0.05],rrbs_adult_RNAseq$locations,sep="_")[rrbs_adult_RNAseq$adj.P.Val<0.05]))
length(unique(paste(rrbs_age_RNAseq$chr[rrbs_age_RNAseq$adj.P.Val<0.05],rrbs_age_RNAseq$locations,sep="_")[rrbs_age_RNAseq$adj.P.Val<0.05]))

n_DM_RNAseq <- length(c(unique(paste(rrbs_1dph_RNAseq$chr[rrbs_1dph_RNAseq$adj.P.Val<0.05],rrbs_1dph_RNAseq$locations,sep="_")[rrbs_1dph_RNAseq$adj.P.Val<0.05]),
                      unique(paste(rrbs_20dph_RNAseq$chr[rrbs_20dph_RNAseq$adj.P.Val<0.05],rrbs_20dph_RNAseq$locations,sep="_")[rrbs_20dph_RNAseq$adj.P.Val<0.05]),
                      unique(paste(rrbs_65dph_RNAseq$chr[rrbs_65dph_RNAseq$adj.P.Val<0.05],rrbs_65dph_RNAseq$locations,sep="_")[rrbs_65dph_RNAseq$adj.P.Val<0.05]),
                      unique(paste(rrbs_adult_RNAseq$chr[rrbs_adult_RNAseq$adj.P.Val<0.05],rrbs_adult_RNAseq$locations,sep="_")[rrbs_adult_RNAseq$adj.P.Val<0.05]),
                      unique(paste(rrbs_age_RNAseq$chr[rrbs_age_RNAseq$adj.P.Val<0.05],rrbs_age_RNAseq$locations,sep="_")[rrbs_age_RNAseq$adj.P.Val<0.05])))

# Substract one of this result since there is an unique NA in this set!
length(unique(unlist(c(rrbs_1dph_RNAseq[,grepl("gene",colnames(rrbs_1dph_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_1dph_RNAseq))],
                rrbs_20dph_RNAseq[,grepl("gene",colnames(rrbs_20dph_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_20dph_RNAseq))],
                rrbs_65dph_RNAseq[,grepl("gene",colnames(rrbs_65dph_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_65dph_RNAseq))],
                rrbs_adult_RNAseq[,grepl("gene",colnames(rrbs_adult_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_adult_RNAseq))],
                rrbs_age_RNAseq[,grepl("gene",colnames(rrbs_age_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_age_RNAseq))]))))
sum(is.na(unique(unlist(c(rrbs_1dph_RNAseq[,grepl("gene",colnames(rrbs_1dph_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_1dph_RNAseq))],
                     rrbs_20dph_RNAseq[,grepl("gene",colnames(rrbs_20dph_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_20dph_RNAseq))],
                     rrbs_65dph_RNAseq[,grepl("gene",colnames(rrbs_65dph_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_65dph_RNAseq))],
                     rrbs_adult_RNAseq[,grepl("gene",colnames(rrbs_adult_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_adult_RNAseq))],
                     rrbs_age_RNAseq[,grepl("gene",colnames(rrbs_age_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_age_RNAseq))])))))
length(unique(unlist(c(rrbs_1dph_RNAseq[,grepl("gene",colnames(rrbs_1dph_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_1dph_RNAseq))],
                       rrbs_20dph_RNAseq[,grepl("gene",colnames(rrbs_20dph_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_20dph_RNAseq))],
                       rrbs_65dph_RNAseq[,grepl("gene",colnames(rrbs_65dph_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_65dph_RNAseq))],
                       rrbs_adult_RNAseq[,grepl("gene",colnames(rrbs_adult_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_adult_RNAseq))]))))


## CpH
print("CpH")
print("Full")
length(unique(paste(rrbs_nonCpG_1dph_full$chr[rrbs_nonCpG_1dph_full$adj.P.Val<0.05],rrbs_nonCpG_1dph_full$locations,sep="_")[rrbs_nonCpG_1dph_full$adj.P.Val<0.05]))
length(unique(paste(rrbs_nonCpG_20dph_full$chr[rrbs_nonCpG_20dph_full$adj.P.Val<0.05],rrbs_nonCpG_20dph_full$locations,sep="_")[rrbs_nonCpG_20dph_full$adj.P.Val<0.05]))
length(unique(paste(rrbs_nonCpG_65dph_full$chr[rrbs_nonCpG_65dph_full$adj.P.Val<0.05],rrbs_nonCpG_65dph_full$locations,sep="_")[rrbs_nonCpG_65dph_full$adj.P.Val<0.05]))
length(unique(paste(rrbs_nonCpG_adult_full$chr[rrbs_nonCpG_adult_full$adj.P.Val<0.05],rrbs_nonCpG_adult_full$locations,sep="_")[rrbs_nonCpG_adult_full$adj.P.Val<0.05]))
length(unique(paste(rrbs_nonCpG_age_full$chr[rrbs_nonCpG_age_full$adj.P.Val<0.05],rrbs_nonCpG_age_full$locations,sep="_")[rrbs_nonCpG_age_full$adj.P.Val<0.05]))

n_DM_full_nonCpG <- length(c(unique(paste(rrbs_nonCpG_1dph_full$chr[rrbs_nonCpG_1dph_full$adj.P.Val<0.05],rrbs_nonCpG_1dph_full$locations,sep="_")[rrbs_nonCpG_1dph_full$adj.P.Val<0.05]),
         unique(paste(rrbs_nonCpG_20dph_full$chr[rrbs_nonCpG_20dph_full$adj.P.Val<0.05],rrbs_nonCpG_20dph_full$locations,sep="_")[rrbs_nonCpG_20dph_full$adj.P.Val<0.05]),
         unique(paste(rrbs_nonCpG_65dph_full$chr[rrbs_nonCpG_65dph_full$adj.P.Val<0.05],rrbs_nonCpG_65dph_full$locations,sep="_")[rrbs_nonCpG_65dph_full$adj.P.Val<0.05]),
         unique(paste(rrbs_nonCpG_adult_full$chr[rrbs_nonCpG_adult_full$adj.P.Val<0.05],rrbs_nonCpG_adult_full$locations,sep="_")[rrbs_nonCpG_adult_full$adj.P.Val<0.05]),
         unique(paste(rrbs_nonCpG_age_full$chr[rrbs_nonCpG_age_full$adj.P.Val<0.05],rrbs_nonCpG_age_full$locations,sep="_")[rrbs_nonCpG_age_full$adj.P.Val<0.05])))

print("RNAseq")
length(unique(paste(rrbs_nonCpG_1dph_RNAseq$chr[rrbs_nonCpG_1dph_RNAseq$adj.P.Val<0.05],rrbs_nonCpG_1dph_RNAseq$locations,sep="_")[rrbs_nonCpG_1dph_RNAseq$adj.P.Val<0.05]))
length(unique(paste(rrbs_nonCpG_20dph_RNAseq$chr[rrbs_nonCpG_20dph_RNAseq$adj.P.Val<0.05],rrbs_nonCpG_20dph_RNAseq$locations,sep="_")[rrbs_nonCpG_20dph_RNAseq$adj.P.Val<0.05]))
length(unique(paste(rrbs_nonCpG_65dph_RNAseq$chr[rrbs_nonCpG_65dph_RNAseq$adj.P.Val<0.05],rrbs_nonCpG_65dph_RNAseq$locations,sep="_")[rrbs_nonCpG_65dph_RNAseq$adj.P.Val<0.05]))
length(unique(paste(rrbs_nonCpG_adult_RNAseq$chr[rrbs_nonCpG_adult_RNAseq$adj.P.Val<0.05],rrbs_nonCpG_adult_RNAseq$locations,sep="_")[rrbs_nonCpG_adult_RNAseq$adj.P.Val<0.05]))
length(unique(paste(rrbs_nonCpG_age_RNAseq$chr[rrbs_nonCpG_age_RNAseq$adj.P.Val<0.05],rrbs_nonCpG_age_RNAseq$locations,sep="_")[rrbs_nonCpG_age_RNAseq$adj.P.Val<0.05]))

n_DM_RNAseq_nonCpG <- length(c(unique(paste(rrbs_nonCpG_1dph_RNAseq$chr[rrbs_nonCpG_1dph_RNAseq$adj.P.Val<0.05],rrbs_nonCpG_1dph_RNAseq$locations,sep="_")[rrbs_nonCpG_1dph_RNAseq$adj.P.Val<0.05]),
                             unique(paste(rrbs_nonCpG_20dph_RNAseq$chr[rrbs_nonCpG_20dph_RNAseq$adj.P.Val<0.05],rrbs_nonCpG_20dph_RNAseq$locations,sep="_")[rrbs_nonCpG_20dph_RNAseq$adj.P.Val<0.05]),
                             unique(paste(rrbs_nonCpG_65dph_RNAseq$chr[rrbs_nonCpG_65dph_RNAseq$adj.P.Val<0.05],rrbs_nonCpG_65dph_RNAseq$locations,sep="_")[rrbs_nonCpG_65dph_RNAseq$adj.P.Val<0.05]),
                             unique(paste(rrbs_nonCpG_adult_RNAseq$chr[rrbs_nonCpG_adult_RNAseq$adj.P.Val<0.05],rrbs_nonCpG_adult_RNAseq$locations,sep="_")[rrbs_nonCpG_adult_RNAseq$adj.P.Val<0.05]),
                             unique(paste(rrbs_nonCpG_age_RNAseq$chr[rrbs_nonCpG_age_RNAseq$adj.P.Val<0.05],rrbs_nonCpG_age_RNAseq$locations,sep="_")[rrbs_nonCpG_age_RNAseq$adj.P.Val<0.05])))



chisq.test(cbind(c(n_DM_full,nrow(rrbs_1dph_full)),c(n_DM_full_nonCpG,nrow(rrbs_nonCpG_1dph_full))))
chisq.test(cbind(c(n_DM_full,nrow(rrbs_1dph_full)),c(n_DM_full_nonCpG,nrow(rrbs_nonCpG_1dph_full))))$expected


## Comparison amount of genes
# not sign
genes_RRBS <- unique(c(unlist(rrbs_1dph_full[,grepl("gene",colnames(rrbs_1dph_full)) & !grepl("gene_symbol",colnames(rrbs_1dph_full))]),
                       unlist(rrbs_20dph_full[,grepl("gene",colnames(rrbs_20dph_full)) & !grepl("gene_symbol",colnames(rrbs_20dph_full))]),
                       unlist(rrbs_65dph_full[,grepl("gene",colnames(rrbs_65dph_full)) & !grepl("gene_symbol",colnames(rrbs_65dph_full))]),
                       unlist(rrbs_adult_full[,grepl("gene",colnames(rrbs_adult_full)) & !grepl("gene_symbol",colnames(rrbs_adult_full))]),
                       unlist(rrbs_age_full[,grepl("gene",colnames(rrbs_age_full)) & !grepl("gene_symbol",colnames(rrbs_age_full))])))
genes_RRBS <- genes_RRBS[!is.na(genes_RRBS)]
genes_RRBS_nonCpG <- unique(c(unlist(rrbs_nonCpG_1dph_full[,grepl("gene",colnames(rrbs_nonCpG_1dph_full)) & !grepl("gene_symbol",colnames(rrbs_nonCpG_1dph_full))]),
                              unlist(rrbs_nonCpG_20dph_full[,grepl("gene",colnames(rrbs_nonCpG_20dph_full)) & !grepl("gene_symbol",colnames(rrbs_nonCpG_20dph_full))]),
                              unlist(rrbs_nonCpG_65dph_full[,grepl("gene",colnames(rrbs_nonCpG_65dph_full)) & !grepl("gene_symbol",colnames(rrbs_nonCpG_65dph_full))]),
                              unlist(rrbs_nonCpG_adult_full[,grepl("gene",colnames(rrbs_nonCpG_adult_full)) & !grepl("gene_symbol",colnames(rrbs_nonCpG_adult_full))]),
                              unlist(rrbs_nonCpG_age_full[,grepl("gene",colnames(rrbs_nonCpG_age_full)) & !grepl("gene_symbol",colnames(rrbs_nonCpG_age_full))])))
genes_RRBS_nonCpG <- genes_RRBS_nonCpG[!is.na(genes_RRBS_nonCpG)]
length(genes_RRBS)
sum(grepl("ENSTGUG",genes_RRBS))
length(genes_RRBS_nonCpG)
sum(grepl("ENSTGUG",genes_RRBS_nonCpG))

# sign
genes_RRBS_sign_1dph <- unique(unlist(rrbs_1dph_full[rrbs_1dph_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_1dph_full)) & !grepl("gene_symbol",colnames(rrbs_1dph_full))]))
genes_RRBS_sign_20dph <- unique(unlist(rrbs_20dph_full[rrbs_20dph_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_20dph_full)) & !grepl("gene_symbol",colnames(rrbs_20dph_full))]))
genes_RRBS_sign_65dph <- unique(unlist(rrbs_65dph_full[rrbs_65dph_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_65dph_full)) & !grepl("gene_symbol",colnames(rrbs_65dph_full))]))
genes_RRBS_sign_adult <- unique(unlist(rrbs_adult_full[rrbs_adult_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_adult_full)) & !grepl("gene_symbol",colnames(rrbs_adult_full))]))
genes_RRBS_sign_age <- unique(unlist(rrbs_age_full[rrbs_age_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_age_full)) & !grepl("gene_symbol",colnames(rrbs_age_full))]))

genes_RRBS_sign_1dph <- genes_RRBS_sign_1dph[!is.na(genes_RRBS_sign_1dph)]
genes_RRBS_sign_20dph <- genes_RRBS_sign_20dph[!is.na(genes_RRBS_sign_20dph)]
genes_RRBS_sign_65dph <- genes_RRBS_sign_65dph[!is.na(genes_RRBS_sign_65dph)]
genes_RRBS_sign_adult <- genes_RRBS_sign_adult[!is.na(genes_RRBS_sign_adult)]
genes_RRBS_sign_age <- genes_RRBS_sign_age[!is.na(genes_RRBS_sign_age)]

length(genes_RRBS_sign_1dph)
length(genes_RRBS_sign_20dph)
length(genes_RRBS_sign_65dph)
length(genes_RRBS_sign_adult)
length(genes_RRBS_sign_age)
length(unique(c(genes_RRBS_sign_1dph,genes_RRBS_sign_20dph,genes_RRBS_sign_65dph,genes_RRBS_sign_adult,genes_RRBS_sign_age)))


genes_RRBS_nonCpG_sign_1dph <- unique(unlist(rrbs_nonCpG_1dph_full[rrbs_nonCpG_1dph_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_nonCpG_1dph_full)) & !grepl("gene_symbol",colnames(rrbs_nonCpG_1dph_full))]))
genes_RRBS_nonCpG_sign_20dph <- unique(unlist(rrbs_nonCpG_20dph_full[rrbs_nonCpG_20dph_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_nonCpG_20dph_full)) & !grepl("gene_symbol",colnames(rrbs_nonCpG_20dph_full))]))
genes_RRBS_nonCpG_sign_65dph <- unique(unlist(rrbs_nonCpG_65dph_full[rrbs_nonCpG_65dph_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_nonCpG_65dph_full)) & !grepl("gene_symbol",colnames(rrbs_nonCpG_65dph_full))]))
genes_RRBS_nonCpG_sign_adult <- unique(unlist(rrbs_nonCpG_adult_full[rrbs_nonCpG_adult_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_nonCpG_adult_full)) & !grepl("gene_symbol",colnames(rrbs_nonCpG_adult_full))]))
genes_RRBS_nonCpG_sign_age <- unique(unlist(rrbs_nonCpG_age_full[rrbs_nonCpG_age_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_nonCpG_age_full)) & !grepl("gene_symbol",colnames(rrbs_nonCpG_age_full))]))

genes_RRBS_nonCpG_sign_1dph <- genes_RRBS_nonCpG_sign_1dph[!is.na(genes_RRBS_nonCpG_sign_1dph)]
genes_RRBS_nonCpG_sign_20dph <- genes_RRBS_nonCpG_sign_20dph[!is.na(genes_RRBS_nonCpG_sign_20dph)]
genes_RRBS_nonCpG_sign_65dph <- genes_RRBS_nonCpG_sign_65dph[!is.na(genes_RRBS_nonCpG_sign_65dph)]
genes_RRBS_nonCpG_sign_adult <- genes_RRBS_nonCpG_sign_adult[!is.na(genes_RRBS_nonCpG_sign_adult)]
genes_RRBS_nonCpG_sign_age <- genes_RRBS_nonCpG_sign_age[!is.na(genes_RRBS_nonCpG_sign_age)]

length(genes_RRBS_nonCpG_sign_1dph)
length(genes_RRBS_nonCpG_sign_20dph)
length(genes_RRBS_nonCpG_sign_65dph)
length(genes_RRBS_nonCpG_sign_adult)
length(genes_RRBS_nonCpG_sign_age)
length(unique(c(genes_RRBS_nonCpG_sign_1dph,genes_RRBS_nonCpG_sign_20dph,genes_RRBS_nonCpG_sign_65dph,genes_RRBS_nonCpG_sign_adult,genes_RRBS_nonCpG_sign_age)))




print("1dph")
print(sum(rrbs_1dph_full$adj.P.Val<0.05))
print(sum(rrbs_1dph_RNAseq$adj.P.Val<0.05))
print("20dph")
print(sum(rrbs_20dph_full$adj.P.Val<0.05))
print(sum(rrbs_20dph_RNAseq$adj.P.Val<0.05))
print("65dph")
print(sum(rrbs_65dph_full$adj.P.Val<0.05))
print(sum(rrbs_65dph_RNAseq$adj.P.Val<0.05))
print("adult")
print(sum(rrbs_adult_full$adj.P.Val<0.05))
print(sum(rrbs_adult_RNAseq$adj.P.Val<0.05))
print("age")
print(sum(rrbs_age_full$adj.P.Val<0.05))
print(sum(rrbs_age_RNAseq$adj.P.Val<0.05))


#######################################
## General correlation Meth vs Exprs ##
#######################################

## Additional data

# RNAseq
rnadata<-read.table("Plastoscine.txt",header=T,sep="\t",quote="")
rna_annot <- read.table("Plastoscine_annot.txt",header=F,sep="\t",quote="")
row.names(rnadata) <- unlist(rna_annot)

y <- DGEList(rnadata)
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep,  ,keep.lib.sizes=FALSE]
summary(keep)

logcpm <- rowMeans(cpm(y,log=T))
rna_annot <- rna_annot$V1[keep]

## CpG

# RRBS
load(file="allrrbs.Rda") # allrrbs.Rda can be retrieved using RRBS_analyis.R script

mincov <- 5

rrbs_CpG <-filterBySharedRegions(object=allrrbs,perc.samples=1,minCov=mincov)
rrbs.rel_CpG <- rawToRel(rrbs_CpG)

locs_CpG<-rowRanges(rrbs_CpG)

if (file.exists("Genes_RRBS.Rda")){
  load("Genes_RRBS.Rda")
  dim(genes)
} else {
  ensembl=useMart(biomart="ensembl",host="jan2019.archive.ensembl.org")
  # ensembl=useMart("ensembl",host="http://sep2019.archive.ensembl.org")
  ensembl = useDataset("tguttata_gene_ensembl",mart=ensembl)
  genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','transcription_start_site'),mart = ensembl)
  colnames(genes) <- c("gene","gene_symbol","chr","start","end","TSS")
  dim(genes) 
  genes <- genes[!duplicated(genes$gene),]
  dim(genes)
  save(genes,file="Genes_RRBS.Rda")
}

sense <- abs(genes$start-genes$TSS)<abs(genes$end-genes$TSS)

new_start <- cbind(c(genes$start[sense]-1999,genes$end[!sense]-499),
                   c(genes$start[sense]-1999,genes$start[!sense]+1))
new_start <- apply(new_start,1,max)
new_end <- cbind(c(genes$start[sense]+499,genes$end[!sense]+1999),
                 c(genes$end[sense]-1,genes$end[!sense]+1999))
new_end <- apply(new_end,1,min)
prom_gr <- GRanges(data.frame(chr=c(genes$chr[sense],genes$chr[!sense]),
                              start=new_start,
                              end=new_end,
                              gene=c(genes$gene[sense],genes$gene[!sense]),
                              gene_symbol=c(genes$gene_symbol[sense],genes$gene_symbol[!sense])))
overlap_p <- subsetByOverlaps(locs_CpG,prom_gr)
promoter_CpG <- rownames(as.data.frame(locs_CpG))%in%rownames(as.data.frame(overlap_p))

new_start <- cbind(c(genes$start[sense]+500,genes$start[!sense]),
                   c(genes$end[sense],genes$start[!sense]))
new_start <- apply(new_start,1,min)
new_end <- cbind(c(genes$end[sense],genes$end[!sense]-500),
                 c(genes$end[sense],genes$start[!sense]))
new_end <- apply(new_end,1,max)
genic_gr <- GRanges(data.frame(chr=c(genes$chr[sense],genes$chr[!sense]),
                               start=new_start,
                               end=new_end,
                               gene=c(genes$gene[sense],genes$gene[!sense]),
                               gene_symbol=c(genes$gene_symbol[sense],genes$gene_symbol[!sense])))
overlap_g <- subsetByOverlaps(locs_CpG,genic_gr)
genic_CpG <- rownames(as.data.frame(locs_CpG))%in%rownames(as.data.frame(overlap_g))


genic_CpG[genic_CpG & promoter_CpG] <- FALSE

overlaps_genic <- findOverlaps(locs_CpG,genic_gr)
overlaps_prom <- findOverlaps(locs_CpG,prom_gr)

all_overlap <- c(as.character(genic_gr[subjectHits(overlaps_genic)]$gene),as.character(prom_gr[subjectHits(overlaps_prom)]$gene))
all_overlap <- cbind(all_overlap,c(as.character(genic_gr[subjectHits(overlaps_genic)]$gene_symbol),as.character(prom_gr[subjectHits(overlaps_prom)]$gene_symbol)))
all_overlap <- cbind(all_overlap,c(locs_CpG[queryHits(overlaps_genic)]@ranges@NAMES,locs_CpG[queryHits(overlaps_prom)]@ranges@NAMES))
colnames(all_overlap) <- c("gene","gene_symbol","rowname")


rn <- as.matrix(locs_CpG@ranges@NAMES,ncol=1)
colnames(rn) <- "rowname"

merger <- merge(rn,all_overlap,all=T)
merger <- merger[sort(as.numeric(as.character(merger$rowname)),index.return=T)$ix,]

merger$gene <- as.character(merger$gene)
merger$gene_symbol <- as.character(merger$gene_symbol)

gene_CpG <- as.character(merger$gene[!duplicated(merger$rowname)])
gene_symbol_CpG <- as.character(merger$gene_symbol[!duplicated(merger$rowname)])

max_n <- max(table(merger$rowname))
print(max_n)
gene_CpG <- cbind(gene_CpG,matrix(rep(NA,length(gene_CpG)*(max_n-1)),ncol=max_n-1))
gene_symbol_CpG <- cbind(gene_symbol_CpG,matrix(rep(NA,length(gene_symbol_CpG)*(max_n-1)),ncol=max_n-1))
for (i in 2:max_n){
  n_times <- merger[merger$rowname%in%names(table(merger$rowname)[table(merger$rowname)==i]),]
  if(nrow(n_times)==0){
    print(paste("no locus with ",i,"-fold annotation",sep=""))
    next
  }
  if (i==2 | nrow(n_times)==i){
    gene_CpG[locs_CpG@ranges@NAMES%in%unique(n_times$rowname),2:i] <- n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]
    gene_symbol_CpG[locs_CpG@ranges@NAMES%in%unique(n_times$rowname),2:i] <- n_times$gene_symbol[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]
  } else {
    gene_CpG[locs_CpG@ranges@NAMES%in%unique(n_times$rowname),2:i] <- t(matrix(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]],nrow=i-1))
    gene_symbol_CpG[locs_CpG@ranges@NAMES%in%unique(n_times$rowname),2:i] <- t(matrix(n_times$gene_symbol[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]],nrow=i-1))
  }
  
}

minmethreads <- 6

inclocsM2_CpG <- rowSums(methReads(rrbs_CpG))>=minmethreads

locs_CpG <- locs_CpG[inclocsM2_CpG]
rrbsredB_CpG <- rrbs_CpG[inclocsM2_CpG,]

genic_CpG <- genic_CpG[inclocsM2_CpG]
promoter_CpG <- promoter_CpG[inclocsM2_CpG]

gene_CpG <- gene_CpG[inclocsM2_CpG,]
gene_symbol_CpG <- gene_symbol_CpG[inclocsM2_CpG,]

rrbsredB_CpG_meth <- methLevel(rawToRel(rrbsredB_CpG))

# Correlation
uniq_annot_RNAseq <- rowSums(is.na(rrbs_age_RNAseq[,grepl("gene_alt",colnames(rrbs_age_RNAseq))]))==10 

rrbs_age_RNAseq_uniq <- rrbs_age_RNAseq[uniq_annot_RNAseq,]
rrbs_age_RNAseq_prom <- rrbs_age_RNAseq[uniq_annot_RNAseq & rrbs_age_RNAseq$func=="promoter",]
rrbs_age_RNAseq_genic <- rrbs_age_RNAseq[uniq_annot_RNAseq & rrbs_age_RNAseq$func=="genic",]

colnumbers <- which(grepl("gene",colnames(rrbs_age_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_age_RNAseq)))

rnaseq_logcpm <- NULL
rrbs_meth <- NULL
genes_uniq <- unique(rrbs_age_RNAseq_uniq$gene)
print(length(genes_uniq))
for(i in 1:length(genes_uniq)){
  if (i %% 100 == 0){print(i)}
  rnaseq_logcpm <- c(rnaseq_logcpm,logcpm[rna_annot==genes_uniq[i]])
  rrbs_meth <- c(rrbs_meth,mean(rrbsredB_CpG_meth[(gene_CpG[,1]==genes_uniq[i]) & !is.na(gene_CpG[,1]),]))
}
cor.test(rrbs_meth,rnaseq_logcpm)

rnaseq_logcpm_prom <- NULL
rrbs_meth_prom <- NULL
genes_uniq_prom <- unique(rrbs_age_RNAseq_prom$gene)
print(length(genes_uniq_prom))
for(i in 1:length(genes_uniq_prom)){
  if (i %% 100 == 0){print(i)}
  rnaseq_logcpm_prom  <- c(rnaseq_logcpm_prom,logcpm[rna_annot==genes_uniq_prom[i]])
  rrbs_meth_prom  <- c(rrbs_meth_prom,mean(rrbsredB_CpG_meth[(gene_CpG[,1]==genes_uniq_prom[i]) & !is.na(gene_CpG[,1]) & promoter_CpG,]))
}
cor.test(rrbs_meth_prom,rnaseq_logcpm_prom)


rnaseq_logcpm_genic <- NULL
rrbs_meth_genic <- NULL
genes_uniq_genic <- unique(rrbs_age_RNAseq_genic$gene)
print(length(genes_uniq_genic))
for(i in 1:length(genes_uniq_genic)){
  if (i %% 100 == 0){print(i)}
  rnaseq_logcpm_genic  <- c(rnaseq_logcpm_genic,logcpm[rna_annot==genes_uniq_genic[i]])
  rrbs_meth_genic  <- c(rrbs_meth_genic,mean(rrbsredB_CpG_meth[(gene_CpG[,1]==genes_uniq_genic[i]) & !is.na(gene_CpG[,1]) & genic_CpG,]))
}
cor.test(rrbs_meth_genic,rnaseq_logcpm_genic)

dat_prom_CpG_summ <- data.frame(Expression=rnaseq_logcpm_prom,
                                Methylation=rrbs_meth_prom,
                                stringsAsFactors = F)
lm_prom_CpG_summ <-  lm(Expression~Methylation,dat_prom_CpG_summ)
dat_genic_CpG_summ <- data.frame(Expression=rnaseq_logcpm_genic,
                                Methylation=rrbs_meth_genic,
                                stringsAsFactors = F)
lm_genic_CpG_summ <-  lm(Expression~Methylation,dat_genic_CpG_summ)
tiff("summ_scatter.tiff",height=4,width=8,units="in",res=150)
par(mfrow=c(1,2))
plot(dat_prom_CpG_summ$Methylation,dat_prom_CpG_summ$Expression,main="Promoter",pch=20)
abline(a=summary(lm_prom_CpG_summ)$coefficients[1,1],b=summary(lm_prom_CpG_summ)$coefficients[2,1],col="red")
plot(dat_genic_CpG_summ$Methylation,dat_genic_CpG_summ$Expression,main="Genic",pch=20)
abline(a=summary(lm_genic_CpG_summ)$coefficients[1,1],b=summary(lm_genic_CpG_summ)$coefficients[2,1],col="red")
dev.off()

## nonCpG

# RRBS

# Load in RRBS annotation
load("rrbs_nonCpG.Rda")
rrbs_nonCpG <- rrbs
rownames(rrbs_nonCpG) <- 1:nrow(rrbs_nonCpG) 
rrbs.rel_nonCpG  <- rawToRel(rrbs_nonCpG)

locs_nonCpG <-rowRanges(rrbs_nonCpG)

if (file.exists("Genes_RRBS.Rda")){
  load("Genes_RRBS.Rda")
  dim(genes)
} else {
  ensembl=useMart(biomart="ensembl",host="jan2019.archive.ensembl.org")
  # ensembl=useMart("ensembl",host="http://sep2019.archive.ensembl.org")
  ensembl = useDataset("tguttata_gene_ensembl",mart=ensembl)
  genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','transcription_start_site'),mart = ensembl)
  colnames(genes) <- c("gene","gene_symbol","chr","start","end","TSS")
  dim(genes) 
  genes <- genes[!duplicated(genes$gene),]
  dim(genes)
  save(genes,file="Genes_RRBS.Rda")
}
sense <- abs(genes$start-genes$TSS)<abs(genes$end-genes$TSS)

new_start <- cbind(c(genes$start[sense]-1999,genes$end[!sense]-499),
                   c(genes$start[sense]-1999,genes$start[!sense]+1))
new_start <- apply(new_start,1,max)
new_end <- cbind(c(genes$start[sense]+499,genes$end[!sense]+1999),
                 c(genes$end[sense]-1,genes$end[!sense]+1999))
new_end <- apply(new_end,1,min)
prom_gr <- GRanges(data.frame(chr=c(genes$chr[sense],genes$chr[!sense]),
                              start=new_start,
                              end=new_end,
                              gene=c(genes$gene[sense],genes$gene[!sense]),
                              gene_symbol=c(genes$gene_symbol[sense],genes$gene_symbol[!sense])))
overlap_p <- subsetByOverlaps(locs_nonCpG,prom_gr)
promoter_nonCpG <- rownames(as.data.frame(locs_nonCpG))%in%rownames(as.data.frame(overlap_p))

new_start <- cbind(c(genes$start[sense]+500,genes$start[!sense]),
                   c(genes$end[sense],genes$start[!sense]))
new_start <- apply(new_start,1,min)
new_end <- cbind(c(genes$end[sense],genes$end[!sense]-500),
                 c(genes$end[sense],genes$start[!sense]))
new_end <- apply(new_end,1,max)
genic_gr <- GRanges(data.frame(chr=c(genes$chr[sense],genes$chr[!sense]),
                               start=new_start,
                               end=new_end,
                               gene=c(genes$gene[sense],genes$gene[!sense]),
                               gene_symbol=c(genes$gene_symbol[sense],genes$gene_symbol[!sense])))
overlap_g <- subsetByOverlaps(locs_nonCpG,genic_gr)
genic_nonCpG <- rownames(as.data.frame(locs_nonCpG))%in%rownames(as.data.frame(overlap_g))


genic_nonCpG[genic_nonCpG & promoter_nonCpG] <- FALSE

overlaps_genic <- findOverlaps(locs_nonCpG,genic_gr)
overlaps_prom <- findOverlaps(locs_nonCpG,prom_gr)

all_overlap <- c(as.character(genic_gr[subjectHits(overlaps_genic)]$gene),as.character(prom_gr[subjectHits(overlaps_prom)]$gene))
all_overlap <- cbind(all_overlap,c(as.character(genic_gr[subjectHits(overlaps_genic)]$gene_symbol),as.character(prom_gr[subjectHits(overlaps_prom)]$gene_symbol)))
all_overlap <- cbind(all_overlap,c(locs_nonCpG[queryHits(overlaps_genic)]@ranges@NAMES,locs_nonCpG[queryHits(overlaps_prom)]@ranges@NAMES))
colnames(all_overlap) <- c("gene","gene_symbol","rowname")


rn <- as.matrix(locs_nonCpG@ranges@NAMES,ncol=1)
colnames(rn) <- "rowname"

merger <- merge(rn,all_overlap,all=T)
merger <- merger[sort(as.numeric(as.character(merger$rowname)),index.return=T)$ix,]

merger$gene <- as.character(merger$gene)
merger$gene_symbol <- as.character(merger$gene_symbol)

gene_nonCpG <- as.character(merger$gene[!duplicated(merger$rowname)])
gene_symbol_nonCpG <- as.character(merger$gene_symbol[!duplicated(merger$rowname)])

max_n <- max(table(merger$rowname))
print(max_n)
gene_nonCpG <- cbind(gene_nonCpG,matrix(rep(NA,length(gene_nonCpG)*(max_n-1)),ncol=max_n-1))
gene_symbol_nonCpG <- cbind(gene_symbol_nonCpG,matrix(rep(NA,length(gene_symbol_nonCpG)*(max_n-1)),ncol=max_n-1))
for (i in 2:max_n){
  n_times <- merger[merger$rowname%in%names(table(merger$rowname)[table(merger$rowname)==i]),]
  if(nrow(n_times)==0){
    print(paste("no locus with ",i,"-fold annotation",sep=""))
    next
  }
  if (i==2 | nrow(n_times)==i){
    gene_nonCpG[locs_nonCpG@ranges@NAMES%in%unique(n_times$rowname),2:i] <- n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]
    gene_symbol_nonCpG[locs_nonCpG@ranges@NAMES%in%unique(n_times$rowname),2:i] <- n_times$gene_symbol[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]
  } else {
    gene_nonCpG[locs_nonCpG@ranges@NAMES%in%unique(n_times$rowname),2:i] <- t(matrix(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]],nrow=i-1))
    gene_symbol_nonCpG[locs_nonCpG@ranges@NAMES%in%unique(n_times$rowname),2:i] <- t(matrix(n_times$gene_symbol[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]],nrow=i-1))
  }
  
}

minmethreads <- 6

inclocsM2_nonCpG <- rowSums(methReads(rrbs_nonCpG))>=minmethreads

locs_nonCpG <- locs_nonCpG[inclocsM2_nonCpG]
rrbsredB_nonCpG <- rrbs_nonCpG[inclocsM2_nonCpG,]

genic_nonCpG <- genic_nonCpG[inclocsM2_nonCpG]
promoter_nonCpG <- promoter_nonCpG[inclocsM2_nonCpG]

gene_nonCpG <- gene_nonCpG[inclocsM2_nonCpG,]
gene_symbol_nonCpG <- gene_symbol_nonCpG[inclocsM2_nonCpG,]

rrbsredB_nonCpG_meth <- methLevel(rawToRel(rrbsredB_nonCpG))

# Correlation
uniq_nonCpG_annot_RNAseq <- rowSums(is.na(rrbs_nonCpG_age_RNAseq[,grepl("gene_alt",colnames(rrbs_age_RNAseq))]))==10 

rrbs_nonCpG_age_RNAseq_uniq <- rrbs_nonCpG_age_RNAseq[uniq_nonCpG_annot_RNAseq,]
rrbs_nonCpG_age_RNAseq_prom <- rrbs_nonCpG_age_RNAseq[uniq_nonCpG_annot_RNAseq & rrbs_nonCpG_age_RNAseq$func=="promoter",]
rrbs_nonCpG_age_RNAseq_genic <- rrbs_nonCpG_age_RNAseq[uniq_nonCpG_annot_RNAseq & rrbs_nonCpG_age_RNAseq$func=="genic",]


rnaseq_logcpm_nonCpG <- NULL
rrbs_meth_nonCpG <- NULL
genes_uniq_nonCpG <- unique(rrbs_nonCpG_age_RNAseq_uniq$gene)
for(i in 1:length(genes_uniq_nonCpG)){
  if (i %% 100 == 0){print(i)}
  rnaseq_logcpm_nonCpG  <- c(rnaseq_logcpm_nonCpG,logcpm[rna_annot==genes_uniq_nonCpG[i]])
  rrbs_meth_nonCpG  <- c(rrbs_meth_nonCpG,mean(rrbsredB_nonCpG_meth[(gene_nonCpG[,1]==genes_uniq_nonCpG[i]) & !is.na(gene_nonCpG[,1]),]))
}
cor.test(rrbs_meth_nonCpG,rnaseq_logcpm_nonCpG)


rnaseq_logcpm_nonCpG_prom <- NULL
rrbs_meth_nonCpG_prom <- NULL
genes_uniq_nonCpG_prom <- unique(rrbs_nonCpG_age_RNAseq_prom$gene)
for(i in 1:length(genes_uniq_nonCpG_prom)){
  if (i %% 100 == 0){print(i)}
  rnaseq_logcpm_nonCpG_prom  <- c(rnaseq_logcpm_nonCpG_prom,logcpm[rna_annot==genes_uniq_nonCpG_prom[i]])
  rrbs_meth_nonCpG_prom  <- c(rrbs_meth_nonCpG_prom,mean(rrbsredB_nonCpG_meth[(gene_nonCpG[,1]==genes_uniq_nonCpG_prom[i]) & !is.na(gene_nonCpG[,1]) & promoter_nonCpG,]))
}
cor.test(rrbs_meth_nonCpG_prom,rnaseq_logcpm_nonCpG_prom)


rnaseq_logcpm_nonCpG_genic <- NULL
rrbs_meth_nonCpG_genic <- NULL
genes_uniq_nonCpG_genic <- unique(rrbs_nonCpG_age_RNAseq_genic$gene)
for(i in 1:length(genes_uniq_nonCpG_genic)){
  if (i %% 100 == 0){print(i)}
  rnaseq_logcpm_nonCpG_genic  <- c(rnaseq_logcpm_nonCpG_genic,logcpm[rna_annot==genes_uniq_nonCpG_genic[i]])
  rrbs_meth_nonCpG_genic  <- c(rrbs_meth_nonCpG_genic,mean(rrbsredB_nonCpG_meth[(gene_nonCpG[,1]==genes_uniq_nonCpG_genic[i]) & !is.na(gene_nonCpG[,1]) & genic_nonCpG,]))
}
cor.test(rrbs_meth_nonCpG_genic,rnaseq_logcpm_nonCpG_genic)



##################
## Mixed Models ##
##################


## CpG

rrbsredB.rel_CpG <- rawToRel(rrbsredB_CpG)
rrbsredB.M_CpG <- rrbsredB.rel_CpG
ct <- 0.01
for (i in 1:ncol(rrbsredB.rel_CpG)){
  ML<-methLevel(rrbsredB.rel_CpG)[,i]
  methLevel(rrbsredB.M_CpG)[,i]<-log2((ML+ct)/(1-ML+ct))
}

logcpm_sample <- cpm(y,log=T)

order_col_CpG <- sort(as.numeric(gsub("O","",gsub("B","",colnames(rrbsredB.M_CpG)))),index.return=T)$ix
head(rrbsredB.M_CpG[,order_col_CpG])

# Model 1 Prom
dat_lme_1_prom_CpG <- NULL
genes_uniq_prom <- unique(rrbs_age_RNAseq_prom$gene)
print(length(genes_uniq_prom))
for(i in 1:length(genes_uniq_prom)){
  if (i %% 100 == 0){print(i)}
  if (sum((gene_CpG[,1]==genes_uniq_prom[i]) & !is.na(gene_CpG[,1]) & promoter_CpG)>1){
    dat_lme_1_prom_CpG <- rbind(dat_lme_1_prom_CpG,
                                cbind(logcpm_sample[rna_annot==genes_uniq_prom[i]],
                                      colMeans(methLevel(rrbsredB.M_CpG)[(gene_CpG[,1]==genes_uniq_prom[i]) & !is.na(gene_CpG[,1]) & promoter_CpG,order_col_CpG]),
                                      rep(genes_uniq_prom[i],length(colnames(rrbsredB.M_CpG))),
                                      colnames(rrbsredB.M_CpG)[order_col_CpG]))
  } else {
    dat_lme_1_prom_CpG <- rbind(dat_lme_1_prom_CpG,
                                cbind(logcpm_sample[rna_annot==genes_uniq_prom[i]],  
                                      methLevel(rrbsredB.M_CpG)[(gene_CpG[,1]==genes_uniq_prom[i]) & !is.na(gene_CpG[,1]) & promoter_CpG,order_col_CpG],
                                      rep(genes_uniq_prom[i],length(colnames(rrbsredB.M_CpG))),
                                      colnames(rrbsredB.M_CpG)[order_col_CpG]))
  }
}
dat_lme_1_prom_CpG <- data.frame(dat_lme_1_prom_CpG,stringsAsFactors = F)
colnames(dat_lme_1_prom_CpG) <- c("Expression","Methylation","Gene","Sample")
dat_lme_1_prom_CpG$Expression <- as.numeric(dat_lme_1_prom_CpG$Expression)  
dat_lme_1_prom_CpG$Methylation <- as.numeric(dat_lme_1_prom_CpG$Methylation) 
dat_lme_1_prom_CpG$Sample <- as.factor(dat_lme_1_prom_CpG$Sample) 
dat_lme_1_prom_CpG$Gene <- as.factor(dat_lme_1_prom_CpG$Gene) 
dim(dat_lme_1_prom_CpG)
head(dat_lme_1_prom_CpG)
str(dat_lme_1_prom_CpG)
write.table(dat_lme_1_prom_CpG,file="dat_lme_1_prom_CpG.txt",col.names = T,row.names = F,quote=F,sep="\t")

mixed_lmer_1_prom_CpG <- lmerTest::lmer(Expression ~ Methylation + (1|Sample) + (1+ Methylation|Gene), data = dat_lme_1_prom_CpG)
summary(mixed_lmer_1_prom_CpG)

# dat_lme_1_prom_CpG_scaled <- dat_lme_1_prom_CpG
# dat_lme_1_prom_CpG_scaled[,1:2] <- scale(dat_lme_1_prom_CpG_scaled[, 1:2], center = T, scale = T)
# 
# mixed_lmer_1_prom_CpG_scaled <- lmerTest::lmer(Expression ~ Methylation + (1|Sample) + (1+ Methylation|Gene), data = dat_lme_1_prom_CpG_scaled)
# summary(mixed_lmer_1_prom_CpG_scaled)


# Model 1 Genic
dat_lme_1_genic_CpG <- NULL
genes_uniq_genic <- unique(rrbs_age_RNAseq_genic$gene)
print(length(genes_uniq_genic))
for(i in 1:length(genes_uniq_genic)){
  if (i %% 100 == 0){print(i)}
  if (sum((gene_CpG[,1]==genes_uniq_genic[i]) & !is.na(gene_CpG[,1]) & genic_CpG)>1){
    dat_lme_1_genic_CpG <- rbind(dat_lme_1_genic_CpG,
                                 cbind(logcpm_sample[rna_annot==genes_uniq_genic[i]],
                                      colMeans(methLevel(rrbsredB.M_CpG)[(gene_CpG[,1]==genes_uniq_genic[i]) & !is.na(gene_CpG[,1]) & genic_CpG,order_col_CpG]),
                                      rep(genes_uniq_genic[i],length(colnames(rrbsredB.M_CpG))),
                                      colnames(rrbsredB.M_CpG)[order_col_CpG]))
  } else {
    dat_lme_1_genic_CpG <- rbind(dat_lme_1_genic_CpG,
                                 cbind(logcpm_sample[rna_annot==genes_uniq_genic[i]],
                                       methLevel(rrbsredB.M_CpG)[(gene_CpG[,1]==genes_uniq_genic[i]) & !is.na(gene_CpG[,1]) & genic_CpG,order_col_CpG],
                                       rep(genes_uniq_genic[i],length(colnames(rrbsredB.M_CpG))),
                                       colnames(rrbsredB.M_CpG)[order_col_CpG]))
  }
}
dat_lme_1_genic_CpG <- data.frame(dat_lme_1_genic_CpG,stringsAsFactors = F)
colnames(dat_lme_1_genic_CpG) <- c("Expression","Methylation","Gene","Sample")
dat_lme_1_genic_CpG$Expression <- as.numeric(dat_lme_1_genic_CpG$Expression)  
dat_lme_1_genic_CpG$Methylation <- as.numeric(dat_lme_1_genic_CpG$Methylation) 
dat_lme_1_genic_CpG$Sample <- as.factor(dat_lme_1_genic_CpG$Sample) 
dat_lme_1_genic_CpG$Gene <- as.factor(dat_lme_1_genic_CpG$Gene) 
dim(dat_lme_1_genic_CpG)
head(dat_lme_1_genic_CpG)
str(dat_lme_1_genic_CpG)
write.table(dat_lme_1_genic_CpG,file="dat_lme_1_genic_CpG.txt",col.names = T,row.names = F,quote=F,sep="\t")

mixed_lmer_1_genic_CpG <- lmerTest::lmer(Expression ~ Methylation + (1|Sample) + (1+ Methylation|Gene), data = dat_lme_1_genic_CpG)
summary(mixed_lmer_1_genic_CpG)

# dat_lme_1_genic_CpG_scaled <- dat_lme_1_genic_CpG
# dat_lme_1_genic_CpG_scaled[,1:2] <- scale(dat_lme_1_genic_CpG_scaled[, 1:2], center = T, scale = T)
# 
# mixed_lmer_1_genic_CpG_scaled <- lmerTest::lmer(Expression ~ Methylation + (1|Sample) + (1+ Methylation|Gene), data = dat_lme_1_genic_CpG_scaled)
# summary(mixed_lmer_1_genic_CpG_scaled)


# Plots 

library(readbitmap)
library(grid)
library(gridExtra)
library(ggplot2)
tiff("mixed_lmer_prom_qc_A.tiff",height=4,width=4,units="in",res=75)
qqnorm(resid(mixed_lmer_1_prom_CpG_scaled))
qqline(resid(mixed_lmer_1_prom_CpG_scaled))
dev.off()
tiff("mixed_lmer_prom_qc_B.tiff",height=4,width=4,units="in",res=75)
plot(mixed_lmer_1_prom_CpG_scaled)
dev.off()
ll <- list("mixed_lmer_prom_qc_A.tiff","mixed_lmer_prom_qc_B.tiff")
plots <- lapply(ll,function(x){
  img <- as.raster(read.bitmap(paste(x,sep="")))
  rasterGrob(img, interpolate = FALSE)
})
ggsave("mixed_lmer_prom_qc.tiff",width=8, height=4, units="in",
       marrangeGrob(grobs = plots, nrow=1, ncol=2,top=NULL))

lm_prom_CpG <- lm(Expression~Methylation,dat_lme_1_prom_CpG)
lm_genic_CpG <- lm(Expression~Methylation,dat_lme_1_genic_CpG)
tiff("mixed_lmer_scatter.tiff",height=4,width=8,units="in",res=150)
par(mfrow=c(1,2))
plot(dat_lme_1_prom_CpG$Methylation,dat_lme_1_prom_CpG$Expression,main="Promoter",pch=20)
abline(a=summary(lm_prom_CpG)$coefficients[1,1],b=summary(lm_prom_CpG)$coefficients[2,1],col="red")
plot(dat_lme_1_genic_CpG$Methylation,dat_lme_1_genic_CpG$Expression,main="Genic",pch=20)
abline(a=summary(lm_genic_CpG)$coefficients[1,1],b=summary(lm_genic_CpG)$coefficients[2,1],col="red")
dev.off()

dat_lme_1_prom_CpG_mean <- NULL
for (i in 1:length(unique(dat_lme_1_prom_CpG$Gene))){
  dat_lme_1_prom_CpG_mean <- rbind(dat_lme_1_prom_CpG_mean,
                                   c(mean(dat_lme_1_prom_CpG$Methylation[dat_lme_1_prom_CpG$Gene==unique(dat_lme_1_prom_CpG$Gene)[i]]),
                                     mean(dat_lme_1_prom_CpG$Expression[dat_lme_1_prom_CpG$Gene==unique(dat_lme_1_prom_CpG$Gene)[i]])))
}
cor.test(dat_lme_1_prom_CpG_mean[,1],dat_lme_1_prom_CpG_mean[,2],method="pearson")
dat_lme_1_genic_CpG_mean <- NULL
for (i in 1:length(unique(dat_lme_1_genic_CpG$Gene))){
  dat_lme_1_genic_CpG_mean <- rbind(dat_lme_1_genic_CpG_mean,
                                    c(mean(dat_lme_1_genic_CpG$Methylation[dat_lme_1_genic_CpG$Gene==unique(dat_lme_1_genic_CpG$Gene)[i]]),
                                      mean(dat_lme_1_genic_CpG$Expression[dat_lme_1_genic_CpG$Gene==unique(dat_lme_1_genic_CpG$Gene)[i]])))
}
cor.test(dat_lme_1_genic_CpG_mean[,1],dat_lme_1_genic_CpG_mean[,2],method="pearson")


## nonCpG


rrbsredB.rel_nonCpG <- rawToRel(rrbsredB_nonCpG)
rrbsredB.M_nonCpG <- rrbsredB.rel_nonCpG
ct <- 0.01
for (i in 1:ncol(rrbsredB.rel_nonCpG)){
  ML<-methLevel(rrbsredB.rel_nonCpG)[,i]
  methLevel(rrbsredB.M_nonCpG)[,i]<-log2((ML+ct)/(1-ML+ct))
}

logcpm_sample <- cpm(y,log=T)

order_col_nonCpG <- sort(as.numeric(gsub("O","",gsub("B","",colnames(rrbsredB.M_nonCpG)))),index.return=T)$ix
head(rrbsredB.M_nonCpG[,order_col_nonCpG])

# Model 1 Prom
dat_lme_1_prom_nonCpG <- NULL
genes_uniq_prom <- unique(rrbs_nonCpG_age_RNAseq_prom$gene)
print(length(genes_uniq_prom))
for(i in 1:length(genes_uniq_prom)){
  if (i %% 100 == 0){print(i)}
  if (sum((gene_nonCpG[,1]==genes_uniq_prom[i]) & !is.na(gene_nonCpG[,1]) & promoter_nonCpG)>1){
    dat_lme_1_prom_nonCpG <- rbind(dat_lme_1_prom_nonCpG,
                                cbind(logcpm_sample[rna_annot==genes_uniq_prom[i]],
                                      colMeans(methLevel(rrbsredB.M_nonCpG)[(gene_nonCpG[,1]==genes_uniq_prom[i]) & !is.na(gene_nonCpG[,1]) & promoter_nonCpG,order_col_nonCpG]),
                                      rep(genes_uniq_prom[i],length(colnames(rrbsredB.M_nonCpG))),
                                      colnames(rrbsredB.M_nonCpG)[order_col_nonCpG]))
  } else {
    dat_lme_1_prom_nonCpG <- rbind(dat_lme_1_prom_nonCpG,
                                cbind(logcpm_sample[rna_annot==genes_uniq_prom[i]],  
                                      methLevel(rrbsredB.M_nonCpG)[(gene_nonCpG[,1]==genes_uniq_prom[i]) & !is.na(gene_nonCpG[,1]) & promoter_nonCpG,order_col_nonCpG],
                                      rep(genes_uniq_prom[i],length(colnames(rrbsredB.M_nonCpG))),
                                      colnames(rrbsredB.M_nonCpG)[order_col_nonCpG]))
  }
}
dat_lme_1_prom_nonCpG <- data.frame(dat_lme_1_prom_nonCpG,stringsAsFactors = F)
colnames(dat_lme_1_prom_nonCpG) <- c("Expression","Methylation","Gene","Sample")
dat_lme_1_prom_nonCpG$Expression <- as.numeric(dat_lme_1_prom_nonCpG$Expression)  
dat_lme_1_prom_nonCpG$Methylation <- as.numeric(dat_lme_1_prom_nonCpG$Methylation) 
dat_lme_1_prom_nonCpG$Sample <- as.factor(dat_lme_1_prom_nonCpG$Sample) 
dat_lme_1_prom_nonCpG$Gene <- as.factor(dat_lme_1_prom_nonCpG$Gene) 
dim(dat_lme_1_prom_nonCpG)
head(dat_lme_1_prom_nonCpG)
str(dat_lme_1_prom_nonCpG)
write.table(dat_lme_1_prom_nonCpG,file="dat_lme_1_prom_nonCpG.txt",col.names = T,row.names = F,quote=F,sep="\t")

mixed_lmer_1_prom_nonCpG <- lmerTest::lmer(Expression ~ Methylation + (1|Sample) + (1+ Methylation|Gene), data = dat_lme_1_prom_nonCpG)
summary(mixed_lmer_1_prom_nonCpG)


# Model 1 Genic
dat_lme_1_genic_nonCpG <- NULL
genes_uniq_genic <- unique(rrbs_nonCpG_age_RNAseq_genic$gene)
print(length(genes_uniq_genic))
for(i in 1:length(genes_uniq_genic)){
  if (i %% 100 == 0){print(i)}
  if (sum((gene_nonCpG[,1]==genes_uniq_genic[i]) & !is.na(gene_nonCpG[,1]) & genic_nonCpG)>1){
    dat_lme_1_genic_nonCpG <- rbind(dat_lme_1_genic_nonCpG,
                                 cbind(logcpm_sample[rna_annot==genes_uniq_genic[i]],
                                       colMeans(methLevel(rrbsredB.M_nonCpG)[(gene_nonCpG[,1]==genes_uniq_genic[i]) & !is.na(gene_nonCpG[,1]) & genic_nonCpG,order_col_nonCpG]),
                                       rep(genes_uniq_genic[i],length(colnames(rrbsredB.M_nonCpG))),
                                       colnames(rrbsredB.M_nonCpG)[order_col_nonCpG]))
  } else {
    dat_lme_1_genic_nonCpG <- rbind(dat_lme_1_genic_nonCpG,
                                 cbind(logcpm_sample[rna_annot==genes_uniq_genic[i]],
                                       methLevel(rrbsredB.M_nonCpG)[(gene_nonCpG[,1]==genes_uniq_genic[i]) & !is.na(gene_nonCpG[,1]) & genic_nonCpG,order_col_nonCpG],
                                       rep(genes_uniq_genic[i],length(colnames(rrbsredB.M_nonCpG))),
                                       colnames(rrbsredB.M_nonCpG)[order_col_nonCpG]))
  }
}
dat_lme_1_genic_nonCpG <- data.frame(dat_lme_1_genic_nonCpG,stringsAsFactors = F)
colnames(dat_lme_1_genic_nonCpG) <- c("Expression","Methylation","Gene","Sample")
dat_lme_1_genic_nonCpG$Expression <- as.numeric(dat_lme_1_genic_nonCpG$Expression)  
dat_lme_1_genic_nonCpG$Methylation <- as.numeric(dat_lme_1_genic_nonCpG$Methylation) 
dat_lme_1_genic_nonCpG$Sample <- as.factor(dat_lme_1_genic_nonCpG$Sample) 
dat_lme_1_genic_nonCpG$Gene <- as.factor(dat_lme_1_genic_nonCpG$Gene) 
dim(dat_lme_1_genic_nonCpG)
head(dat_lme_1_genic_nonCpG)
str(dat_lme_1_genic_nonCpG)
write.table(dat_lme_1_genic_nonCpG,file="dat_lme_1_genic_nonCpG.txt",col.names = T,row.names = F,quote=F,sep="\t")

mixed_lmer_1_genic_nonCpG <- lmerTest::lmer(Expression ~ Methylation + (1|Sample) + (1+ Methylation|Gene), data = dat_lme_1_genic_nonCpG)
summary(mixed_lmer_1_genic_nonCpG)







# Model 2 Prom
dat_lme_2_prom_CpG <- NULL
genes_uniq <- unique(rrbs_age_RNAseq_prom$gene)
print(length(genes_uniq))
for(i in 1:length(genes_uniq)){
  if (i %% 100 == 0){print(i)}
  if (sum(gene_CpG[,1]==genes_uniq[i] & !is.na(gene_CpG[,1]))>1){
    for (j in 1:sum(gene_CpG[,1]==genes_uniq[i] & !is.na(gene_CpG[,1]))){
      wh <- which(gene_CpG[,1]==genes_uniq[i] & !is.na(gene_CpG[,1]))[j]
      dat_lme_2_prom_CpG <- rbind(dat_lme_2_prom_CpG,
                                  cbind(logcpm_sample[rna_annot==genes_uniq[i]],
                                        methLevel(rrbsredB.M_CpG)[wh,order_col_CpG],
                                        rep(genes_uniq[i],length(colnames(rrbsredB.M_CpG))),
                                        colnames(rrbsredB.M_CpG)[order_col_CpG],
                                        rep(rownames(rrbsredB.M_CpG)[wh],length(colnames(rrbsredB.M_CpG)))))
      
    }
  } else {
    dat_lme_2_prom_CpG <- rbind(dat_lme_2_prom_CpG,
                                cbind(logcpm_sample[rna_annot==genes_uniq[i]],
                                      methLevel(rrbsredB.M_CpG)[gene_CpG[,1]==genes_uniq[i] & !is.na(gene_CpG[,1]),order_col_CpG],
                                      rep(genes_uniq[i],length(colnames(rrbsredB.M_CpG))),
                                      colnames(rrbsredB.M_CpG)[order_col_CpG],
                                      rep(rownames(rrbsredB.M_CpG)[gene_CpG[,1]==genes_uniq[i] & !is.na(gene_CpG[,1])],length(colnames(rrbsredB.M_CpG)))))
  }
}
dat_lme_2_prom_CpG <- data.frame(dat_lme_2_prom_CpG,stringsAsFactors = F)
colnames(dat_lme_2_prom_CpG) <- c("Expression","Methylation","Gene","Sample","CpG")
dat_lme_2_prom_CpG$Expression <- as.numeric(dat_lme_2_prom_CpG$Expression)  
dat_lme_2_prom_CpG$Methylation <- as.numeric(dat_lme_2_prom_CpG$Methylation) 
dat_lme_2_prom_CpG$Sample <- as.factor(dat_lme_2_prom_CpG$Sample) 
dat_lme_2_prom_CpG$Gene <- as.factor(dat_lme_2_prom_CpG$Gene) 
dat_lme_2_prom_CpG$CpG <- as.factor(dat_lme_2_prom_CpG$CpG) 
dim(dat_lme_2_prom_CpG)
head(dat_lme_2_prom_CpG)
str(dat_lme_2_prom_CpG)

mixed_lmer_2_prom_CpG <- lmerTest::lmer(Expression ~ Methylation + (1|Sample) + (1+Methylation|Gene/CpG), data = dat_lme_2_prom_CpG)
summary(mixed_lmer_2_prom_CpG)


######################################################
# Difference in DM between DE genes and not DE genes #
######################################################

## All
genes_DE <- fulldata_output$ensembl.gene.id[(fulldata_output$Pvalue.age<0.05) & !is.na(fulldata_output$Pvalue.age)]
genes_notDE <- fulldata_output$ensembl.gene.id[!(fulldata_output$Pvalue.age<0.05) | is.na(fulldata_output$Pvalue.age)]

# genes_notDE <- read.table("Plastoscine_annot.txt",header=F,sep="\t",quote="",stringsAsFactors = F)$V1

genes_present <- unique(rrbs_age_full$gene)
genes_present <- genes_present[!is.na(genes_present)]

uniq_annot_full <- rowSums(is.na(rrbs_age_full[,grepl("gene_alt",colnames(rrbs_age_full))]))==10 

genes_DM <- unique(unlist(rrbs_age_full$gene[rrbs_age_full$adj.P.Val<0.05 & uniq_annot_full]))
genes_DM <- genes_DM[!is.na(genes_DM)]
genes_DM <- genes_DM[genes_DM%in%genes_present]

genes_notDM <-c(genes_DE,genes_notDE)[!(c(genes_DE,genes_notDE)%in%genes_DM)]
genes_notDM <- genes_notDM[genes_notDM%in%genes_present]

ct <- rbind(c(sum(genes_DE%in%genes_DM),sum(genes_DE%in%genes_notDM)),
            c(sum(genes_notDE%in%genes_DM),sum(genes_notDE%in%genes_notDM)))
chisq.test(ct)
ct
ct[,1]/(ct[,1]+ct[,2])


## Promoter
genes_present_prom <- unique(rrbs_age_full$gene[(rrbs_age_full$func=="promoter")])
genes_present_prom <- genes_present_prom[!is.na(genes_present_prom)]

genes_DM_prom <- genes_DM[genes_DM%in%genes_present_prom]
genes_notDM_prom <- genes_notDM[genes_notDM%in%genes_present_prom]

ct_prom <- rbind(c(sum(genes_DE%in%genes_DM_prom),sum(genes_DE%in%genes_notDM_prom)),
            c(sum(genes_notDE%in%genes_DM_prom),sum(genes_notDE%in%genes_notDM_prom)))
chisq.test(ct_prom)
ct_prom
ct_prom[,1]/(ct_prom[,1]+ct_prom[,2])

## Genic
genes_present_genic <- unique(rrbs_age_full$gene[ (rrbs_age_full$func=="genic")])
genes_present_genic <- genes_present_genic[!is.na(genes_present_genic)]
genes_present_genic <- genes_present_genic[!(genes_present_genic%in%genes_present_prom)]

genes_DM_genic <- genes_DM[genes_DM%in%genes_present_genic]
genes_notDM_genic <- genes_notDM[(genes_notDM%in%genes_present_genic)]

ct_genic <- rbind(c(sum(genes_DE%in%genes_DM_genic),sum(genes_DE%in%genes_notDM_genic)),
                 c(sum(genes_notDE%in%genes_DM_genic),sum(genes_notDE%in%genes_notDM_genic)))
chisq.test(ct_genic)
ct_genic
ct_genic[,1]/(ct_genic[,1]+ct_genic[,2])


##########################
## Transcription Factor ##
##########################


## Annotation
############

# annot <- read.table("MATRIX.txt",stringsAsFactors = F,fill=T)
# print(head(annot))
# species <- read.table("MATRIX_SPECIES.txt",stringsAsFactors = F,fill=T)
# print(head(species))
# species <- species[!duplicated(species$V1),]
# annot2 <- read.table("MATRIX_ANNOTATION.txt",stringsAsFactors = F,sep="\t",quote='')
# print(head(annot2))
# print(dim(annot2))
# print(str(annot2))
# class <- annot2[grepl("class",annot2$V2),c(1,3)]
# print(head(class))
# class <- class[!duplicated(class$V1),]
# family <- annot2[grepl("family",annot2$V2),c(1,3)]
# print(head(family))
# family <- family[!duplicated(family$V1),]
# medline <- annot2[grepl("medline",annot2$V2),c(1,3)]
# print(head(medline))
# medline <- medline[!duplicated(medline$V1),]
# print(sum(duplicated(annot$V1)))
# annot_new <- as.data.frame(cbind(rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot))),stringsAsFactors=F)
# colnames(annot_new) <- c("ID","M_ID","Gene_Symbol","Species","Family","Class","Medline")
# for (i in 1:nrow(annot)){
#   annot_new$ID[i] <- annot[i,1]
#   annot_new$M_ID[i] <- annot[i,3]
#   annot_new$Gene_Symbol[i] <- annot[i,5]
#   if (!identical(species$V2[species$V1==annot$V1[i]],character(0))){
#     annot_new$Species[i] <- species$V2[species$V1==annot$V1[i]]
#   }
#   if (!identical(family$V3[family$V1==annot$V1[i]],character(0))){
#     annot_new$Family[i] <- family$V3[family$V1==annot$V1[i]]
#   }
#   if (!identical(class$V3[class$V1==annot$V1[i]],character(0))){
#     annot_new$Class[i] <- class$V3[class$V1==annot$V1[i]]
#   }
#   if (!identical(medline$V3[medline$V1==annot$V1[i]],character(0))){
#     annot_new$Medline[i] <-medline$V3[medline$V1==annot$V1[i]]
#   }
# }
# annot <- annot_new

annot <- read.table("JASPAR_annot.txt",header=T,sep="\t",stringsAsFactors=F)
head(annot)

annot$ID <- unlist(strsplit(annot$ID,"\\.[0-9]"))

annot_human <- annot[annot$Species=="Homo sapiens",]
head(annot_human)

## Annotation
############

genes_TF <- fulldata_output$ensembl.gene.id[(fulldata_output$gene.symbol%in%annot_human$Name) | (fulldata_output$uniprot.annot%in%annot_human$Name)]
genes_notTF <- fulldata_output$ensembl.gene.id[!((fulldata_output$gene.symbol%in%annot_human$Name) | (fulldata_output$uniprot.annot%in%annot_human$Name))]
 
ct_TF <- rbind(c(sum(genes_TF%in%genes_DM),sum(genes_TF%in%genes_notDM)),
            c(sum(genes_notTF%in%genes_DM),sum(genes_notTF%in%genes_notDM)))
chisq.test(ct_TF)
ct_TF
ct_TF[,1]/(ct_TF[,1]+ct_TF[,2])

fulldata_output_jaspar <- fulldata_output[fulldata_output$ensembl.gene.id%in%genes_TF,]

length(unique(fulldata_output_jaspar$uniprot.annot[(fulldata_output_jaspar$ensembl.gene.id%in%genes_DM) & (fulldata_output_jaspar$ensembl.gene.id%in%genes_DE)],
            fulldata_output_jaspar$gene.symbol[(fulldata_output_jaspar$ensembl.gene.id%in%genes_DM) & (fulldata_output_jaspar$ensembl.gene.id%in%genes_DE)]))
sort(unique(fulldata_output_jaspar$uniprot.annot[fulldata_output_jaspar$ensembl.gene.id%in%genes_DM],fulldata_output_jaspar$gene.symbol[fulldata_output_jaspar$ensembl.gene.id%in%genes_DM]))


##############
## BedGraph ##
##############

makeBedGraph <- function(dataobj,cn,trackname,description,col=NULL,range=NULL,startview=NULL,isSign=NULL){
  if (!is.null(isSign)){
    if (!is.null(range)){
      range[1] <- -range[2] 
    }
    dataobj[!isSign,cn[4]] <- -dataobj[!isSign,cn[4]]
  }
  
  if (is.null(startview)){
    write.table("browser position chrZ:56,283,356-56,343,168",file=paste(trackname,".BedGraph",sep=""),col.names = F, quote = F, row.names = F)
  } else {
    write.table(paste("browser position ",startview,sep=""),file=paste(trackname,".BedGraph",sep=""),col.names = F, quote = F, row.names = F)
  }
  
  write.table("browser hide all",file=paste(trackname,".BedGraph",sep=""),col.names = F, quote = F, row.names = F, append = T)
  
  if (is.null(range)){
    header = paste('track type=bedGraph name=',trackname,' description=',description,
                   ' visibility=full color=',paste(col2rgb(col[1]),collapse=","),' altColor=',paste(col2rgb(col[2]),collapse=","),' gridDefault=on',sep="")
  } else {
    header = paste('track type=bedGraph name=',trackname,' description=',description,
                   ' visibility=full color=',paste(col2rgb(col[1]),collapse=","),' altColor=',paste(col2rgb(col[2]),collapse=","),
                   ' gridDefault=on viewLimits=',range[1],":",range[2]," autoScale=off",sep="")
  }
  
  write.table(header,file=paste(trackname,".BedGraph",sep=""),col.names = F, quote = F, row.names = F, append = T)
  
  write.table(cbind(paste("chr",dataobj[,cn[1]],sep=""),dataobj[,cn[2:4]]),file=paste(trackname,".BedGraph",sep=""),
              col.names = F, quote = F, row.names = F,sep="\t", append = T)
  
}

rrbs_age_prom_sign <- rrbs_age_RNAseq[(rrbs_age_RNAseq$func=="promoter"),]
rrbs_age_prom_sign <- rrbs_age_prom_sign[!grepl("random",rrbs_age_prom_sign$chr),]
rrbs_age_prom_sign <- rrbs_age_prom_sign[!grepl("Un",rrbs_age_prom_sign$chr),]
rrbs_age_prom_sign <- rrbs_age_prom_sign[!grepl("MT",rrbs_age_prom_sign$chr),]
makeBedGraph(rrbs_age_prom_sign,c(1,2,2,12),"RRBS_age_1dph","Meth_perc",col=c("green","red"),isSign = (rrbs_age_prom_sign$adj.P.Val<0.05), range=c(0,0.3))
makeBedGraph(rrbs_age_prom_sign,c(1,2,2,14),"RRBS_age_older","Meth_perc",col=c("blue","brown"),isSign = (rrbs_age_prom_sign$adj.P.Val<0.05), range=c(0,0.3))

# Comet plot
rrbs_age_prom_ARX <- rrbs_age_prom_sign[rrbs_age_prom_sign$gene_symbol=="ARX",]
rrbs_age_prom_ARX <- rrbs_age_prom_ARX[sort(rrbs_age_prom_ARX$locations,index.return=T)$ix,]
dim(rrbs_age_prom_sign_ARX)
head(rrbs_age_prom_ARX[,1:15])

rrbs_age_prom_sign_ARX <- rrbs_age_prom_ARX[rrbs_age_prom_ARX$adj.P.Val<0.05,]
dim(rrbs_age_prom_sign_ARX)
head(rrbs_age_prom_ARX[,1:15])

# mat_infofile <- cbind(rownames(rrbs_age_prom_ARX),rrbs_age_prom_ARX$chr,
#                       rrbs_age_prom_ARX$locations,rrbs_age_prom_ARX$adj.P.Val)
# colnames(mat_infofile) <- c("TargetID","CHR","MAPINFO","Pval")
# write.table(mat_infofile,"ARX_infofile.txt",col.names = T,row.names = F,quote = F,sep="\t")
# 
# library("coMET")
# library("rtracklayer")
# extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
# configfile <- "ARX_configfile.txt"
# myinfofile <- "ARX_infofile.txt"
# myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
# mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")
# chrom <- "chr1"
# start <- 12462436
# end <- 12466436
# gen <- "hg19"
# strand <- "*"
# BROWSER.SESSION="UCSC"
# mySession <- browserSession(BROWSER.SESSION)
# genome(mySession) <- gen
# genetrack <-genes_ENSEMBL(gen,chrom,start,end,showId=TRUE)
# # snptrack <- snpBiomart_ENSEMBL(gen,chrom, start, end,
# #                                dataset="hsapiens_snp_som",showId=FALSE,
# #                                title="SNPs Ensembl")
# cpgIstrack <- cpgIslands_UCSC(gen,chrom,start,end)
# prombedFilePath <- file.path(extdata, "/RoadMap/regions_prom_E063.bed")
# promRMtrackE063<- DNaseI_RoadMap(gen,chrom,start, end, prombedFilePath,
#                                  featureDisplay='promotor', type_stacking="squish")
# bedFilePath <- file.path(extdata, "RoadMap/E063_15_coreMarks_mnemonics.bed")
# chromHMM_RoadMapAllE063 <- chromHMM_RoadMap(gen,chrom,start, end,
#                                             bedFilePath, featureDisplay = "all",
#                                             colorcase='roadmap15' )
# listgviz <- list(genetrack,cpgIstrack,promRMtrackE063,chromHMM_RoadMapAllE063)
# comet(config.file=configfile, mydata.file=myinfofile, mydata.type="file",
#       cormatrix.file=mycorrelation, cormatrix.type="listfile",
#       mydata.large.file=myexpressfile, mydata.large.type="listfile",
#       tracks.gviz=listgviz, verbose=TRUE, print.image=FALSE)
# 
