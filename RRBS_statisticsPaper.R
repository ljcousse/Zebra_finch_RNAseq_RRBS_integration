#############################
#############################
##                         ##
## RRBS_statistics_Paper.R ##
##                         ##
#############################
#############################

library("edgeR")
library("BiSeq")


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

rrbs<-filterBySharedRegions(object=allrrbs,perc.samples=1,minCov=mincov)
rrbs.rel <- rawToRel(rrbs)

locs<-rowRanges(rrbs)

ensembl=useMart("ensembl")
ensembl = useDataset("tguttata_gene_ensembl",mart=ensembl)
genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','transcription_start_site'),mart = ensembl)
colnames(genes) <- c("gene","gene_symbol","chr","start","end","TSS")
dim(genes) 
genes <- genes[!duplicated(genes$gene),]
dim(genes)
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
overlap_p <- subsetByOverlaps(locs,prom_gr)
promoter <- rownames(as.data.frame(locs))%in%rownames(as.data.frame(overlap_p))

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
overlap_g <- subsetByOverlaps(locs,genic_gr)
genic <- rownames(as.data.frame(locs))%in%rownames(as.data.frame(overlap_g))


genic[genic & promoter] <- FALSE

overlaps_genic <- findOverlaps(locs,genic_gr)
overlaps_prom <- findOverlaps(locs,prom_gr)

all_overlap <- c(as.character(genic_gr[subjectHits(overlaps_genic)]$gene),as.character(prom_gr[subjectHits(overlaps_prom)]$gene))
all_overlap <- cbind(all_overlap,c(as.character(genic_gr[subjectHits(overlaps_genic)]$gene_symbol),as.character(prom_gr[subjectHits(overlaps_prom)]$gene_symbol)))
all_overlap <- cbind(all_overlap,c(locs[queryHits(overlaps_genic)]@ranges@NAMES,locs[queryHits(overlaps_prom)]@ranges@NAMES))
colnames(all_overlap) <- c("gene","gene_symbol","rowname")


rn <- as.matrix(locs@ranges@NAMES,ncol=1)
colnames(rn) <- "rowname"

merger <- merge(rn,all_overlap,all=T)
merger <- merger[sort(as.numeric(as.character(merger$rowname)),index.return=T)$ix,]

merger$gene <- as.character(merger$gene)
merger$gene_symbol <- as.character(merger$gene_symbol)

gene <- as.character(merger$gene[!duplicated(merger$rowname)])
gene_symbol <- as.character(merger$gene_symbol[!duplicated(merger$rowname)])

max_n <- max(table(merger$rowname))
print(max_n)
gene <- cbind(gene,matrix(rep(NA,length(gene)*(max_n-1)),ncol=max_n-1))
gene_symbol <- cbind(gene_symbol,matrix(rep(NA,length(gene_symbol)*(max_n-1)),ncol=max_n-1))
for (i in 2:max_n){
  n_times <- merger[merger$rowname%in%names(table(merger$rowname)[table(merger$rowname)==i]),]
  if(nrow(n_times)==0){
    print(paste("no locus with ",i,"-fold annotation",sep=""))
    next
  }
  if (i==2 | nrow(n_times)==i){
    gene[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]
    gene_symbol[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- n_times$gene_symbol[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]
  } else {
    gene[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- t(matrix(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]],nrow=i-1))
    gene_symbol[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- t(matrix(n_times$gene_symbol[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]],nrow=i-1))
  }
  
}

minmethreads <- 6

inclocsM2 <- rowSums(methReads(rrbs))>=minmethreads

locs <- locs[inclocsM2]
rrbsredB <- rrbs[inclocsM2,]

genic <- genic[inclocsM2]
promoter <- promoter[inclocsM2]

gene <- gene[inclocsM2,]
gene_symbol <- gene_symbol[inclocsM2,]

rrbsredB_meth <- methLevel(rawToRel(rrbsredB))

# Correlation
uniq_annot_RNAseq <- rowSums(is.na(rrbs_age_RNAseq[,grepl("gene_alt",colnames(rrbs_age_RNAseq))]))==10 

rrbs_age_RNAseq_uniq <- rrbs_age_RNAseq[uniq_annot_RNAseq,]
rrbs_age_RNAseq_prom <- rrbs_age_RNAseq[uniq_annot_RNAseq & rrbs_age_RNAseq$func=="promoter",]
rrbs_age_RNAseq_genic <- rrbs_age_RNAseq[uniq_annot_RNAseq & rrbs_age_RNAseq$func=="genic",]

colnumbers <- which(grepl("gene",colnames(rrbs_age_RNAseq)) & !grepl("gene_symbol",colnames(rrbs_age_RNAseq)))

rnaseq_logcpm <- NULL
rrbs_meth <- NULL
genes_uniq <- unique(rrbs_age_RNAseq_uniq$gene)
for(i in 1:length(genes_uniq)){
  if (i %% 100 == 0){print(i)}
  rnaseq_logcpm <- c(rnaseq_logcpm,logcpm[rna_annot==genes_uniq[i]])
  rrbs_meth <- c(rrbs_meth,mean(rrbsredB_meth[gene[,1]==genes_uniq[i] & !is.na(gene[,1]),]))
}
cor.test(rrbs_meth,rnaseq_logcpm)

rnaseq_logcpm_prom <- NULL
rrbs_meth_prom <- NULL
genes_uniq_prom <- unique(rrbs_age_RNAseq_prom$gene)
for(i in 1:length(genes_uniq_prom)){
  if (i %% 100 == 0){print(i)}
  rnaseq_logcpm_prom  <- c(rnaseq_logcpm_prom,logcpm[rna_annot==genes_uniq_prom[i]])
  rrbs_meth_prom  <- c(rrbs_meth_prom,mean(rrbsredB_meth[gene[,1]==genes_uniq_prom[i] & !is.na(gene[,1]),]))
}
cor.test(rrbs_meth_prom,rnaseq_logcpm_prom)

rnaseq_logcpm_genic <- NULL
rrbs_meth_genic <- NULL
genes_uniq_genic <- unique(rrbs_age_RNAseq_genic$gene)
for(i in 1:length(genes_uniq_genic)){
  if (i %% 100 == 0){print(i)}
  rnaseq_logcpm_genic  <- c(rnaseq_logcpm_genic,logcpm[rna_annot==genes_uniq_genic[i]])
  rrbs_meth_genic  <- c(rrbs_meth_genic,mean(rrbsredB_meth[gene[,1]==genes_uniq_genic[i] & !is.na(gene[,1]),]))
}
cor.test(rrbs_meth_genic,rnaseq_logcpm_genic)


## nonCpG

# RRBS

# Load in RRBS annotation
load("rrbs_nonCpG.Rda")
rownames(rrbs) <- 1:nrow(rrbs) 
rrbs.rel <- rawToRel(rrbs)

locs<-rowRanges(rrbs)

ensembl=useMart("ensembl")
ensembl = useDataset("tguttata_gene_ensembl",mart=ensembl)
genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','transcription_start_site'),mart = ensembl)
colnames(genes) <- c("gene","gene_symbol","chr","start","end","TSS")
dim(genes) 
genes <- genes[!duplicated(genes$gene),]
dim(genes)
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
overlap_p <- subsetByOverlaps(locs,prom_gr)
promoter <- rownames(as.data.frame(locs))%in%rownames(as.data.frame(overlap_p))

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
overlap_g <- subsetByOverlaps(locs,genic_gr)
genic <- rownames(as.data.frame(locs))%in%rownames(as.data.frame(overlap_g))


genic[genic & promoter] <- FALSE

overlaps_genic <- findOverlaps(locs,genic_gr)
overlaps_prom <- findOverlaps(locs,prom_gr)

all_overlap <- c(as.character(genic_gr[subjectHits(overlaps_genic)]$gene),as.character(prom_gr[subjectHits(overlaps_prom)]$gene))
all_overlap <- cbind(all_overlap,c(as.character(genic_gr[subjectHits(overlaps_genic)]$gene_symbol),as.character(prom_gr[subjectHits(overlaps_prom)]$gene_symbol)))
all_overlap <- cbind(all_overlap,c(locs[queryHits(overlaps_genic)]@ranges@NAMES,locs[queryHits(overlaps_prom)]@ranges@NAMES))
colnames(all_overlap) <- c("gene","gene_symbol","rowname")


rn <- as.matrix(locs@ranges@NAMES,ncol=1)
colnames(rn) <- "rowname"

merger <- merge(rn,all_overlap,all=T)
merger <- merger[sort(as.numeric(as.character(merger$rowname)),index.return=T)$ix,]

merger$gene <- as.character(merger$gene)
merger$gene_symbol <- as.character(merger$gene_symbol)

gene <- as.character(merger$gene[!duplicated(merger$rowname)])
gene_symbol <- as.character(merger$gene_symbol[!duplicated(merger$rowname)])

max_n <- max(table(merger$rowname))
print(max_n)
gene <- cbind(gene,matrix(rep(NA,length(gene)*(max_n-1)),ncol=max_n-1))
gene_symbol <- cbind(gene_symbol,matrix(rep(NA,length(gene_symbol)*(max_n-1)),ncol=max_n-1))
for (i in 2:max_n){
  n_times <- merger[merger$rowname%in%names(table(merger$rowname)[table(merger$rowname)==i]),]
  if(nrow(n_times)==0){
    print(paste("no locus with ",i,"-fold annotation",sep=""))
    next
  }
  if (i==2 | nrow(n_times)==i){
    gene[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]
    gene_symbol[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- n_times$gene_symbol[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]
  } else {
    gene[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- t(matrix(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]],nrow=i-1))
    gene_symbol[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- t(matrix(n_times$gene_symbol[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]],nrow=i-1))
  }
  
}

minmethreads <- 6

inclocsM2 <- rowSums(methReads(rrbs))>=minmethreads

locs <- locs[inclocsM2]
rrbsredB <- rrbs[inclocsM2,]

genic <- genic[inclocsM2]
promoter <- promoter[inclocsM2]

gene <- gene[inclocsM2,]
gene_symbol <- gene_symbol[inclocsM2,]

rrbsredB_meth <- methLevel(rawToRel(rrbsredB))

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
  rrbs_meth_nonCpG  <- c(rrbs_meth_nonCpG,mean(rrbsredB_meth[gene[,1]==genes_uniq_nonCpG[i] & !is.na(gene[,1]),]))
}
cor.test(rrbs_meth_nonCpG,rnaseq_logcpm_nonCpG)


rnaseq_logcpm_nonCpG_prom <- NULL
rrbs_meth_nonCpG_prom <- NULL
genes_uniq_nonCpG_prom <- unique(rrbs_nonCpG_age_RNAseq_prom$gene)
for(i in 1:length(genes_uniq_nonCpG_prom)){
  if (i %% 100 == 0){print(i)}
  rnaseq_logcpm_nonCpG_prom  <- c(rnaseq_logcpm_nonCpG_prom,logcpm[rna_annot==genes_uniq_nonCpG_prom[i]])
  rrbs_meth_nonCpG_prom  <- c(rrbs_meth_nonCpG_prom,mean(rrbsredB_meth[gene[,1]==genes_uniq_nonCpG_prom[i] & !is.na(gene[,1]),]))
}
cor.test(rrbs_meth_nonCpG_prom,rnaseq_logcpm_nonCpG_prom)


rnaseq_logcpm_nonCpG_genic <- NULL
rrbs_meth_nonCpG_genic <- NULL
genes_uniq_nonCpG_genic <- unique(rrbs_nonCpG_age_RNAseq_genic$gene)
for(i in 1:length(genes_uniq_nonCpG_genic)){
  if (i %% 100 == 0){print(i)}
  rnaseq_logcpm_nonCpG_genic  <- c(rnaseq_logcpm_nonCpG_genic,logcpm[rna_annot==genes_uniq_nonCpG_genic[i]])
  rrbs_meth_nonCpG_genic  <- c(rrbs_meth_nonCpG_genic,mean(rrbsredB_meth[gene[,1]==genes_uniq_nonCpG_genic[i] & !is.na(gene[,1]),]))
}
cor.test(rrbs_meth_nonCpG_genic,rnaseq_logcpm_nonCpG_genic)



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

## Annotation
############

genes_TF <- fulldata_output$ensembl.gene.id[(fulldata_output$gene.symbol%in%annot_human$Gene_Symbol) | (fulldata_output$uniprot.annot%in%annot_human$Gene_Symbol)]
genes_notTF <- fulldata_output$ensembl.gene.id[!((fulldata_output$gene.symbol%in%annot_human$Gene_Symbol) | (fulldata_output$uniprot.annot%in%annot_human$Gene_Symbol))]
 
ct_TF <- rbind(c(sum(genes_TF%in%genes_DM),sum(genes_TF%in%genes_notDM)),
            c(sum(genes_notTF%in%genes_DM),sum(genes_notTF%in%genes_notDM)))
chisq.test(ct_TF)
ct_TF
ct_TF[,1]/(ct_TF[,1]+ct_TF[,2])


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


