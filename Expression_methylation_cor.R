##################################
##################################
##                              ##
## Expression_methylation_cor.R ##
##                              ##
##################################
##################################


## Load Packages
#####################

library(edgeR)
library(BiSeq)
library(biomaRt)


## Import data 
##############

# Load in RRBS raw data
load(file="allrrbs.Rda") # allrrbs.Rda can be retrieved using RRBS_analyis.R script
locs<-rowRanges(allrrbs)
allrrbs.rel <- rawToRel(allrrbs)
head(locs@ranges@start)

# Load in RNA  raw data
rnadata<-read.table("Plastoscine.txt",header=T,sep="\t",quote="")
rna_annot <- read.table("Plastoscine_annot.txt",header=F,sep="\t",quote="")
row.names(rnadata) <- unlist(rna_annot)

# Load in RRBS annotation
annot<-c("Pool","Sex","Age","8","F","20","12","M","20","14","F","65","18","M","65","1","F","1","7","F","20","10","M","20","24","M","adult","2","F","1","15","F","65","19","F","adult","22","M","adult","9","F","20","11","M","20","16","M","65","17","M","65","4","M","1","5","M","1","13","F","65","20","F","adult","3","F","1","6","M","1","21","F","adult","23","M","adult")
annot<-matrix(annot,nrow=25,ncol=3,byrow=T)
colnames(annot)<-annot[1,]
annot<-annot[-1,]
rownames(annot)<-annot[,1]
annot <- data.frame(annot)
annot$Pool <- as.numeric(as.character(annot$Pool))

# Get order of RRBS data
ix <- sort(annot$Pool,index.return=T)$ix
annot[ix,]

# Get significant differentially methylated loci for the contrast age
rrbs <- read.table(file="RRBS_age_RNAseq.txt",sep="\t",header=T)
rrbs$gene <- as.character(rrbs$gene)
print(dim(rrbs))
rrbs <- rrbs[rrbs$adj.P.Val<0.05,]
print(dim(rrbs))
rrbs <- rrbs[as.character(rrbs$func)=="promoter",]
print(dim(rrbs))
print(str(rrbs))
print(head(rrbs))


## Get methylation percentages for significant differentially methylated loci 
##############

locs_df <- as.data.frame(locs)
locs_pos <- paste(locs_df$seqnames,locs_df$start,sep="_")
wh <- NULL
teller <- 1
for (i in 1:nrow(rrbs)){
  wh <- c(wh,which(locs_pos%in%paste(rrbs$chr[i],rrbs$locations[i],sep="_")))
  if (i >= 50*teller){
    print(paste(i,"done out of",dim(rrbs)[1]))
    teller <- teller + 1
  }
}

rrbs_raw <- methLevel(allrrbs.rel[wh,])
ix <- sort(as.numeric(unlist(strsplit(colnames(rrbs_raw),"[BO]"))[unlist(strsplit(colnames(rrbs_raw),"[BO]"))!=""]),index.return=T)$ix
rrbs_raw <- rrbs_raw[,ix]


## Get logcpm values
##############

y <-DGEList(rnadata)
y <- calcNormFactors(y)
logcpm <- cpm(y,log = T)


## Get raw (expression) counts for significant differentially methylated loci
##############

logcpm_sub <- NULL
for (i in 1:length(rrbs$gene)){
  # print(rnadata[rownames(rnadata)==rrbs$gene[i],])
  logcpm_sub <- rbind(logcpm_sub ,logcpm[rownames(logcpm)==rrbs$gene[i],])
}


## Correlation methylation-expression in function of position in promotor
########################################################################

# Get TSS and extra annotation
if (!file.exists("Genes.Rda")){
  ensembl=useMart("ensembl")
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

# Get order RRBS 
ix <- sort(as.integer(unlist(strsplit(colnames(rrbs_raw),"[A-Z]"))[unlist(strsplit(colnames(rrbs_raw),"[A-Z]"))!=""]),index.return=T)$ix

distanceFromTSS <- function(t){
  if(genes$start_position[genes$ensembl_gene_id==rrbs$gene[t]]==genes$transcription_start_site[genes$ensembl_gene_id==rrbs$gene[t]]){
    allTSS <- genes$transcription_start_site[genes$ensembl_gene_id%in%rrbs[t,grep("gene",colnames(rrbs))][!is.na(rrbs[t,grep("gene",colnames(rrbs))])]]
    wh <- which(abs(rrbs$locations[t]-allTSS)==min(abs(rrbs$locations[t]-allTSS)))
    return(rrbs$locations[t]-allTSS[wh])
  } else {
    allTSS <- genes$transcription_start_site[genes$ensembl_gene_id%in%rrbs[t,grep("gene",colnames(rrbs))][!is.na(rrbs[t,grep("gene",colnames(rrbs))])]]
    wh <- which(abs(-(rrbs$locations[t]-allTSS))==min(abs(-(rrbs$locations[t]-allTSS))))
    return(-(rrbs$locations[t]-allTSS[wh]))
  }
}

# Get methylation percentages and expression counts
pos <- NULL
meth_perc <- NULL
express <- NULL
teller <- 0
for(t in 1:length(rrbs$gene_symbol)){
  pos <- c(pos,distanceFromTSS(t))
}

# Calculate correlation per position in the promotor
pos_plot_raw <- NULL
cor_20 <- NULL
cor_65 <- NULL
cor_adult <- NULL
cor_all <- NULL
cor_raw <- NULL
pvals <- NULL
for(i in -2000:800){
  if (sum(pos==i)!=0){
    pos_plot_raw <- c(pos_plot_raw,i)
    # cor_raw <- c(cor_raw,cor.test(c(rrbs_raw[pos==i,]),c(logcpm_sub[pos==i,]),method="spearman")$estimate)
    # pvals <- c(pvals,cor.test(c(rrbs_raw[pos==i,]),c(logcpm_sub[pos==i,]),method="spearman")$p.value)
    cor_raw <- c(cor_raw,cor.test(c(rrbs_raw[pos==i,]),c(logcpm_sub[pos==i,]))$estimate)
    pvals <- c(pvals,cor.test(c(rrbs_raw[pos==i,]),c(logcpm_sub[pos==i,]))$p.value)
  }

}

pos_plot_raw <- pos_plot_raw[pvals*length(pvals)<0.05]
cor_raw <- cor_raw[pvals*length(pvals)<0.05]


## Plot
###########

jpeg(paste('Correlation_pos_raw.jpg',sep=''),width = 600, height = 600, units = "px")
pos_plot_raw <- pos_plot_raw[!is.na(cor_raw)]
cor_raw <- cor_raw[!is.na(cor_raw)]
df_raw <- as.data.frame(cbind(pos_plot_raw,cor_raw))
colnames(df_raw) <- c("pos","cor")
plot(df_raw$pos,df_raw$cor,main="Raw",ylab=c("R"),xlab=c("Pos relative to TSS"),pch=20)
loc<-loess(cor~pos,df_raw,span=0.55)
sm_loc <- predict(loc) 
lines(df_raw$pos,sm_loc,col="red")
abline(h=0)
dev.off()

####
## CH*

# Load in RRBS raw data
load(file="allrrbs.Rda") # allrrbs.Rda can be retrieved using RRBS_analyis.R script
locs<-rowRanges(allrrbs)
allrrbs.rel <- rawToRel(allrrbs)
head(locs@ranges@start)


# Get significant differentially methylated loci for the contrast age
rrbs_nonCpG <- read.table(file="RRBS_age_RNAseq_nonCpG.txt",sep="\t",header=T)
rrbs_nonCpG$gene <- as.character(rrbs_nonCpG$gene)
print(dim(rrbs_nonCpG))
rrbs_nonCpG <- rrbs_nonCpG[rrbs_nonCpG$adj.P.Val<0.05,]
print(dim(rrbs_nonCpG))
rrbs_nonCpG <- rrbs_nonCpG[as.character(rrbs_nonCpG$func)=="promoter",]
print(dim(rrbs_nonCpG))
print(str(rrbs_nonCpG))
print(head(rrbs_nonCpG))


## Get methylation percentages for significant differentially methylated loci 
##############

locs_df <- as.data.frame(locs)
locs_pos <- paste(locs_df$seqnames,locs_df$start,sep="_")
wh <- NULL
teller <- 1
for (i in 1:nrow(rrbs_nonCpG)){
  wh <- c(wh,which(locs_pos%in%paste(rrbs_nonCpG$chr[i],rrbs_nonCpG$locations[i],sep="_")))
  if (i >= 50*teller){
    print(paste(i,"done out of",dim(rrbs_nonCpG)[1]))
    teller <- teller + 1
  }
}

rrbs_nonCpG_raw <- methLevel(allrrbs_nonCpG.rel[wh,])
ix <- sort(as.numeric(unlist(strsplit(colnames(rrbs_nonCpG_raw),"[BO]"))[unlist(strsplit(colnames(rrbs_nonCpG_raw),"[BO]"))!=""]),index.return=T)$ix
rrbs_nonCpG_raw <- rrbs_nonCpG_raw[,ix]


## Get logcpm values
##############

y <-DGEList(rnadata)
y <- calcNormFactors(y)
logcpm <- cpm(y,log = T)


## Get raw (expression) counts for significant differentially methylated loci
##############

logcpm_sub <- NULL
for (i in 1:length(rrbs_nonCpG$gene)){
  # print(rnadata[rownames(rnadata)==rrbs_nonCpG$gene[i],])
  logcpm_sub <- rbind(logcpm_sub ,logcpm[rownames(logcpm)==rrbs_nonCpG$gene[i],])
}


## Correlation methylation-expression in function of position in promotor
########################################################################

# Get TSS and extra annotation
if (!file.exists("Genes.Rda")){
  ensembl=useMart("ensembl")
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

# Get order rrbs_nonCpG 
ix <- sort(as.integer(unlist(strsplit(colnames(rrbs_nonCpG_raw),"[A-Z]"))[unlist(strsplit(colnames(rrbs_nonCpG_raw),"[A-Z]"))!=""]),index.return=T)$ix

distanceFromTSS <- function(t){
  if(genes$start_position[genes$ensembl_gene_id==rrbs_nonCpG$gene[t]]==genes$transcription_start_site[genes$ensembl_gene_id==rrbs_nonCpG$gene[t]]){
    allTSS <- genes$transcription_start_site[genes$ensembl_gene_id%in%rrbs_nonCpG[t,grep("gene",colnames(rrbs_nonCpG))][!is.na(rrbs_nonCpG[t,grep("gene",colnames(rrbs_nonCpG))])]]
    wh <- which(abs(rrbs_nonCpG$locations[t]-allTSS)==min(abs(rrbs_nonCpG$locations[t]-allTSS)))
    return(rrbs_nonCpG$locations[t]-allTSS[wh])
  } else {
    allTSS <- genes$transcription_start_site[genes$ensembl_gene_id%in%rrbs_nonCpG[t,grep("gene",colnames(rrbs_nonCpG))][!is.na(rrbs_nonCpG[t,grep("gene",colnames(rrbs_nonCpG))])]]
    wh <- which(abs(-(rrbs_nonCpG$locations[t]-allTSS))==min(abs(-(rrbs_nonCpG$locations[t]-allTSS))))
    return(-(rrbs_nonCpG$locations[t]-allTSS[wh]))
  }
}

# Get methylation percentages and expression counts
pos <- NULL
meth_perc <- NULL
express <- NULL
teller <- 0
for(t in 1:length(rrbs_nonCpG$gene_symbol)){
  pos <- c(pos,distanceFromTSS(t))
}

# Calculate correlation per position in the promotor
pos_plot_raw <- NULL
cor_20 <- NULL
cor_65 <- NULL
cor_adult <- NULL
cor_all <- NULL
cor_raw <- NULL
pvals <- NULL
for(i in -2000:800){
  if (sum(pos==i)!=0){
    pos_plot_raw <- c(pos_plot_raw,i)
    # cor_raw <- c(cor_raw,cor.test(c(rrbs_nonCpG_raw[pos==i,]),c(logcpm_sub[pos==i,]),method="spearman")$estimate)
    # pvals <- c(pvals,cor.test(c(rrbs_nonCpG_raw[pos==i,]),c(logcpm_sub[pos==i,]),method="spearman")$p.value)
    cor_raw <- c(cor_raw,cor.test(c(rrbs_nonCpG_raw[pos==i,]),c(logcpm_sub[pos==i,]))$estimate)
    pvals <- c(pvals,cor.test(c(rrbs_nonCpG_raw[pos==i,]),c(logcpm_sub[pos==i,]))$p.value)
  }
  
}

pos_plot_raw <- pos_plot_raw[pvals*length(pvals)<0.05]
cor_raw <- cor_raw[pvals*length(pvals)<0.05]


## Plot
###########

jpeg(paste('Correlation_pos_raw.jpg',sep=''),width = 600, height = 600, units = "px")
pos_plot_raw <- pos_plot_raw[!is.na(cor_raw)]
cor_raw <- cor_raw[!is.na(cor_raw)]
df_raw <- as.data.frame(cbind(pos_plot_raw,cor_raw))
colnames(df_raw) <- c("pos","cor")
plot(df_raw$pos,df_raw$cor,main="Raw",ylab=c("R"),xlab=c("Pos relative to TSS"),pch=20)
loc<-loess(cor~pos,df_raw,span=0.55)
sm_loc <- predict(loc) 
lines(df_raw$pos,sm_loc,col="red")
abline(h=0)
dev.off()

