####################
####################
##                ##
## RRBS_analyse.R ##
##                ##
####################
####################


## Parameter settings 
#####################

home_dir <- getwd()
data_dir <- "/data/sandras/BiSeq/"
options(warn=2L)
mincov<-5 
minmethreads<-6
ct<-0.01 #constant value for calculation M-values from methylation degrees


## Load Packages
#####################

library(BiSeq)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(limma)
library(biomaRt)
library(UniProt.ws)
library(RColorBrewer)
library(gplots)

sink(file='RRBS_analyse.log', type = c("output", "message"))


## Extra annotation Uniprot 
###########################

setwd(home_dir)

if (file.exists("uniprot.Rda")){
  load(file="uniprot.Rda")
} else {
  up <- UniProt.ws(taxId=59729)
  # availableUniprotSpecies(pattern="guttata")
  egs = keys(up, "ENSEMBL")
  columns <- c("ENSEMBL","GENES")
  kt <- "ENSEMBL"
  uniprot <- NULL
  for (i in 1:ceiling(length(egs)/300)){
    if (i != ceiling(length(egs)/300)){
      uniprot <- rbind(uniprot,select(up, egs[((i-1)*300+1):(i*300)], columns, kt))
    } else {
      uniprot <- rbind(uniprot,select(up, egs[((i-1)*300+1):length(egs)], columns, kt))
    }
    print(paste(i*300,"done out of",length(egs)))
  }
  uniprot[is.na(uniprot)] <- ""
  save(uniprot,file="uniprot.Rda")
}


## Import data 
##############

setwd(data_dir)

# annotation
annot<-c("Pool","Sex","Age","8","F","20","12","M","20","14","F","65","18","M","65","1","F","1","7","F","20","10","M","20","24","M","adult","2","F","1","15","F","65","19","F","adult","22","M","adult","9","F","20","11","M","20","16","M","65","17","M","65","4","M","1","5","M","1","13","F","65","20","F","adult","3","F","1","6","M","1","21","F","adult","23","M","adult")
annot<-matrix(annot,nrow=25,ncol=3,byrow=T)
colnames(annot)<-annot[1,]
annot<-annot[-1,]
rownames(annot)<-annot[,1]

poolfiles<-dir()
poolfiles<-poolfiles[grepl("\\.B",poolfiles)]
# Mix up: labels switched
poolfiles[grepl("\\.B16",poolfiles)] <- dir()[grepl("\\.O16",dir())]
pids<-unlist(strsplit(substr(poolfiles,29,31),split="_"))
annot<-cbind(annot[substr(pids,2,3),],substr(pids,1,1))
annot<-cbind(annot,apply(cbind(annot[,"Sex"],annot[,"Age"]),1,paste,collapse="_"))
rownames(annot)<-pids
colnames(annot)[4:5]<-c("Ox_vs_RRBS","Fenolabel")

annotN<-DataFrame(annot)

# read RRBS data given colData object
allrrbs<-readBismark(poolfiles,colData=annotN)
setwd(home_dir)
save(allrrbs,file="allrrbs.Rda")

locs<-rowRanges(allrrbs)

setwd(home_dir)


## Data filtering 
##############

#Filtering based on shared regions (should be shared with coverage >= mincov in all samples)
rrbs<-filterBySharedRegions(object=allrrbs,perc.samples=1,minCov=mincov)
rrbs.rel <- rawToRel(rrbs)

print("After filtering on coverage")
print(dim(rrbs))


## Functional annotation per CpG
##############

locs<-rowRanges(rrbs)

if (file.exists("Genes_RRBS.Rda")){
  load("Genes_RRBS.Rda")
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

genelist_allgenes_RRBS <- unique(unlist(gene))[-1]
head(genelist_allgenes_RRBS)
length(genelist_allgenes_RRBS)
write.table(genelist_allgenes_RRBS,"GeneLists/ListGenes_allRRBS.txt",col.names = F,row.names = F,quote = F, sep="\n")

## Boxplot Methylation percentage
################################

f1 <- paste0(annot[,3],"dph")
f2 <- annot[,2]
f1 <- as.factor(f1)
f2 <- as.factor(f2)

data_lm <- NULL
for(i in 1:ncol(methLevel(rrbs.rel))){
  data_lm <- rbind(data_lm,cbind(mean(methLevel(rrbs.rel)[,i]),f1[i],f2[i]))
  
}
data_lm <- as.data.frame(data_lm,stringsAsFactors=F)
colnames(data_lm) <- c("y","f1","f2")
data_lm$y <- as.numeric(data_lm$y)
data_lm$f1 <- as.factor(data_lm$f1)
data_lm$f2 <- as.factor(data_lm$f2)

data_boxplot <- data_lm
print(data_boxplot)

labels <- c("1dph,F","1dph,M","20dph,F","20dph,M","65dph,F","65dph,M","adult,F","adult,M")
labels <- labels[c(1,3,5,7,2,4,6,8)]
print(labels)
print(levels(f1))
print(levels(f2))

print(paste(rep(as.character(levels(f1)),2),rep(as.character(levels(f2)),4),sep=","))
jpeg("Methylation_plot.jpg")
boxplot(y~f1+f2,data_boxplot,at=c(1,3,5,7,2,4,6,8),las=2,
        names=labels[])
dev.off()


## Linear model differences methylation percentage
#################################################

lm_diffmeth <- lm(y~f1*f2,data=data_lm)

print(lm_diffmeth)
print(summary(lm_diffmeth))

anova_diffmeth <- aov(y~f1,data=data_lm)

print(anova_diffmeth)
print(summary(anova_diffmeth))
print(TukeyHSD(anova_diffmeth))


## Density plot
##########

jpeg("Density_Methperc.jpg",width=1000,height=1000)
max_y <- 0
for (i in 1:ncol(rrbs.rel)){
  max_y <- max(max_y,density(methLevel(rrbs.rel)[,i])$y)
}
for (i in 1:ncol(rrbs.rel)){
  if (i==1){
    plot(density(methLevel(rrbs.rel)[,i]),col=i,ylim=c(0,max_y))
  } else {
    lines(density(methLevel(rrbs.rel)[,i]),col=i)
  }
}
legend("topright",colnames(rrbs),pch=19,col=1:ncol(rrbs.rel))
dev.off()


## data filtering 2 and processing
################################

inclocsM2 <- rowSums(methReads(rrbs))>=minmethreads

locs <- locs[inclocsM2]
rrbsredB <- rrbs[inclocsM2,]

genic <- genic[inclocsM2]
promoter <- promoter[inclocsM2]

gene <- gene[inclocsM2,]
gene_symbol <- gene_symbol[inclocsM2,]

rrbsredB.rel<-rawToRel(rrbsredB)
rrbsredB.M<-rrbsredB.rel 


for (i in 1:ncol(rrbsredB.rel)){
  ML<-methLevel(rrbsredB.rel)[,i]
  methLevel(rrbsredB.M)[,i]<-log2((ML+ct)/(1-ML+ct))
}

varB <- apply(methLevel(rrbsredB.M),1,var)
meanB <- apply(methLevel(rrbsredB.M),1,mean)
l <- lowess(meanB,varB, f = 0.5)
jpeg('VarianceTrend.jpg')
plot(meanB,varB)
lines(l, col = "red")

jpeg("Density_Methperc_Mval.jpg",width=1000,height=1000)
max_y <- 0
for (i in 1:ncol(rrbs.rel)){
  max_y <- max(max_y,density(methLevel(rrbsredB.M)[,i])$y)
}
for (i in 1:ncol(rrbs.rel)){
 if (i==1){
    plot(density(methLevel(rrbsredB.M)[,i]),col=i,ylim=c(0,max_y))
  } else {
    lines(density(methLevel(rrbsredB.M)[,i]),col=i)
  }
}
legend("topright",colnames(rrbs),pch=19,col=1:ncol(rrbs.rel))
dev.off()


## Heatmap 
##########

select = order(apply(methLevel(rrbsredB.M),1,var),decreasing=TRUE)[1:10000]

hmcol = colorRampPalette(brewer.pal(9, "RdBu"))(100)

print(head(methLevel(rrbsredB.M)))
m_matrix <- matrix(as.matrix(methLevel(rrbsredB.M)),ncol=24)[select,]
colnames(m_matrix) <- paste(f1,f2,sep=",")
print(paste(f1,f2,sep=","))
print(annot)
print(head(m_matrix))
print(dim(m_matrix))
print(str(m_matrix))
cor_t <- 1- abs(cor(t(m_matrix)))
cor <- 1- abs(cor(m_matrix))
distancet <- as.dist(cor_t)
distance <- as.dist(cor)
hclust_complete <- hclust(distancet, method = "complete")
dend_tcomplete <- as.dendrogram(hclust_complete)
hclus_complete <- hclust(distance, method = "complete")
dend_complete <- as.dendrogram(hclus_complete)
jpeg("Heatmap_methylation.jpg")
heatmap.2(m_matrix, col = hmcol,Rowv=dend_tcomplete, Colv=dend_complete, trace="non", margin=c(5,1),labRow = "")
dev.off()

jpeg("Clustering_RRBS.jpg")
par(xpd=TRUE)
par(mar=c(2.5,2.5,2.5,11))
Bclust<-plot(hclus_complete)
dev.off()


## Statistical testing 
########################

f1 <- paste0(annot[,3],"dph")
f2 <- annot[,2]
f1 <- as.factor(f1)
f2 <- as.factor(f2)

c_names <- rep("",length(f1))
for (i in 1:length(f1)){
  c_names[i] <- paste(f1[i],f2[i],i,sep=".")
}
colnames(rrbsredB.rel) <- c_names

GOI_1dph <- read.table(file="GeneLists/ListGenes_1dph.txt",header=F,stringsAsFactors =F)$V1
GOI_20dph <- read.table(file="GeneLists/ListGenes_20dph.txt",header=F,stringsAsFactors =F)$V1
GOI_65dph <- read.table(file="GeneLists/ListGenes_65dph.txt",header=F,stringsAsFactors =F)$V1
GOI_adult <- read.table(file="GeneLists/ListGenes_adult.txt",header=F,stringsAsFactors =F)$V1
GOI_age <- read.table(file="GeneLists/ListGenes_age.txt",header=F,stringsAsFactors =F)$V1


# General #

design<-model.matrix(~f1+f2+f1:f2)

L<- matrix(0,ncol=5,nrow=8)
rownames(L)<-colnames(design)
colnames(L)<-c("1dph","20dph","65dph","adult","age")
# 1dph 
L[5,1] <- 1
# 20dph 
L[5,2] <- 1
L[6,2] <- 1
# 65dph
L[5,3] <- 1
L[7,3] <- 1
# adult
L[5,4] <- 1
L[8,4] <- 1
# age
L[2,5] <- 1/3
L[3,5] <- 1/3
L[4,5] <- 1/3
L[6,5] <- 1/6
L[7,5] <- 1/6
L[8,5] <- 1/6

print(L)

fit <- lmFit(methLevel(rrbsredB.M),design)


# 1dph #

lmfit.cont <- contrasts.fit(fit, L[,colnames(L)=="1dph"])
lmfit.cont.ebayes <- eBayes(lmfit.cont)

Bres_1dph <- topTable(lmfit.cont.ebayes ,n=nrow(methLevel(rrbsredB.M)))

Bres_1dph <- Bres_1dph[sort(as.numeric(rownames(Bres_1dph)),index.return=T)$ix,]

B_obs_stat_1dph <- cbind(apply(methLevel(rrbsredB.rel)[,f1=="1dph" & f2=="F"],1,mean),
                    apply(methLevel(rrbsredB.rel)[,f1=="1dph" & f2=="F"],1,sd),
                    apply(methLevel(rrbsredB.rel)[,f1=="1dph" & f2=="M"],1,mean),
                    apply(methLevel(rrbsredB.rel)[,f1=="1dph" & f2=="M"],1,sd))
                

if (sum(rownames(rrbsredB)==rownames(Bres_1dph))!=nrow(Bres_1dph)){
  stop("Annotation will be faulty!")
}
if (sum(rownames(rrbsredB)==rownames(B_obs_stat_1dph))!=nrow(Bres_1dph)){
  stop("Annotation will be faulty!")
}

func <- rep("intergenic",length(locs))
func[genic] <- "genic"
func[promoter] <- "promoter"

cn_ori <- colnames(Bres_1dph)

Bres_1dph_complete <- cbind(as.data.frame(locs)[1:2],gene[,1],gene_symbol[,1],func,Bres_1dph,B_obs_stat_1dph,gene[,2:ncol(gene)],gene_symbol[,2:ncol(gene_symbol)])
print(head(Bres_1dph_complete))
print(cn_ori)
colnames(Bres_1dph_complete) <- c("chr","locations","gene","gene_symbol","func",cn_ori,
                                  "mean_1dph,F","sd_1dph,F","mean_1dph,M","sd_1dph,M",
                                  paste("gene_alt",1:(ncol(gene)-1),sep="_"),
                                  paste("gene_symbol_alt",1:(ncol(gene)-1),sep="_"))

am_gene_alt <- sum(grepl("gene",colnames(Bres_1dph_complete)) & !grepl("gene_symbol",colnames(Bres_1dph_complete)))
inRNAseq <- rep(F,nrow(Bres_1dph_complete))
for (i in 1:am_gene_alt){
  if (i == 1){
    print(sum(inRNAseq))
    inRNAseq <- inRNAseq | (Bres_1dph_complete[,colnames(Bres_1dph_complete)=="gene"]%in%GOI_1dph)
  } else {
    print(sum(inRNAseq))
    inRNAseq <- inRNAseq | (Bres_1dph_complete[,colnames(Bres_1dph_complete)==paste("gene_alt",i-1,sep="_")]%in%GOI_1dph)
  }
}
print(sum(inRNAseq))

Bres_1dph_RNAseq <- Bres_1dph_complete[inRNAseq,]

print(head(Bres_1dph_RNAseq))
print(length(Bres_1dph_RNAseq$P.Value))
Bres_1dph_RNAseq$adj.P.Val <- p.adjust(Bres_1dph_RNAseq$P.Value,method="BH",n=length(Bres_1dph_RNAseq$P.Value))

Bres_1dph_complete <- Bres_1dph_complete[sort(Bres_1dph_complete$P.Value,index.return=T)$ix,]
Bres_1dph_RNAseq <- Bres_1dph_RNAseq[sort(Bres_1dph_RNAseq$P.Value,index.return=T)$ix,]

write.table(Bres_1dph_complete,file="RRBS_1dph_full.txt",sep="\t",col.names = T,row.names = F, quote = F)
write.table(Bres_1dph_RNAseq,file="RRBS_1dph_RNAseq.txt",sep="\t",col.names = T,row.names = F, quote = F)


# 20dph #

lmfit.cont <- contrasts.fit(fit, L[,colnames(L)=="20dph"])
lmfit.cont.ebayes <- eBayes(lmfit.cont)

Bres_20dph <- topTable(lmfit.cont.ebayes ,n=nrow(methLevel(rrbsredB.M)))

Bres_20dph <- Bres_20dph[sort(as.numeric(rownames(Bres_20dph)),index.return=T)$ix,]

B_obs_stat_20dph <- cbind(apply(methLevel(rrbsredB.rel)[,f1=="20dph" & f2=="F"],1,mean),
                          apply(methLevel(rrbsredB.rel)[,f1=="20dph" & f2=="F"],1,sd),
                          apply(methLevel(rrbsredB.rel)[,f1=="20dph" & f2=="M"],1,mean),
                          apply(methLevel(rrbsredB.rel)[,f1=="20dph" & f2=="M"],1,sd))


if (sum(rownames(rrbsredB)==rownames(Bres_20dph))!=nrow(Bres_20dph)){
  stop("Annotation will be faulty!")
}
if (sum(rownames(rrbsredB)==rownames(B_obs_stat_20dph))!=nrow(Bres_20dph)){
  stop("Annotation will be faulty!")
}

func <- rep("intergenic",length(locs))
func[genic] <- "genic"
func[promoter] <- "promoter"

cn_ori <- colnames(Bres_20dph)

Bres_20dph_complete <- cbind(as.data.frame(locs)[1:2],gene[,1],gene_symbol[,1],func,Bres_20dph,B_obs_stat_20dph,gene[,2:ncol(gene)],gene_symbol[,2:ncol(gene_symbol)])
head(Bres_20dph_complete)
colnames(Bres_20dph_complete) <- c("chr","locations","gene","gene_symbol","func",cn_ori,
                                  "mean_20dph,F","sd_20dph,F","mean_20dph,M","sd_20dph,M",
                                  paste("gene_alt",1:(ncol(gene)-1),sep="_"),
                                  paste("gene_symbol_alt",1:(ncol(gene)-1),sep="_"))

am_gene_alt <- sum(grepl("gene",colnames(Bres_20dph_complete)) & !grepl("gene_symbol",colnames(Bres_20dph_complete)))
inRNAseq <- rep(F,nrow(Bres_20dph_complete))
for (i in 1:am_gene_alt){
  if (i == 1){
    inRNAseq <- inRNAseq | (Bres_20dph_complete[,colnames(Bres_20dph_complete)=="gene"]%in%GOI_20dph)
  } else {
    inRNAseq <- inRNAseq | (Bres_20dph_complete[,colnames(Bres_20dph_complete)==paste("gene_alt",i-1,sep="_")]%in%GOI_20dph)
  }
}

Bres_20dph_RNAseq <- Bres_20dph_complete[inRNAseq,]
Bres_20dph_RNAseq$adj.P.Val <- p.adjust(Bres_20dph_RNAseq$P.Value,method="BH",length(Bres_20dph_RNAseq$P.Value))

Bres_20dph_complete <- Bres_20dph_complete[sort(Bres_20dph_complete$P.Value,index.return=T)$ix,]
Bres_20dph_RNAseq <- Bres_20dph_RNAseq[sort(Bres_20dph_RNAseq$P.Value,index.return=T)$ix,]

write.table(Bres_20dph_complete,file="RRBS_20dph_full.txt",sep="\t",col.names = T,row.names = F, quote = F)
write.table(Bres_20dph_RNAseq,file="RRBS_20dph_RNAseq.txt",sep="\t",col.names = T,row.names = F, quote = F)


# 65dph #

lmfit.cont <- contrasts.fit(fit, L[,colnames(L)=="65dph"])
lmfit.cont.ebayes <- eBayes(lmfit.cont)

Bres_65dph <- topTable(lmfit.cont.ebayes ,n=nrow(methLevel(rrbsredB.M)))

Bres_65dph <- Bres_65dph[sort(as.numeric(rownames(Bres_65dph)),index.return=T)$ix,]

B_obs_stat_65dph <- cbind(apply(methLevel(rrbsredB.rel)[,f1=="65dph" & f2=="F"],1,mean),
                          apply(methLevel(rrbsredB.rel)[,f1=="65dph" & f2=="F"],1,sd),
                          apply(methLevel(rrbsredB.rel)[,f1=="65dph" & f2=="M"],1,mean),
                          apply(methLevel(rrbsredB.rel)[,f1=="65dph" & f2=="M"],1,sd))


if (sum(rownames(rrbsredB)==rownames(Bres_65dph))!=nrow(Bres_65dph)){
  stop("Annotation will be faulty!")
}
if (sum(rownames(rrbsredB)==rownames(B_obs_stat_65dph))!=nrow(Bres_65dph)){
  stop("Annotation will be faulty!")
}

func <- rep("intergenic",length(locs))
func[genic] <- "genic"
func[promoter] <- "promoter"

cn_ori <- colnames(Bres_65dph)

Bres_65dph_complete <- cbind(as.data.frame(locs)[1:2],gene[,1],gene_symbol[,1],func,Bres_65dph,B_obs_stat_65dph,gene[,2:ncol(gene)],gene_symbol[,2:ncol(gene_symbol)])
head(Bres_65dph_complete)
colnames(Bres_65dph_complete) <- c("chr","locations","gene","gene_symbol","func",cn_ori,
                                   "mean_65dph,F","sd_65dph,F","mean_65dph,M","sd_65dph,M",
                                   paste("gene_alt",1:(ncol(gene)-1),sep="_"),
                                   paste("gene_symbol_alt",1:(ncol(gene)-1),sep="_"))

am_gene_alt <- sum(grepl("gene",colnames(Bres_65dph_complete)) & !grepl("gene_symbol",colnames(Bres_65dph_complete)))
inRNAseq <- rep(F,nrow(Bres_65dph_complete))
for (i in 1:am_gene_alt){
  if (i == 1){
    inRNAseq <- inRNAseq | (Bres_65dph_complete[,colnames(Bres_65dph_complete)=="gene"]%in%GOI_65dph)
  } else {
    inRNAseq <- inRNAseq | (Bres_65dph_complete[,colnames(Bres_65dph_complete)==paste("gene_alt",i-1,sep="_")]%in%GOI_65dph)
  }
}

Bres_65dph_RNAseq <- Bres_65dph_complete[inRNAseq,]
Bres_65dph_RNAseq$adj.P.Val <- p.adjust(Bres_65dph_RNAseq$P.Value,method="BH",length(Bres_65dph_RNAseq$P.Value))

Bres_65dph_complete <- Bres_65dph_complete[sort(Bres_65dph_complete$P.Value,index.return=T)$ix,]
Bres_65dph_RNAseq <- Bres_65dph_RNAseq[sort(Bres_65dph_RNAseq$P.Value,index.return=T)$ix,]

write.table(Bres_65dph_complete,file="RRBS_65dph_full.txt",sep="\t",col.names = T,row.names = F, quote = F)
write.table(Bres_65dph_RNAseq,file="RRBS_65dph_RNAseq.txt",sep="\t",col.names = T,row.names = F, quote = F)


# adult #

lmfit.cont <- contrasts.fit(fit, L[,colnames(L)=="adult"])
lmfit.cont.ebayes <- eBayes(lmfit.cont)

Bres_adult <- topTable(lmfit.cont.ebayes ,n=nrow(methLevel(rrbsredB.M)))

Bres_adult <- Bres_adult[sort(as.numeric(rownames(Bres_adult)),index.return=T)$ix,]

B_obs_stat_adult <- cbind(apply(methLevel(rrbsredB.rel)[,f1=="adultdph" & f2=="F"],1,mean),
                          apply(methLevel(rrbsredB.rel)[,f1=="adultdph" & f2=="F"],1,sd),
                          apply(methLevel(rrbsredB.rel)[,f1=="adultdph" & f2=="M"],1,mean),
                          apply(methLevel(rrbsredB.rel)[,f1=="adultdph" & f2=="M"],1,sd))


if (sum(rownames(rrbsredB)==rownames(Bres_adult))!=nrow(Bres_adult)){
  stop("Annotation will be faulty!")
}
if (sum(rownames(rrbsredB)==rownames(B_obs_stat_adult))!=nrow(Bres_adult)){
  stop("Annotation will be faulty!")
}

func <- rep("intergenic",length(locs))
func[genic] <- "genic"
func[promoter] <- "promoter"

cn_ori <- colnames(Bres_adult)

Bres_adult_complete <- cbind(as.data.frame(locs)[1:2],gene[,1],gene_symbol[,1],func,Bres_adult,B_obs_stat_adult,gene[,2:ncol(gene)],gene_symbol[,2:ncol(gene_symbol)])
head(Bres_adult_complete)
colnames(Bres_adult_complete) <- c("chr","locations","gene","gene_symbol","func",cn_ori,
                                   "mean_adult,F","sd_adult,F","mean_adult,M","sd_adult,M",
                                   paste("gene_alt",1:(ncol(gene)-1),sep="_"),
                                   paste("gene_symbol_alt",1:(ncol(gene)-1),sep="_"))

am_gene_alt <- sum(grepl("gene",colnames(Bres_adult_complete)) & !grepl("gene_symbol",colnames(Bres_adult_complete)))
inRNAseq <- rep(F,nrow(Bres_adult_complete))
for (i in 1:am_gene_alt){
  if (i == 1){
    inRNAseq <- inRNAseq | (Bres_adult_complete[,colnames(Bres_adult_complete)=="gene"]%in%GOI_adult)
  } else {
    inRNAseq <- inRNAseq | (Bres_adult_complete[,colnames(Bres_adult_complete)==paste("gene_alt",i-1,sep="_")]%in%GOI_adult)
  }
}

Bres_adult_RNAseq <- Bres_adult_complete[inRNAseq,]
Bres_adult_RNAseq$adj.P.Val <- p.adjust(Bres_adult_RNAseq$P.Value,method="BH",length(Bres_adult_RNAseq$P.Value))

Bres_adult_complete <- Bres_adult_complete[sort(Bres_adult_complete$P.Value,index.return=T)$ix,]
Bres_adult_RNAseq <- Bres_adult_RNAseq[sort(Bres_adult_RNAseq$P.Value,index.return=T)$ix,]

write.table(Bres_adult_complete,file="RRBS_adult_full.txt",sep="\t",col.names = T,row.names = F, quote = F)
write.table(Bres_adult_RNAseq,file="RRBS_adult_RNAseq.txt",sep="\t",col.names = T,row.names = F, quote = F)


# age #

lmfit.cont <- contrasts.fit(fit, L[,colnames(L)=="age"])
lmfit.cont.ebayes <- eBayes(lmfit.cont)

Bres_age <- topTable(lmfit.cont.ebayes ,n=nrow(methLevel(rrbsredB.M)))

Bres_age <- Bres_age[sort(as.numeric(rownames(Bres_age)),index.return=T)$ix,]

B_obs_stat_age <- cbind(apply(methLevel(rrbsredB.rel)[,f1=="1dph"],1,mean),
                          apply(methLevel(rrbsredB.rel)[,f1=="1dph"],1,sd),
                          apply(methLevel(rrbsredB.rel)[,f1!="1dph"],1,mean),
                          apply(methLevel(rrbsredB.rel)[,f1!="1dph"],1,sd))


if (sum(rownames(rrbsredB)==rownames(Bres_age))!=nrow(Bres_age)){
  stop("Annotation will be faulty!")
}
if (sum(rownames(rrbsredB)==rownames(B_obs_stat_age))!=nrow(Bres_age)){
  stop("Annotation will be faulty!")
}

func <- rep("intergenic",length(locs))
func[genic] <- "genic"
func[promoter] <- "promoter"

cn_ori <- colnames(Bres_age)

Bres_age_complete <- cbind(as.data.frame(locs)[1:2],gene[,1],gene_symbol[,1],func,Bres_age,B_obs_stat_age,gene[,2:ncol(gene)],gene_symbol[,2:ncol(gene_symbol)])
head(Bres_age_complete)
colnames(Bres_age_complete) <- c("chr","locations","gene","gene_symbol","func",cn_ori,
                                   "mean_1dph","sd_1dph","mean_older","sd_older",
                                   paste("gene_alt",1:(ncol(gene)-1),sep="_"),
                                   paste("gene_symbol_alt",1:(ncol(gene)-1),sep="_"))

am_gene_alt <- sum(grepl("gene",colnames(Bres_age_complete)) & !grepl("gene_symbol",colnames(Bres_age_complete)))
inRNAseq <- rep(F,nrow(Bres_age_complete))
for (i in 1:am_gene_alt){
  if (i == 1){
    inRNAseq <- inRNAseq | (Bres_age_complete[,colnames(Bres_age_complete)=="gene"]%in%GOI_age)
  } else {
    inRNAseq <- inRNAseq | (Bres_age_complete[,colnames(Bres_age_complete)==paste("gene_alt",i-1,sep="_")]%in%GOI_age)
  }
}

Bres_age_RNAseq <- Bres_age_complete[inRNAseq,]
Bres_age_RNAseq$adj.P.Val <- p.adjust(Bres_age_RNAseq$P.Value,method="BH",length(Bres_age_RNAseq$P.Value))

Bres_age_complete <- Bres_age_complete[sort(Bres_age_complete$P.Value,index.return=T)$ix,]
Bres_age_RNAseq <- Bres_age_RNAseq[sort(Bres_age_RNAseq$P.Value,index.return=T)$ix,]

write.table(Bres_age_complete,file="RRBS_age_full.txt",sep="\t",col.names = T,row.names = F, quote = F)
write.table(Bres_age_RNAseq,file="RRBS_age_RNAseq.txt",sep="\t",col.names = T,row.names = F, quote = F)
