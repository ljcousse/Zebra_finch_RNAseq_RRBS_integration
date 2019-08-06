####################
####################
##                ##
## DosageEffect.R ##
##                ##
####################
####################



## Load libraries
################

library(biomaRt)
library(edgeR)
library(BiSeq)
library(xlsx)


# Function for calculation of percentage of genes per bin
###############

CalcPercentage <- function(ratio,am_bins,cte){
  borders <- seq(min(ratio)-cte,max(ratio)+cte,((max(ratio)+cte)-(min(ratio)-cte))/am_bins)
  x <- NULL
  y <- NULL
  for(i in 1:(am_bins-1)){
    x <- c(x,(borders[i]+borders[i+1])/2)
    y <- c(y,sum((ratio>=borders[i]) & (ratio<borders[i+1]))*100/length(ratio))
  }
  results <- data.frame(x=x,y=y,stringsAsFactors = F)
  lo <- loess(y~x,results,span=0.05)
  results <- data.frame(x=results$x,y=lo$fitted)
  results$y[results$y < 0] <- 0
  return(results)
}


## Load Data & annotation
###############

setwd("C:/Users/Louis/OneDrive/Louis/AAP/Zebra_Finch_new/")
rnadata <- read.table("Plastoscine.txt",header=T,sep="\t",quote="")
rna_annot <- read.table("Plastoscine_annot.txt",header=F,sep="\t",quote="",stringsAsFactors = F)$V1

out2excel <- rnadata
rownames(out2excel) <- rna_annot
write.xlsx(out2excel,"RNAseq_Raw.xlsx")

# Load in ensembl annotation
ensembl=useMart("ensembl")
ensembl = useDataset("tguttata_gene_ensembl",mart=ensembl)
genes <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id','chromosome_name'),mart = ensembl)
genes <- genes[!duplicated(genes$ensembl_gene_id),]

# Sort genes (same order as rna_annot)
genes <- genes[sort(genes$ensembl_gene_id,index.return=T)$ix,]

# Check if they are identical
sum(genes$ensembl_gene_id==rna_annot)

# Make DGElist object
y <- DGEList(rnadata)
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep,  ,keep.lib.sizes=FALSE]
summary(keep)

# Perform filtering 
rna_annot <- rna_annot[keep]
chrom_names <- genes$chromosome_name[keep]


## Ratio Calculation (Log transformation)
########################################

# Make log-transformed cpm counts
log_cpm <- cpm(y,log=F)

# Get annotation for RNAseq
f1 <- "1dph,1dph,1dph,1dph,1dph,1dph,20dph,20dph,20dph,20dph,20dph,20dph,65dph,65dph,65dph,65dph,65dph,65dph,adult,adult,adult,adult,adult,adult"
f2 <- "F,F,F,M,M,M,F,F,F,M,M,M,F,F,F,M,M,M,F,F,F,M,M,M"
f1 <- unlist(strsplit(f1,","))
f2 <- unlist(strsplit(f2,","))

# Divide log counts in female and male part
log_cpm_F <- log_cpm[,f2=="F"]
log_cpm_M <- log_cpm[,f2=="M"]

# Make Z-related vs Autosomal annotation
genes_ZvsA <- rep("",length(chrom_names))
genes_ZvsA[grepl("Z",chrom_names)] <- "Z"
genes_ZvsA[!grepl("Z",chrom_names)] <- "A"


## test without random annotation 
########################################

# Calculate logcpm means for both Z and A per age and sex
means_Z_norand <- NULL
means_A_norand <- NULL
for (t in unique(paste(f1,f2,sep=","))){
  means_Z_norand <- cbind(means_Z_norand,colMeans(log_cpm[genes_ZvsA=="Z" & (!grepl("random",chrom_names)),paste(f1,f2,sep=",")==t]))
  means_A_norand <- cbind(means_A_norand,colMeans(log_cpm[genes_ZvsA=="A"& (!grepl("random",chrom_names)),paste(f1,f2,sep=",")==t]))
}
colnames(means_Z_norand) <- unique(paste(f1,f2,sep=","))
colnames(means_A_norand) <- unique(paste(f1,f2,sep=","))


# Calculate M:F ratio's for both Z and A per age
print("M:F on Z-chromosome")
MF_Z_norand <- NULL
for (i in 1:4){
  MF_Z_norand <- c(MF_Z_norand,colMeans(means_Z_norand)[i*2]/colMeans(means_Z_norand)[(i-1)*2+1])
}
print("M:F on autosomes")
MF_A_norand<- NULL
for (i in 1:4){
  MF_A_norand<- c(MF_A_norand,colMeans(means_A_norand)[i*2]/colMeans(means_A_norand)[(i-1)*2+1])
}


## test without random annotation 
########################################

# Calculate logcpm means for both Z and A per age and sex
means_Z <- NULL
means_A <- NULL
for (t in unique(paste(f1,f2,sep=","))){
  means_Z <- cbind(means_Z,colMeans(log_cpm[genes_ZvsA=="Z",paste(f1,f2,sep=",")==t]))
  means_A <- cbind(means_A,colMeans(log_cpm[genes_ZvsA=="A",paste(f1,f2,sep=",")==t]))
}
colnames(means_Z) <- unique(paste(f1,f2,sep=","))
colnames(means_A) <- unique(paste(f1,f2,sep=","))

rna_annot_Z <- as.character(rna_annot[genes_ZvsA=="Z"])
rna_annot_A <- as.character(rna_annot[genes_ZvsA=="A"])

# Calculate M:F ratio's for both Z and A per age
print("M:F on Z-chromosome")
MF_Z <- NULL
for (i in 1:4){
  MF_Z <- c(MF_Z,colMeans(means_Z)[i*2]/colMeans(means_Z)[(i-1)*2+1])
}
print("M:F on autosomes")
MF_A<- NULL
for (i in 1:4){
  MF_A <- c(MF_A,colMeans(means_A)[i*2]/colMeans(means_A)[(i-1)*2+1])
}

# Check for statistical differences between Z and A
MF_all <- data.frame(cbind(c(MF_Z,MF_A),rep(unique(f1),2),c(rep("Z",4),rep("A",4))))
colnames(MF_all) <- c("MF","f1","f2")
MF_all$MF <- as.numeric(as.character(MF_all$MF))
anova_MF <- aov(MF~f1,MF_all)
summary(anova_MF)

# Check assumptions
ks.test(anova_MF$residuals,pnorm)
qqnorm(anova_MF$residuals)
qqline(anova_MF$residuals)

# Calculate Z:A ratios for all ages and sexes
print("Z:A")
ZA <- NULL
for (t in 1:8){
  ZA <- c(ZA,colMeans(means_Z)[t]/colMeans(means_A)[t])
}

# Check for statistical differences between M and F on all ages
ZA_all <- data.frame(cbind(c(means_Z/means_A),f1,f2))
colnames(ZA_all) <- c("ZA","f1","f2")
ZA_all$ZA <- as.numeric(as.character(ZA_all$ZA))
anova_ZA <- aov(ZA~f1*f2,ZA_all)
summary(anova_ZA)
TukeyHSD(anova_ZA)

# Check assumptions
ks.test(anova_ZA$residuals,pnorm)
qqnorm(anova_ZA$residuals)
qqline(anova_ZA$residuals)

# Calculate Z:A ratio's for M and F seperatly
ZA_F <- means_Z[,c(1,3,5,7)]/means_A[,c(1,3,5,7)]
ZA_M <- means_Z[,c(2,4,6,8)]/means_A[,c(2,4,6,8)]

# Plot ZA ratios
tiff("ZAratios.tiff",width=2000, height=2000,units="px",type="windows",res = 300)
plot(1:12,ZA_F,pch=17,xlab="",ylab="Z:A ratio",xaxt="n",ylim=c(min(c(ZA_F,ZA_M))*0.90,max(c(ZA_F,ZA_M))*1.15))
points(1:12,ZA_M,pch=18)
segments(0.7, 0.44, 3.3, 0.44, xpd=TRUE)
segments(3.7, 0.44, 6.3, 0.44, xpd=TRUE)
segments(6.7, 0.44, 9.3, 0.44, xpd=TRUE)
segments(9.7, 0.44, 12.3, 0.44, xpd=TRUE)
legend("topright",pch=c(17,18),c("Females","Males"))
text(x=c(2,5,8,11),y=rep(0.4,4),labels=c("1dph","20dph","65dph","adult"),las=2,xpd=T)
dev.off()

# Plot boxplot of logcpm values for both Z and A
jpeg("Boxplot_expression_log.jpg",width = 1000,height = 500)
par(mar=c(7,4.1,4.1,2.1))
boxplot(rowMeans(log_cpm)~genes_ZvsA,las=2)
dev.off()

# Statistical differences between both?
summary(aov(rowMeans(log_cpm)~genes_ZvsA))
summary(aov(rowMeans(log_cpm_F)~genes_ZvsA))
summary(aov(rowMeans(log_cpm_M)~genes_ZvsA))


## Plot dosage effect + genelist for stat analysis
###################################

# Calculate ratio's for Autosomal genes
ratio_A <- NULL
logcpm_A_M <- log_cpm[genes_ZvsA=="A",f2=="M"]
logcpm_A_F <- log_cpm[genes_ZvsA=="A",f2=="F"]
cte <- 0.1
for (t in 1:nrow(logcpm_A_F)){
  ratio_A <- c(ratio_A,(mean(logcpm_A_M[t,])+cte)/(mean(logcpm_A_F[t,])+cte))
}
ratio_A <- ratio_A[!is.na(ratio_A)]

# Calculate ratio's for Z-related genes
ratio_Z <- NULL
logcpm_Z_M <- log_cpm[genes_ZvsA=="Z",f2=="M"]
logcpm_Z_F <- log_cpm[genes_ZvsA=="Z",f2=="F"]
cte <- 0.1
for (t in 1:nrow(logcpm_Z_F)){
  ratio_Z <- c(ratio_Z,(mean(logcpm_Z_M[t,])+cte)/(mean(logcpm_Z_F[t,])+cte))
}
ratio_Z <- ratio_Z[!is.na(ratio_Z)]

# Calculate percentages plot
per_Z <- CalcPercentage(ratio_Z,100,0.1)
per_A <- CalcPercentage(ratio_A,100,0.1)

# Make percentages plots
jpeg("dosage_compensation_percentages.jpg")
plot(per_A$x,per_A$y,type="l",main="Dosage Compensation",ylab="Percentage of genes",xlab="Z:A ratio",xlim=c(0,3))
lines(per_Z$x,per_Z$y,col="red")
dev.off()

# Density plots
den_A <- density(ratio_A)
den_Z <- density(ratio_Z)
jpeg("dosage_compensation.jpg")
plot(den_A,xlim=c(0,3),main="Dosage compensation",xlab="",type="l")
lines(den_Z,col="red")
dev.off()


## Make Gene lists for further analysis
###################################

# Gene lists dosage vs no dosage
genes_A <- genes[genes_ZvsA=="A",]
genes_Z <- genes[genes_ZvsA=="Z",]

# Calculate 95% interval for Z-related ratios
quantile((logcpm_Z_M)/(logcpm_Z_F),probs=c(0.05,0.95),na.rm=T)

# Calculate 95% interval for Autosomal ratios: borders will be used as cutoffs for "dosage compensated" genes
cutoffs <- quantile((logcpm_A_M+cte)/(logcpm_A_F+cte),probs=c(0.05,0.95,0.975),na.rm=T)
print(cutoffs)

# Separation autosomal ratio's per age
ratio_A_1dph <- NULL
logcpm_A_M_1dph <- log_cpm[genes_ZvsA=="A",f2=="M" & f1=="1dph"]
logcpm_A_F_1dph <- log_cpm[genes_ZvsA=="A",f2=="F" & f1=="1dph"]
ratio_A_20dph <- NULL
logcpm_A_M_20dph <- log_cpm[genes_ZvsA=="A",f2=="M" & f1=="20dph"]
logcpm_A_F_20dph <- log_cpm[genes_ZvsA=="A",f2=="F" & f1=="20dph"]
ratio_A_65dph <- NULL
logcpm_A_M_65dph <- log_cpm[genes_ZvsA=="A",f2=="M" & f1=="65dph"]
logcpm_A_F_65dph <- log_cpm[genes_ZvsA=="A",f2=="F" & f1=="65dph"]
ratio_A_adult <- NULL
logcpm_A_M_adult <- log_cpm[genes_ZvsA=="A",f2=="M" & f1=="adult"]
logcpm_A_F_adult <- log_cpm[genes_ZvsA=="A",f2=="F" & f1=="adult"]
cte <- 0.1
for (t in 1:nrow(logcpm_A_F_1dph)){
  ratio_A_1dph <- c(ratio_A_1dph,(mean(logcpm_A_M_1dph[t,])+cte)/(mean(logcpm_A_F_1dph[t,])+cte))
  ratio_A_20dph <- c(ratio_A_20dph,(mean(logcpm_A_M_20dph[t,])+cte)/(mean(logcpm_A_F_20dph[t,])+cte))
  ratio_A_65dph <- c(ratio_A_65dph,(mean(logcpm_A_M_65dph[t,])+cte)/(mean(logcpm_A_F_65dph[t,])+cte))
  ratio_A_adult <- c(ratio_A_adult,(mean(logcpm_A_M_adult[t,])+cte)/(mean(logcpm_A_F_adult[t,])+cte))
}

# Separation z-related ratio's per age
ratio_Z_1dph <- NULL
logcpm_Z_M_1dph <- log_cpm[genes_ZvsA=="Z",f2=="M" & f1=="1dph"]
logcpm_Z_F_1dph <- log_cpm[genes_ZvsA=="Z",f2=="F" & f1=="1dph"]
ratio_Z_20dph <- NULL
logcpm_Z_M_20dph <- log_cpm[genes_ZvsA=="Z",f2=="M" & f1=="20dph"]
logcpm_Z_F_20dph <- log_cpm[genes_ZvsA=="Z",f2=="F" & f1=="20dph"]
ratio_Z_65dph <- NULL
logcpm_Z_M_65dph <- log_cpm[genes_ZvsA=="Z",f2=="M" & f1=="65dph"]
logcpm_Z_F_65dph <- log_cpm[genes_ZvsA=="Z",f2=="F" & f1=="65dph"]
ratio_Z_adult <- NULL
logcpm_Z_M_adult <- log_cpm[genes_ZvsA=="Z",f2=="M" & f1=="adult"]
logcpm_Z_F_adult <- log_cpm[genes_ZvsA=="Z",f2=="F" & f1=="adult"]
cte <- 0.1
for (t in 1:nrow(logcpm_Z_F_1dph)){
  ratio_Z_1dph <- c(ratio_Z_1dph,(mean(logcpm_Z_M_1dph[t,])+cte)/(mean(logcpm_Z_F_1dph[t,])+cte))
  ratio_Z_20dph <- c(ratio_Z_20dph,(mean(logcpm_Z_M_20dph[t,])+cte)/(mean(logcpm_Z_F_20dph[t,])+cte))
  ratio_Z_65dph <- c(ratio_Z_65dph,(mean(logcpm_Z_M_65dph[t,])+cte)/(mean(logcpm_Z_F_65dph[t,])+cte))
  ratio_Z_adult <- c(ratio_Z_adult,(mean(logcpm_Z_M_adult[t,])+cte)/(mean(logcpm_Z_F_adult[t,])+cte))
}


# Make percentage plots per age
tiff("dosage_compensation_diffage_percentage.tiff",width=2000, height=2000,units="px",type="windows",res = 300)
par(mfrow=c(2,2))
per_Z <- CalcPercentage(ratio_Z_1dph,1000,0.1)
per_A <- CalcPercentage(ratio_A_1dph,1000,0.1)
plot(per_A$x,per_A$y,type="l",main="1dph",ylab="Percentage of genes",xlab="M:F ratio",xlim=c(0,3),ylim=c(0.1,5.5))
lines(per_Z$x,per_Z$y,col="red")
legend("topright",c(paste("A:",round(MF_A[1],digits=2)),paste("Z:",round(MF_Z[1],digits=2))))
per_Z <- CalcPercentage(ratio_Z_20dph,1000,0.1)
per_A <- CalcPercentage(ratio_A_20dph,1000,0.1)
plot(per_A$x,per_A$y,type="l",main="20dph",ylab="Percentage of genes",xlab="M:F ratio",xlim=c(0,3),ylim=c(0.1,5.5))
lines(per_Z$x,per_Z$y,col="red")
legend("topright",c(paste("A:",round(MF_A[2],digits=2)),paste("Z:",round(MF_Z[2],digits=2))))
per_Z <- CalcPercentage(ratio_Z_65dph,1000,0.1)
per_A <- CalcPercentage(ratio_A_65dph,1000,0.1)
plot(per_A$x,per_A$y,type="l",main="65dph",ylab="Percentage of genes",xlab="M:F ratio",xlim=c(0,3),ylim=c(0.1,5.5))
lines(per_Z$x,per_Z$y,col="red")
legend("topright",c(paste("A:",round(MF_A[3],digits=2)),paste("Z:",round(MF_Z[3],digits=2))))
per_Z <- CalcPercentage(ratio_Z_adult,1000,0.1)
per_A <- CalcPercentage(ratio_A_adult,1000,0.1)
plot(per_A$x,per_A$y,type="l",main="adult",ylab="Percentage of genes",xlab="M:F ratio",xlim=c(0,3),ylim=c(0.1,5.5))
lines(per_Z$x,per_Z$y,col="red")
legend("topright",c(paste("A:",round(MF_A[4],digits=2)),paste("Z:",round(MF_Z[4],digits=2))))
par(mfrow=c(1,1))
dev.off()

# Make density plots per age
den_A_1dph <- density(ratio_A_1dph)
den_A_20dph <- density(ratio_A_20dph)
den_A_65dph <- density(ratio_A_65dph)
den_A_adult <- density(ratio_A_adult)
den_Z_1dph <- density(ratio_Z_1dph)
den_Z_20dph <- density(ratio_Z_20dph)
den_Z_65dph <- density(ratio_Z_65dph)
den_Z_adult <- density(ratio_Z_adult)
jpeg("dosage_compensation_diffage.jpg")
par(mfrow=c(2,2))
plot(den_A_1dph,xlim=c(0,3),main="20dph",xlab="",type="l")
lines(den_Z_1dph,col="red")
plot(den_A_20dph,xlim=c(0,3),main="20dph",xlab="",type="l")
lines(den_Z_20dph,col="red")
plot(den_A_65dph,xlim=c(0,3),main="65dph",xlab="",type="l")
lines(den_Z_65dph,col="red")
plot(den_A_adult,xlim=c(0,3),main="adult",xlab="",type="l")
lines(den_Z_adult,col="red")
par(mfrow=c(1,1))
dev.off()

# Make genelists "diffdosage": meaning not differentially no dosage compensation on 1dph but dosage compensation for adults
genelist_dosage <- c(rna_annot_Z[(ratio_Z_1dph>cutoffs[2]) & (ratio_Z_adult<cutoffs[2])])
write.table(genelist_dosage,file="GeneLists/ListGenes_dosage.txt",col.names = F,row.names = F,quote = F)
genelist_nodosage <- c(rna_annot_Z[(ratio_Z_1dph>cutoffs[2]) & (ratio_Z_adult>cutoffs[2])])
write.table(genelist_nodosage,file="GeneLists/ListGenes_nodosage.txt",col.names = F,row.names = F,quote = F)
genelist_dosage_bis <- c(rna_annot_Z[(ratio_Z_1dph<cutoffs[2]) & (ratio_Z_adult<cutoffs[2])])
write.table(genelist_dosage_bis,file="GeneLists/ListGenes_dosage_bis.txt",col.names = F,row.names = F,quote = F)
genelist_nodosage_bis <- c(rna_annot_Z[(ratio_Z_1dph<cutoffs[2]) & (ratio_Z_adult>cutoffs[2])])
length(rna_annot_Z)


cbind(ratio_Z_1dph[rna_annot_Z%in%genelist_dosage],ratio_Z_adult[rna_annot_Z%in%genelist_dosage])
cbind(ratio_Z_1dph[rna_annot_Z%in%genelist_nodosage],ratio_Z_adult[rna_annot_Z%in%genelist_nodosage])


test <- rbind(cbind(rowMeans(logcpm[rna_annot%in%genelist_dosage,]),"dosage"),cbind(rowMeans(logcpm[rna_annot%in%genelist_nodosage,]),"no_dosage"),cbind(rowMeans(logcpm[rna_annot%in%genelist_dosage_bis,]),"dosage_bis"))
test <- data.frame(test,stringsAsFactors = F)
colnames(test) <- c("logcpm","category")
test$category <- as.factor(test$category)
test$logcpm <- as.numeric(test$logcpm)
boxplot(logcpm~category,data=test)

anova_dosage <- aov(logcpm~category,data=test)
summary(anova_dosage)
TukeyHSD(anova_dosage)

## Methylation Comparison between dosage - no dosage
#######################################"


# Load RRBS data
load(file="allrrbs.Rda") # allrrbs.Rda can be retrieved using RRBS_analyis.R script

# Set minimal coverages
mincov <- 5

# Load in RRBS annotation
annot<-c("Pool","Sex","Age","8","F","20","12","M","20","14","F","65","18","M","65","1","F","1","7","F","20","10","M","20","24","M","adult","2","F","1","15","F","65","19","F","adult","22","M","adult","9","F","20","11","M","20","16","M","65","17","M","65","4","M","1","5","M","1","13","F","65","20","F","adult","3","F","1","6","M","1","21","F","adult","23","M","adult")
annot<-matrix(annot,nrow=25,ncol=3,byrow=T)
colnames(annot)<-annot[1,]
annot<-annot[-1,]
rownames(annot)<-annot[,1]
pids<-colnames(allrrbs)
annot<-cbind(annot[substr(pids,2,3),],substr(pids,1,1))
annot<-cbind(annot,apply(cbind(annot[,"Sex"],annot[,"Age"]),1,paste,collapse="_"))
rownames(annot)<-pids
colnames(annot)[4:5]<-c("Ox_vs_RRBS","Fenolabel")


# Annotation RRBS
#Filtering based on shared regions (should be shared with coverage >= mincov in all samples)
rrbs<-filterBySharedRegions(object=allrrbs,perc.samples=1,minCov=mincov)
rrbs.rel <- rawToRel(rrbs)

locs<-rowRanges(rrbs)
# print("locs")
# locations <- locs@ranges@start
# print(head(locations))
# print(str(locs))


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
  # print(head(n_times))
  # print(dim(n_times))
  # print(head(as.character(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]])))
  # print(length(as.character(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]])))
  # print(head(c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]))
  # print(length(c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]))
  # print(head(-c(c(1:(nrow(n_times)/i))*i-(i-1))))
  # print(tail(-c(c(1:(nrow(n_times)/i))*i-(i-1))))
  # print(length(-c(c(1:(nrow(n_times)/i))*i-(i-1))))
  
  if (i==2 | nrow(n_times)==i){
    # print(length(as.character(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]])))
    # print(length(gene[locs@ranges@NAMES%in%unique(n_times$rowname),2:i]))
    gene[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]
    gene_symbol[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- n_times$gene_symbol[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]
  } else {
    # print(head(t(matrix(as.character(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]),nrow=i-1))))
    
    # print(dim(t(matrix(as.character(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]]),nrow=i-1))))
    # print(dim(gene[locs@ranges@NAMES%in%unique(n_times$rowname),2:i]))
    gene[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- t(matrix(n_times$gene[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]],nrow=i-1))
    gene_symbol[locs@ranges@NAMES%in%unique(n_times$rowname),2:i] <- t(matrix(n_times$gene_symbol[c(1:((nrow(n_times)/i)*i))[-c(c(1:(nrow(n_times)/i))*i-(i-1))]],nrow=i-1))
  }
  
}


# Filtering on minimal methylated reads
minmethreads <- 6

inclocsM2 <- rowSums(methReads(rrbs))>=minmethreads

locs <- locs[inclocsM2]
rrbsredB <- rrbs[inclocsM2,]

genic <- genic[inclocsM2]
promoter <- promoter[inclocsM2]

gene <- gene[inclocsM2,]
gene_symbol <- gene_symbol[inclocsM2,]


# RRBS_annotation
f1 <- paste0(annot[,3],"dph")
f2 <- annot[,2]
f1 <- as.factor(f1)
f2 <- as.factor(f2)


out2excel <- cbind(as.data.frame(locs),gene[,1],methLevel(rawToRel(rrbsredB))[,sort(as.numeric(annot[,1]),index.return=T)$ix])
out2excel <- out2excel[,c(-4,-5)]
colnames(out2excel)[1] <- "chr"
colnames(out2excel)[4] <- "gene"
write.table(out2excel,"RRBS_Raw.txt",col.names = T,row.names = F,quote = F, sep="\t")


## Testje
RRBS_1dph_full <- read.table("RRBS_1dph_full.txt",header=T,stringsAsFactors = F,sep="\t")  
RRBS_20dph_full <- read.table("RRBS_20dph_full.txt",header=T,stringsAsFactors = F,sep="\t")  
RRBS_65dph_full <- read.table("RRBS_65dph_full.txt",header=T,stringsAsFactors = F,sep="\t")  
RRBS_adult_full <- read.table("RRBS_adult_full.txt",header=T,stringsAsFactors = F,sep="\t")  


goi <- RRBS_adult_full$gene_symbol[RRBS_adult_full$gene[RRBS_adult_full$P.Value<0.05]%in%genelist_dosage]
goi <- unique(goi[!(is.na(goi)) & (goi!="")])

for (i in 1:length(goi)){
  if(sum((RRBS_adult_full$gene_symbol==goi[i]) & (!is.na(RRBS_adult_full$gene_symbol))& (RRBS_adult_full$func=="promoter"))==0){
    next
  }
  # plot(RRBS_adult_full$locations[(RRBS_adult_full$gene_symbol==goi[i]) & (!is.na(RRBS_adult_full$gene_symbol)) & (RRBS_adult_full$func=="promoter")],
  #      RRBS_adult_full$logFC[(RRBS_adult_full$gene_symbol==goi[i]) & (!is.na(RRBS_adult_full$gene_symbol)) & (RRBS_adult_full$func=="promoter")],
  #      main=goi[i],xlab='pos',ylab="logFC",ylim=c(-lim,lim))
  # points(RRBS_1dph_full$locations[(RRBS_1dph_full$gene_symbol==goi[i]) & (!is.na(RRBS_1dph_full$gene_symbol)) & (RRBS_1dph_full$func=="promoter")],
  #        RRBS_1dph_full$logFC[(RRBS_1dph_full$gene_symbol==goi[i]) & (!is.na(RRBS_1dph_full$gene_symbol)) & (RRBS_1dph_full$func=="promoter")],
  #        main=goi[i],xlab='pos',ylab="logFC",ylim=c(-lim,lim),col="blue")
  ix_adult <- sort(as.numeric(RRBS_adult_full$locations[(RRBS_adult_full$gene_symbol==goi[i]) & (!is.na(RRBS_adult_full$gene_symbol)) & (RRBS_adult_full$func=="promoter")]),index.return=T)$ix
  ix_1dph <- sort(as.numeric(RRBS_adult_full$locations[(RRBS_1dph_full$gene_symbol==goi[i]) & (!is.na(RRBS_1dph_full$gene_symbol)) & (RRBS_1dph_full$func=="promoter")]),index.return=T)$ix
  lim <- 1.2*max(abs(RRBS_adult_full$logFC[(RRBS_adult_full$gene_symbol==goi[i]) & (!is.na(RRBS_adult_full$gene_symbol)) & (RRBS_adult_full$func=="promoter")][ix_adult] -RRBS_1dph_full$logFC[(RRBS_1dph_full$gene_symbol==goi[i]) & (!is.na(RRBS_1dph_full$gene_symbol)) & (RRBS_1dph_full$func=="promoter")][ix_1dph]))
  plot(RRBS_adult_full$locations[(RRBS_adult_full$gene_symbol==goi[i]) & (!is.na(RRBS_adult_full$gene_symbol)) & (RRBS_adult_full$func=="promoter")][ix_adult],
       RRBS_adult_full$logFC[(RRBS_adult_full$gene_symbol==goi[i]) & (!is.na(RRBS_adult_full$gene_symbol)) & (RRBS_adult_full$func=="promoter")][ix_adult] -RRBS_1dph_full$logFC[(RRBS_1dph_full$gene_symbol==goi[i]) & (!is.na(RRBS_1dph_full$gene_symbol)) & (RRBS_1dph_full$func=="promoter")][ix_1dph],
       main=goi[i],xlab='pos',ylab="logFC",ylim=c(-lim,lim))
  
  abline(h=0,col="red")
  Sys.sleep(3)
}

RRBS_adult_full[RRBS_adult_full$gene_symbol=="CHAC1" & !is.na(RRBS_adult_full$gene_symbol),]

counter <- 1
meth_colmeans_dosage <- NULL
make_submat <- function(dat){
  res <- cbind(dat,as.character(f1),as.character(f2),i,1:24)
  return(data.frame(res))
}
for (i in 1:length(genelist_dosage)){
  if (sum((gene[,1]==genelist_dosage[i]) & (!is.na(gene[,1])) & promoter)==0){
    next
  } else if (sum((gene[,1]==genelist_dosage[i]) & (!is.na(gene[,1])) & promoter)==1){
    meth_colmeans_dosage <- rbind(meth_colmeans_dosage,cbind(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_dosage[i]) & (!is.na(gene[,1])) & promoter,],as.character(f1),as.character(f2),i,1:24))
    points(1:24,methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_dosage[i]) & (!is.na(gene[,1])),][sort(as.numeric(annot[,1]),index.return=T)$ix])
  } else {
    # meth_colmeans_dosage <- rbind(meth_colmeans_dosage,as.matrix(do.call(rbind,apply(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_dosage[i]) & (!is.na(gene[,1])) & promoter,],1,make_submat))))
    meth_colmeans_dosage <- rbind(meth_colmeans_dosage,cbind(colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_dosage[i]) & (!is.na(gene[,1])) & promoter,]),as.character(f1),as.character(f2),i,1:24))
    if (counter == 1){
      plot(1:24,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_dosage[i]) & (!is.na(gene[,1])) & promoter,sort(as.numeric(annot[,1]),index.return=T)$ix]),ylim=c(0,1))
      counter <- counter +1
    } else {
      points(1:24,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_dosage[i]) & (!is.na(gene[,1])) & promoter,sort(as.numeric(annot[,1]),index.return=T)$ix]))
    }
  }
}
print(dim(meth_colmeans_dosage))
meth_colmeans_dosage <- as.data.frame(meth_colmeans_dosage)
colnames(meth_colmeans_dosage) <- c("meth","f1","f2","gene","individual")
meth_colmeans_dosage$meth <- as.numeric(as.character(meth_colmeans_dosage$meth))
anova_dosage <- aov(meth~gene+f1*f2,data=meth_colmeans_dosage)
summary(anova_dosage)
TukeyHSD(anova_dosage)$'f1:f2'[TukeyHSD(anova_dosage)$'f1:f2'[,4]<0.05,]
car::leveneTest(meth~f1*f2,data=meth_colmeans_dosage)
par(mfrow=c(2,2))
plot(anova_dosage)
par(mfrow=c(1,1))

counter <- 1
meth_colmeans_nodosage <- NULL
for (i in 1:length(genelist_nodosage)){
  if (sum((gene[,1]==genelist_nodosage[i]) & (!is.na(gene[,1])))==0){
    next
  } else if (sum((gene[,1]==genelist_nodosage[i]) & (!is.na(gene[,1])))==1){
    meth_colmeans_nodosage <- rbind(meth_colmeans_nodosage,cbind(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_nodosage[i]) & (!is.na(gene[,1])),],as.character(f1),as.character(f2),i))
    points(1:24,methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_nodosage[i]) & (!is.na(gene[,1])),sort(as.numeric(annot[,1]),index.return=T)$ix])
  } else {
    meth_colmeans_nodosage <- rbind(meth_colmeans_nodosage,cbind(colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_nodosage[i]) & (!is.na(gene[,1])),]),as.character(f1),as.character(f2),i))
    if (counter == 1){
      plot(1:24,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_nodosage[i]) & (!is.na(gene[,1])),sort(as.numeric(annot[,1]),index.return=T)$ix]),ylim=c(0,1))
      counter <- counter +1
    } else {
      points(1:24,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_nodosage[i]) & (!is.na(gene[,1])),sort(as.numeric(annot[,1]),index.return=T)$ix]))
    }
  }
}
meth_colmeans_nodosage <- as.data.frame(meth_colmeans_nodosage)
colnames(meth_colmeans_nodosage) <- c("meth","f1","f2","gene")
meth_colmeans_nodosage$meth <- as.numeric(as.character(meth_colmeans_nodosage$meth))
anova_nodosage <- aov(meth~gene+f1*f2,data=meth_colmeans_nodosage)
summary(anova_nodosage)
TukeyHSD(anova_nodosage)$'f1:f2'[TukeyHSD(anova_nodosage)$'f1:f2'[,4]<0.05,]
car::leveneTest(meth~f1*f2,data=meth_colmeans_nodosage)
par(mfrow=c(2,2))
plot(anova_nodosage)


counter <- 1
meth_colmeans_nodosage_bis <- NULL
for (i in 1:length(genelist_nodosage_bis)){
  if (sum((gene[,1]==genelist_nodosage_bis[i]) & (!is.na(gene[,1])))==0){
    next
  } else if (sum((gene[,1]==genelist_nodosage_bis[i]) & (!is.na(gene[,1])))==1){
    meth_colmeans_nodosage_bis <- rbind(meth_colmeans_nodosage_bis,cbind(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_nodosage_bis[i]) & (!is.na(gene[,1])),],as.character(f1),as.character(f2),i))
    points(1:24,methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_nodosage_bis[i]) & (!is.na(gene[,1])),sort(as.numeric(annot[,1]),index.return=T)$ix])
  } else {
    meth_colmeans_nodosage_bis <- rbind(meth_colmeans_nodosage_bis,cbind(colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_nodosage_bis[i]) & (!is.na(gene[,1])),]),as.character(f1),as.character(f2),i))
    if (counter == 1){
      plot(1:24,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_nodosage_bis[i]) & (!is.na(gene[,1])),sort(as.numeric(annot[,1]),index.return=T)$ix]),ylim=c(0,1))
      counter <- counter +1
    } else {
      points(1:24,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_nodosage_bis[i]) & (!is.na(gene[,1])),sort(as.numeric(annot[,1]),index.return=T)$ix]))
    }
  }
}
meth_colmeans_nodosage_bis <- as.data.frame(meth_colmeans_nodosage_bis)
colnames(meth_colmeans_nodosage_bis) <- c("meth","f1","f2","gene")
meth_colmeans_nodosage_bis$meth <- as.numeric(as.character(meth_colmeans_nodosage_bis$meth))
anova_nodosage_bis <- aov(meth~gene+f1*f2,data=meth_colmeans_nodosage_bis)
summary(anova_nodosage_bis)
TukeyHSD(anova_nodosage_bis)$'f1:f2'[TukeyHSD(anova_nodosage_bis)$'f1:f2'[,4]<0.05,]
car::leveneTest(meth~f1*f2,data=meth_colmeans_nodosage_bis)
par(mfrow=c(2,2))
plot(anova_nodosage_bis)

rownames(TukeyHSD(anova_nodosage)$'f1:f2'[TukeyHSD(anova_nodosage)$'f1:f2'[,4]<0.05,])[!(rownames(TukeyHSD(anova_nodosage)$'f1:f2'[TukeyHSD(anova_nodosage)$'f1:f2'[,4]<0.05,])%in%rownames(TukeyHSD(anova_dosage)$'f1:f2'[TukeyHSD(anova_dosage)$'f1:f2'[,4]<0.05,]))]
rownames(TukeyHSD(anova_dosage)$'f1:f2'[TukeyHSD(anova_dosage)$'f1:f2'[,4]<0.05,])[!(rownames(TukeyHSD(anova_dosage)$'f1:f2'[TukeyHSD(anova_dosage)$'f1:f2'[,4]<0.05,])%in%rownames(TukeyHSD(anova_nodosage)$'f1:f2'[TukeyHSD(anova_nodosage)$'f1:f2'[,4]<0.05,]))]


counter <- 1
meth_colmeans_diffdosage <- NULL
for (i in 1:length(genelist_diffdosage)){
  if (sum((gene[,1]==genelist_diffdosage[i]) & (!is.na(gene[,1])))==0){
    next
  } else if (sum((gene[,1]==genelist_diffdosage[i]) & (!is.na(gene[,1])))==1){
    meth_colmeans_diffdosage <- rbind(meth_colmeans_diffdosage,cbind(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_diffdosage[i]) & (!is.na(gene[,1])),],f1,f2,i))
    points(1:24,methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_diffdosage[i]) & (!is.na(gene[,1])),])
  } else {
    meth_colmeans_diffdosage <- rbind(meth_colmeans_diffdosage,cbind(colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_diffdosage[i]) & (!is.na(gene[,1])),]),f1,f2,i))
    if (counter == 1){
      plot(1:24,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_diffdosage[i]) & (!is.na(gene[,1])),]),ylim=c(0,1))
      counter <- counter +1
    } else {
      points(1:24,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist_diffdosage[i]) & (!is.na(gene[,1])),]))
    }
  }
}
meth_colmeans_diffdosage <- as.data.frame(meth_colmeans_diffdosage)
colnames(meth_colmeans_diffdosage) <- c("meth","f1","f2","gene")
meth_colmeans_diffdosage$meth <- as.numeric(as.character(meth_colmeans_diffdosage$meth))
anova_diffdosage <- aov(meth~gene+f1*f2,data=meth_colmeans_diffdosage)
summary(anova_diffdosage)
TukeyHSD(anova_diffdosage)$'f1:f2'[TukeyHSD(anova_diffdosage)$'f1:f2'[,4]<0.05,]
car::leveneTest(meth~f1*f2,data=meth_colmeans_diffdosage)
par(mfrow=c(2,2))
plot(anova_diffdosage)

# DIFF METHYLATION
#####################

meth_analysis <- NULL

# ALL #
#######

type <- "both"

# Make subsets
meth_dosage_M <- methLevel(rrbs.rel)[gene[,1]%in%genelist_dosage,f2=="M"]
meth_dosage_F <- methLevel(rrbs.rel)[gene[,1]%in%genelist_dosage,f2=="F"]
meth_nodosage_M <- methLevel(rrbs.rel)[gene[,1]%in%genelist_nodosage,f2=="M"]
meth_nodosage_F <- methLevel(rrbs.rel)[gene[,1]%in%genelist_nodosage,f2=="F"]

# Male
sex <- "M"

# Make dataframe for boxplot and anova
dat_boxplot <- as.data.frame(cbind(c(rowMeans(meth_dosage_M),rowMeans(meth_nodosage_M))
                                   ,c(rep("D",nrow(meth_dosage_M)),rep("ND",nrow(meth_nodosage_M)))))
colnames(dat_boxplot) <- c("y","DC")
dat_boxplot$y <- as.numeric(as.character(dat_boxplot$y))

# Make boxplot
boxplot(y~DC,dat_boxplot)

# Anova
summary(aov(y~DC,dat_boxplot))
print(mean(dat_boxplot$y[dat_boxplot$DC=="D"]))
print(mean(dat_boxplot$y[dat_boxplot$DC=="ND"]))

# Add results to result object
meth_analysis <- rbind(meth_analysis,c(type,sex,mean(dat_boxplot$y[dat_boxplot$DC=="D"]),
                                       mean(dat_boxplot$y[dat_boxplot$DC=="ND"]),
                                       summary(aov(y~DC,dat_boxplot))[[1]][["Pr(>F)"]][1]))

# Female
sex <- "F"

# Make dataframe for boxplot and anova
dat_boxplot <- as.data.frame(cbind(c(rowMeans(meth_dosage_F),rowMeans(meth_nodosage_F)),
                                   c(rep("D",nrow(meth_dosage_F)),rep("ND",nrow(meth_nodosage_F)))))
colnames(dat_boxplot) <- c("y","DC")
dat_boxplot$y <- as.numeric(as.character(dat_boxplot$y))

# Make boxplot
boxplot(y~DC,dat_boxplot)

# Anova
summary(aov(y~DC,dat_boxplot))
mean(dat_boxplot$y[dat_boxplot$DC=="D"])
mean(dat_boxplot$y[dat_boxplot$DC=="ND"])

# Add results to result object
meth_analysis <- rbind(meth_analysis,c(type,sex,mean(dat_boxplot$y[dat_boxplot$DC=="D"]),
                                       mean(dat_boxplot$y[dat_boxplot$DC=="ND"]),
                                       summary(aov(y~DC,dat_boxplot))[[1]][["Pr(>F)"]][1]))



# PROMOTER #
############

type <- "promotor"

# Make subsets
meth_dosage_M <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_dosage) & promoter),f2=="M"]
meth_dosage_F <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_dosage) & promoter),f2=="F"]
meth_nodosage_M <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_nodosage) & promoter),f2=="M"]
meth_nodosage_F <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_nodosage) & promoter),f2=="F"]

# Male
sex <- "M"

# Make dataframe for boxplot and anova
dat_boxplot <- as.data.frame(cbind(c(rowMeans(meth_dosage_M),rowMeans(meth_nodosage_M))
                                   ,c(rep("D",nrow(meth_dosage_M)),rep("ND",nrow(meth_nodosage_M)))))
colnames(dat_boxplot) <- c("y","DC")
dat_boxplot$y <- as.numeric(as.character(dat_boxplot$y))

# Make boxplot
boxplot(y~DC,dat_boxplot)

# Anova
summary(aov(y~DC,dat_boxplot))
print(mean(dat_boxplot$y[dat_boxplot$DC=="D"]))
print(mean(dat_boxplot$y[dat_boxplot$DC=="ND"]))

# Add results to result object
meth_analysis <- rbind(meth_analysis,c(type,sex,mean(dat_boxplot$y[dat_boxplot$DC=="D"]),
                                       mean(dat_boxplot$y[dat_boxplot$DC=="ND"]),
                                       summary(aov(y~DC,dat_boxplot))[[1]][["Pr(>F)"]][1]))

# Female
sex <- "F"

# Make dataframe for boxplot and anova
dat_boxplot <- as.data.frame(cbind(c(rowMeans(meth_dosage_F),rowMeans(meth_nodosage_F)),
                                   c(rep("D",nrow(meth_dosage_F)),rep("ND",nrow(meth_nodosage_F)))))
colnames(dat_boxplot) <- c("y","DC")
dat_boxplot$y <- as.numeric(as.character(dat_boxplot$y))

# Make boxplot
boxplot(y~DC,dat_boxplot)

# Anova
summary(aov(y~DC,dat_boxplot))
mean(dat_boxplot$y[dat_boxplot$DC=="D"])
mean(dat_boxplot$y[dat_boxplot$DC=="ND"])

# Add results to result object
meth_analysis <- rbind(meth_analysis,c(type,sex,mean(dat_boxplot$y[dat_boxplot$DC=="D"]),
                                       mean(dat_boxplot$y[dat_boxplot$DC=="ND"]),
                                       summary(aov(y~DC,dat_boxplot))[[1]][["Pr(>F)"]][1]))


# GENIC #
#########

type <- "genic"

# Make subsets
meth_dosage_M <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_dosage) & genic),f2=="M"]
meth_dosage_F <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_dosage) & genic),f2=="F"]
meth_nodosage_M <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_nodosage) & genic),f2=="M"]
meth_nodosage_F <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_nodosage) & genic),f2=="F"]

# Male
sex <- "M"

# Make dataframe for boxplot and anova
dat_boxplot <- as.data.frame(cbind(c(rowMeans(meth_dosage_M),rowMeans(meth_nodosage_M))
                                   ,c(rep("D",nrow(meth_dosage_M)),rep("ND",nrow(meth_nodosage_M)))))
colnames(dat_boxplot) <- c("y","DC")
dat_boxplot$y <- as.numeric(as.character(dat_boxplot$y))

# Make boxplot
boxplot(y~DC,dat_boxplot)

# Anova
summary(aov(y~DC,dat_boxplot))
print(mean(dat_boxplot$y[dat_boxplot$DC=="D"]))
print(mean(dat_boxplot$y[dat_boxplot$DC=="ND"]))

# Add results to result object
meth_analysis <- rbind(meth_analysis,c(type,sex,mean(dat_boxplot$y[dat_boxplot$DC=="D"]),
                                       mean(dat_boxplot$y[dat_boxplot$DC=="ND"]),
                                       summary(aov(y~DC,dat_boxplot))[[1]][["Pr(>F)"]][1]))

# Female
sex <- "F"

# Make dataframe for boxplot and anova
dat_boxplot <- as.data.frame(cbind(c(rowMeans(meth_dosage_F),rowMeans(meth_nodosage_F)),
                                   c(rep("D",nrow(meth_dosage_F)),rep("ND",nrow(meth_nodosage_F)))))
colnames(dat_boxplot) <- c("y","DC")
dat_boxplot$y <- as.numeric(as.character(dat_boxplot$y))

# Make boxplot
boxplot(y~DC,dat_boxplot)

# Anova
summary(aov(y~DC,dat_boxplot))
mean(dat_boxplot$y[dat_boxplot$DC=="D"])
mean(dat_boxplot$y[dat_boxplot$DC=="ND"])

# Add results to result object
meth_analysis <- rbind(meth_analysis,c(type,sex,mean(dat_boxplot$y[dat_boxplot$DC=="D"]),
                                       mean(dat_boxplot$y[dat_boxplot$DC=="ND"]),
                                       summary(aov(y~DC,dat_boxplot))[[1]][["Pr(>F)"]][1]))


# DIFF DOSAGE COMPENSATION 
############################

# ALL #
#######

type <- "both"

# Make subsets
meth_dosage_M <- methLevel(rrbs.rel)[gene[,1]%in%genelist_diffdosage,f2=="M" & f1=="adultdph"]
meth_dosage_F <- methLevel(rrbs.rel)[gene[,1]%in%genelist_diffdosage,f2=="F" & f1=="adultdph"]
meth_nodosage_M <- methLevel(rrbs.rel)[gene[,1]%in%genelist_diffdosage,f2=="M" & f1=="1dph"]
meth_nodosage_F <- methLevel(rrbs.rel)[gene[,1]%in%genelist_diffdosage,f2=="F" & f1=="1dph"]

# Male
sex <- "M"

# Make dataframe for boxplot and anova
dat_boxplot <- as.data.frame(cbind(c(rowMeans(meth_dosage_M),rowMeans(meth_nodosage_M))
                                   ,c(rep("D",nrow(meth_dosage_M)),rep("ND",nrow(meth_nodosage_M)))))
colnames(dat_boxplot) <- c("y","DC")
dat_boxplot$y <- as.numeric(as.character(dat_boxplot$y))

# Make boxplot
boxplot(y~DC,dat_boxplot)

# Anova
summary(aov(y~DC,dat_boxplot))
print(mean(dat_boxplot$y[dat_boxplot$DC=="D"]))
print(mean(dat_boxplot$y[dat_boxplot$DC=="ND"]))

# Add results to result object
meth_analysis <- rbind(meth_analysis,c(type,sex,mean(dat_boxplot$y[dat_boxplot$DC=="D"]),
                                       mean(dat_boxplot$y[dat_boxplot$DC=="ND"]),
                                       summary(aov(y~DC,dat_boxplot))[[1]][["Pr(>F)"]][1]))

# Female
sex <- "F"

# Make dataframe for boxplot and anova
dat_boxplot <- as.data.frame(cbind(c(rowMeans(meth_dosage_F),rowMeans(meth_nodosage_F)),
                                   c(rep("D",nrow(meth_dosage_F)),rep("ND",nrow(meth_nodosage_F)))))
colnames(dat_boxplot) <- c("y","DC")
dat_boxplot$y <- as.numeric(as.character(dat_boxplot$y))

# Make boxplot
boxplot(y~DC,dat_boxplot)

# Anova
summary(aov(y~DC,dat_boxplot))
mean(dat_boxplot$y[dat_boxplot$DC=="D"])
mean(dat_boxplot$y[dat_boxplot$DC=="ND"])

# Add results to result object
meth_analysis <- rbind(meth_analysis,c(type,sex,mean(dat_boxplot$y[dat_boxplot$DC=="D"]),
                                       mean(dat_boxplot$y[dat_boxplot$DC=="ND"]),
                                       summary(aov(y~DC,dat_boxplot))[[1]][["Pr(>F)"]][1]))


# PROMOTER #
############

type <- "promotor"

# Make subsets
meth_dosage_M <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_diffdosage) & promoter),f2=="M" & f1=="adultdph"]
meth_dosage_F <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_diffdosage) & promoter),f2=="F" & f1=="adultdph"]
meth_nodosage_M <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_diffdosage) & promoter),f2=="M" & f1=="1dph"]
meth_nodosage_F <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_diffdosage) & promoter),f2=="F" & f1=="1dph"]

# Male
sex <- "M"

# Make dataframe for boxplot and anova
dat_boxplot <- as.data.frame(cbind(c(rowMeans(meth_dosage_M),rowMeans(meth_nodosage_M))
                                   ,c(rep("D",nrow(meth_dosage_M)),rep("ND",nrow(meth_nodosage_M)))))
colnames(dat_boxplot) <- c("y","DC")
dat_boxplot$y <- as.numeric(as.character(dat_boxplot$y))

# Make boxplot
boxplot(y~DC,dat_boxplot)

# Anova
summary(aov(y~DC,dat_boxplot))
print(mean(dat_boxplot$y[dat_boxplot$DC=="D"]))
print(mean(dat_boxplot$y[dat_boxplot$DC=="ND"]))

# Add results to result object
meth_analysis <- rbind(meth_analysis,c(type,sex,mean(dat_boxplot$y[dat_boxplot$DC=="D"]),
                                       mean(dat_boxplot$y[dat_boxplot$DC=="ND"]),
                                       summary(aov(y~DC,dat_boxplot))[[1]][["Pr(>F)"]][1]))

# Female
sex <- "F"

# Make dataframe for boxplot and anova
dat_boxplot <- as.data.frame(cbind(c(rowMeans(meth_dosage_F),rowMeans(meth_nodosage_F)),
                                   c(rep("D",nrow(meth_dosage_F)),rep("ND",nrow(meth_nodosage_F)))))
colnames(dat_boxplot) <- c("y","DC")
dat_boxplot$y <- as.numeric(as.character(dat_boxplot$y))

# Make boxplot
boxplot(y~DC,dat_boxplot)

# Anova
summary(aov(y~DC,dat_boxplot))
mean(dat_boxplot$y[dat_boxplot$DC=="D"])
mean(dat_boxplot$y[dat_boxplot$DC=="ND"])

# Add results to result object
meth_analysis <- rbind(meth_analysis,c(type,sex,mean(dat_boxplot$y[dat_boxplot$DC=="D"]),
                                       mean(dat_boxplot$y[dat_boxplot$DC=="ND"]),
                                       summary(aov(y~DC,dat_boxplot))[[1]][["Pr(>F)"]][1]))


# GENIC #
#########

type <- "genic"

# Make subsets
meth_dosage_M <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_diffdosage) & genic),f2=="M" & f1=="adultdph"]
meth_dosage_F <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_diffdosage) & genic),f2=="F" & f1=="adultdph"]
meth_nodosage_M <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_diffdosage) & genic),f2=="M" & f1=="1dph"]
meth_nodosage_F <- methLevel(rrbs.rel)[((gene[,1]%in%genelist_diffdosage) & genic),f2=="F" & f1=="1dph"]

# Male
sex <- "M"

# Make dataframe for boxplot and anova
dat_boxplot <- as.data.frame(cbind(c(rowMeans(meth_dosage_M),rowMeans(meth_nodosage_M))
                                   ,c(rep("D",nrow(meth_dosage_M)),rep("ND",nrow(meth_nodosage_M)))))
colnames(dat_boxplot) <- c("y","DC")
dat_boxplot$y <- as.numeric(as.character(dat_boxplot$y))

# Make boxplot
boxplot(y~DC,dat_boxplot)

# Anova
summary(aov(y~DC,dat_boxplot))
print(mean(dat_boxplot$y[dat_boxplot$DC=="D"]))
print(mean(dat_boxplot$y[dat_boxplot$DC=="ND"]))

# Add results to result object
meth_analysis <- rbind(meth_analysis,c(type,sex,mean(dat_boxplot$y[dat_boxplot$DC=="D"]),
                                       mean(dat_boxplot$y[dat_boxplot$DC=="ND"]),
                                       summary(aov(y~DC,dat_boxplot))[[1]][["Pr(>F)"]][1]))

# Female
sex <- "F"

# Make dataframe for boxplot and anova
dat_boxplot <- as.data.frame(cbind(c(rowMeans(meth_dosage_F),rowMeans(meth_nodosage_F)),
                                   c(rep("D",nrow(meth_dosage_F)),rep("ND",nrow(meth_nodosage_F)))))
colnames(dat_boxplot) <- c("y","DC")
dat_boxplot$y <- as.numeric(as.character(dat_boxplot$y))

# Make boxplot
boxplot(y~DC,dat_boxplot)

# Anova
summary(aov(y~DC,dat_boxplot))
mean(dat_boxplot$y[dat_boxplot$DC=="D"])
mean(dat_boxplot$y[dat_boxplot$DC=="ND"])

# Add results to result object
meth_analysis <- rbind(meth_analysis,c(type,sex,mean(dat_boxplot$y[dat_boxplot$DC=="D"]),
                                       mean(dat_boxplot$y[dat_boxplot$DC=="ND"]),
                                       summary(aov(y~DC,dat_boxplot))[[1]][["Pr(>F)"]][1]))


# Write results object
colnames(meth_analysis) <- c("CpGs in","sex","meanD","meanND","pvalue")
write.table(meth_analysis,file="dosage_effect_meth.txt",col.names = T, row.names = F, sep="\t", quote = F)


#####???

fulldata_output <- read.table("FullData_output.txt",header=T,sep="\t")

getZorA <- function(genelist){
  # print(head(genelist))
  ix <- rank(as.character(genelist))
  genelist_sorted <- sort(as.character(genelist))
  chr_sorted <- genes$chromosome_name[genes$ensembl_gene_id%in%genelist_sorted]
  # print(head(genelist_sorted[ix]))
  chr <- chr_sorted[ix]
  return(chr)
}

## logFC are M/F
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

length(genes_male)
length(genes_female)
length(genes_male_A)
length(genes_female_A)
length(genes_male_Z)
length(genes_female_Z)

genes_male_Z <- genes_male_Z[!(genes_male_Z%in%genes_female_Z)]

write.table(genes_male_Z,"GeneLists/ListGenes_MaleOver.txt",col.names = F,row.names = F)
write.table(genes_female_Z,"GeneLists/ListGenes_FemaleOver.txt",col.names = F,row.names = F)


## NEWNEWNEWNEW
dosage <- read.table("/home/louisc/Zebra_Finch_new/GeneLists/ListGenes_dosage.txt",header=F,stringsAsFactors=F)$V1
dosage_bis <- read.table("/home/louisc/Zebra_Finch_new/GeneLists/ListGenes_dosage_bis.txt",header=F,stringsAsFactors=F)$V1
nodosage <- read.table("/home/louisc/Zebra_Finch_new/GeneLists/ListGenes_nodosage.txt",header=F,stringsAsFactors=F)$V1

fulldata_output <- read.table("FullData_output.txt",header=T,sep="\t")

annot <- read.table("JASPAR_annot.txt",header=T,sep="\t",stringsAsFactors=F)
head(annot)

rates <- cbind(log(fulldata_output$mean.20dph.F/fulldata_output$mean.1dph.F),
               log(fulldata_output$mean.65dph.F/fulldata_output$mean.20dph.F),
               log(fulldata_output$mean.adult.F/fulldata_output$mean.65dph.F),
               log(fulldata_output$mean.20dph.M/fulldata_output$mean.1dph.M),
               log(fulldata_output$mean.65dph.M/fulldata_output$mean.20dph.M),
               log(fulldata_output$mean.adult.M/fulldata_output$mean.65dph.M))

coi <- c("dosage","dosage_bis","nodosage")

for (i in 1:length(coi)){
  rates_part <- rates[fulldata_output$ensembl.gene.id%in%get(coi[i]),]
  jpeg(paste("test_",coi[i],".jpg",sep=""))
  plot(1:3,rep(0,3),ylim=c(-0.5,0.5),"l",col="black",xaxt="n")
  axis(side=1, at=1:3, labels = c("log(20dph/1dph)","log(65dph/20dph)","log(adult/65dph)"))
  lines(1:3,colMeans(rates_part)[1:3],col="blue")
  lines(1:3,colMeans(rates_part)[4:6],col="red")
  legend("topright",c("Female","Male"),col=c("blue","red"),pch=16)
}

 

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

annot_y <- rna_annot[keep,]

titles <- c("Acquired Dosage Compensation", "Dosage Compensation From Birth", "No Dosage Compensation")
par(mfrow=c(1,3))
for (i in 1:length(coi)){
  means <- cbind(mean(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="1dph" & f2=="F"]+0.01))),
                 mean(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="20dph" & f2=="F"]+0.01))),
                 mean(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="65dph" & f2=="F"]+0.01))),
                 mean(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="adult" & f2=="F"]+0.01))),
                 mean(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="1dph" & f2=="M"]+0.01))),
                 mean(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="20dph" & f2=="M"]+0.01))),
                 mean(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="65dph" & f2=="M"]+0.01))),
                 mean(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="adult" & f2=="M"]+0.01))))
  
  SEM <- cbind(sd(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="1dph" & f2=="F"]+0.01))),
               sd(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="20dph" & f2=="F"]+0.01))),
               sd(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="65dph" & f2=="F"]+0.01))),
               sd(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="adult" & f2=="F"]+0.01))),
               sd(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="1dph" & f2=="M"]+0.01))),
               sd(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="20dph" & f2=="M"]+0.01))),
               sd(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="65dph" & f2=="M"]+0.01))),
               sd(colMeans(log(y$counts[annot_y%in%get(coi[i]),f1=="adult" & f2=="M"]+0.01))))
  
  rates_part <- rates[fulldata_output$ensembl.gene.id%in%get(coi[i]),]
  
  tiff(paste("DosageCompensation_panel",i,".tiff",sep=""))
  plot(1:4,rep(0,4),ylim=c(2,7),xaxt="n",type="n",xlab="Age",ylab="Log(expression)",main=titles[i])
  axis(side=1, at=1:4, labels = c("1dph","20dph","65dph","adult"))
  lines(1:4,means[1:4],col="blue")
  arrows(1:4,means[1:4]-SEM[1:4],1:4,means[1:4]+SEM[1:4],length=0.05,angle=90,code=3,col="blue")
  lines(1:4,means[5:8],col="red")
  arrows(1:4,means[5:8]-SEM[5:8],1:4,means[5:8]+SEM[5:8],length=0.05,angle=90,code=3,col="red")
  legend("topright",c("Female","Male"),col=c("blue","red"),pch=16)      
  dev.off()
}

dosage_up <- dosage[rates_part[,3]>0]
dosage_down <- dosage[rates_part[,3]<0]
coi <- c("dosage_up","dosage_down")

TF_genes <- read.table("test.txt",header=F,stringsAsFactors = F)$V1
old <- which(f1=="65dph" | f1=="adultdph")
for (j in 1:length(coi)){
  print(coi[j])
  if (sum(gene%in%get(coi[j]))==0){
    print("No CpG probes available")
    next
    }
  counter <- 1
  meth_colmeans_nodosage <- NULL
  meth_bis <- NULL
  for (i in 1:length(get(coi[j]))){
    if (sum((gene[,1]==get(coi[j])[i]) & (!is.na(gene[,1])))==0){
      next
    } else if (sum((gene[,1]==get(coi[j])[i]) & (!is.na(gene[,1])))==1){
      meth_bis <- rbind(meth_bis,methLevel(rawToRel(rrbsredB))[(gene[,1]==get(coi[j])[i]) & (!is.na(gene[,1])),old])
      meth_colmeans_nodosage <- rbind(meth_colmeans_nodosage,cbind(methLevel(rawToRel(rrbsredB))[(gene[,1]==get(coi[j])[i]) & (!is.na(gene[,1])),old],as.character(f1)[old],as.character(f2)[old],gene[i,1]))
      points(1:12,methLevel(rawToRel(rrbsredB))[(gene[,1]==get(coi[j])[i]) & (!is.na(gene[,1])),sort(as.numeric(annot[,1]),index.return=T)$ix[old]])
    } else {
      meth_bis <- rbind(meth_bis,methLevel(rawToRel(rrbsredB))[(gene[,1]==get(coi[j])[i]) & (!is.na(gene[,1])),old])
      meth_colmeans_nodosage <- rbind(meth_colmeans_nodosage,cbind(colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==get(coi[j])[i]) & (!is.na(gene[,1])),old]),as.character(f1)[old],as.character(f2)[old],gene[i,1]))
      if (counter == 1){
        plot(1:12,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==get(coi[j])[i]) & (!is.na(gene[,1])),sort(as.numeric(annot[,1]),index.return=T)$ix[old]]),ylim=c(0,1))
        counter <- counter +1
      } else {
        points(1:12,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==get(coi[j])[i]) & (!is.na(gene[,1])),sort(as.numeric(annot[,1]),index.return=T)$ix[old]]))
      }
    }
  }
  meth_colmeans_nodosage <- as.data.frame(meth_colmeans_nodosage)
  colnames(meth_colmeans_nodosage) <- c("meth","f1","f2","gene")
  meth_colmeans_nodosage$meth <- as.numeric(as.character(meth_colmeans_nodosage$meth))
  meth_colmeans_nodosage$TF <- meth_colmeans_nodosage$gene%in%TF_genes
  meth_colmeans_nodosage <- meth_colmeans_nodosage[meth_colmeans_nodosage$TF,]
  anova_nodosage <- aov(meth~f1*f2,data=meth_colmeans_nodosage)
  print(summary(anova_nodosage))
  print(TukeyHSD(anova_nodosage)[2:4])
  car::leveneTest(meth~f1*f2,data=meth_colmeans_nodosage)
  par(mfrow=c(2,2))
  plot(anova_nodosage)
  par(mfrow=c(1,1))
}


dosage_meanlogFC <- NULL
for (i in 1:length(dosage)){
  if (nrow(RRBS_adult_full[RRBS_adult_full$gene%in%dosage[i] & RRBS_adult_full$func=="promoter",])==0){
    dosage_meanlogFC <- c(dosage_meanlogFC,0)
  } else if (nrow(RRBS_adult_full[RRBS_adult_full$gene%in%dosage[i] & RRBS_adult_full$func=="promoter",])==1){
    dosage_meanlogFC <- c(dosage_meanlogFC,RRBS_adult_full[RRBS_adult_full$gene%in%dosage[i] & RRBS_adult_full$func=="promoter",]$logFC)
  } else {
    dosage_meanlogFC <- c(dosage_meanlogFC,mean(RRBS_adult_full[RRBS_adult_full$gene%in%dosage[i] & RRBS_adult_full$func=="promoter",]$logFC))
  }
  
}

sum(rates_part[,3]>0 & dosage%in%TF_genes & dosage_meanlogFC>0)
sum(rates_part[,3]>0 & dosage%in%TF_genes & dosage_meanlogFC<0)
sum(rates_part[,3]>0 & dosage%in%TF_genes & dosage_meanlogFC>0)
sum(rates_part[,3]>0 & dosage_meanlogFC<0)
