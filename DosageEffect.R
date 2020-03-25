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

working_dir <- getwd()

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

rnadata <- read.table("Plastoscine.txt",header=T,sep="\t",quote="")
rna_annot <- read.table("Plastoscine_annot.txt",header=F,sep="\t",quote="",stringsAsFactors = F)$V1

out2excel <- rnadata
rownames(out2excel) <- rna_annot
write.xlsx(out2excel,"RNAseq_Raw.xlsx")

# Load in ensembl annotation
if (!file.exists("Genes.Rda")){
  ensembl=useMart("ensembl",host="http://sep2019.archive.ensembl.org")
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

MF_A_norand
MF_Z_norand


## test with random annotation 
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

MF_A
MF_Z


## Plot dosage effect 
###################################

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
per_Z <- CalcPercentage(ratio_Z_1dph,1000,0.1)
per_A <- CalcPercentage(ratio_A_1dph,1000,0.1)
tiff(paste("MoverFratios_panel_1.tiff",sep=""))
plot(per_A$x,per_A$y,type="l",main="1dph",ylab="Percentage of genes",xlab="M:F ratio",xlim=c(0,3),ylim=c(0.1,5.5))
lines(per_Z$x,per_Z$y,col="red")
legend("topright",c(paste("A:",round(MF_A[1],digits=2)),paste("Z:",round(MF_Z[1],digits=2))))
dev.off()
per_Z <- CalcPercentage(ratio_Z_20dph,1000,0.1)
per_A <- CalcPercentage(ratio_A_20dph,1000,0.1)
tiff(paste("MoverFratios_panel_2.tiff",sep=""))
plot(per_A$x,per_A$y,type="l",main="20dph",ylab="Percentage of genes",xlab="M:F ratio",xlim=c(0,3),ylim=c(0.1,5.5))
lines(per_Z$x,per_Z$y,col="red")
legend("topright",c(paste("A:",round(MF_A[2],digits=2)),paste("Z:",round(MF_Z[2],digits=2))))
dev.off()
per_Z <- CalcPercentage(ratio_Z_65dph,1000,0.1)
per_A <- CalcPercentage(ratio_A_65dph,1000,0.1)
tiff(paste("MoverFratios_panel_3.tiff",sep=""))
plot(per_A$x,per_A$y,type="l",main="65dph",ylab="Percentage of genes",xlab="M:F ratio",xlim=c(0,3),ylim=c(0.1,5.5))
lines(per_Z$x,per_Z$y,col="red")
legend("topright",c(paste("A:",round(MF_A[3],digits=2)),paste("Z:",round(MF_Z[3],digits=2))))
dev.off()
per_Z <- CalcPercentage(ratio_Z_adult,1000,0.1)
per_A <- CalcPercentage(ratio_A_adult,1000,0.1)
tiff(paste("MoverFratios_panel_4.tiff",sep=""))
plot(per_A$x,per_A$y,type="l",main="adult",ylab="Percentage of genes",xlab="M:F ratio",xlim=c(0,3),ylim=c(0.1,5.5))
lines(per_Z$x,per_Z$y,col="red")
legend("topright",c(paste("A:",round(MF_A[4],digits=2)),paste("Z:",round(MF_Z[4],digits=2))))
dev.off()


## Make Gene lists for further analysis
###################################

cte <- 0.1

# Calculate ratio's for Autosomal genes
logcpm_A_M <- c(rowMeans(log_cpm[genes_ZvsA=="A",f2=="M" & f1=="1dph"]),
                rowMeans(log_cpm[genes_ZvsA=="A",f2=="M" & f1=="20dph"]),
                rowMeans(log_cpm[genes_ZvsA=="A",f2=="M" & f1=="65dph"]),
                rowMeans(log_cpm[genes_ZvsA=="A",f2=="M" & f1=="adult"]))
logcpm_A_F <- c(rowMeans(log_cpm[genes_ZvsA=="A",f2=="F" & f1=="1dph"]),
                rowMeans(log_cpm[genes_ZvsA=="A",f2=="F" & f1=="20dph"]),
                rowMeans(log_cpm[genes_ZvsA=="A",f2=="F" & f1=="65dph"]),
                rowMeans(log_cpm[genes_ZvsA=="A",f2=="F" & f1=="adult"]))
logcpm_Z_M <- c(rowMeans(log_cpm[genes_ZvsA=="Z",f2=="M" & f1=="1dph"]),
                rowMeans(log_cpm[genes_ZvsA=="Z",f2=="M" & f1=="20dph"]),
                rowMeans(log_cpm[genes_ZvsA=="Z",f2=="M" & f1=="65dph"]),
                rowMeans(log_cpm[genes_ZvsA=="Z",f2=="M" & f1=="adult"]))
logcpm_Z_F <- c(rowMeans(log_cpm[genes_ZvsA=="Z",f2=="F" & f1=="1dph"]),
                rowMeans(log_cpm[genes_ZvsA=="Z",f2=="F" & f1=="20dph"]),
                rowMeans(log_cpm[genes_ZvsA=="Z",f2=="F" & f1=="65dph"]),
                rowMeans(log_cpm[genes_ZvsA=="Z",f2=="F" & f1=="adult"]))

# Calculate 95% interval for Z-related ratios
quantile((logcpm_Z_M)/(logcpm_Z_F),probs=c(0.025,0.975),na.rm=T)

# Calculate 95% interval for Autosomal ratios: borders will be used as cutoffs for "dosage compensated" genes
cutoffs <- quantile((logcpm_A_M+cte)/(logcpm_A_F+cte),probs=c(0.05,0.95,0.975),na.rm=T)
print(cutoffs)

# Make genelists "diffdosage": meaning not differentially no dosage compensation on 1dph but dosage compensation for adults
genelist_dosage <- c(rna_annot_Z[(ratio_Z_1dph>cutoffs[2]) & (ratio_Z_adult<cutoffs[2])])
write.table(genelist_dosage,file="GeneLists/ListGenes_dosage.txt",col.names = F,row.names = F,quote = F)
length(genelist_dosage)
genelist_nodosage <- c(rna_annot_Z[(ratio_Z_1dph>cutoffs[2]) & (ratio_Z_adult>cutoffs[2])])
write.table(genelist_nodosage,file="GeneLists/ListGenes_nodosage.txt",col.names = F,row.names = F,quote = F)
length(genelist_nodosage)
genelist_dosage_bis <- c(rna_annot_Z[(ratio_Z_1dph<cutoffs[2]) & (ratio_Z_adult<cutoffs[2])])
write.table(genelist_dosage_bis,file="GeneLists/ListGenes_dosage_bis.txt",col.names = F,row.names = F,quote = F)
length(genelist_dosage_bis)
length(rna_annot_Z)


## Higher expression for dosage compensated genes
#################################

# Additional data
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

# Plot
coi <- c("dosage","dosage_bis","nodosage")
titles <- c("Acquired Dosage Compensation", "Dosage Compensation From Birth", "No Dosage Compensation")
par(mfrow=c(1,3))
for (i in 1:length(coi)){
  goi <- get(paste("genelist",coi[i],sep="_"))
  means <- cbind(mean(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="1dph" & f2=="F"])),
                 mean(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="20dph" & f2=="F"])),
                 mean(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="65dph" & f2=="F"])),
                 mean(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="adult" & f2=="F"])),
                 mean(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="1dph" & f2=="M"])),
                 mean(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="20dph" & f2=="M"])),
                 mean(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="65dph" & f2=="M"])),
                 mean(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="adult" & f2=="M"])))
  
  SEM <- cbind(sd(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="1dph" & f2=="F"])),
                 sd(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="20dph" & f2=="F"])),
                 sd(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="65dph" & f2=="F"])),
                 sd(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="adult" & f2=="F"])),
                 sd(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="1dph" & f2=="M"])),
                 sd(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="20dph" & f2=="M"])),
                 sd(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="65dph" & f2=="M"])),
                 sd(colMeans(cpm(y,log=T)[annot_y%in%goi,f1=="adult" & f2=="M"])))
                 
                 # mean(colMeans(log(y$counts[annot_y%in%goi,f1=="20dph" & f2=="F"]+0.01))),
                 # mean(colMeans(log(y$counts[annot_y%in%goi,f1=="65dph" & f2=="F"]+0.01))),
                 # mean(colMeans(log(y$counts[annot_y%in%goi,f1=="adult" & f2=="F"]+0.01))),
                 # mean(colMeans(log(y$counts[annot_y%in%goi,f1=="1dph" & f2=="M"]+0.01))),
                 # mean(colMeans(log(y$counts[annot_y%in%goi,f1=="20dph" & f2=="M"]+0.01))),
                 # mean(colMeans(log(y$counts[annot_y%in%goi,f1=="65dph" & f2=="M"]+0.01))),
                 # mean(colMeans(log(y$counts[annot_y%in%goi,f1=="adult" & f2=="M"]+0.01))))
  
  # SEM <- cbind(sd(colMeans(log(y$counts[annot_y%in%goi,f1=="1dph" & f2=="F"]+0.01))),
  #              sd(colMeans(log(y$counts[annot_y%in%goi,f1=="20dph" & f2=="F"]+0.01))),
  #              sd(colMeans(log(y$counts[annot_y%in%goi,f1=="65dph" & f2=="F"]+0.01))),
  #              sd(colMeans(log(y$counts[annot_y%in%goi,f1=="adult" & f2=="F"]+0.01))),
  #              sd(colMeans(log(y$counts[annot_y%in%goi,f1=="1dph" & f2=="M"]+0.01))),
  #              sd(colMeans(log(y$counts[annot_y%in%goi,f1=="20dph" & f2=="M"]+0.01))),
  #              sd(colMeans(log(y$counts[annot_y%in%goi,f1=="65dph" & f2=="M"]+0.01))),
  #              sd(colMeans(log(y$counts[annot_y%in%goi,f1=="adult" & f2=="M"]+0.01))))
  
  tiff(paste("DosageCompensation_panel",i,".tiff",sep=""))
  plot(1:4,rep(0,4),ylim=c(2.5,5.5),xaxt="n",type="n",xlab="Age",ylab="Log(expression)",main=titles[i])
  axis(side=1, at=1:4, labels = c("1dph","20dph","65dph","adult"))
  lines(1:4,means[1:4],col="blue")
  arrows(1:4,means[1:4]-SEM[1:4],1:4,means[1:4]+SEM[1:4],length=0.05,angle=90,code=3,col="blue")
  lines(1:4,means[5:8],col="red")
  arrows(1:4,means[5:8]-SEM[5:8],1:4,means[5:8]+SEM[5:8],length=0.05,angle=90,code=3,col="red")
  legend("topright",c("Female","Male"),col=c("blue","red"),pch=16)      
  dev.off()
}

# Statistical test

logcpm <- rowMeans(cpm(y,log=T))
rna_annot <- rna_annot$V1[keep]

mat <- NULL
mat <- rbind(mat,cbind(logcpm[rna_annot%in%genelist_dosage],rep("dosage",sum(rna_annot%in%genelist_dosage))))
mat <- rbind(mat,cbind(logcpm[rna_annot%in%genelist_nodosage],rep("nodosage",sum(rna_annot%in%genelist_nodosage))))
mat <- rbind(mat,cbind(logcpm[rna_annot%in%genelist_dosage_bis],rep("dosage_bis",sum(rna_annot%in%genelist_dosage_bis))))

mat <- data.frame(mat,stringsAsFactors = T)
colnames(mat) <- c("logcpm","group")
mat$logcpm <- as.numeric(as.character(mat$logcpm))
print(summary(aov(logcpm~group,mat)))
print(TukeyHSD(aov(logcpm~group,mat)))


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

if (!file.exists("Genes_RRBS.Rda")){
  ensembl=useMart("ensembl",host="http://sep2019.archive.ensembl.org")
  ensembl = useDataset("tguttata_gene_ensembl",mart=ensembl)
  genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position','transcription_start_site'),mart = ensembl)
  colnames(genes) <- c("gene","gene_symbol","chr","start","end","TSS")
  print(dim(genes))
  genes <- genes[!duplicated(genes$ensembl_gene_id),]
  print(dim(genes))
  save(genes,file="Genes.Rda")
} else {
  load(file="Genes_RRBS.Rda")
  genes <- genes[!duplicated(genes$ensembl_gene_id),]
  print(dim(genes))
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

make_submat <- function(dat){
  res <- cbind(dat,as.character(f1),as.character(f2),i,1:24)
  return(data.frame(res))
}


## Promoter

for (t in 1:length(coi)){
  print(coi[t])

  counter <- 1
  meth_colmeans <- NULL
  genelist <- get(paste("genelist",coi[t],sep="_"))
  
  for (i in 1:length(genelist)){
    if (sum((gene[,1]==genelist[i]) & (!is.na(gene[,1])) & promoter)==0){
      next
    } else if (sum((gene[,1]==genelist[i]) & (!is.na(gene[,1])) & promoter)==1){
      meth_colmeans <- rbind(meth_colmeans,cbind(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])) & promoter,],as.character(f1),as.character(f2),i,1:24))
      if (counter == 1){
        jpeg(paste(coi[t],"_methylation.jpg",sep=""))
        plot(1:24,methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])),][sort(as.numeric(annot[,1]),index.return=T)$ix],ylim=c(0,1))
        counter <- counter +1
      } else {
        points(1:24,methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])),][sort(as.numeric(annot[,1]),index.return=T)$ix])
      }
    } else {
      # meth_colmeans <- rbind(meth_colmeans,as.matrix(do.call(rbind,apply(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])) & promoter,],1,make_submat))))
      meth_colmeans <- rbind(meth_colmeans,cbind(colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])) & promoter,]),as.character(f1),as.character(f2),i,1:24))
      if (counter == 1){
        jpeg(paste(coi[t],"_methylation.jpg",sep=""))
        plot(1:24,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])) & promoter,sort(as.numeric(annot[,1]),index.return=T)$ix]),ylim=c(0,1))
        counter <- counter +1
      } else {
        points(1:24,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])) & promoter,sort(as.numeric(annot[,1]),index.return=T)$ix]))
      }
    }
  }
  
  if(is.null(nrow(meth_colmeans))){
    print("No datapoints available")
  } else {
    dev.off()
    meth_colmeans <- as.data.frame(meth_colmeans)
    colnames(meth_colmeans) <- c("meth","f1","f2","gene","individual")
    meth_colmeans$meth <- as.numeric(as.character(meth_colmeans$meth))
    print(mean(meth_colmeans$meth))
    print(median(meth_colmeans$meth))
    anova_dosage <- aov(meth~gene+f1*f2,data=meth_colmeans)
    print(anova_dosage)
    print(summary(anova_dosage))
    print(TukeyHSD(anova_dosage)$'f1')
    print(car::leveneTest(meth~f1*f2,data=meth_colmeans))
    par(mfrow=c(2,2))
    jpeg(paste("anova_",coi[t],".jpg",sep=""))
    plot(anova_dosage)
    dev.off()
    par(mfrow=c(1,1))
  }
}


## Genic

for (t in 1:length(coi)){
  print(coi[t])
  
  counter <- 1
  meth_colmeans <- NULL
  genelist <- get(paste("genelist",coi[t],sep="_"))
  
  for (i in 1:length(genelist)){
    if (sum((gene[,1]==genelist[i]) & (!is.na(gene[,1])) & genic)==0){
      next
    } else if (sum((gene[,1]==genelist[i]) & (!is.na(gene[,1])) & genic)==1){
      meth_colmeans <- rbind(meth_colmeans,cbind(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])) & genic,],as.character(f1),as.character(f2),i,1:24))
      if (counter == 1){
        jpeg(paste(coi[t],"_methylation.jpg",sep=""))
        plot(1:24,methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])),][sort(as.numeric(annot[,1]),index.return=T)$ix],ylim=c(0,1))
        counter <- counter +1
      } else {
        points(1:24,methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])),][sort(as.numeric(annot[,1]),index.return=T)$ix])
      }
    } else {
      # meth_colmeans <- rbind(meth_colmeans,as.matrix(do.call(rbind,apply(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])) & genic,],1,make_submat))))
      meth_colmeans <- rbind(meth_colmeans,cbind(colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])) & genic,]),as.character(f1),as.character(f2),i,1:24))
      if (counter == 1){
        jpeg(paste(coi[t],"_methylation.jpg",sep=""))
        plot(1:24,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])) & genic,sort(as.numeric(annot[,1]),index.return=T)$ix]),ylim=c(0,1))
        counter <- counter +1
      } else {
        points(1:24,colMeans(methLevel(rawToRel(rrbsredB))[(gene[,1]==genelist[i]) & (!is.na(gene[,1])) & genic,sort(as.numeric(annot[,1]),index.return=T)$ix]))
      }
    }
  }
  
  if(is.null(nrow(meth_colmeans))){
    print("No datapoints available")
  } else {
    dev.off()
    meth_colmeans <- as.data.frame(meth_colmeans)
    colnames(meth_colmeans) <- c("meth","f1","f2","gene","individual")
    meth_colmeans$meth <- as.numeric(as.character(meth_colmeans$meth))
    print(mean(meth_colmeans$meth))
    print(median(meth_colmeans$meth))
    anova_dosage <- aov(meth~gene+f1*f2,data=meth_colmeans)
    print(anova_dosage)
    print(summary(anova_dosage))
    print(TukeyHSD(anova_dosage)$'f1')
    print(car::leveneTest(meth~f1*f2,data=meth_colmeans))
    par(mfrow=c(2,2))
    jpeg(paste("anova_",coi[t],".jpg",sep=""))
    plot(anova_dosage)
    dev.off()
    par(mfrow=c(1,1))
  }
}


