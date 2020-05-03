############################
############################
##                        ##
## RRBS_analyse_hydroxy.R ##
##                        ##
############################
############################


## Parameter settings 
#####################

home_dir <- getwd()
data_dir <- "/data/sandras/BiSeq"
options(warn=2L)
mincov<-5 
minmethreads<-6
ct<-0.01 #constant value for calculation M-values from methylation degrees
quartzFonts(sans=c("Arial","Arial Black","Arial Narrow","Arial Rounded MT Bold"))

## load packages
################

library(BiSeq)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(limma)
library(biomaRt)
library(UniProt.ws)
library(RColorBrewer)
library(gplots)

sink(file='RRBS_analyse_hydroxy.log', type = c("output", "message"))

setwd(data_dir)


## Import data 
##############

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

annot<-c("Pool","Sex","Age","8","F","20","12","M","20","14","F","65","18","M","65","1","F","1","7","F","20","10","M","20","24","M","adult","2","F","1","15","F","65","19","F","adult","22","M","adult","9","F","20","11","M","20","16","M","65","17","M","65","4","M","1","5","M","1","13","F","65","20","F","adult","3","F","1","6","M","1","21","F","adult","23","M","adult")
annot<-matrix(annot,nrow=25,ncol=3,byrow=T)
colnames(annot)<-annot[1,]
annot<-annot[-1,]
rownames(annot)<-annot[,1]

# read RRBS data given colData object
allrrbs_B<-readBismark(poolfiles,colData=annotN)

poolfiles<-dir()
poolfiles<-poolfiles[grepl("\\.O",poolfiles)]
# Mix up: labels switched
poolfiles[grepl("\\.O16",poolfiles)] <- dir()[grepl("\\.B16",dir())]
pids<-unlist(strsplit(substr(poolfiles,29,31),split="_"))
annot<-cbind(annot[substr(pids,2,3),],substr(pids,1,1))
annot<-cbind(annot,apply(cbind(annot[,"Sex"],annot[,"Age"]),1,paste,collapse="_"))
rownames(annot)<-pids
colnames(annot)[4:5]<-c("Ox_vs_RRBS","Fenolabel")

annotN<-DataFrame(annot)

# read RRBS data given colData object
allrrbs_O <-readBismark(poolfiles,colData=annotN)

setwd(home_dir)
locs_B<-rowRanges(allrrbs_B)
locs_O<-rowRanges(allrrbs_O)

print(head(locs_B))
print(head(locs_O))

setwd(home_dir)

## Data filtering 
##################

#Filtering based on shared regions (should be shared with coverage >= mincov in all samples)
print(dim(allrrbs_B))
rrbs_B<-filterBySharedRegions(object=allrrbs_B,perc.samples=1,minCov=mincov)
print(dim(rrbs_B))

inclocsM2<-names(which(rowSums(methReads(rrbs_B))>=minmethreads))
rrbs_B<-rrbs_B[inclocsM2,]
print(dim(rrbs_B))

locs_B <- rowRanges(rrbs_B)
locs_O <- rowRanges(allrrbs_O)

print(head(paste(locs_O@seqnames,locs_O@ranges@start,sep=",")))
print(head(paste(locs_B@seqnames,locs_B@ranges@start,sep=",")))
rrbs_O <- allrrbs_O[paste(locs_O@seqnames,locs_O@ranges@start,sep=",")%in%paste(locs_B@seqnames,locs_B@ranges@start,sep=","),]
print(dim(rrbs_B))
print(dim(rrbs_O))

rrbs_B.rel <- rawToRel(rrbs_B)
rrbs_O.rel <- rawToRel(rrbs_O)

print(head(methLevel(rrbs_B.rel)))
print(head(methLevel(rrbs_O.rel)))
print("NA")
for (t in 1:ncol(methLevel(rrbs_O.rel))){
  print(paste(colnames(methLevel(rrbs_O.rel))[t],":",sum(is.na(methLevel(rrbs_O.rel)[,t])),"out of",
              length(methLevel(rrbs_O.rel)[,14]),"are NA's! (",sum(is.na(methLevel(rrbs_O.rel)[,t]))/length(methLevel(rrbs_O.rel)[,14]),"%)"))
}

diff <- methLevel(rrbs_B.rel)-methLevel(rrbs_O.rel)


## Plot Hydroxylation distribution
##################

breaks <- seq(-1.05,1.05,0.1)


tiff(file="RRBS_hydroxy.tiff",units="in", width=4, height=4,res=300,family = "sans")
par(mar=c(4,4,2,2))
hist_diff <- hist(diff,breaks=breaks,xlab="Hydroxymethylation",cex.lab=0.6,cex.axis=0.5,cex.main=0.6)
dev.off()
print(hist_diff$breaks)
print(hist_diff$counts)


## Boxplot Methylation percentage
################################

meth_levels <- methLevel(rrbs_B.rel)-methLevel(rrbs_O.rel)
save(meth_levels,file="methLevel_hydroxy.Rda")

f1 <- paste0(annot[,3],"dph")
f2 <- annot[,2]
f1 <- as.factor(f1)
f2 <- as.factor(f2)
print(f1)
print(f2)

data_lm <- NULL
for(i in 1:ncol(methLevel(rrbs_B.rel))){
  data_lm <- rbind(data_lm,cbind(mean(methLevel(rrbs_B.rel)[,i]-methLevel(rrbs_O.rel)[,i],na.rm=TRUE),f1[i],f2[i]))
  
}
print(data_lm)
data_lm <- data_lm[!(rownames(data_lm)%in%c("O21")),]
print(data_lm)
data_lm <- as.data.frame(data_lm,stringsAsFactors=F)

colnames(data_lm) <- c("y","f1","f2")
data_lm$f3 <- 0

print(data_lm)
for (i in 1:nrow(data_lm)){
  if (data_lm$f1[i] <= 2){
    data_lm$f3[i] <- 1
  } else {
    data_lm$f3[i] <- 2
  }
}
data_lm$y <- as.numeric(data_lm$y)
data_lm$f1 <- as.factor(data_lm$f1)
data_lm$f2 <- as.factor(data_lm$f2)
data_lm$f3 <- as.factor(data_lm$f3)

data_boxplot <- data_lm[,1:3]

naampjes <- c("1dph,F","1dph,M","20dph,F","20dph,M","65dph,F","65dph,M","adult,F","adult,M")
naampjes <- naampjes[c(1,3,5,7,2,4,6,8)]
print(naampjes)
print(levels(f1))
print(levels(f2))

print(paste(rep(as.character(levels(f1)),2),rep(as.character(levels(f2)),4),sep=","))
tiff("Methylation_hydroxy_plot.tiff",width = 1000, height = 1000, units = "px")
boxplot(y~f1+f2,data_boxplot,at=c(1,3,5,7,2,4,6,8),las=2,
        names=naampjes[])
dev.off()
save(data_boxplot,file="hydroxy_means.Rda")


## Linear model differences methylation percentage
######################

lm_diffmeth <- lm(y~f1*f2,data=data_lm)

print(lm_diffmeth)
print(summary(lm_diffmeth))

anova_diffmeth <- aov(y~f1*f2,data=data_lm)
print(anova_diffmeth)
print(summary(anova_diffmeth))
anova_diffmeth <- aov(y~f1+f2,data=data_lm)
print(anova_diffmeth)
print(summary(anova_diffmeth))
anova_diffmeth <- aov(y~f1,data=data_lm)
print(anova_diffmeth)
print(summary(anova_diffmeth))
print(TukeyHSD(anova_diffmeth))