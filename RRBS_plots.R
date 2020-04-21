####################
####################
##                ##
## DosageEffect.R ##
##                ##
####################
####################


## Load libraries
################

library(BiSeq)
quartzFonts(sans=c("Arial","Arial Black","Arial Narrow","Arial Rounded MT Bold"))


## Parameter settings 
#####################

home_dir <- getwd()
data_dir <- "/data/sandras/BiSeq"
options(warn=2L)
mincov<-5 
minmethreads<-6
ct<-0.01 

## Load Data
################

load("allrrbs.Rda")
rrbs_CpG <-filterBySharedRegions(object=allrrbs,perc.samples=1,minCov=mincov)
rrbs_CpG.rel <- rawToRel(rrbs_CpG)

load("methLevel_hydroxy.Rda")

load("rrbs_nonCpG.Rda")
rownames(rrbs) <- 1:nrow(rrbs) 
rrbs_CpH <- rrbs
rrbs_CpH.rel <- rawToRel(rrbs_CpH)


## CpG
##############

# annotation
setwd(data_dir)

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

setwd(home_dir)

f1_CpG <- paste0(annot[,3],"dph")
f2_CpG <- annot[,2]
f1_CpG <- as.factor(f1_CpG)
f2_CpG <- as.factor(f2_CpG)

data_CpG <- NULL
for(i in 1:ncol(methLevel(rrbs_CpG.rel))){
  data_CpG <- rbind(data_CpG,cbind(mean(methLevel(rrbs_CpG.rel)[,i]),f1_CpG[i],f2_CpG[i]))
  
}
data_CpG <- as.data.frame(data_CpG,stringsAsFactors=T)
colnames(data_CpG) <- c("y","f1","f2")
data_CpG$y <- as.numeric(data_CpG$y)
data_CpG$f1 <- as.factor(data_CpG$f1)
data_CpG$f2 <- as.factor(data_CpG$f2)


labels <- c("1dph,F","1dph,M","20dph,F","20dph,M","65dph,F","65dph,M","adult,F","adult,M")
labels <- labels[c(1,3,5,7,2,4,6,8)]

tiff("Methylation_plot.tiff",units="in", width=4, height=4,res=300,family = "sans")
par(mar=c(4,4,2,2))
boxplot(y~f1+f2,data_CpG,at=c(1,3,5,7,2,4,6,8),las=2,
        names=labels,col=c(rep("gray",4),rep("white",4)),main="CpG methylation",ylab= "Methylation (%)",xlab="",cex.lab=0.8,cex.axis=0.7,cex.main=0.8)
legend("topleft",c("female","male"),pch=21,col="black",pt.bg=c("gray","white"),cex=0.7)
dev.off()


## Hydroxy
##############

f1_hydroxy <- f1_CpG
f2_hydroxy <- f2_CpG

data_hydroxy <- NULL
for(i in 1:ncol(meth_levels)){
  data_hydroxy <- rbind(data_hydroxy,cbind(mean(meth_levels[,i],na.rm=T),f1_hydroxy[i],f2_hydroxy[i]))
  
}
data_hydroxy <- as.data.frame(data_hydroxy,stringsAsFactors=T)
colnames(data_hydroxy) <- c("y","f1","f2")
data_hydroxy$y <- as.numeric(data_hydroxy$y)
data_hydroxy$f1 <- as.factor(data_hydroxy$f1)
data_hydroxy$f2 <- as.factor(data_hydroxy$f2)

data_hydroxy <- data_hydroxy[!(rownames(data_hydroxy)%in%c("B21")),]

labels <- c("1dph,F","1dph,M","20dph,F","20dph,M","65dph,F","65dph,M","adult,F","adult,M")
labels <- labels[c(1,3,5,7,2,4,6,8)]

tiff("Methylation_hydroxy_plot.tiff",units="in", width=4, height=4,res=300,family = "sans")
boxplot(y~f1+f2,data_hydroxy,at=c(1,3,5,7,2,4,6,8),las=2,
        names=labels,col=c(rep("gray",4),rep("white",4)),main="Hydroxymethylation",ylab= "Methylation (%)",xlab="",cex.lab=0.8,cex.axis=0.7,cex.main=0.8)
legend("topleft",c("female","male"),pch=21,col="black",pt.bg=c("gray","white"),cex=0.7)
dev.off()


## CpH
##############

# annotation
annot <- rbind(annot[-which(rownames(annot)=="O16"),],annot[which(rownames(annot)=="O16"),])

f1_CpH <- paste0(annot[,3],"dph")
f2_CpH <- annot[,2]
f1_CpH <- as.factor(f1_CpH)
f2_CpH <- as.factor(f2_CpH)


data_CpH <- NULL
for(i in 1:ncol(methLevel(rrbs_CpH.rel))){
  data_CpH <- rbind(data_CpH,cbind(mean(methLevel(rrbs_CpH.rel)[,i]),f1_CpH[i],f2_CpH[i]))
  
}
data_CpH <- as.data.frame(data_CpH,stringsAsFactors=T)
colnames(data_CpH) <- c("y","f1","f2")
data_CpH$y <- as.numeric(data_CpH$y)
data_CpH$f1 <- as.factor(data_CpH$f1)
data_CpH$f2 <- as.factor(data_CpH$f2)


labels <- c("1dph,F","1dph,M","20dph,F","20dph,M","65dph,F","65dph,M","adult,F","adult,M")
labels <- labels[c(1,3,5,7,2,4,6,8)]

tiff("Methylation_nonCpG_plot.tiff",units="in", width=4, height=4,res=300,family = "sans")
boxplot(y~f1+f2,data_CpH,at=c(1,3,5,7,2,4,6,8),las=2,
        names=labels,col=c(rep("gray",4),rep("white",4)),main="nonCpG methylation",ylab= "Methylation (%)",xlab="",cex.lab=0.8,cex.axis=0.7,cex.main=0.8)
legend("topleft",c("female","male"),pch=21,col="black",pt.bg=c("gray","white"),cex=0.7)
dev.off()

data_barplot <- cbind(c(mean(data_CpG$y[data_CpG$f1==1 & data_CpG$f2==1]),mean(data_hydroxy$y[data_hydroxy$f1==1 & data_hydroxy$f2==1]),mean(data_CpH$y[data_CpH$f1==1 & data_CpH$f2==1])),
                      c(mean(data_CpG$y[data_CpG$f1==1 & data_CpG$f2==2]),mean(data_hydroxy$y[data_hydroxy$f1==1 & data_hydroxy$f2==2]),mean(data_CpH$y[data_CpH$f1==1 & data_CpH$f2==2])),
                      c(mean(data_CpG$y[data_CpG$f1==2 & data_CpG$f2==1]),mean(data_hydroxy$y[data_hydroxy$f1==2 & data_hydroxy$f2==1]),mean(data_CpH$y[data_CpH$f1==2 & data_CpH$f2==1])),
                      c(mean(data_CpG$y[data_CpG$f1==2 & data_CpG$f2==2]),mean(data_hydroxy$y[data_hydroxy$f1==2 & data_hydroxy$f2==2]),mean(data_CpH$y[data_CpH$f1==2 & data_CpH$f2==2])),
                      c(mean(data_CpG$y[data_CpG$f1==3 & data_CpG$f2==1]),mean(data_hydroxy$y[data_hydroxy$f1==3 & data_hydroxy$f2==1]),mean(data_CpH$y[data_CpH$f1==3 & data_CpH$f2==1])),
                      c(mean(data_CpG$y[data_CpG$f1==3 & data_CpG$f2==2]),mean(data_hydroxy$y[data_hydroxy$f1==3 & data_hydroxy$f2==2]),mean(data_CpH$y[data_CpH$f1==3 & data_CpH$f2==2])),
                      c(mean(data_CpG$y[data_CpG$f1==4 & data_CpG$f2==1]),mean(data_hydroxy$y[data_hydroxy$f1==4 & data_hydroxy$f2==1]),mean(data_CpH$y[data_CpH$f1==4 & data_CpH$f2==1])),
                      c(mean(data_CpG$y[data_CpG$f1==4 & data_CpG$f2==2]),mean(data_hydroxy$y[data_hydroxy$f1==4 & data_hydroxy$f2==2]),mean(data_CpH$y[data_CpH$f1==4 & data_CpH$f2==2])))
tiff("Comparison_RRBS.tiff",units="in", width=4, height=4,res=300,family = "sans")
barplot(data_barplot,beside=T,main="Comparison",names=labels,las=2,cex=0.7,
        ylab= "Methylation (%)",xlab="",cex.lab=0.8,cex.axis=0.7,cex.main=0.8)
legend("topleft",c("CpG","hydroxy","nonCpG"),pt.bg=c("darkgrey","grey","white"),pch=21,col="black",cex=0.7)
dev.off()