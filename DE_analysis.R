###################
###################
##               ##
## DE_analysis.R ##
##               ##
###################
###################



## Load libraries
################

# library("xlsx")
library("UniProt.ws")
library("stringr")
library("edgeR")
library("tidyverse")
library("stageR")
library("biomaRt")


## Extra annotation Uniprot
###########################

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


## Ensembl annotation
############

ensembl=useMart("ensembl",host="http://jan2019.archive.ensembl.org",mirror="uswest")
ensembl = useDataset("tguttata_gene_ensembl",mart=ensembl)
genes <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id','chromosome_name'),mart = ensembl)


## Load Data & annotation
###############

rnadata <- read.table("Plastoscine.txt",header=T,sep="\t",quote="")
rna_annot <- read.table("Plastoscine_annot.txt",header=F,sep="\t",quote="")
row.names(rnadata) <- unlist(rna_annot)

f1 <- "1dph,1dph,1dph,1dph,1dph,1dph,20dph,20dph,20dph,20dph,20dph,20dph,65dph,65dph,65dph,65dph,65dph,65dph,adult,adult,adult,adult,adult,adult"
f2 <- "F,F,F,M,M,M,F,F,F,M,M,M,F,F,F,M,M,M,F,F,F,M,M,M"
f1 <- unlist(strsplit(f1,","))
f2 <- unlist(strsplit(f2,","))


## Filtering
###############

y <- DGEList(rnadata)
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep,  ,keep.lib.sizes=FALSE]
summary(keep)


## Normalisation & design
###############

# Make DGElist objects for differential methylation analysis
y <- calcNormFactors(y)

design <- model.matrix(~f1*f2)
  
y <- estimateDisp(y, design, robust=TRUE)  
jpeg("MDSplot.jpg")
plotMDS(y)
dev.off()
jpeg("BCV.jpg")
plotBCV(y)
dev.off()


## DE analysis (omnibus test)
###############

fit = glmFit(y,design)
jpeg("GOF.jpg")
gof(fit,plot=TRUE)
dev.off()

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

L
omnibus<-glmLRT(fit,contrast =L)
summary(dtOmnibus<-decideTests(omnibus))


## DE analysis (Post-hoc testing)
###############

# Calculation adjusted alpha
sigTests <- summary(dtOmnibus<-decideTests(omnibus))
alphaS1 <- (0.05 / sum(sigTests)) * (sigTests[2])  
alphaS1

# Stagewise testing
dph1=glmLRT(fit,contrast=L[,1])
summary(dtdph1<- decideTestsDGE(dph1))
dph20=glmLRT(fit,contrast=L[,2])
summary(dtdph20<- decideTestsDGE(dph20))
dph65=glmLRT(fit,contrast=L[,3])
summary(dtdph65<- decideTestsDGE(dph65))
dphad=glmLRT(fit,contrast=L[,4])
summary(dtdphad<- decideTestsDGE(dphad))
dphage=glmLRT(fit,contrast=L[,5])
summary(dtdphage<- decideTestsDGE(dphage))

# Stagewise confirmation
pConfirmation<-cbind(dph1$table$PValue,dph20$table$PValue,dph65$table$PValue,dphad$table$PValue,dphage$table$PValue)
fcConfirmation<-cbind(dph1$table$logFC,dph20$table$logFC,dph65$table$logFC,dphad$table$logFC,dphage$table$logFC)
colnames(pConfirmation)<- colnames(fcConfirmation) <- colnames(L)
rownames(pConfirmation) <- rownames(fcConfirmation) <- rownames(dph1)

library(stageR)
pScreen <- omnibus$table$PValue
names(pScreen) <-rownames(omnibus)
stageRObj <- stageR(pScreen=pScreen,pConfirmation=pConfirmation, pScreenAdjusted=FALSE)
stageRObj <- stageWiseAdjustment(object=stageRObj, method="holm", alpha=0.05)
stageWiseRes <- getAdjustedPValues(stageRObj,onlySignificantGenes = TRUE,order = TRUE)
stageWiseRes_all <- getAdjustedPValues(stageRObj,onlySignificantGenes = FALSE,order = TRUE)

colSums(stageWiseRes<0.05)


## Write results & gene lists
###############

res_1dph <- dph1$table[sort(dph1$table$PValue,index.return=T)$ix,]
res_20dph <- dph20$table[sort(dph20$table$PValue,index.return=T)$ix,]
res_65dph <- dph65$table[sort(dph65$table$PValue,index.return=T)$ix,]
res_adult <- dphad$table[sort(dphad$table$PValue,index.return=T)$ix,]
res_age <- dphage$table[sort(dphage$table$PValue,index.return=T)$ix,]
res <- list(res_1dph,res_20dph,res_65dph,res_adult,res_age)

list_all <- rownames(stageWiseRes[stageWiseRes[,1]<0.05,])
list_1dph <- rownames(stageWiseRes[stageWiseRes[,2]<0.05,])
list_20dph <- rownames(stageWiseRes[stageWiseRes[,3]<0.05,])
list_65dph <- rownames(stageWiseRes[stageWiseRes[,4]<0.05,])
list_adult <- rownames(stageWiseRes[stageWiseRes[,5]<0.05,])
list_age <- rownames(stageWiseRes[stageWiseRes[,6]<0.05,])
list_nonsig <- rna_annot$V1[keep][!(rna_annot$V1[keep]%in%unique(c(list_1dph,list_20dph,list_65dph,list_adult,list_age)))]

write.table(list_all,file="GeneLists/ListGenes_all.txt",col.names = F,row.names = F,quote = F,sep="\t")
write.table(list_1dph,file="GeneLists/ListGenes_1dph.txt",col.names = F,row.names = F,quote = F,sep="\t")
write.table(list_20dph,file="GeneLists/ListGenes_20dph.txt",col.names = F,row.names = F,quote = F,sep="\t")
write.table(list_65dph,file="GeneLists/ListGenes_65dph.txt",col.names = F,row.names = F,quote = F,sep="\t")
write.table(list_adult,file="GeneLists/ListGenes_adult.txt",col.names = F,row.names = F,quote = F,sep="\t")
write.table(list_age,file="GeneLists/ListGenes_age.txt",col.names = F,row.names = F,quote = F,sep="\t")
write.table(list_nonsig,file="GeneLists/ListGenes_nonsig.txt",col.names = F,row.names = F,quote = F,sep="\t")

coi <- c("1dph","20dph","65dph","adult","age")
for (i in 1:length(coi)){
  sorted_genes <- rownames(stageWiseRes[sort(stageWiseRes[,1+i],index.return=T)$ix,])
  for (j in 1:20){
    dat_boxplot <- data.frame(countspm=cpm(y)[sorted_genes[j],],f1=as.factor(f1),f2=as.factor(f2))
    jpeg(paste("boxplots/",coi[i],"_",j,"_",sorted_genes[j],".jpg",sep=""))
    boxplot(countspm ~ f1*f2,data=dat_boxplot, names=unique(paste(f1,f2,sep="_")),at=c(1,3,5,7,2,4,6,8),las=2)
    dev.off()
  }
} 

out2excel <- cbind(rowMeans(rnadata[keep,f1=="1dph" & f2=="F"]),apply(rnadata[keep,f1=="1dph" & f2=="F"],1,sd),
                   rowMeans(rnadata[keep,f1=="1dph" & f2=="M"]),apply(rnadata[keep,f1=="1dph" & f2=="M"],1,sd),
                   rowMeans(rnadata[keep,f1=="20dph" & f2=="F"]),apply(rnadata[keep,f1=="20dph" & f2=="F"],1,sd),
                   rowMeans(rnadata[keep,f1=="20dph" & f2=="M"]),apply(rnadata[keep,f1=="20dph" & f2=="M"],1,sd),
                   rowMeans(rnadata[keep,f1=="65dph" & f2=="F"]),apply(rnadata[keep,f1=="65dph" & f2=="F"],1,sd),
                   rowMeans(rnadata[keep,f1=="65dph" & f2=="M"]),apply(rnadata[keep,f1=="65dph" & f2=="M"],1,sd),
                   rowMeans(rnadata[keep,f1=="adult" & f2=="F"]),apply(rnadata[keep,f1=="adult" & f2=="F"],1,sd),
                   rowMeans(rnadata[keep,f1=="adult" & f2=="M"]),apply(rnadata[keep,f1=="adult" & f2=="M"],1,sd))


genes_new <- rep("",dim(out2excel)[1])
uniprot_new <- rep("",dim(out2excel)[1])
for (i in 1:nrow(out2excel)){
  if (rownames(out2excel)[i]%in%genes[,2]){
    genes_new[i] <- genes[which(genes[,2]==rownames(out2excel)[i]),1]
    if (rownames(out2excel)[i]%in%uniprot[,1]){
      uniprot_new[i] <- uniprot[which(uniprot[,1]==rownames(out2excel)[i]),2]
    }
  }
}

ix <- sort(rownames(stageWiseRes_all),index.return=T)$ix

out2excel <- cbind(rownames(out2excel),genes_new,uniprot_new,stageWiseRes_all[ix,1],
                   dph1$table$logFC,stageWiseRes_all[ix,2],
                   dph20$table$logFC,stageWiseRes_all[ix,3],
                   dph65$table$logFC,stageWiseRes_all[ix,4],
                   dphad$table$logFC,stageWiseRes_all[ix,5],
                   dphage$table$logFC,stageWiseRes_all[ix,6],
                   out2excel)

colnames(out2excel) <- c("ensembl.gene.id","gene.symbol","uniprot.annot","Pvalue.omnibus",
                         paste("logFC",coi[1],sep="."),paste("Pvalue",coi[1],sep="."),
                         paste("logFC",coi[2],sep="."),paste("Pvalue",coi[2],sep="."),
                         paste("logFC",coi[3],sep="."),paste("Pvalue",coi[3],sep="."),
                         paste("logFC",coi[4],sep="."),paste("Pvalue",coi[4],sep="."),
                         paste("logFC",coi[5],sep="."),paste("Pvalue",coi[5],sep="."),
                         "mean.1dph,F","sd.1dph,F","mean.1dph,M","sd.1dph,M",
                         "mean.20dph,F","sd.20dph,F","mean.20dph,M","sd.20dph,M",
                         "mean.65dph,F","sd.65dph,F","mean.65dph,M","sd.65dph,M",
                         "mean.adult,F","sd.adult,F","mean.adult,M","sd.adult,M")


out2excel <- out2excel[sort(as.numeric(out2excel[,colnames(out2excel)=="Pvalue.omnibus"]),index.return=T)$ix,]
out2excel <- data.frame(out2excel,stringsAsFactors = F)

for (i in 4:ncol(out2excel)){
  out2excel[,i] <- as.numeric(out2excel[,i])
}

# write.xlsx(out2excel,file="FullData_output.xlsx",row.names = F)
write.table(out2excel,file="FullData_output.txt",col.names = T, row.names = F, quote = F, sep="\t")
