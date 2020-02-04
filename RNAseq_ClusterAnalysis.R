##############################
##############################
##                          ##
## RNAseq_ClusterAnalysis.R ##
##                          ##
##############################
##############################


## Load libraries
################

library(edgeR)
library(ggplot2)


## Load Data & annotation
###############

setwd("C:/Users/Louis/OneDrive/Louis/AAP/Zebra_Finch_new/")
rnadata<-read.table("Plastoscine.txt",header=T,sep="\t",quote="")
rna_annot <- read.table("Plastoscine_annot.txt",header=F,sep="\t",quote="")$V1

# Make DGElist object
y <- DGEList(rnadata)
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep,  ,keep.lib.sizes=FALSE]
summary(keep)

# Perform filtering 
rna_annot <- rna_annot[keep]

# Calculate logcpm counts
logcpm <- cpm(y,log=T)

# Annotation RNAseq
f1 <- "1dph,1dph,1dph,1dph,1dph,1dph,20dph,20dph,20dph,20dph,20dph,20dph,65dph,65dph,65dph,65dph,65dph,65dph,adult,adult,adult,adult,adult,adult"
f2 <- "F,F,F,M,M,M,F,F,F,M,M,M,F,F,F,M,M,M,F,F,F,M,M,M"
f1 <- unlist(strsplit(f1,","))
f2 <- unlist(strsplit(f2,","))


## Mean & sd calculations
###############

print(dim(logcpm))
logcpm_rowm <- NULL
for (t in unique(f1)){
  logcpm_rowm <- cbind(logcpm_rowm,rowMeans(logcpm[,f1==t]))
}
print(dim(logcpm_rowm))

# Calculate standard deviation per gene
sd_row <- sqrt(apply(logcpm,1,var))
# elevated <- NULL
# for(i in 1:nrow(logcpm_rowm)){
#   elevated <- rbind(elevated,logcpm_rowm[i,2:4] > logcpm_rowm[i,1]+sd_row)
# }


## Definition of borders
###############

am_sd <- 1.65
pnorm(am_sd)
elevated <- logcpm_rowm[,2:4] > logcpm_rowm[,1]+am_sd*sd_row
equal <- (logcpm_rowm[,2:4] < logcpm_rowm[,1]+1*sd_row) & (logcpm_rowm[,2:4] > logcpm_rowm[,1]-1*sd_row)
lower <- logcpm_rowm[,2:4] < logcpm_rowm[,1]-am_sd*sd_row
direction <- list(elevated,equal,lower)
names(direction) <- c("elevated","equal","lower")
head(elevated)
dim(elevated)


## Trend extraction
###############

pattern1 <- c("lower","lower","equal")
pattern2 <- c("elevated","elevated","equal")
pattern3 <- c("elevated","elevated","elevated")
pattern4 <- c("lower","lower","lower")
pattern5 <- c("elevated","equal","equal")
pattern6 <- c("lower","equal","equal")
pattern7 <- c("equal","lower","equal")
pattern8 <- c("equal","elevated","equal")
pattern9 <- c("equal","equal","lower")
pattern10 <- c("equal","equal","elevated")
pattern11 <- c("equal","lower","lower")
pattern12 <- c("equal","elevated","elevated")

# Get amount of patterns
am_patt <- max(as.numeric(unlist(strsplit(ls()[grep("pattern",ls())],"pattern"))[unlist(strsplit(ls()[grep("pattern",ls())],"pattern"))!=""]))
annot_clust <- list()
# jpeg("Trends.jpg")
# par(mfrow=c(3,4))
for (t in 1:am_patt){
  
  # Make mask for pattern
  mask_pat <- NULL
  for (i in 1:nrow(elevated)){
    pattern <- eval(parse(text = paste("pattern",t,sep="")))
    mask_pat <- c(mask_pat,sum(c(eval(parse(text = paste("direction$",pattern[1],sep="")))[i,1],
                                 eval(parse(text = paste("direction$",pattern[2],sep="")))[i,2],
                                 eval(parse(text = paste("direction$",pattern[3],sep="")))[i,3]))==3)
  }
  print(head(mask_pat))
  print(sum(mask_pat))
  print(length(mask_pat))
  
  # Select annotation
  annot_mask <- rna_annot[mask_pat]
  annot_clust[[t]] <- annot_mask
  
  # Write genenames for pattern
  write.table(annot_mask,file=paste("GeneLists/ListGenes_patt_",paste(pattern,collapse="_"),".txt",sep=""), quote = F, row.names = F, col.names = F, sep = "\t")
  
  # Select logcpm values for pattern
  logcpm_mask <- logcpm_rowm[mask_pat,1:4]
  x <- as.factor(unique(f1))
  
  # Plot trends
  png(paste("Trends/trend_patt_",paste(pattern,collapse="_"),".png",sep=""))
  if(nrow(logcpm_mask)>500){
    rest <- nrow(logcpm_mask) %% 500
    times <- (nrow(logcpm_mask)-rest)/500
    for (i in 1:times){
      if (i!=times){
        logcpm_mask_sub <- logcpm_mask[((500*(i-1))+1):(500*i),]
        temp <- data.frame(Age = x, LogCpm = colMeans(logcpm_mask_sub))
        g <- ggplot(data=temp,aes(Age,LogCpm,group=1)) +
          theme_bw() +
          geom_line() +
          ylim(min(logcpm_rowm),max(logcpm_rowm))
        for (i in 1:nrow(logcpm_mask_sub)) {
          #print(temp)
          temp <- data.frame(Age = x, LogCpm = logcpm_mask_sub[i,])
          g <- g + geom_line(data = temp, aes(Age,LogCpm),col="grey") 
        }
        print(g)
      } else {
        logcpm_mask_sub <- logcpm_mask[((500*(i-1))+1):((500*(i-1))+rest),]
        temp <- data.frame(Age = x, LogCpm = colMeans(logcpm_mask_sub))
        g <- ggplot(data=temp,aes(Age,LogCpm,group=1)) +
          theme_bw() +
          geom_line()
        for (i in 1:nrow(logcpm_mask_sub)) {
          #print(temp)
          temp <- data.frame(Age = x, LogCpm = logcpm_mask_sub[i,])
          g <- g + geom_line(data = temp, aes(Age,LogCpm),col="grey")
        }
        textdf <- data.frame(Age = 2.5, LogCpm = min(logcpm_mask)+(max(logcpm_mask)-min(logcpm_mask))*0.10)
        temp <- data.frame(Age = x, LogCpm = colMeans(logcpm_mask))
        g <- g + geom_label(data=textdf,aes(Age,LogCpm),label=paste(paste(toupper(substr(gsub("elevated","higher",pattern),1,1)),collapse=""),": ",nrow(logcpm_mask),sep=""),size = 6,label.size = 0) + 
          geom_line(data = temp, aes(Age,LogCpm),col="blue")
        print(g)
      }
    }
  } else { 
    temp <- data.frame(Age = x, LogCpm = colMeans(logcpm_mask))
    g <- ggplot(data=temp,aes(Age,LogCpm,group=1)) +
      theme_bw() +
      geom_line()
    for (i in 1:nrow(logcpm_mask)) {
      #print(temp)
      temp <- data.frame(Age = x, LogCpm = logcpm_mask[i,])
      g <- g +  geom_line(data = temp, aes(Age,LogCpm),col="grey")
    }
    temp <- data.frame(Age = x, LogCpm = colMeans(logcpm_mask))
    textdf <- data.frame(Age = 2.5, LogCpm = min(logcpm_mask)+(max(logcpm_mask)-min(logcpm_mask))*0.10)
    g <- g + geom_label(data=textdf,aes(Age,LogCpm),label=paste(paste(toupper(substr(gsub("elevated","higher",pattern),1,1)),collapse=""),": ",nrow(logcpm_mask),sep=""),size = 6,label.size = 0) + 
      geom_line(data = temp, aes(Age,LogCpm),col="blue")
    print(g)

  }
  dev.off()
}

library(grid)
library(png)
library(ggplot2)
library(gridExtra)

unlink("Trends_overall.png")
plots <- lapply(ll <- list.files("Trends/",patt='*[.]png'),function(x){
  img <- as.raster(readPNG(paste("Trends/",x,sep="")))
  rasterGrob(img, interpolate = FALSE)
})

ggsave("Trends_overall.png",width=20, height=15, units="cm",
       marrangeGrob(grobs = plots, nrow=3, ncol=4,top=NULL))



## Trend metrics
###############
n_clusters <- 0
for (t in 1:am_patt){
  pattern <- eval(parse(text = paste("pattern",t,sep="")))
  print(pattern)
  gene_list_patt <-read.table(file=paste("GeneLists/ListGenes_patt_",paste(pattern,collapse="_"),".txt",sep=""), header = F, sep="\t")$V1
  n_clusters <- n_clusters + length(gene_list_patt)
}
print(n_clusters)
print(n_clusters/13891)



####################
#####################


# VOOR WEGIO


#################"
###################


## Check for differential methylation
###############

RRBS_age_full <- read.table("RRBS_age_full.txt",header=T,stringsAsFactors = F,sep="\t")
genes_RRBS_age <- as.character(unique(unlist(RRBS_age_full[RRBS_age_full$P.Value<0.05,grepl("gene",colnames(RRBS_age_full)) & !grepl("gene_symbol",colnames(RRBS_age_full))]))[-1])
for (t in 1:am_patt){
  print(paste(pattern,collapse="_"))
  pattern <- eval(parse(text = paste("pattern",t,sep="")))
  genes_cluster <- read.table(paste("GeneLists/ListGenes_patt_",paste(pattern,collapse="_"),".txt",sep=""),header=F,stringsAsFactors = F)$V1
  print(length(genes_cluster))
  # print(sum(genes_cluster%in%genes_RRBS_age))
  # print(sum(genes_cluster%in%genes_RRBS_age)/length(genes_RRBS_age))
  mat <- cbind(c(sum(genes_cluster%in%genes_RRBS_age),length(genes_cluster)),
               c(length(genes_RRBS_age),10574))
  print(mat)
  print(chisq.test(mat)$p.value)
}




####
## Sex ##


# Calculate mean per gene
print(dim(logcpm))
logcpm_rowm <- NULL
for (t in unique(f1)){
  logcpm_rowm <- cbind(logcpm_rowm,rowMeans(logcpm[,f2=="M" & f1==t])-rowMeans(logcpm[,f2=="F" & f1==t]))
}
print(dim(logcpm_rowm))

declining <- (logcpm_rowm[,1]>logcpm_rowm[,2]) & (logcpm_rowm[,2]>logcpm_rowm[,3]) & (logcpm_rowm[,3]>logcpm_rowm[,4])

for (i in which(declining)){
  if (i == which(declining)[1]){
    plot(1:4,logcpm_rowm[i,])
  } else {
    lines(1:4,logcpm_rowm)
  }
}

fulldata_output_sorted <- fulldata_output[sort(as.character(fulldata_output$ensembl.gene.id),index.return=T)$ix,]

sort(unique(genes$hgnc_symbol[genes$ensembl_gene_id%in%rna_annot[declining]]))

ascending <- (logcpm_rowm[,4]>logcpm_rowm[,3]) & (logcpm_rowm[,3]>logcpm_rowm[,2]) & (logcpm_rowm[,2]>logcpm_rowm[,1])

sort(unique(genes$hgnc_symbol[genes$ensembl_gene_id%in%rna_annot[ascending]]))

declining_plus <- declining & ((fulldata_output_sorted$Pvalue.1dph<0.05) & !is.na(fulldata_output_sorted$Pvalue.1dph))
sort(unique(genes$hgnc_symbol[genes$ensembl_gene_id%in%rna_annot[declining_plus]]))

# Calculate standard deviation per gene
sd_row <- sqrt(apply(logcpm,1,var))
# elevated <- NULL
# for(i in 1:nrow(logcpm_rowm)){
#   elevated <- rbind(elevated,logcpm_rowm[i,2:4] > logcpm_rowm[i,1]+sd_row)
# }

# Set cutoff for elevated/lower or equal
am_sd <- 1.65
pnorm(am_sd)
elevated <- logcpm_rowm > am_sd*sd_row
equal <- (logcpm_rowm < sd_row) & (logcpm_rowm > -sd_row)
lower <- logcpm_rowm < -am_sd*sd_row
direction <- list(elevated,equal,lower)
names(direction) <- c("elevated","equal","lower")
head(elevated)
dim(elevated)

# Define patterns
rm(list=ls()[grep("pattern",ls())])
pattern1 <- c("equal","elevated","elevated","equal")
pattern2 <- c("equal","equal","elevated","equal")
pattern3 <- c("equal","equal","equal","elevated")
pattern4 <- c("equal","elevated","equal","equal")
# pattern2 <- c("elevated","elevated","equal")
# pattern3 <- c("elevated","elevated","elevated")
# pattern4 <- c("lower","lower","lower")
# pattern5 <- c("elevated","equal","equal")
# pattern6 <- c("lower","equal","equal")
# pattern7 <- c("equal","lower","equal")
# pattern8 <- c("equal","elevated","equal")
# pattern9 <- c("equal","equal","lower")
# pattern10 <- c("equal","equal","elevated")
# pattern11 <- c("equal","lower","lower")
# pattern12 <- c("equal","elevated","elevated")

am_patt <- max(as.numeric(unlist(strsplit(ls()[grep("pattern",ls())],"pattern"))[unlist(strsplit(ls()[grep("pattern",ls())],"pattern"))!=""]))
annot_clust <- list()
for (t in 1:am_patt){
  
  # Make mask for pattern
  mask_pat <- NULL
  for (i in 1:nrow(elevated)){
    pattern <- eval(parse(text = paste("pattern",t,sep="")))
    mask_pat <- c(mask_pat,sum(c(eval(parse(text = paste("direction$",pattern[1],sep="")))[i,1],
                                 eval(parse(text = paste("direction$",pattern[2],sep="")))[i,2],
                                 eval(parse(text = paste("direction$",pattern[3],sep="")))[i,3],
                                 eval(parse(text = paste("direction$",pattern[4],sep="")))[i,4]))==4)
  }
  print(head(mask_pat))
  print(sum(mask_pat))
  print(length(mask_pat))
  if(sum(mask_pat)==0){
    next
  }
  
  # Select annotation
  annot_mask <- rna_annot[mask_pat]
  annot_clust[[t]] <- annot_mask
  
  # Write genenames for pattern
  write.table(annot_mask,file=paste("Trends_sex/genelist_patt",paste(pattern,collapse="_"),".txt",sep=""), quote = F, row.names = F, col.names = F, sep = "\t")
  
  # Select logcpm values for pattern
  logcpm_mask <- logcpm_rowm[mask_pat,1:4]
  x <- as.factor(unique(f1))
  
  # Plot trends
  jpeg(paste("Trends_sex/trend_patt_",paste(pattern,collapse="_"),".jpg",sep=""))
  if(nrow(logcpm_mask)>500){
    rest <- nrow(logcpm_mask) %% 500
    times <- (nrow(logcpm_mask)-rest)/500
    for (i in 1:times){
      if (i!=times){
        logcpm_mask_sub <- logcpm_mask[((500*(i-1))+1):(500*i),]
        temp <- data.frame(Age = x, LogCpm = colMeans(logcpm_mask_sub))
        g <- ggplot(data=temp,aes(Age,LogCpm,group=1)) +
          theme_bw() +
          geom_line()
        for (i in 1:nrow(logcpm_mask_sub)) {
          #print(temp)
          temp <- data.frame(Age = x, LogCpm = logcpm_mask_sub[i,])
          g <- g + geom_line(data = temp, aes(Age,LogCpm),col="grey")
        }
        print(g)
      } else {
        logcpm_mask_sub <- logcpm_mask[((500*(i-1))+1):((500*(i-1))+rest),]
        temp <- data.frame(Age = x, LogCpm = colMeans(logcpm_mask_sub))
        g <- ggplot(data=temp,aes(Age,LogCpm,group=1)) +
          theme_bw() +
          geom_line()
        for (i in 1:nrow(logcpm_mask_sub)) {
          #print(temp)
          temp <- data.frame(Age = x, LogCpm = logcpm_mask_sub[i,])
          g <- g + geom_line(data = temp, aes(Age,LogCpm),col="grey")
        }
        temp <- data.frame(Age = x, LogCpm = colMeans(logcpm_mask))
        g <- g + geom_line(data = temp, aes(Age,LogCpm),col="blue")
        print(g)
      }
    }
  } else { 
    temp <- data.frame(Age = x, LogCpm = colMeans(logcpm_mask))
    g <- ggplot(data=temp,aes(Age,LogCpm,group=1)) +
      theme_bw() +
      geom_line()
    for (i in 1:nrow(logcpm_mask)) {
      #print(temp)
      temp <- data.frame(Age = x, LogCpm = logcpm_mask[i,])
      g <- g + geom_line(data = temp, aes(Age,LogCpm),col="grey")
    }
    temp <- data.frame(Age = x, LogCpm = colMeans(logcpm_mask))
    g <- g + geom_line(data = temp, aes(Age,LogCpm),col="blue")
    print(g)
    
  }
  dev.off()
}



# Load in Expression and Jaspar analysis
DE_results <- read.table("FullData_output.txt",header=T,sep="\t",stringsAsFactors = F)
JASPAR_results_1dph <- read.table("JASPAR_TF_1dph_comp_h.txt",header=T,stringsAsFactors = F,sep="\t")
JASPAR_results_20dph <- read.table("JASPAR_TF_20dph_comp_h.txt",header=T,stringsAsFactors = F,sep="\t")
JASPAR_results_65dph <- read.table("JASPAR_TF_65dph_comp_h.txt",header=T,stringsAsFactors = F,sep="\t")
JASPAR_results_adult <- read.table("JASPAR_TF_adult_comp_h.txt",header=T,stringsAsFactors = F,sep="\t")
JASPAR_results_age <- read.table("JASPAR_TF_age_comp_h.txt",header=T,stringsAsFactors = F,sep="\t")

genes_JASPAR_1dph <- DE_results$ensembl.gene.id[DE_results$gene.symbol%in%JASPAR_results_1dph$TF]
DE_results[(DE_results$ensembl.gene.id%in%genes_JASPAR_1dph) & (DE_results$Pvalue.1dph<0.05 & !is.na(DE_results$Pvalue.1dph)),]

genes_JASPAR_20dph <- DE_results$ensembl.gene.id[DE_results$gene.symbol%in%JASPAR_results_20dph$TF]
DE_results[(DE_results$ensembl.gene.id%in%genes_JASPAR_20dph) & (DE_results$Pvalue.20dph<0.05 & !is.na(DE_results$Pvalue.20dph)),]

genes_JASPAR_65dph <- DE_results$ensembl.gene.id[DE_results$gene.symbol%in%JASPAR_results_65dph$TF]
DE_results[(DE_results$ensembl.gene.id%in%genes_JASPAR_65dph) & (DE_results$Pvalue.65dph<0.05 & !is.na(DE_results$Pvalue.65dph)),]

genes_JASPAR_adult <- DE_results$ensembl.gene.id[DE_results$gene.symbol%in%JASPAR_results_adult$TF]
DE_results[(DE_results$ensembl.gene.id%in%genes_JASPAR_adult) & (DE_results$Pvalue.adult<0.05 & !is.na(DE_results$Pvalue.adult)),]

genes_JASPAR_age <- DE_results$ensembl.gene.id[DE_results$gene.symbol%in%JASPAR_results_age$TF]
DE_results[(DE_results$ensembl.gene.id%in%genes_JASPAR_age) & (DE_results$Pvalue.age<0.05 & !is.na(DE_results$Pvalue.age)),]

for (i in 1:length(annot_clust)){
  pattern <- eval(parse(text = paste("pattern",i,sep="")))
  if (pattern[1]!="equal"){
    print(paste(pattern,collapse="-"))
    print(sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.20dph<0.05 & !is.na(DE_results$Pvalue.20dph)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]]) & (DE_results$ensembl.gene.id%in%genes_JASPAR_20dph)]))))
  }
}


genes_DE_inCluster_1dph <- list()
for (i in 1:length(annot_clust)){
  print(sum(annot_clust[[i]]%in%DE_results$ensembl.gene.id[DE_results$Pvalue.1dph<0.05 & !is.na(DE_results$Pvalue.1dph)]))
  genes_DE_inCluster_1dph[[i]] <- sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.1dph<0.05 & !is.na(DE_results$Pvalue.1dph)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]])])))
}
genes_DE_inCluster_20dph <- list()
for (i in 1:length(annot_clust)){
  print(sum(annot_clust[[i]]%in%DE_results$ensembl.gene.id[DE_results$Pvalue.20dph<0.05 & !is.na(DE_results$Pvalue.20dph)]))
  genes_DE_inCluster_20dph[[i]] <- sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.20dph<0.05 & !is.na(DE_results$Pvalue.20dph)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]])])))
}
genes_DE_inCluster_65dph <- list()
for (i in 1:length(annot_clust)){
  print(sum(annot_clust[[i]]%in%DE_results$ensembl.gene.id[DE_results$Pvalue.65dph<0.05 & !is.na(DE_results$Pvalue.65dph)]))
  genes_DE_inCluster_65dph[[i]] <- sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.65dph<0.05 & !is.na(DE_results$Pvalue.65dph)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]])])))
}
genes_DE_inCluster_adult <- list()
for (i in 1:length(annot_clust)){
  print(sum(annot_clust[[i]]%in%DE_results$ensembl.gene.id[DE_results$Pvalue.adult<0.05 & !is.na(DE_results$Pvalue.adult)]))
  genes_DE_inCluster_adult[[i]] <- sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.adult<0.05 & !is.na(DE_results$Pvalue.adult)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]])])))
}

compareGenes <- function(genes_in,genes_out){
  all_genes <- sort(unique(unlist(genes_in,genes_out)))
  for (i in 1:length(genes_in)){
    if (i == 1){
      mask_in <- all_genes%in%genes_in[[i]]
    } else {
      mask_in <- mask_in & (all_genes%in%genes_in[[i]])
    }
  }
  for (j in 1:length(genes_out)){
    if (j == 1){
      mask_out <- all_genes%in%genes_out[[j]]
    } else {
      mask_out <- mask_out | (all_genes%in%genes_out[[j]])
    }
  }
  mask <- mask_in & !(mask_out)
  return(all_genes[mask])
}

## Only in 1dph
for (i in 1:length(genes_DE_inCluster_1dph)){
  print(paste("pattern",i,":",paste(eval(parse(text = paste("pattern",i,sep=""))),collapse="-")))
  print(compareGenes(list(genes_DE_inCluster_1dph[[i]]),
                     list(genes_DE_inCluster_65dph[[i]],genes_DE_inCluster_20dph[[i]],genes_DE_inCluster_adult[[i]])))
}

## Only in 20dph
for (i in 1:length(genes_DE_inCluster_20dph)){
  print(paste("pattern",i,":",paste(eval(parse(text = paste("pattern",i,sep=""))),collapse="-")))
  print(compareGenes(list(genes_DE_inCluster_20dph[[i]]),
                     list(genes_DE_inCluster_65dph[[i]],genes_DE_inCluster_1dph[[i]],genes_DE_inCluster_adult[[i]])))
}

## Only in 65dph
for (i in 1:length(genes_DE_inCluster_65dph)){
  print(paste("pattern",i,":",paste(eval(parse(text = paste("pattern",i,sep=""))),collapse="-")))
  print(compareGenes(list(genes_DE_inCluster_65dph[[i]]),
                     list(genes_DE_inCluster_20dph[[i]],genes_DE_inCluster_1dph[[i]],genes_DE_inCluster_adult[[i]])))
}

## Only in adult
for (i in 1:length(genes_DE_inCluster_adult)){
  print(paste("pattern",i,":",paste(eval(parse(text = paste("pattern",i,sep=""))),collapse="-")))
  print(compareGenes(list(genes_DE_inCluster_adult[[i]]),
                     list(genes_DE_inCluster_65dph[[i]],genes_DE_inCluster_1dph[[i]],genes_DE_inCluster_20dph[[i]])))
}

## Stale

for (i in 1:length(genes_DE_inCluster_1dph)){
  print(compareGenes(list(genes_DE_inCluster_1dph[[i]],genes_DE_inCluster_20dph[[i]],genes_DE_inCluster_65dph[[i]]),
                     list(genes_DE_inCluster_adult[[i]])))
}

for (i in 1:length(genes_DE_inCluster_1dph)){
  print(compareGenes(list(genes_DE_inCluster_20dph[[i]],genes_DE_inCluster_65dph[[i]]),
                     list(genes_DE_inCluster_1dph[[i]],genes_DE_inCluster_adult[[i]])))
}




for (i in 1:length(annot_clust)){
  print(sum(annot_clust[[i]]%in%DE_results$ensembl.gene.id[DE_results$Pvalue.age<0.05 & !is.na(DE_results$Pvalue.age)]))
  print(sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.age<0.05 & !is.na(DE_results$Pvalue.age)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]]) & (DE_results$ensembl.gene.id%in%JASPAR_results_age$TF[JASPAR_results_age$FDR<0.05])]))))
}


for (i in 1:length(annot_clust)){
  print(sum(annot_clust[[i]]%in%DE_results$ensembl.gene.id[DE_results$Pvalue.1dph<0.05 & !is.na(DE_results$Pvalue.1dph)]))
  print(sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.1dph<0.05 & !is.na(DE_results$Pvalue.1dph)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]]) & (DE_results$ensembl.gene.id%in%JASPAR_results_1dph$TF[JASPAR_results_1dph$FDR<0.05])]))))
}

for (i in 1:length(annot_clust)){
  print(sum(annot_clust[[i]]%in%DE_results$ensembl.gene.id[DE_results$Pvalue.20dph<0.05 & !is.na(DE_results$Pvalue.20dph)]))
  print(sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.20dph<0.05 & !is.na(DE_results$Pvalue.20dph)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]]) & (DE_results$ensembl.gene.id%in%JASPAR_results_20dph$TF[JASPAR_results_20dph$FDR<0.05])]))))
}

for (i in 1:length(annot_clust)){
  print(sum(annot_clust[[i]]%in%DE_results$ensembl.gene.id[DE_results$Pvalue.65dph<0.05 & !is.na(DE_results$Pvalue.65dph)]))
  print(sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.65dph<0.05 & !is.na(DE_results$Pvalue.65dph)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]]) & (DE_results$ensembl.gene.id%in%JASPAR_results_65dph$TF[JASPAR_results_65dph$FDR<0.05])]))))
}

for (i in 1:length(annot_clust)){
  print(sum(annot_clust[[i]]%in%DE_results$ensembl.gene.id[DE_results$Pvalue.adult<0.05 & !is.na(DE_results$Pvalue.adult)]))
  print(sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.adult<0.05 & !is.na(DE_results$Pvalue.adult)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]]) & (DE_results$ensembl.gene.id%in%JASPAR_results_adult$TF[JASPAR_results_adult$FDR<0.05])]))))
}

for (i in 1:length(annot_clust)){
  print(sum(annot_clust[[i]]%in%DE_results$ensembl.gene.id[DE_results$Pvalue.65dph<0.05 & !is.na(DE_results$Pvalue.65dph)]))
  print(sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.65dph<0.05 & !is.na(DE_results$Pvalue.65dph)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]])]))))
}

for (i in 1:length(annot_clust)){
  print(sum(annot_clust[[i]]%in%DE_results$ensembl.gene.id[DE_results$Pvalue.65dph<0.05 & !is.na(DE_results$Pvalue.65dph)]))
  print(sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.65dph<0.05 & !is.na(DE_results$Pvalue.65dph)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]])]))))
}

for (i in 1:length(annot_clust)){
  print(sum(annot_clust[[i]]%in%DE_results$ensembl.gene.id[DE_results$Pvalue.adult<0.05 & !is.na(DE_results$Pvalue.adult)]))
  print(sort(unique(as.character(DE_results$gene.symbol[(DE_results$Pvalue.adult<0.05 & !is.na(DE_results$Pvalue.adult)) & (DE_results$ensembl.gene.id%in%annot_clust[[i]])]))))
}

# plot(mean(logcpm_mask), type="n", xlab="Age", ylab="logcpm",ylim=c(min(logcpm_mask),max(logcpm_mask)))
# 
# for (i in 1:nrow(logcpm_mask)){
#   lines(x,logcpm_mask[i,],col="grey")
# }
# lines(x,colMeans(logcpm_mask),col="blue")


JASPAR_results_1dph$TF[JASPAR_results_1dph$FDR<0.05][JASPAR_results_1dph$TF[JASPAR_results_1dph$FDR<0.05]%in%DE_results$gene.symbol[(DE_results$Pvalue.1dph<0.05 & !is.na(DE_results$Pvalue.1dph))]]
JASPAR_results_20dph$TF[JASPAR_results_20dph$FDR<0.05][JASPAR_results_20dph$TF[JASPAR_results_20dph$FDR<0.05]%in%DE_results$gene.symbol[(DE_results$Pvalue.20dph<0.05 & !is.na(DE_results$Pvalue.20dph))]]
JASPAR_results_65dph$TF[JASPAR_results_65dph$FDR<0.05][JASPAR_results_65dph$TF[JASPAR_results_65dph$FDR<0.05]%in%DE_results$gene.symbol[(DE_results$Pvalue.65dph<0.05 & !is.na(DE_results$Pvalue.65dph))]]
JASPAR_results_adult$TF[JASPAR_results_adult$FDR<0.05][JASPAR_results_adult$TF[JASPAR_results_adult$FDR<0.05]%in%DE_results$gene.symbol[(DE_results$Pvalue.adult<0.05 & !is.na(DE_results$Pvalue.adult))]]
JASPAR_results_age$TF[JASPAR_results_age$FDR<0.05][JASPAR_results_age$TF[JASPAR_results_age$FDR<0.05]%in%DE_results$gene.symbol[(DE_results$Pvalue.age<0.05 & !is.na(DE_results$Pvalue.age))]]


