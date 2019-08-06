###########################
###########################
##                       ##
## RRBS_analyse_nonCpG.R ##
##                       ##
###########################
###########################


setwd("C:/Users/Louis/Onedrive/Louis/AAP/Zebra_Finch_new")

##################
## Load in data ##
##################

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

genes_DE <- fulldata_output$ensembl.gene.id[(fulldata_output$Pvalue.age<0.05) & !is.na(fulldata_output$Pvalue.age)]
genes_notDE <- fulldata_output$ensembl.gene.id[!(fulldata_output$Pvalue.age<0.05) | is.na(fulldata_output$Pvalue.age)]
# genes_DM <- unlist(rrbs_age_full[rrbs_age_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_age_full)) & !grepl("gene_symbol",colnames(rrbs_age_full))])
genes_DM <- unlist(rrbs_age_full[rrbs_age_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_age_full)) & !grepl("gene_symbol",colnames(rrbs_age_full))])
genes_DM <- genes_DM[!is.na(genes_DM)]

# genes_DM_nonCpG <- unlist(rrbs_nonCpG_age_full[rrbs_nonCpG_age_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_nonCpG_age_full)) & !grepl("gene_symbol",colnames(rrbs_nonCpG_age_full))])
genes_DM_nonCpG  <- unlist(rrbs_nonCpG_age_full[rrbs_nonCpG_age_full$adj.P.Val<0.05,grepl("gene",colnames(rrbs_nonCpG_age_full)) & !grepl("gene_symbol",colnames(rrbs_nonCpG_age_full))])
genes_DM_nonCpG <- genes_DM_nonCpG[!is.na(genes_DM_nonCpG)]

ct <- rbind(c(sum(genes_DE%in%genes_DM),sum(!(genes_DE%in%genes_DM))),
            c(sum(genes_notDE%in%genes_DM),sum(!(genes_notDE%in%genes_DM))))
chisq.test(ct)
ct
ct[,1]/(ct[,1]+ct[,2])

genes_TF <- fulldata_output_jaspar$ensembl.gene.id[(fulldata_output_jaspar$Pvalue.age<0.05) & !is.na(fulldata_output_jaspar$Pvalue.age)]
genes_notTF <- genes_DE[!(genes_DE%in%genes_TF)]

ct <- rbind(c(sum(genes_notTF%in%genes_DM),sum(!(genes_notTF%in%genes_DM))),
            c(sum(genes_TF%in%genes_DM),sum(!(genes_TF%in%genes_DM))))
chisq.test(ct)
ct
ct[,1]/(ct[,1]+ct[,2])

sort(unique(fulldata_output_jaspar$uniprot.annot[(fulldata_output_jaspar$ensembl.gene.id%in%genes_DM) & (fulldata_output_jaspar$ensembl.gene.id%in%genes_DE)],
            fulldata_output_jaspar$gene.symbol[(fulldata_output_jaspar$ensembl.gene.id%in%genes_DM) & (fulldata_output_jaspar$ensembl.gene.id%in%genes_DE)]))
sort(unique(fulldata_output_jaspar$uniprot.annot[fulldata_output_jaspar$ensembl.gene.id%in%genes_DM_nonCpG],fulldata_output_jaspar$gene.symbol[fulldata_output_jaspar$ensembl.gene.id%in%genes_DM_nonCpG]))


## check matches ESR1
sort(unique(unlist(rrbs_age_RNAseq[paste(rrbs_age_RNAseq$chr,rrbs_age_RNAseq$locations,sep="_")%in%matches_ESR1$Gene,grepl("gene_symbol",colnames(rrbs_age_RNAseq))])))


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
makeBedGraph(rrbs_age_prom_sign,c(1,2,2,12),"RRBS_age_Female","Meth_perc",col=c("green","red"),isSign = (rrbs_age_prom_sign$adj.P.Val<0.05), range=c(0,0.3))
makeBedGraph(rrbs_age_prom_sign,c(1,2,2,14),"RRBS_age_Male","Meth_perc",col=c("blue","brown"),isSign = (rrbs_age_prom_sign$adj.P.Val<0.05), range=c(0,0.3))


head(sort(table(rrbs_age_RNAseq$gene[(rrbs_age_RNAseq$func=="promoter") & (rrbs_age_RNAseq$adj.P.Val<0.05)]),decreasing = T))
head(sort(table(rrbs_age_RNAseq$gene[(rrbs_age_RNAseq$func=="promoter")]),decreasing = T))

rrbs_adult_prom_sign <- rrbs_adult_RNAseq[(rrbs_adult_RNAseq$func=="promoter"),]
rrbs_adult_prom_sign <- rrbs_adult_prom_sign[!grepl("random",rrbs_adult_prom_sign$chr),]
rrbs_adult_prom_sign <- rrbs_adult_prom_sign[!grepl("Un",rrbs_adult_prom_sign$chr),]
rrbs_adult_prom_sign <- rrbs_adult_prom_sign[!grepl("MT",rrbs_adult_prom_sign$chr),]
makeBedGraph(rrbs_adult_prom_sign,c(1,2,2,12),"RRBS_adult_Female","Meth_perc",col=c("green","red"),isSign = (rrbs_adult_prom_sign$adj.P.Val<0.05), range=c(0,0.3))
makeBedGraph(rrbs_adult_prom_sign,c(1,2,2,14),"RRBS_adult_Male","Meth_perc",col=c("blue","brown"),isSign = (rrbs_adult_prom_sign$adj.P.Val<0.05), range=c(0,0.3))

head(sort(table(rrbs_adult_RNAseq$gene[(rrbs_adult_RNAseq$func=="promoter") & (rrbs_adult_RNAseq$adj.P.Val<0.05)]),decreasing = T))
head(sort(table(rrbs_adult_RNAseq$gene[(rrbs_adult_RNAseq$func=="promoter")]),decreasing = T))


