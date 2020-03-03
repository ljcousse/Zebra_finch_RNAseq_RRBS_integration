###################################
# ExtractSeq_TFsearching_RNAseq.R #
###################################

working_dir <- getwd()

# Loading Contrast of Interest
#################################

args <- commandArgs(TRUE)
if(grepl(",",args)){
  args <- unlist(strsplit(args,","))
}
print(args)

# Annotation extraction
#################################

annot <- read.table("JASPAR_annot.txt",header=T,sep="\t",stringsAsFactors=F)
head(annot)

annot$ID <- unlist(strsplit(annot$ID,"\\.[0-9]"))


# List extraction 
##################

fimo_all <- read.table(paste(working_dir,"/fimo_RRBS_out/fimo.txt",sep=""),header=T,stringsAsFactors = F, sep="\t")
# fimo_all <- fimo_all[fimo_all$Adj.Pval<0.05,]

fimo_all$ID <- unlist(strsplit(fimo_all$ID,"\\.[0-9]"))

chrom_pos_in <- read.table(paste("RRBS_List.txt",sep=""),header=F,stringsAsFactors = F)$V1
fimo <- fimo_all[fimo_all$chrom_pos%in%chrom_pos_in,]
fimo_ns <- fimo_all[!(fimo_all$chrom_pos%in%chrom_pos_in),]
s <- length(chrom_pos_in)
ns <- 0
fasta_f <- read.table(paste("seq_RRBS_TFsearching.fa",sep=""),stringsAsFactors = F,sep="\t")$V1
fasta_f <- fasta_f[grep(">",fasta_f)]
ns <- ns + length(fasta_f)
ns <- ns - s
print(head(fimo))
print(head(fimo_ns))
print(s)
print(ns)


if(sum(unique(fimo$chrom_pos)%in%unique(fimo_ns$chrom_pos))){
  stop("Bad gene lists! Overlapping")
}


## signifant loci

ID_chrom_pos <- paste(fimo$ID,fimo$chrom_pos,sep="_")
fimo_sub <- fimo[!duplicated(ID_chrom_pos),]
counts <- table(fimo_sub$ID)
counts <- as.matrix(cbind(names(counts),counts))


tf <- NULL
counts_torm <- NULL
for (t in 1:nrow(counts)){
  if (sum(annot$ID==counts[t,1])==0){
    print(paste(counts[t,1],"is not human"))
    counts_torm <- c(counts_torm,t)
    next
  } 
  if (sum(annot$ID==counts[t,1])>1){
    if (sum(toupper(annot$Name[annot$ID==counts[t,1]])==toupper(annot$Name[annot$ID==counts[t,1]][1]))==sum(annot$ID==counts[t,1])){
      tf <- c(tf,annot$Name[annot$ID==counts[t,1]][1])
    } else {
      tf <- c(tf,paste(annot$Name[annot$ID==counts[t,1]],collapse=","))
    }
  } else {
    tf <- c(tf,annot$Name[annot$ID==counts[t,1]])
  }
}

if(!is.null(counts_torm)){
  counts <- counts[-counts_torm,]
}
counts <- as.data.frame(counts,stringsAsFactors=F)
colnames(counts) <- c("Matrix_ID","Counts")

print(dim(counts))
print(length(tf))
counts$TF <- tf

counts <- counts[c("Matrix_ID","Counts","TF")]
counts$Counts <- as.integer(counts$Counts)
print(head(counts))
print(str(counts))


# Non significant loci
######################

ns_ID_chrom_pos <- paste(fimo_ns$ID,fimo_ns$chrom_pos,sep="_")
fimo_ns_sub <- fimo_ns[!duplicated(ns_ID_chrom_pos),]
counts_ns <- table(fimo_ns_sub$ID)
counts_ns <- as.matrix(cbind(names(counts_ns),counts_ns))

tf_ns <- NULL
counts_ns_torm <- NULL
for (t in 1:nrow(counts_ns)){
  if (sum(annot$ID==counts_ns[t,1])==0){
    print(paste(counts_ns[t,1],"is not human"))
    counts_ns_torm <- c(counts_ns_torm,t)
    next
  } 
  if (sum(annot$ID==counts_ns[t,1])>1){
    if (sum(toupper(annot$Name[annot$ID==counts_ns[t,1]])==toupper(annot$Name[annot$ID==counts_ns[t,1]][1]))==sum(annot$ID==counts_ns[t,1])){
      tf_ns <- c(tf_ns,annot$Name[annot$ID==counts_ns[t,1]][1])
    } else {
      tf_ns <- c(tf_ns,paste(annot$Name[annot$ID==counts_ns[t,1]],collapse=","))
    }
  } else {
    tf_ns <- c(tf_ns,annot$Name[annot$ID==counts_ns[t,1]])
  }
}

if(!is.null(counts_ns_torm)){
  counts_ns <- counts_ns[-counts_ns_torm,]
}
counts_ns <- as.data.frame(counts_ns,stringsAsFactors=F)
colnames(counts_ns) <- c("Matrix_ID","Counts")


print(dim(counts_ns))
print(length(tf_ns))
counts_ns$TF <- tf_ns

counts_ns <- counts_ns[c("Matrix_ID","Counts","TF")]
counts_ns$Counts <- as.integer(counts_ns$Counts)
print(head(counts_ns))
print(str(counts_ns))

write.table(counts,file="JASPAR_counts/JASPAR_TF_counts_RRBS.txt",quote = F,row.names = F,col.names = T)
write.table(counts_ns,file="JASPAR_counts/JASPAR_TF_counts_RRBS_ns.txt",quote = F,row.names = F,col.names = T)



# Comparison
################


TF_comp_h <- as.data.frame(cbind(rep("",length(unique(c(counts$TF,counts_ns$TF)))),
                               rep(0,length(unique(c(counts$TF,counts_ns$TF)))),
                               rep(0,length(unique(c(counts$TF,counts_ns$TF))))),stringsAsFactors=F)
                               # rep("",length(unique(c(counts$TF,counts_ns$TF))))),

for (t in 1:length(unique(c(counts$TF,counts_ns$TF)))){
  TF_comp_h[t,1] <- unique(c(counts$TF,counts_ns$TF))[t]
  
  if (unique(c(counts$TF,counts_ns$TF))[t]%in%counts$TF){
    TF_comp_h[t,2] <- counts$Counts[counts$TF==unique(c(counts$TF,counts_ns$TF))[t]]
    # TF_comp_h[t,4] <- counts$Seq[counts$TF==unique(c(counts$TF,counts_ns$TF))[t]]
  }
  if (unique(c(counts$TF,counts_ns$TF))[t]%in%counts_ns$TF){
    TF_comp_h[t,3] <- counts_ns$Counts[counts_ns$TF==unique(c(counts$TF,counts_ns$TF))[t]]
  }
}
# colnames(TF_comp_h) <- c("TF","Counts","Counts_ns","Seq")
colnames(TF_comp_h) <- c("TF","Counts","Counts_ns")
TF_comp_h$Counts <- as.integer(TF_comp_h$Counts)
TF_comp_h$Counts_ns <- as.integer(TF_comp_h$Counts_ns)
print(head(TF_comp_h))
print(str(TF_comp_h))

print(s)
print(ns)
# binom exact test
TF_comp_h$class <- ""
TF_comp_h$family <- ""
TF_comp_h$p_value <- 1
TF_comp_h$direction <- ""
for (i in 1:nrow(TF_comp_h)){
  mat <- rbind(c(TF_comp_h$Counts[i],s-TF_comp_h$Counts[i]),c(TF_comp_h$Counts_ns[i],ns-TF_comp_h$Counts_ns[i]))
  # print(mat)
  TF_comp_h$p_value[i] <- chisq.test(mat)$p.value
  TF_comp_h$class[i] <- annot$Class[annot$Name==TF_comp_h$TF[i]][1]
  TF_comp_h$family[i] <- annot$Family[annot$Name==TF_comp_h$TF[i]][1]
  # print(str(chisq.test(mat,s/(s+ns))))
  if(chisq.test(mat)$observed[1] > chisq.test(mat)$expected[1]){
    TF_comp_h$direction[i] <- "enriched"
  } else {
    TF_comp_h$direction[i] <- "depleted"
  }
  
}

TF_comp_h <- TF_comp_h[sort(TF_comp_h$p_value,index.return=T)$ix,]
TF_comp_h$FDR <- p.adjust(TF_comp_h$p_value,"BH",length(TF_comp_h$p_value))
colnames(TF_comp_h) <- c("TF",paste("Counts (n=",s,")",sep=""),paste("Counts_ns (n=",ns,")",sep=""),colnames(TF_comp_h)[-c(1:3)])
print(head(TF_comp_h))
logFC <- log2((TF_comp_h[,2]/s)/(TF_comp_h[,3]/ns))
logFC[is.infinite(logFC) | is.na(logFC)] <- 0
TF_comp_h$logFC <- logFC
TF_comp_h$FC <- 2^TF_comp_h$logFC
TF_comp_h$filter <- (abs(TF_comp_h$logFC)>1) & (TF_comp_h$FDR<0.10)


write.table(TF_comp_h,file="JASPAR_output/JASPAR_TF_RRBS_comp_h.txt",quote = F,row.names = F,col.names = T,sep="\t")



