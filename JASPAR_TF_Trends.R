###################################
# ExtractSeq_TFsearching_RNAseq.R #
###################################

# Loading Contrast of Interest
#################################

setwd("/home/louisc/Zebra_Finch_new/Trends/")
coi <- list.files()
coi <- coi[grep("genelist",coi)]
coi <- unlist(strsplit(coi,".txt"))
coi <- unlist(strsplit(coi,"genelist_patt"))[unlist(strsplit(coi,"genelist_patt"))!=""]

s <- NULL
ns <- NULL
genes_sign <- NULL
for (i in 1:length(coi)){
  fasta_f <- read.table(paste("seq_",coi[i],"_TFsearching.fa",sep=""),stringsAsFactors = F,sep="\t")$V1
  fasta_f <- fasta_f[grep(">",fasta_f)]
  s <- c(s,length(fasta_f))
  fasta_f <- read.table(paste("seq_",coi[i],"_ns_TFsearching.fa",sep=""),stringsAsFactors = F,sep="\t")$V1
  fasta_f <- fasta_f[grep(">",fasta_f)]
  ns <- c(ns,length(fasta_f))
}


# Annotation extraction
#################################

setwd("/home/louisc/Zebra_Finch_new/")
annot <- read.table("MATRIX.txt",stringsAsFactors = F,fill=T)
print(head(annot))
species <- read.table("MATRIX_SPECIES.txt",stringsAsFactors = F,fill=T)
print(head(species))
species <- species[!duplicated(species$V1),]
annot2 <- read.table("MATRIX_ANNOTATION.txt",stringsAsFactors = F,sep="\t",quote='')
print(head(annot2))
print(dim(annot2))
print(str(annot2))
class <- annot2[grepl("class",annot2$V2),c(1,3)]
print(head(class))
class <- class[!duplicated(class$V1),]
family <- annot2[grepl("family",annot2$V2),c(1,3)]
print(head(family))
family <- family[!duplicated(family$V1),]
medline <- annot2[grepl("medline",annot2$V2),c(1,3)]
print(head(medline))
medline <- medline[!duplicated(medline$V1),]
print(sum(duplicated(annot$V1)))
annot_new <- as.data.frame(cbind(rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot)),rep("",nrow(annot))),stringsAsFactors=F)
colnames(annot_new) <- c("ID","M_ID","Gene_Symbol","Species","Family","Class","Medline")
for (i in 1:nrow(annot)){
  annot_new$ID[i] <- annot[i,1]
  annot_new$M_ID[i] <- annot[i,3]
  annot_new$Gene_Symbol[i] <- annot[i,5]
  if (!identical(species$V2[species$V1==annot$V1[i]],character(0))){
    annot_new$Species[i] <- species$V2[species$V1==annot$V1[i]]
  }
  if (!identical(family$V3[family$V1==annot$V1[i]],character(0))){
    annot_new$Family[i] <- family$V3[family$V1==annot$V1[i]]
  }
  if (!identical(class$V3[class$V1==annot$V1[i]],character(0))){
    annot_new$Class[i] <- class$V3[class$V1==annot$V1[i]]
  }
  if (!identical(medline$V3[medline$V1==annot$V1[i]],character(0))){
    annot_new$Medline[i] <-medline$V3[medline$V1==annot$V1[i]]
  }
}
annot <- annot_new
print(head(annot))



# List extraction 
##################

setwd("/home/louisc/Zebra_Finch_new/Trends/")
for (j in 1:length(coi)){
  fimo <- read.table(paste("fimo_",coi[j],"_out/fimo.txt",sep=""),stringsAsFactors = F)
  print(head(fimo))

  fimo <- fimo[fimo$V8<0.05,]
  fimo$V1 <- unlist(strsplit(fimo$V1,"\\.[0-9]"))
  counts <- NULL
  matches <- list()
  for (t in 1:length(unique(fimo$V1))){
    matches[[t]] <- fimo[fimo$V1==unique(fimo$V1)[t],c(2,3,4,9)]
    colnames(matches[[t]]) <- c("Gene","start","stop","seq")
    sub <- fimo[fimo$V1==unique(fimo$V1)[t],]
    sub <- sub[!duplicated(sub$V2),]
    counts <- rbind(counts,c(as.character(unique(fimo$V1)[t]),nrow(sub),as.character(fimo$V9[!(duplicated(fimo$V1))][t])))
  }
  matches <- matches[sort(as.integer(counts[,2]),index.return=T,decreasing=T)$ix]
  counts <- counts[sort(as.integer(counts[,2]),index.return=T,decreasing=T)$ix,]
  
  # counts[,2] <- as.integer(counts[,2])/1003
  
  tf <- NULL
  for (t in 1:nrow(counts)){
    if (sum(annot$M_ID==counts[t,1])>1){
      if (sum(toupper(annot$Gene_Symbol[annot$M_ID==counts[t,1]])==toupper(annot$Gene_Symbol[annot$M_ID==counts[t,1]][1]))==sum(annot$M_ID==counts[t,1])){
        tf <- c(tf,annot$Gene_Symbol[annot$M_ID==counts[t,1]][1])
      } else {
        tf <- c(tf,paste(annot$Gene_Symbol[annot$M_ID==counts[t,1]],collapse=","))
      }
    } else {
      tf <- c(tf,annot$Gene_Symbol[annot$M_ID==counts[t,1]])
    }
    
  }
  
  counts <- as.data.frame(counts,stringsAsFactors=F)
  colnames(counts) <- c("Matrix_ID","Counts","Seq")
  
  print(dim(counts))
  print(length(tf))
  counts$TF <- tf
  
  names(matches) <- tf
  
  save(matches,file=paste("Matches_",coi[j],".Rda",sep=""))
  
  counts <- counts[c("Matrix_ID","Counts","TF","Seq")]
  counts$Counts <- as.integer(counts$Counts)
  print(head(counts))
  print(str(counts))
  
  write.table(counts,file=paste("JASPAR_TF_counts_",coi[j],".txt",sep=""),quote = F,row.names = F,col.names = T)
  
  # Non significant loci
  ######################
  
  
  fimo_ns <- read.table(paste("fimo_",coi[j],"_ns_out/fimo.txt",sep=""),stringsAsFactors = F)
  print(head(fimo_ns))
  
  fimo_ns <- fimo_ns[fimo_ns$V8<0.05,]
  fimo_ns$V1 <- unlist(strsplit(fimo_ns$V1,"\\.[0-9]"))
  counts_ns <- NULL
  matches_ns <- list()
  for (t in 1:length(unique(fimo_ns$V1))){
    matches_ns[[t]] <- fimo_ns[fimo_ns$V1==unique(fimo_ns$V1)[t],c(2,3,4,9)]
    colnames(matches_ns[[t]]) <- c("Gene","start","stop","seq")
    sub_ns <- fimo_ns[fimo_ns$V1==unique(fimo_ns$V1)[t],]
    sub_ns <- sub_ns[!duplicated(sub_ns$V2),]
    counts_ns <- rbind(counts_ns,c(as.character(unique(fimo_ns$V1)[t]),nrow(sub_ns)))
  }
  matches <- matches_ns[sort(as.integer(counts_ns[,2]),index.return=T,decreasing=T)$ix]
  counts_ns <- counts_ns[sort(as.integer(counts_ns[,2]),index.return=T,decreasing=T)$ix,]
  
  # counts[,2] <- as.integer(counts[,2])/1003
  
  tf_ns <- NULL
  prev_len <- 0
  for (t in 1:nrow(counts_ns)){
    if (sum(annot$M_ID==counts_ns[t,1])>1){
      if (sum(toupper(annot$Gene_Symbol[annot$M_ID==counts_ns[t,1]])==toupper(annot$Gene_Symbol[annot$M_ID==counts_ns[t,1]][1]))==sum(annot$M_ID==counts_ns[t,1])){
        tf_ns <- c(tf_ns,annot$Gene_Symbol[annot$M_ID==counts_ns[t,1]][1])
      } else {
        tf_ns <- c(tf_ns,paste(annot$Gene_Symbol[annot$M_ID==counts_ns[t,1]],collapse=","))
      }
    } else {
      tf_ns <- c(tf_ns,annot$Gene_Symbol[annot$M_ID==counts_ns[t,1]])
    }
    # print(tf_ns)
    if (length(tf_ns)>=(prev_len+2)){
      stop(paste("error in loop",t))
    }
    prev_len <- length(tf_ns)
  }
  
  
  counts_ns <- as.data.frame(counts_ns,stringsAsFactors=F)
  colnames(counts_ns) <- c("Matrix_ID","Counts")
  
  print(dim(counts_ns))
  print(length(tf_ns))
  counts_ns$TF <- tf_ns
  
  names(matches_ns) <- tf_ns
  if (j == 1){
    save(matches_ns,file="Matches_RNAseq_ns.Rda")
  }
  
  
  counts_ns <- counts_ns[c("Matrix_ID","Counts","TF")]
  counts_ns$Counts <- as.integer(counts_ns$Counts)
  print(head(counts_ns))
  print(str(counts_ns))
  
  write.table(counts_ns,file=paste("JASPAR_TF_counts_RNAseq_ns.txt",sep=""),quote = F,row.names = F,col.names = T)
  
  
  # Comparison
  ################
  
  # print(counts)
  # print(counts_ns)
  print(counts[duplicated(counts$TF),])
  if(sum(duplicated(counts$TF))>0){
    counts$Counts[duplicated(counts$TF,fromLast=T)] <- counts$Counts[duplicated(counts$TF,fromLast=T)] + counts$Counts[duplicated(counts$TF)]
    counts <- counts[-which(duplicated(counts$TF)),]  
    print(dim(counts))
  }
  
  print(counts_ns[duplicated(counts_ns$TF),])
  if(sum(duplicated(counts_ns$TF))>0){
    counts_ns$Counts[duplicated(counts_ns$TF,fromLast=T)] <- counts_ns$Counts[duplicated(counts_ns$TF,fromLast=T)] + counts_ns$Counts[duplicated(counts_ns$TF)]
    counts_ns <- counts_ns[-which(duplicated(counts_ns$TF)),]  
    print(dim(counts_ns))
  }
  print(counts$TF[!(counts$TF%in%counts_ns$TF)])
  print(counts_ns$TF[!(counts_ns$TF%in%counts$TF)])
  
  TF_comp <- as.data.frame(cbind(rep("",length(unique(c(counts$TF,counts_ns$TF)))),
                                 rep(0,length(unique(c(counts$TF,counts_ns$TF)))),
                                 rep(0,length(unique(c(counts$TF,counts_ns$TF)))),
                                 rep("",length(unique(c(counts$TF,counts_ns$TF))))),stringsAsFactors=F)
  
  for (t in 1:length(unique(c(counts$TF,counts_ns$TF)))){
    TF_comp[t,1] <- unique(c(counts$TF,counts_ns$TF))[t]
    
    if (unique(c(counts$TF,counts_ns$TF))[t]%in%counts$TF){
      TF_comp[t,2] <- counts$Counts[counts$TF==unique(c(counts$TF,counts_ns$TF))[t]]
      TF_comp[t,4] <- counts$Seq[counts$TF==unique(c(counts$TF,counts_ns$TF))[t]]
    }
    if (unique(c(counts$TF,counts_ns$TF))[t]%in%counts_ns$TF){
      TF_comp[t,3] <- counts_ns$Counts[counts_ns$TF==unique(c(counts$TF,counts_ns$TF))[t]]
    }
  }
  colnames(TF_comp) <- c("TF","Counts","Counts_ns","Seq")
  TF_comp$Counts <- as.integer(TF_comp$Counts)
  TF_comp$Counts_ns <- as.integer(TF_comp$Counts_ns)
  print(head(TF_comp))
  print(str(TF_comp))
  
  # Select human TF
  TF_comp_h <- NULL 
  for (i in 1:nrow(TF_comp)){
    print(annot$Species[annot$Gene_Symbol==TF_comp$TF[i]])
    if (9606%in%annot$Species[annot$Gene_Symbol==TF_comp$TF[i]]){
      TF_comp_h <- rbind(TF_comp_h,TF_comp[i,])
    }
    
  }

  print(TF_comp_h)

  
  # binom exact test
  TF_comp$class <- ""
  TF_comp$family <- ""
  TF_comp$medline <- ""
  TF_comp$p_value <- 1
  TF_comp$direction <- ""
  for (i in 1:nrow(TF_comp)){
    mat <- matrix(c(TF_comp$Counts[i],s[j]-TF_comp$Counts[i],TF_comp$Counts_ns[i],ns[j]-TF_comp$Counts_ns[i]),nrow=2)
    TF_comp$p_value[i] <- chisq.test(mat,s[j]/(s[j]+ns[j]))$p.value
    TF_comp$class[i] <- annot$Class[annot$Gene_Symbol==TF_comp$TF[i]][1]
    TF_comp$family[i] <- annot$Family[annot$Gene_Symbol==TF_comp$TF[i]][1]
    TF_comp$medline[i] <- annot$Medline[annot$Gene_Symbol==TF_comp$TF[i]][1]
    # print(str(chisq.test(mat,s/(s+ns))))
    if(chisq.test(mat,s[j]/(s[j]+ns[j]))$observed[1] > chisq.test(mat,s[j]/(s[j]+ns[j]))$expected[1]){
      TF_comp$direction[i] <- "enriched"
    } else {
      TF_comp$direction[i] <- "depleted"
    }
    
  }
  
  TF_comp <- TF_comp[sort(TF_comp$p_value,index.return=T)$ix,]
  TF_comp$FDR <- p.adjust(TF_comp$p_value,"BH",length(TF_comp$p_value))
  print(head(TF_comp))


  write.table(TF_comp,file=paste("JASPAR_TF_",coi[j],"_comp.txt",sep=""),quote = F,row.names = F,col.names = T,sep="\t")
  
  # binom exact test: only human
  TF_comp_h$class <- ""
  TF_comp_h$family <- ""
  TF_comp_h$medline <- ""
  TF_comp_h$p_value <- 1
  TF_comp_h$direction <- ""
  for (i in 1:nrow(TF_comp_h)){
    mat <- matrix(c(TF_comp_h$Counts[i],s[j]-TF_comp_h$Counts[i],TF_comp_h$Counts_ns[i],ns[j]-TF_comp_h$Counts_ns[i]),nrow=2)
    TF_comp_h$p_value[i] <- chisq.test(mat,s[j]/(s[j]+ns[j]))$p.value
    TF_comp_h$class[i] <- annot$Class[annot$Gene_Symbol==TF_comp_h$TF[i]][1]
    TF_comp_h$family[i] <- annot$Family[annot$Gene_Symbol==TF_comp_h$TF[i]][1]
    TF_comp_h$medline[i] <- annot$Medline[annot$Gene_Symbol==TF_comp_h$TF[i]][1]
    # print(str(chisq.test(mat,s/(s+ns))))
    if(chisq.test(mat,s[j]/(s[j]+ns[j]))$observed[1] > chisq.test(mat,s[j]/(s[j]+ns[j]))$expected[1]){
      TF_comp_h$direction[i] <- "enriched"
    } else {
      TF_comp_h$direction[i] <- "depleted"
    }
    
  }
  
  TF_comp_h <- TF_comp_h[sort(TF_comp_h$p_value,index.return=T)$ix,]
  TF_comp_h$FDR <- p.adjust(TF_comp_h$p_value,"BH",length(TF_comp_h$p_value))
  print(head(TF_comp_h))
  
  colnames(TF_comp_h) <- c("TF",paste("Counts (n=",s[j],")",sep=""),paste("Counts_ns (n=",ns[j],")",sep=""),colnames(TF_comp_h)[-c(1:3)])
  write.table(TF_comp_h,file=paste("JASPAR_TF_",coi[j],"_comp_h.txt",sep=""),quote = F,row.names = F,col.names = T,sep="\t")

  TF_comp_h_filtered <- TF_comp_h[(TF_comp_h[,2]+TF_comp_h[,3])>3,]
  logFC <- log2((TF_comp_h_filtered[,3]/ns[j])/(TF_comp_h_filtered[,2]/s[j]))
  logFC[is.infinite(logFC) | is.na(logFC)] <- 0
  TF_comp_h_filtered$logFC <- logFC
  TF_comp_h_filtered <- TF_comp_h_filtered[abs(TF_comp_h_filtered$logFC)>0.5,]
  
  write.table(TF_comp_h_filtered,file=paste("JASPAR_TF_",coi[j],"_comp_h_filtered.txt",sep=""),quote = F,row.names = F,col.names = T,sep="\t")
  
}
