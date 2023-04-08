###################################
# ExtractSeq_TFsearching_RNAseq.R #
###################################


# Loading libraries
####################

library("biomaRt")
library("parallel")
options(warn=1)

working_dir <- getwd()
data_dir <- "/share/louisc/Porthos_backup/louisc/Project_Plastoscine/"


# Loading Contrast of Interest
#################################

coi <- c("1dph","20dph","65dph","adult","age")
  
# Loading Chromosomes and defining functions
##############################################

setwd(data_dir)

chromosome <- c(1, "1A", "1B", 2, 3, 4, "4A", 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 'LGE22', 'Z', 'LG2', 'LG5', 'MT')

chr <- list()
for (i in 1:length(chromosome)){
  chr[[i]] <- scan(paste("Taeniopygia_guttata.taeGut3.2.4.dna.chromosome.",chromosome[i],".fa",sep=""),what="character",skip=1)
  
  if(sum(nchar(chr[[i]])==60)!=(length(chr[[i]])-1)){
    print(sum(nchar(chr[[i]])==60))
    print(length(chr[[i]]))
    stop("For chromosome",chromosome[i],"all contigs are not of equal length")
  }
}
print(length(chr))

getPromSeq <- function(pos,chrom,sense){
  # print("new")
  chrom_standard <- c(1, "1A", "1B", 2, 3, 4, "4A", 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 'LGE22', 'Z', 'LG2', 'LG5', 'MT')
  if (!(sum(chromosome%in%chrom_standard)==36 & length(chromosome)==36)){
    stop("Vector chromosome should be defined as c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 'LGE22', 'Z','LG2', 'LG5', 'MT') for this function")
  }
  # print(pos)
  # print(chrom)
  # print(which(chromosome==chrom))
  chr_sp <- chr[[which(chromosome==chrom)]]
  # print(length(chr_sp))
  contig_len <- nchar(chr[[which(chromosome==chrom)]][1])
  if (sense){
    start <- pos - 2000
    if (start < 0){
      start <- 1
    }
    end <- pos + 500
  } else {
    start <- pos - 500
    if (start < 0){
      start <- 1
    }
    end <- pos + 2000
  }
  
  # print(start)
  # print(end)
  base_st <- start %% contig_len
  base_end <- end %% contig_len
  
  
  if (base_st == 0){
    contig_st <- ((start-base_st)/contig_len)
    base_st <- contig_len
  } else {
    contig_st <- ((start-base_st)/contig_len)+1
  }
  if (base_end == 0){
    contig_end <- ((end-base_end)/contig_len)
    base_end <- contig_len
  } else {
    contig_end <- ((end-base_end)/contig_len)+1
  }
  if (contig_end > length(chr_sp)){
    contig_end <- length(chr_sp)
    base_end <- nchar(chr_sp[contig_end])
  }
  # print(base_st)
  # print(base_end)
  # print(contig_st)
  # print(contig_end)
  # print(chr_sp[(contig_st-1):(contig_end+1)])
  if(contig_st == contig_end){
    seq <- paste(unlist(strsplit(chr_sp[contig_st],""))[base_st:base_end],collapse="")
  } else {
    seq_st <- unlist(strsplit(chr_sp[contig_st],""))[base_st:contig_len]
    seq_end <- unlist(strsplit(chr_sp[contig_end],""))[1:base_end]
    if (contig_st == (contig_end - 1)){
      seq <- paste(paste(seq_st,collapse=""),paste(seq_end,collapse=""),sep="")
    } else {
      seq <- paste(paste(seq_st,collapse=""),paste(chr_sp[(contig_st+1):(contig_end-1)],collapse=""),paste(seq_end,collapse=""),sep="")
    }
  }
  # print(seq)
  if (grepl("NANANANANANANANANANA",seq)){
    print(chrom)
    print(start)
    print(end)
    print(contig_st)
    print(contig_end)
    print(length(chr_sp))
  }
  return(seq)
}

getSeq <- function(pos,chrom){
  
  chrom_standard <- c(1, "1A", "1B", 2, 3, 4, "4A", 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 'LGE22', 'Z', 'LG2', 'LG5', 'MT')
  if (!(sum(chromosome%in%chrom_standard)==33 & length(chromosome)==33)){
    stop("Vector chromosome should be defined as c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 'LGE22', 'Z', 'LG2', 'LG5', 'MT') for this function")
  }
  # print(chrom)
  # print(which(chromosome==chrom))
  chr_sp <- chr[[which(chromosome==chrom)]]
  contig_len <- nchar(chr[[which(chromosome==chrom)]][1])
  start <- pos - 50
  end <- pos + 50
  # print(start)
  # print(end)
  base_st <- start %% contig_len
  base_end <- end %% contig_len
  
  if (base_st == 0){
    contig_st <- ((start-base_st)/contig_len)
    base_st <- contig_len
  } else {
    contig_st <- ((start-base_st)/contig_len)+1
  }
  if (base_end == 0){
    contig_end <- ((end-base_end)/contig_len)
    base_end <- contig_len
  } else {
    contig_end <- ((end-base_end)/contig_len)+1
  }
  # print(base_st)
  # print(base_end)
  # print(contig_st)
  # print(contig_end)
  # print(chr_sp[(contig_st-1):(contig_end+1)])
  if(contig_st == contig_end){
    seq <- paste(unlist(strsplit(chr_sp[contig_st],""))[base_st:base_end],collapse="")
  } else {
    seq_st <- unlist(strsplit(chr_sp[contig_st],""))[base_st:contig_len]
    seq_end <- unlist(strsplit(chr_sp[contig_end],""))[1:base_end]
    if (contig_st == (contig_end - 1)){
      seq <- paste(paste(seq_st,collapse=""),paste(seq_end,collapse=""),sep="")
    } else {
      seq <- paste(paste(seq_st,collapse=""),paste(chr_sp[(contig_st+1):(contig_end-1)],collapse=""),paste(seq_end,collapse=""),sep="")
    }
  }
  # print(seq)
  return(seq)
}

# Seq extraction 
##################

setwd(working_dir)

# print(listMarts(host="www.ensembl.org"))
if (file.exists("Genes.Rda")){
  load("Genes.Rda")
  genes <- genes[!duplicated(genes$ensembl_gene_id),]
} else {
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org")
  ensembl <- useDataset("tguttata_gene_ensembl",mart=ensembl)
  genes <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id','chromosome_name','start_position','end_position','transcription_start_site'),mart = ensembl)
  print(head(genes))
  save(genes,file="Genes.Rda")
  genes <- genes[!duplicated(genes$ensembl_gene_id),]
}

myfun <- function(j){
  # Get significant differentially expressed gene for the contrast of interest
  file_in <- paste("ListGenes_",coi[j],".txt",sep="")
  print(file_in)
  genes_coi <- read.table(file_in,header=F,stringsAsFactors=F)$V1
  print(head(genes_coi))
  print(length(genes_coi))
  
  
  coi_new <- NULL
  for (i in 1:length(genes_coi)){
    coi_new <- rbind(coi_new,genes[genes$ensembl_gene_id==genes_coi[i],])
  }
  
  
  print(head(coi_new))
  print(dim(coi_new))
  coi_new <- coi_new[!grepl("random",coi_new$chromosome_name),]
  coi_new <- coi_new[!grepl("Un",coi_new$chromosome_name),]
  print(dim(coi_new))
  seq_coi <- NULL
  for (i in 1:nrow(coi_new)){
    seq_coi  <- c(seq_coi,getPromSeq(coi_new$transcription_start_site[i],
                                     as.character(coi_new$chromosome_name[i]),
                                     abs(coi_new$transcription_start_site[i]-coi_new$start_position[i])<abs(coi_new$transcription_start_site[i]-coi_new$end_position[i])))
  }
  seq_coi  <- toupper(seq_coi)
  
  
  for (i in 1:length(seq_coi)){
    if (!grepl("NANANANANANANANANANA",seq_coi[i])){
      # print(i)
      if (i == 1){
        write.table(paste(">",coi_new$ensembl_gene_id[i]),file=paste("seq_",coi[j],"_TFsearching.fa",sep=""),col.names = F,row.names = F, quote = F)
        write.table(seq_coi[i],file=paste("seq_",coi[j],"_TFsearching.fa",sep=""),col.names = F,row.names = F, quote = F, append = T)
      } else {
        write.table(paste(">",coi_new$ensembl_gene_id[i]),file=paste("seq_",coi[j],"_TFsearching.fa",sep=""),col.names = F,row.names = F, quote = F, append = T)
        write.table(seq_coi[i],file=paste("seq_",coi[j],"_TFsearching.fa",sep=""),col.names = F,row.names = F,quote = F, append = T)
      }
    } else {
      print(seq_coi[i])
      warning("Non-sensical sequences!")
      next()
    }
  }
  return(genes_coi)
}
genes_sign <- mclapply(1:length(coi),myfun,mc.cores=length(coi),mc.preschedule=T)
genes_sign <- unique(unlist(genes_sign))

all_genes <- c(as.matrix(read.table("Plastoscine_annot.txt",header=F,stringsAsFactors=F)))
# print(length(all_genes))
# print(tail(all_genes))
all_genes <- all_genes[!grepl('__',all_genes)]
# print(length(all_genes))
# print(tail(all_genes))

coi_ns <- all_genes[!(all_genes%in%genes_sign)]
print(head(coi_ns))
print(length(coi_ns))
print(str(coi_ns))

coi_ns_new <- NULL
for (i in 1:length(coi_ns)){
  coi_ns_new <- rbind(coi_ns_new,genes[genes$ensembl_gene_id==coi_ns[i],])
}

print(head(coi_ns_new))

coi_ns_new <- coi_ns_new[!grepl("random",coi_ns_new$chromosome_name),]
coi_ns_new <- coi_ns_new[!grepl("Un",coi_ns_new$chromosome_name),]

seq_coi_ns <- NULL
for (i in 1:nrow(coi_ns_new)){
  seq_coi_ns  <- c(seq_coi_ns,getPromSeq(coi_ns_new$transcription_start_site[i],
                                         as.character(coi_ns_new$chromosome_name[i]),
                                         abs(coi_ns_new$transcription_start_site[i]-coi_ns_new$start_position[i])<abs(coi_ns_new$transcription_start_site[i]-coi_ns_new$end_position[i])))
}
seq_coi_ns  <- toupper(seq_coi_ns)



for (i in 1:length(seq_coi_ns)){
  if (!grepl("NANANANANANANANANANA",seq_coi_ns[i])){
    # print(i)
    if (i == 1){
      write.table(paste(">",coi_ns_new$ensembl_gene_id[i]),file=paste("seq_RNAseq_ns_TFsearching.fa",sep=""),col.names = F,row.names = F, quote = F)
      write.table(seq_coi_ns[i],file=paste("seq_RNAseq_ns_TFsearching.fa",sep=""),col.names = F,row.names = F, quote = F, append = T)
    } else {
      write.table(paste(">",coi_ns_new$ensembl_gene_id[i]),file=paste("seq_RNAseq_ns_TFsearching.fa",sep=""),col.names = F,row.names = F, quote = F, append = T)
      write.table(seq_coi_ns[i],file=paste("seq_RNAseq_ns_TFsearching.fa",sep=""),col.names = F,row.names = F,quote = F, append = T)
    }
  } else {
    print(seq_coi_ns[i])
    stop("Non-sensical sequences!")
  }
  
}
