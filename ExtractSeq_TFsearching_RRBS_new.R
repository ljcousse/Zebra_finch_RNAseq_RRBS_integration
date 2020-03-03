###################################
# ExtractSeq_TFsearching_RNAseq.R #
###################################


# Loading libraries
####################

library("biomaRt")
library("parallel")
options(warn=2)

working_dir <- getwd()
data_dir <- "Genomic_data/"

# Loading Contrast of Interest
#################################

coi <- c("age")

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

getSeq <- function(pos,chrom){
  
  chrom_standard <- c(1, "1A", "1B", 2, 3, 4, "4A", 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 'LGE22', 'Z', 'LG2', 'LG5', 'MT')
  if (!(sum(chromosome%in%chrom_standard)==36 & length(chromosome)==36)){
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

for (j in 1:length(coi)){
  rrbs <- read.table(file=paste("RRBS_",coi[j],"_RNAseq.txt",sep=""),header=T,stringsAsFactors = F,sep="\t")
  print(dim(rrbs))
  
  rrbs <- rrbs[!grepl("random",rrbs$chr),]
  rrbs <- rrbs[!grepl("Un",rrbs$chr),]
  
  rrbs <- rrbs[rrbs$func=="promoter",]
  
  print(dim(rrbs))
  
  # Seq extraction
  seq <- NULL
  for (i in 1:length(rrbs$locations)){
    if (nchar(getSeq(rrbs$locations[i],as.character(rrbs$chr[i])))!=101){
      print(rrbs$locations[i])
      print(rrbs$chr[i])
      print(nchar(getSeq(rrbs$locations[i],as.character(rrbs$chr[i]))))
    }
    seq <- c(seq,getSeq(rrbs$locations[i],as.character(rrbs$chr[i])))
  }
  seq <- toupper(seq)
  print(head(seq))

  for (i in 1:length(seq)){
    if (!grepl("NANANANANANANANANANA",seq[i])){
      # print(i)
      if (i == 1){
        write.table(paste(">",paste(rrbs$chr[i],rrbs$locations[i],sep="_")),file="seq_RRBS_TFsearching.fa",col.names = F,row.names = F, quote = F)
        write.table(seq[i],file="seq_RRBS_TFsearching.fa",col.names = F,row.names = F, quote = F, append = T)
      } else {
        write.table(paste(">",paste(rrbs$chr[i],rrbs$locations[i],sep="_")),file="seq_RRBS_TFsearching.fa",col.names = F,row.names = F, quote = F, append = T)
        write.table(seq[i],file="seq_RRBS_TFsearching.fa",col.names = F,row.names = F,quote = F, append = T)
      }
    } else {
      stop("Non-sensical sequences!")
    }
  }
}
