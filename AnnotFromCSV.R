all_annot <- NULL
for (i in 1:4){
  part <- read.csv(paste("JASPAR -     Search profile(s) (",i,").csv",sep=""),row.names = NULL,stringsAsFactors = F)
  part <- part[,-c(1,ncol(part))]
  print(head(part))
  print(dim(part))
  all_annot <- rbind(all_annot,part)
}

write.table(all_annot,"JASPAR_annot.txt",col.names = T, row.names = F,quote = F,sep="\t")