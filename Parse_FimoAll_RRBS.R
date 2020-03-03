library(parallel)


fimo_dirs <- list.dirs("fimo_RRBS_out/")
fimo_dirs <- fimo_dirs[-1]
myfun <- function(n){
  print(paste("fimo_analyse",n,"started"))
  data_part <- read.table(file=paste(fimo_dirs[n],"/fimo.txt",sep=""),header=F,stringsAsFactors=F,sep="\t",comment.char = "")
  if (nrow(data_part)==1){
    return(NULL)
  } else {
    print(head(data_part))
    data_part <- data_part[-1,]
    print(head(data_part))
    return(data_part)
  }
}
results <- mclapply(1:length(fimo_dirs),myfun,mc.cores = 40)
results <- do.call(rbind,results)

colnames(results) <- c("ID","chrom_pos","start","stop","strand","score","Pval","Adj.Pval","Seq")
results$start <- as.numeric(results$start)
results$stop <- as.numeric(results$stop)
results$score <- as.numeric(results$score)
results$Pval <- as.numeric(results$Pval)
results$Adj.Pval <- as.numeric(results$Adj.Pval)

sum(duplicated(paste(results$ID,results$chrom_pos,sep="_")))
results <- results[!duplicated(paste(results$ID,results$chrom_pos,sep="_")),]
print(sum(duplicated(paste(results$ID,results$chrom_pos,sep="_"))))
# results$Adj.Pval <- p.adjust(results$Pval,method = "BH",am_BS)
results <- results[sort(results$Pval,index.return=TRUE)$ix,]
print(head(results))
print(dim(results))

write.table(results,file="fimo_RRBS_out/fimo.txt",col.names = T, row.names = F, sep="\t", quote = F)


