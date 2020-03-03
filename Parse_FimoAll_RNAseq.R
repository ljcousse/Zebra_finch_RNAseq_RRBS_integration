library(parallel)

working_dir <- getwd()

fimo_dirs <- list.dirs("fimo_RNAseq_out/")
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



# myfun <- function(n){
#   print(paste("fimo_analyse",n,"started"))
#   data_part <- read.table(file=paste("seq_RNAseq_TFsearching/fimo_",n,"_out/fimo.txt",sep=""),header=F,stringsAsFactors=F,sep="\t")
#   # n <- mean((data_part$V8[!(duplicated(data_part$V8,fromLast=T))]/data_part$V7[!(duplicated(data_part$V8,fromLast=T))]*which(!(duplicated(data_part$V8,fromLast=T))))[1:100])
#   n <- data_part$V8[1]/data_part$V7[1]
#   
#   return(n)
# }
# results_n <- mclapply(1:40,myfun,mc.cores = 20)
# am_BS <- sum(unlist(results_n))

colnames(results) <- c("ID","gene","start","stop","strand","score","Pval","Adj.Pval","Seq")
results$start <- as.numeric(results$start)
results$stop <- as.numeric(results$stop)
results$score <- as.numeric(results$score)
results$Pval <- as.numeric(results$Pval)
results$Adj.Pval <- as.numeric(results$Adj.Pval)
sum(duplicated(paste(results$ID,results$gene,results$start,results$stop,sep="_")))
results <- results[!duplicated(paste(results$ID,results$gene,results$start,results$stop,sep="_")),]
print(sum(duplicated(paste(results$ID,results$gene,results$start,results$stop,sep="_"))))
# results$Adj.Pval <- p.adjust(results$Pval,method = "BH",am_BS)
results <- results[sort(results$Pval,index.return=TRUE)$ix,]
print(head(results))
print(dim(results))

write.table(results,file="fimo_RNAseq_out/fimo.txt",col.names = T, row.names = F, sep="\t", quote = F)


