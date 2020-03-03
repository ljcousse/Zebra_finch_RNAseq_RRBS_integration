########################
########################
##                    ##
##  Install packages  ##
##                    ##
########################
########################

chooseCRANmirror(ind=8)
list_packages <- c("UniProt.ws","stringr","edgeR","tidyverse","stageR","biomaRt","RCircos","BiSeq",
                   "GenomicRanges","rtracklayer","limma","RColorBrewer","gplots","GenomicFeatures",
                   "xlsx")
for (i in 1:length(list_packages)){
  if (!(list_packages[i] %in% rownames(installed.packages()))) {
    BiocManager::install(list_packages[i],update=F)
  }
}


