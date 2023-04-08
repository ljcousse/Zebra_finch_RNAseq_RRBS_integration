########################
########################
##                    ##
##  Install packages  ##
##                    ##
########################
########################

# Bioconductor
chooseCRANmirror(ind=7)
list_packages <- c("UniProt.ws","stringr","edgeR","tidyverse","stageR","biomaRt","RCircos","BiSeq",
                   "GenomicRanges","rtracklayer","limma","RColorBrewer","gplots","GenomicFeatures",
                   "xlsx")
for (i in 1:length(list_packages)){
  if (!(list_packages[i] %in% rownames(installed.packages()))) {
    BiocManager::install(list_packages[i],update=F)
  }
}

# CRAN
list_packages <- c("car","readbitmap")
for (i in 1:length(list_packages)){
  if (!(list_packages[i] %in% rownames(installed.packages()))) {
    install.packages(list_packages[i],update=F)
  }
}

