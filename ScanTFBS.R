dosage <- read.table("/home/louisc/Zebra_Finch_new/GeneLists/ListGenes_dosage.txt",header=F,stringsAsFactors = F)$V1
dosage_bis <- read.table("/home/louisc/Zebra_Finch_new/GeneLists/ListGenes_dosage_bis.txt",header=F,stringsAsFactors = F)$V1
nodosage <- read.table("/home/louisc/Zebra_Finch_new/GeneLists/ListGenes_nodosage.txt",header=F,stringsAsFactors = F)$V1

fulldata_output <- read.table("/home/louisc/Zebra_Finch_new/FullData_output.txt", header=T, sep="\t", stringsAsFactors = F)

genes_65dph <- fulldata_output$ensembl.gene.id[(fulldata_output$Pvalue.65dph<0.05 & !is.na(fulldata_output$Pvalue.65dph)) &
                                                 (fulldata_output$Pvalue.1dph>0.05 & !is.na(fulldata_output$Pvalue.1dph)) &
                                                 (fulldata_output$Pvalue.20dph>0.05 & !is.na(fulldata_output$Pvalue.20dph)) &
                                                 (fulldata_output$Pvalue.adult>0.05 & !is.na(fulldata_output$Pvalue.adult)) &
                                                 (fulldata_output$Pvalue.age>0.05 & !is.na(fulldata_output$Pvalue.age))]

IDs <- list.files()
IDs <- IDs[grepl("meme",IDs)]
IDs <- unlist(strsplit(IDs,"_"))[c(1:length(IDs))*2]
IDs <- unlist(strsplit(IDs,"\\.meme"))

# FOXP2
system("cat fimo.txt | grep MA0593.1 | awk '{print $2'} > test.txt")
# OTX2
system("cat fimo.txt | grep MA0712.1 | awk '{print $2'} > test.txt")
# FOXO3
system("cat fimo.txt | grep MA0157.1 | awk '{print $2'} > test.txt")


perc_dosage <- cbind(rep(0,length(IDs)),rep(0,length(IDs)))
perc_65dph <- cbind(rep(0,length(IDs)),rep(0,length(IDs)))
for (i in 1:length(IDs)){
  system(paste("cat fimo.txt | grep ",IDs[i]," | awk '{print $2'} > test.txt",sep=""))
  if (system("cat test.txt | wc -l",intern = T)==0){
    next
  }
  if (i %%50==0){print(i)}
  test <- read.table("test.txt",header=F,stringsAsFactors=F)$V1
  
  test <- unique(sort(test))
  length(test)
  
  perc_dosage[i,] <- c(sum(test%in%dosage),length(test))
  perc_65dph[i,] <- c(sum(test%in%genes_65dph),length(test))
}

foxy <- c("ENSTGUG00000002129","ENSTGUG00000004163","ENSTGUG00000004166","ENSTGUG00000005315","ENSTGUG00000005846","ENSTGUG00000007282","ENSTGUG00000007654","ENSTGUG00000007783","ENSTGUG00000012135")
sum(test%in%dosage_bis)
sum(test%in%nodosage)


test <- unique(sort(test))
length(test)

sum(fulldata_output$Pvalue.age[fulldata_output$ensembl.gene.id%in%test & !is.na(fulldata_output$Pvalue.age)]<0.05)

sum(fulldata_output$logFC[fulldata_output$ensembl.gene.id%in%test & !is.na(fulldata_output$Pvalue.age)]>0)

chisq.test(cbind(c(1872,2488),c(10193,13891)))

test_GS <- fulldata_output$gene.symbol[fulldata_output$ensembl.gene.id%in%test]
test_GS <- test_GS[test_GS!=""]
sort(test_GS)




