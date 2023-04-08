################
# Info_paper.R #
################

server <- commandArgs(TRUE)
if (identical(server,character(0))){
  stop("Please give in on which server you are running this script")
} else if (!(server=="midas" | server=="porthos")){
  stop("This script is not supposed to run on this server")
}

# options(warn = 2)

if (server=="porthos"){
  library("xlsx")
  
  print(dir())
  folders <- dir()[grep("42",dir())]
  folders <- folders[!grepl("align",folders)]
  print(folders)
  input <- rep(0,24)
  unique_mapped <- rep(0,24)
  for (i in 1:length(folders)) {
    files <- dir(folders[i])[grep("Log.final.out",dir(folders[i]))]
    for (t in 1:length(files)){
      log <- read.table(paste(folders[i],"/",files[t],sep=""),fill=T,sep="\t",stringsAsFactors = F)
      
      JD <- as.integer(unlist(strsplit(unlist(strsplit(files[t],"_"))[1],"-"))[2])
      print(JD)
      print(log)
      input[JD] <- input[JD] + as.integer(log[5,2])
      unique_mapped[JD] <- unique_mapped[JD] + as.integer(log[8,2])  + as.integer(log[23,2])
    }
  }
  output <- cbind(paste("JD",1:24,sep="-"),input,unique_mapped,unique_mapped/input)
  print(output)
  colnames(output) <- c("Sample","Input","Uniquely mapped","Mapping percentage")
  
  write.table(output,file="RNA_info.txt",row.names = F,quote = F,sep="\t")
  
} else if (server=="midas"){
  setwd("/data/sandras")
  print(dir())
  folders <- dir()[grep("run",dir())]
  trim_report_final <- matrix(rep(0,24*3),ncol=3)
  for (t in 1:length(folders)){
    setwd(paste("/data/sandras/",folders[t],"/Data/Intensities/BaseCalls/",sep=""))
    prefixes <- unique(unlist(strsplit(dir(), "_"))[grep("[B][0-9]+",unlist(strsplit(dir(), "_")))])
    # print(prefixes)
    for (s in 1:length(prefixes)){
      print(prefixes[s])
      fn <- paste(prefixes[s],"_R1.fastq.gz_trimming_report.txt",sep="")
      trim_report <- scan(fn,what="character")
      trim_report <- trim_report[which(trim_report=="STATISTICS"):length(trim_report)]
      trim_report <- trim_report[grep("[0-9][0-9][0-9]",trim_report)]
      print(trim_report)
      print(as.integer(paste(unlist(strsplit(prefixes[s],""))[-1],collapse="")))
      trim_report_final[as.integer(paste(unlist(strsplit(prefixes[s],""))[-1],collapse="")),] <- c(as.integer(trim_report[2]),as.integer(trim_report[5]),as.integer(trim_report[2])-as.integer(trim_report[5]))
      # print(trim_report)
    }
  }
  setwd("/home/louisc/Zebra_Finch")
  print(trim_report_final)
  output <- cbind(paste("JD",1:24,sep="-"),trim_report_final)
  print(output)
  colnames(output) <- c("Sample","Total reads","Trimmed reads","Remaining reads")
  
  write.table(output,file="RRBS_info.txt",row.names = F,quote = F,sep="\t")
}
