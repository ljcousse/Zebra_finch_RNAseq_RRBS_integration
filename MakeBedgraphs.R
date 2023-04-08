makeBedGraph <- function(dataobj,cn,trackname,description,col=NULL,range=NULL,startview=NULL,isSign=NULL){
  if (!is.null(isSign)){
    if (!is.null(range)){
      range[1] <- -range[2] 
    }
    dataobj[!isSign,cn[4]] <- -dataobj[!isSign,cn[4]]
  }
  
  if (is.null(startview)){
    write.table("browser position chrZ:56,283,356-56,343,168",file=paste(trackname,".BedGraph",sep=""),col.names = F, quote = F, row.names = F)
  } else {
    write.table(paste("browser position ",startview,sep=""),file=paste(trackname,".BedGraph",sep=""),col.names = F, quote = F, row.names = F)
  }
  
  write.table("browser hide all",file=paste(trackname,".BedGraph",sep=""),col.names = F, quote = F, row.names = F, append = T)
  
  if (is.null(range)){
    header = paste('track type=bedGraph name=',trackname,' description=',description,
                   ' visibility=full color=',paste(col2rgb(col[1]),collapse=","),' altColor=',paste(col2rgb(col[2]),collapse=","),' gridDefault=on',sep="")
  } else {
    header = paste('track type=bedGraph name=',trackname,' description=',description,
                   ' visibility=full color=',paste(col2rgb(col[1]),collapse=","),' altColor=',paste(col2rgb(col[2]),collapse=","),
                   ' gridDefault=on viewLimits=',range[1],":",range[2]," autoScale=off",sep="")
  }
  
  write.table(header,file=paste(trackname,".BedGraph",sep=""),col.names = F, quote = F, row.names = F, append = T)
  
  write.table(cbind(paste("chr",dataobj[,cn[1]],sep=""),dataobj[,cn[2:4]]),file=paste(trackname,".BedGraph",sep=""),
              col.names = F, quote = F, row.names = F,sep="\t", append = T)
  
}
rrbs_65dph_full <- rrbs_65dph_full[!grepl("random",rrbs_65dph_full$chr),]
rrbs_65dph_full <- rrbs_65dph_full[!grepl("Un",rrbs_65dph_full$chr),]
rrbs_65dph_full <- rrbs_65dph_full[!grepl("MT",rrbs_65dph_full$chr),]
makeBedGraph(rrbs_65dph_full,c(1,2,2,12),"Female_65dph","Meth_perc",col=c("red","black"),range=c(0,0.3))
makeBedGraph(rrbs_65dph_full,c(1,2,2,14),"Male_65dph","Meth_perc",col=c("red","black"),range=c(0,0.3))
