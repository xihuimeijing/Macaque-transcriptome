#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Usage:geneRPKM.barplot.R -i=<rpkmFile>\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\t\tFILE\tInput RPKM tsv format file\n',file=stderr())
  cat('\t-o\t\tSTRING\tOUtput prefix\n',file=stderr())
  cat('\t-h\t\tPrint this help information\n',file=stderr())
  q(save="no")
}
if(length(args)>=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      inFile=arg.split[2]
    }
    if(grepl('^-o=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      prefix=arg.split[2]
    }
    if(grepl('^-h',arg)){usage()}
  }
}else{
  usage()
}
library(ggplot2)
library(reshape2)
library(grid)
if(exists("inFile")){
  data<-read.delim(file=inFile, comment.char = "#",header=F)
}else{
  data<-read.delim(file('stdin'), header = F)
}
colnames(data)<-c("Symbol","Species","Adipose","Brain","Cerebellum","Heart","Kidney","Liver","Lung","Muscle","Spleen","Testis")
data<-melt(data,variable.name = "Tissue",value.name = "RPKM")
outBySpecies=paste(prefix,".BySpecies.svg",sep="")
outByTissue=paste(prefix,".ByTissue.svg",sep="")
p1<-ggplot(data=data,aes(x=Species,y=RPKM,fill=Tissue))+geom_bar(stat="identity",color="black", width=0.7,position=position_dodge(0.7))+theme_bw()+theme(
  axis.title.x = element_blank(),axis.title.y=element_text(size=14),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size=14,colour = "black"),axis.text.y = element_text(colour = "black",size = 14),legend.key.size = unit(0.18,"inches"))+
  scale_fill_manual(values=c('#a6bddb','#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9'))
p2<-ggplot(data=data,aes(x=Tissue,y=RPKM,fill=Species))+geom_bar(stat="identity",color="black", width=0.7,position=position_dodge(0.7))+theme_bw()+theme(
  axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size=14,colour = "black"),axis.title.y=element_text(size=14),axis.title.x = element_blank(),axis.text.y =element_text(colour = "black",size = 14) )+
  scale_fill_manual(values=c("#edf8b1","#7fcdbb","#2c7fb8"))

gl <- lapply(list(p1,p2), ggplotGrob)
ht <- do.call(unit.pmax, lapply(gl, "[[", 'heights'))
gl <- lapply(gl, function(x) {    
  x[['heights']] = ht
  x})
svg(file=outBySpecies,width = 6, height = 3)
grid.newpage(); grid.draw(gl[[1]])
dev.off()
svg(file=outByTissue,width = 6, height = 3)
grid.newpage(); grid.draw(gl[[2]])
dev.off()