#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Usage:barplot.PAmotif.R -i=<file> -names=A,B.. -o=<file.pdf>\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i\tFILE\tThe input file delimited by tab with the first column is PA motif and the following columns are motif number.\n',file=stderr())
  cat('\t-names\tSTRING\tThe names for each bar separated by comma\n',file=stderr())
  cat('\t-combine\tLOGIC\tCombine AUUAAA motif into other\n',file=stderr())
  cat('\t-o\tFILE\tThe output file name.\n',file=stderr())
  q(save="no")  
}
combine=F
if(length(args)>=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      inFile=arg.split[2]
    }
    if(grepl('^-names=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      names=strsplit(arg.split[2],',',fixed = T)[[1]]
    }
    if(grepl('^-combine=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      combine=arg.split[2]
    }
    if(grepl('^-o=',arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      out=arg.split[2]
    }
    if(grepl('^-h', arg)){
      usage()
    }
  }
}else{
  usage()
}
if(exists("inFile")){
  data<-read.delim(file=inFile,comment.char = "#",header=F,na.strings="")
}else{
  data<-read.delim(file('stdin'), header = F,na.strings="")
}
newData<-data[,-1]
rownames(newData)=data[,1]
colnames(newData)=names
x=c(newData[rownames(newData)=="AATAAA",1],sum(newData[,1])-newData[rownames(newData)=="AATAAA",1])
y=c(newData[rownames(newData)=="AATAAA",2],sum(newData[,2])-newData[rownames(newData)=="AATAAA",2])
pvalue=round(fisher.test(cbind(x,y))$p.value,digits = 3)
newData<-as.matrix(newData/colSums(newData))
if(combine == "T"){
  data<-matrix(nrow=3,ncol = ncol(newData))
  colnames(data)<-colnames(newData)
  rownames(data)<-c("AAUAAA","Other","No")
  data[1,]=as.matrix(newData[rownames(newData)=="AATAAA",])
  data[3,]=as.matrix(newData[rownames(newData)=="NA",])
  data[2,]=1-colSums(data[c(1,3),])
  pdf(file=out)
  par(mar=c(5.1,4.1,4.1,8))
  barplot(data,ylab="Fraction",legend.text = T,border = F,col=c("#5B9BD5","#ED7D31","gray"),args.legend = list(x = length(names)+1,y=0.5))
  text(1.5,0.8,paste("p=",pvalue))
  dev.off()
}else{
  data<-matrix(nrow = 4,ncol = ncol(newData))
  colnames(data)<-colnames(newData)
  rownames(data)<-c("AAUAAA","AUUAAA","Other","No")
  data[1,]=as.matrix(newData[rownames(newData)=="AATAAA",])
  data[2,]=as.matrix(newData[rownames(newData)=="ATTAAA",])
  data[4,]=as.matrix(newData[rownames(newData)=="NA",])
  data[3,]=1-colSums(data[c(1,2,4),])
  pdf(file=out)
  par(mar=c(5.1,4.1,4.1,8))
  barplot(data,ylab="Fraction",legend.text = T,border = F,col=c("#5B9BD5","#ED7D31","#FFC000","gray"),args.legend = list(x = length(names)+1,y=0.5))
  text(1.5,0.8,paste("p=",pvalue))
  dev.off()
}