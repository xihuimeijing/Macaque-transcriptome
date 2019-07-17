#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

usage=function(){
  cat('Usage:cor.R -i1=<file1> [-i2=<file2>] -s1=1,5 -s2=6,10 -m=<pearson>\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-i1\tFILE\tThe first input file.\n',file=stderr())
  cat('\t-i2\tFILE\tThe second input file if your input is from two files.File1 and file2 must have the same rows[optional]\n',file=stderr())
  cat('\t-s1\tSTRING\tThe column number separated by comma for the first sample\n',file=stderr())
  cat('\t-s2\tSTRING\tThe column number separated by comma for the second sample\n',file=stderr())
  cat('\t-m\tSTRING\tA character string specifying which correlation coefficient is to be computed("pearson" (default), "kendall", or "spearman")\n',file=stderr())
  cat('\t-o\tFILE\tThe output file name.\n',file=stderr())
  q(save="no")  
}
method="pearson"
if(length(args)>=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(grepl('^-i1=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      inFile1=arg.split[2]
    }
    if(grepl('^-i2=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      inFile2=arg.split[2]
    }
    if(grepl('^-m=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      method=arg.split[2]
    }
    if(grepl('^-s1=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      sample1=arg.split[2]
    }
    if(grepl('^-s2=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      sample2=arg.split[2]
    }
    if(grepl('^-o=', arg)){
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

if(exists("inFile1")){
  data1<-read.delim(file=inFile1,header=F,comment.char='#')
}else{
  data1<-read.delim(file('stdin'), header = F)
}
sample1=as.numeric(strsplit(sample1,",")[[1]])
sample2=as.numeric(strsplit(sample2,",")[[1]])
cor=''
if(exists("inFile2")){
  data2<-read.delim(file=inFile2,header=F)
  for(i in 1:nrow(data2)){
    vector1=as.numeric(data1[i,sample1[1]:sample1[2]])
    vector2=as.numeric(data2[i,sample2[1]:sample2[2]])
    cor[i]=round(cor(vector1,vector2,method=method,use="pairwise.complete.obs"),digits = 4)
  }
}else{
  for(i in 1:nrow(data1)){
    vector1=as.numeric(data1[i,sample1[1]:sample1[2]])
    vector2=as.numeric(data1[i,sample2[1]:sample2[2]])
    cor[i]=round(cor(vector1,vector2,method=method,use="pairwise.complete.obs"),digits = 4)
    }
}
write.table(cor,file = out, sep = "\t",quote = F,row.names = F,col.names = F)