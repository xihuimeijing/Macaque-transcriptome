#!/bin/env Rscript
##### Rscript Rplot.R hCompare/plot.hBrain-hCere.grey hCompare/plot.hBrain-hCere.grey hCompare/plot.hBrain-hCere.red hCere hBrain hBrain-hCere.pdf
args<-commandArgs(T)
datagrey<-read.table(args[1],header=T)
datablue<-read.table(args[2],header=T)
datared<-read.table(args[3],header=T)
x=args[4]
y=args[5]
pdf=args[6]
pdf(pdf)
plot(datagrey$cere,datagrey$brain,xlab = x , ylab = y ,pch=20,col="grey",ylim=c(3,14),xlim=c(3,14),cex=0.2)
par(new=TRUE)
plot (datablue$cere,datablue$brain,pch=20,xlab = "", ylab = "",col="blue",axes = FALSE,ylim=c(3,14),xlim=c(3,14),cex=0.2)
par(new=TRUE)
plot (datared$cere,datared$brain,pch=20,xlab = "", ylab = "",col="red",axes = FALSE,ylim=c(3,14),xlim=c(3,14),cex=0.2)
dev.off()




