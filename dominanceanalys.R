## Script for dominance analysis to calculate relative contribution of CDHE drivers
## Written by Akshay Rajeev

rm(list=ls())
library(dominanceanalysis)
library(stringr)

dir<-"/dir/ISMR_data/"
Files<-list.files(path=dir, pattern='', full.names = TRUE)
as=as.list(Files)


d = NULL
for (i in 1:length(Files)){
  print(i)
  filename<- Files[i]
  print(filename)
  data1 <- read.csv(filename, header=FALSE,sep = ' ') 
  #data2 <- as.matrix(data1)
  lm.1<-lm(V3~.,data1) # v2 is runoff here
  da<-dominanceAnalysis(lm.1)
  #dr.boot<-bootDominanceAnalysis(lm.1,R=1000)
  print(da)
  df=data.frame(Pan=da$contribution.average$r2[1],Tanm=da$contribution.average$r2[2],Dsp=da$contribution.average$r2[3],Wsp=da$contribution.average$r2[4])
  #df=data.frame()
  d=rbind(d, df)
  
  write.table(d,"/outdir/outdominance.txt",row.names = FALSE,sep = ' ',col.names = TRUE)
}