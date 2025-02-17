##这段代码的主要目的是下载和处理两个不同的基因表达数据集，准备数据以供后续分析使用。数据集中的临床信息被提取为 group_list，原始表达值保存在 raw_exprSet 中，并将它们保存为RData文件以备后续使用。
### ---------------
###
### Create: Jianming Zeng
### Date: 
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-07-09  First version
###
### ---------------
setwd("D:/RCODE/GEO/GEO-master/task1-check-specific-genes")
rm(list=ls())
if(F){
  library(GEOquery)
  eSet <- getGEO('GSE3325', destdir=".",#下载数据保存到当前目录
                 AnnotGPL = F,#不获取平台注释信息
                 getGPL = F)#不获取的平台信息
  # Sys.setenv("VROOM_CONNECTION_SIZE"=131072*2)
  save(eSet,file='GSE3325_eSet.Rdata')
}
load('GSE3325_eSet.Rdata')
b = eSet[[1]]
raw_exprSet= exprs(b) #提取这个对象中的原始表达值
phe=pData(b)#提取临床信息
library(stringr)
group_list= str_split(as.character(phe$title),' ',simplify = T)[,1] 
save(raw_exprSet,group_list,
     file='GSE3325_raw_exprSet.Rdata')



rm(list=ls())
if(F){
  library(GEOquery)
  eSet <- getGEO('GSE17708', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE17708_eSet.Rdata')
}
load('GSE17708_eSet.Rdata')
b = eSet[[1]]
raw_exprSet=exprs(b) 
raw_exprSet[1:4,1:4]
phe=pData(b)
phe$title
library(stringr)
group_list= str_split(as.character(phe$title),' ',simplify = T)[,11] 
group_list=paste0('group',group_list,'h')
save(raw_exprSet,group_list,
     file='GSE17708_raw_exprSet.Rdata')


