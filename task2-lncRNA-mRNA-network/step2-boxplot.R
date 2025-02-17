##
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
setwd("D:/RCODE/GEO/GEO-master/task2-lncRNA-mRNA-network")
rm(list=ls())
if(F){
  library(GEOquery)
  eSet <- getGEO('GSE3325', destdir=".",#获取基因
                 AnnotGPL = F,#平台注释信息
                 getGPL = F)#平台信息
  # Sys.setenv("VROOM_CONNECTION_SIZE"=131072*2)
  save(eSet,file='GSE3325_eSet.Rdata')
}
load('GSE3325_eSet.Rdata')
b = eSet[[1]]
raw_exprSet= exprs(b) #样本和探针
phe=pData(b)#提取临床信息
library(stringr)
group_list= str_split(as.character(phe$title),' ',simplify = T)[,1] 
save(raw_exprSet,group_list,
     file='GSE3325_raw_exprSet.Rdata')

rm(list=ls())
load(file='GSE3325_raw_exprSet.Rdata')
library(hgu133plus2.db)
eg2probe=toTable(hgu133plus2SYMBOL)#提取出来探针和基因名
a <- eg2probe[eg2probe$symbol=='TRAF4',]
raw_exprSet[1:4,1:4]
exprSet=log2(raw_exprSet)#缩小表达值
unique(exprSet['211899_s_at',])
dat=data.frame(gene= exprSet['211899_s_at',] ,
               mut= group_list)
head(dat)
if(require('ggpubr')){
  library(ggpubr)
  # google search : ggpubr boxplot add p-value
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  p <- ggboxplot(dat, x = "mut", y = "gene",#箱线图
                 color = "mut", 
                 palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means(method = "t.test")
print(p)
  }

if(require('ggstatsplot')){
  library(ggstatsplot)
  ggbetweenstats(data = dat, #瓶装图
                 x = mut,  y = gene)
}


if(require('ggplot2')){
  library(ggplot2)
  ggplot(dat,aes(x=mut,y=gene))+
    geom_boxplot()+
    theme_bw()
}



