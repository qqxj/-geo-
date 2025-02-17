##检查差异性
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
load(file='GSE3325_raw_exprSet.Rdata')#来自step1中
library(hgu133plus2.db)
eg2probe=toTable(hgu133plus2SYMBOL)#将这个对象中的数据转换为表格或数据框的形式。

#提取指定基因，画图
eg2probe[eg2probe$symbol=='TRAF4',]
raw_exprSet[1:4,1:4]
exprSet=log2(raw_exprSet)    
dat=data.frame(gene= exprSet['211899_s_at',] ,
               mut= group_list) 
head(dat)
if(require('ggpubr')){
  library(ggpubr)
  # google search : ggpubr boxplot add p-value
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  p <- ggboxplot(dat, x = "mut", y = "gene",
                 color = "mut", 
                 palette = "jco",#调色板
                 add = "jitter")#抖动                     
  #  Add p-value
  p + stat_compare_means(method = "t.test")  #添加一个均值比较，使用t检验来比较两组数据的均值差异      
}

if(require('ggstatsplot')){
  library(ggstatsplot)
  ggbetweenstats(data = dat, x = mut,  y = gene)
}


if(require('ggplot2')){
  library(ggplot2)
  ggplot(dat,aes(x=mut,y=gene))+
    geom_boxplot(aes(colour=mut))+  #箱线图 
    theme_bw()#主题
}



