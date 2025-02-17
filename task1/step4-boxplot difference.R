##用箱线图查看数据差异性
### ---------------
rm(list=ls())
options(stringsAsFactors = F)#在数据框中的字符向量不再被自动转换为因子类型

#加载
load(file = '../save/step1-output.Rdata')

#想比较的数据
ids[ids$symbol=='TRAF4',]
dat[1:4,1:4]
dat=log2(dat)  +0.1  
dat=data.frame(gene= dat['211899_s_at',] ,
               mut= group_list) 
head(dat)

if(require('ggpubr')){
  library(ggpubr)
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


