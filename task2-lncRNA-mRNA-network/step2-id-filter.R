##这段代码主要用于对基因表达数据进行预处理、可视化和聚类分析，以便进一步的生物信息学分析和数据挖掘。
### ---------------
###
### Create: Jianming Zeng
### Date: 2018-07-09 20:11:07
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-07-09  First version
###
### ---------------
rm(list=ls())
setwd("D:/RCODE/GEO/GEO-master/task2-lncRNA-mRNA-network")
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
class(phe$title)
group_list= str_split(as.character(phe$title),' ',simplify = T)[,11] 
group_list=paste0('group',group_list,'h')
save(raw_exprSet,group_list,
     file='GSE17708_raw_exprSet.Rdata')

rm(list=ls()) 
load(file='GSE17708_raw_exprSet.Rdata')
exprSet=raw_exprSet
library(hgu133plus2.db)
ids=toTable(hgu133plus2SYMBOL)#探针和基因
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(sort(table(ids$symbol)))
plot(table(sort(table(ids$symbol))))
#将ids和exprset中的基因进行匹配筛出不重复的
table(rownames(exprSet) %in% ids$probe_id)
dim(exprSet)
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)
ids=ids[match(rownames(exprSet),ids$probe_id),]
head(ids)
exprSet[1:5,1:5]
#自定义函数(将exprSet中的行名转换成ids$symbol)
jimmy <- function(exprSet,ids){
  tmp = by(exprSet,
           ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )
  probes = as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  print(dim(exprSet))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]#匹配的第二列的值
  return(exprSet)
}

new_exprSet <- jimmy(exprSet,ids)
save(new_exprSet,group_list,
     file='GSE17708_new_exprSet.Rdata')
load(file='GSE17708_new_exprSet.Rdata')
exprSet=new_exprSet
if(T){
  library(reshape2)
  head(exprSet)
  exprSet_L=melt(exprSet)
  colnames(exprSet_L)=c('probe','sample','value')
  exprSet_L$group=rep(group_list,each=nrow(exprSet))
  head(exprSet_L)
  ### ggplot2
  library(ggplot2)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+
    geom_boxplot()#箱线图
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+
    geom_violin()#琴图
  print(p)
  p=ggplot(exprSet_L,aes(value,fill=group))+
    geom_histogram(bins = 200)+#绘制直方图
    facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+
    geom_density()+#密度图
    facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+
    geom_density()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+
    geom_boxplot()
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
  print(p) 
  ## hclust
  colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
  # Define nodePar
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                  cex = 0.7, col = "blue")
  hc=hclust(dist(t(exprSet)))#层次聚类
  par(mar=c(5,5,5,10))#设置图形格式
  png('hclust.png',res=120)
  plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
  dev.off()
  
  ## PCA
  library(ggfortify)
  df=as.data.frame(t(exprSet))
  df$group=group_list
  png('pca.png',res=120)
  autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,
           colour = 'group')+theme_bw()
  dev.off()
}


