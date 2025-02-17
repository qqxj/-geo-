##用ggpolt2画图（箱线图，瓶图，直方图，树形图）
### ---------------
rm(list=ls())
options(stringsAsFactors = F)#在数据框中的字符向量不再被自动转换为因子类型

exprSet=dat
library(hgu133plus2.db)
ids=toTable(hgu133plus2SYMBOL)   # 提取探针和基因
#看下ids$symbol的频数
length(unique(ids$symbol))
table(ids$symbol)
tail(sort(table(ids$symbol)))    
table(sort(table(ids$symbol)))   
plot(table(sort(table(ids$symbol))))
# %in%判断，接下来删除不匹配的
table(rownames(exprSet) %in% ids$probe_id) 
dim(exprSet)                                  
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)                                    
ids=ids[match(rownames(exprSet),ids$probe_id),] 
head(ids)
exprSet[1:5,1:5]
#把ids$symbol替换掉exprset的行名
jimmy <- function(exprSet,ids){ 
  tmp = by(exprSet,ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )
  #从这两个数据中找出行中最大的数据并赋值给tmp
  probes = as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  print(dim(exprSet))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
  return(exprSet)
}
new_exprSet <- jimmy(exprSet,ids)
save(new_exprSet,group_list,
     file='GSE9770_new_exprSet.Rdata')#原42872
#画图
load(file='GSE9770_new_exprSet.Rdata')
exprSet=new_exprSet
if(T){
  library(reshape2)
  exprSet_L=melt(exprSet)#重塑或重新排列数据集的形状
  str(exprSet_L)
  exprSet_L$value=log2(exprSet_L$value)
  colnames(exprSet_L)=c('probe','sample','value')
  exprSet_L$group=rep(group_list,each=nrow(exprSet))#
  head(exprSet_L)
  ### ggplot2
  library(ggplot2)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+
    geom_boxplot()#箱线图
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+
    geom_violin()#小提琴图
  print(p)
  p=ggplot(exprSet_L,aes(value,fill=group))+
    geom_histogram(bins = 200)+#直方图
    facet_wrap(~sample, nrow = 5)#分成五个小图
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
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")#????T????
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
  print(p) 
  
  colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
  
  # Define nodePar
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                  cex = 0.7, col = "blue")
  #层次聚类分析是一种用于探索数据集中样本之间相似性的方法，通常用于聚类分析和可视化。
  hc=hclust(dist(t(exprSet)))
  #dist计算中样本之间的距离矩阵.
  #hclust将距离矩阵作为输入，执行层次聚类分析。
  par(mar=c(5,5,5,10))#调图形参数
  png('hclust.png')
  plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
  dev.off()
  #导出图片出问题
  
  ## PCA
  library(ggfortify)
  df=as.data.frame(t(exprSet))
  df$group=group_list
  png('pca.png',res=120)
  autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,
           colour = 'group')+theme_bw()
  dev.off()
}


