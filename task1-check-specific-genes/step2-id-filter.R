##这段代码的主要目标是对基因表达数据进行进一步的分析和可视化。使用 ggplot2 包绘制箱线图、小提琴图、直方图和密度图。并根据不同的样本进行分组展示使用层次聚类分析 hclust 对样本进行聚类，并绘制层次聚类树。使用主成分分析（PCA）并绘制 PCA 图。
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
setwd("D:/RCODE/GEO/GEO-master/task1-check-specific-genes")
rm(list=ls()) 
load(file='GSE17708_raw_exprSet.Rdata')

exprSet=raw_exprSet
library(hgu133plus2.db)
ids=toTable(hgu133plus2SYMBOL)   # #将这个对象中的数据转换为表格或数据框的形式
length(unique(ids$symbol)) 
table(ids$symbol)
tail(sort(table(ids$symbol)))    
table(sort(table(ids$symbol)))   #对排序后的频数表进行统计，以查看每个频数出现的次数。这可以帮助你了解数据中有多少个不同的频数。 
plot(table(sort(table(ids$symbol))))

table(rownames(exprSet) %in% ids$probe_id) # %in%判断，在不在
dim(exprSet)                                  
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)                                    
ids=ids[match(rownames(exprSet),ids$probe_id),] 
head(ids)
exprSet[1:5,1:5]

jimmy <- function(exprSet,ids){ 
  tmp = by(exprSet,
           ids$symbol,
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
     file='GSE42872_new_exprSet.Rdata')

load(file='GSE42872_new_exprSet.Rdata')
exprSet=new_exprSet
if(T){
  
  library(reshape2)
  exprSet_L=melt(exprSet)#重塑或重新排列数据集的形状
  #?????ص?һ???????ݵ???״̬??Ȼ???ٸ??ݹ۲?id??ĳ??��???ƣ??????????ݣ?????��??
  colnames(exprSet_L)=c('probe','sample','value')#?᷽??ȡ???????????????? ???? ??????
  exprSet_L$group=rep(group_list,each=nrow(exprSet))#??exprSet_L??group?У???ÿһ????????????Ӧһ?ּ?????չ״̬??
  
  #??ÿһ????????Ӧ??nrow(exprSet)????????????һ?????˵ļ?????չ״̬?ظ?nrow(exprSet)??
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
  ## hclust???????????? ?о???
  colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
  # Define nodePar
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                  cex = 0.7, col = "blue")
  
  hc=hclust(dist(t(exprSet)))
  #dist计算中样本之间的距离矩阵.
  #hclust将距离矩阵作为输入，执行层次聚类分析。
  par(mar=c(5,5,5,10))#调图形参数
  png('hclust.png',res=120)
  plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)#as.dendrogram()??hclust???ɵĶ???ת??Ϊ?????ľ???ͼ
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


