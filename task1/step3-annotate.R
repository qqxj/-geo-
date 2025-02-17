##基因注释和数据差异
### ---------------
rm(list=ls())
options(stringsAsFactors = F)#在数据框中的字符向量不再被自动转换为因子类型


#加载
load(file='../save/GSE3325_raw_exprSet.Rdata')

#将基因注释有两种方法(已知):常用第二种方法
#第一种方法(时间非常慢)：直接下载平台注释
if(F){
  library(GEOquery)
  gpl <- getGEO('GPL570', destdir=".", 
                AnnotGPL = F,    
                getGPL = F)
  colnames(Table(gpl))  
  head(Table(gpl)[,c(1,15)]) ##检查这里，哪一列是注释
  ids=Table(gpl)[,c(1,15)]
  head(ids)
  save(ids,file='ids.Rdata')
}
load(file='../save/ids.Rdata')

#修改symbol
library(stringr) 
ids$s=str_split(probe2gene$gene_assignment,' // ',simplify = T)[,2]

ids=ids[,c(1,3)]
colnames(ids)=c('probe_id','symbol') 

#第二种方法：在csdn网站直接搜索注释包
library(hgu133plus2.db)
ls("package:hgu133plus2.db")#查看
ids=toTable(hgu133plus2SYMBOL)
head(ids)
colnames(ids)=c('probe_id','symbol') 

#将dat中的探针换成注释
if(T){
ids=ids[ids$symbol != '',]
ids=ids[ids$probe_id %in%  rownames(dat),]
ids$probe_id=as.character(ids$probe_id)
dat[1:4,1:4]
dat=dat[ids$probe_id,] 
ids$median=apply(dat,1,median) #apply，1代表行，2代表列
ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$symbol),]#去除重复的gene
dat=dat[ids$probe_id,] 
rownames(dat)=ids$symbol
dat[1:4,1:4]  #保留每个基因ID第一次出现的信息
}

save(dat,ids,group_list,file = '../save/step1-output.Rdata')



