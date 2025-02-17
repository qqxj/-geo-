##下载基因，分组
### ---------------
rm(list=ls())
options(stringsAsFactors = F)#在数据框中的字符向量不再被自动转换为因子类型

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3325
if(F){
  library(GEOquery)
  gset <- getGEO('GSE3325', destdir=".",#下载数据保存到当前目录
                 AnnotGPL = F,#不获取平台注释信息
                 getGPL = F)#不获取的平台信息
  save(gset,file='/home/datahup/syj/GEO/GEO-master/GEO/save/GSE3325_eSet.Rdata')
}
load('../save/GSE3325_eSet.Rdata')

class(gset)
length(gset)
class(gset[[1]])
 a = gset[[1]]
dat= exprs(a) #提取矩阵
boxplot(dat,las=2)
phe=pData(a)#提取临床信息

#用临床数据分组
library(stringr)
#修改
#各种切割格式(可供选择)
group_list= str_split(as.character(phe$title),' ',simplify = T)[,1] 
#group_list=str_split(str_split(pd$title,',',simplify = T)[,2],' ',simplify = T)[,3]
#group_list=paste0('group',group_list,'h')
#group_list= ifelse(grepl('normal',pd$title),'normal','npc')
#group_list=ifelse(grepl('t',as.character(meta$title)),'tumor','normal')
table(group_list)

#将'GSE3325_raw_exprSet.Rdata'改成下载的数据集
save(dat,group_list,
     file='../save/GSE3325_raw_exprSet.Rdata')


