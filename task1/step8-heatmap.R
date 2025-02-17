##注释和分析数据差异
### ---------------
rm(list=ls())
options(stringsAsFactors = F)#在数据框中的字符向量不再被自动转换为因子类型

if(T){
  load(file = 'step1-output.Rdata')
  dat[1:4,1:4]
  table(group_list)
  x=deg$logFC 
  names(x)=rownames(deg) 
  cg=c(names(head(sort(x),100)),
       names(tail(sort(x),100)))
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) 
  n=t(scale(t(dat[cg,])))
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n) 
  pheatmap(n,show_colnames =F,
           show_rownames = F,
           cluster_cols = F, 
           annotation_col=ac,filename = 'heatmap_top200_DEG.png') 
}

write.csv(deg,file = 'deg.csv')
