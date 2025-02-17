##注释和分析数据差异
### ---------------
rm(list=ls())
options(stringsAsFactors = F)#在数据框中的字符向量不再被自动转换为因子类型

#加载
load(file='GSE3325_raw_exprSet.Rdata')

#PCA
table(group_list)
dat[1:4,1:4]
if(T){
dat=t(dat)
dat=as.data.frame(dat)
dat=cbind(dat,group_list) 
library("FactoMineR")
library("factoextra") 
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = dat$group_list, 
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, 
             legend.title = "Groups"
)
ggsave('all_samples_PCA.png')
}

#热图
rm(list = ls()) 
load(file = 'step1-output.Rdata')
dat[1:4,1:4] 
if(t){
cg=names(tail(sort(apply(dat,1,sd)),1000))
library(pheatmap)
pheatmap(dat[cg,],show_colnames =F,show_rownames = F) 
n=t(scale(t(dat[cg,]))) 
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(g=group_list)
rownames(ac)=colnames(n) 
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac,filename = 'heatmap_top1000_sd.png')
}
