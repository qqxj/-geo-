##limma包得出差异数据
### ---------------
rm(list=ls())
options(stringsAsFactors = F)#在数据框中的字符向量不再被自动转换为因子类型

#加载
load(file='step1-output.Rdata')

dat[1:4,1:4] 
table(group_list) 
#通过为每个数据集绘制箱形图，比较数据集中的数据分布
boxplot(dat[1,]~group_list) #按照group_list分组画箱线图
bp=function(g){         
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
bp(dat[1,]) ## 调用上面定义好的函数，避免同样的绘图代码重复多次敲。
bp(dat[2,])
bp(dat[3,])
bp(dat[4,])
bp(dat['CDK9',])
dim(dat)

#limma
library(limma)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
head(design)
exprSet=dat
rownames(design)=colnames(exprSet)
design

contrast.matrix<-makeContrasts("Tumor-Normal",
                               levels = design)
##这个矩阵声明，我们要把 Tumor 组跟 Normal 进行差异分析比较
contrast.matrix 

deg = function(exprSet,design,contrast.matrix){
  fit <- lmFit(exprSet,design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)  
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  head(nrDEG)
  return(nrDEG)
}
deg = deg(exprSet,design,contrast.matrix)
head(deg)
save(dat,group_list,deg,file = 'deg.Rdata')




