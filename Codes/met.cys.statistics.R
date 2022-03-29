# 2021年11月16日
# 计算所有蛋白中met和cys的个数
rm(list=ls())
options(stringsAsFactors = F)
workingdir = "C:/Users/Dell/Documents/R/Oxidative_stress_new"
setwd(workingdir)
memory.limit(400000)
#-------------------------------------------------------------------------------
# 第一个段落： 这段方法是用来读取所有人类蛋白质的序列
#             并且统计某一种氨基酸的个数
#-------------------------------------------------------------------------------
library(stringr)
library(Biostrings) 
s = readDNAStringSet("UP000005640_9606.fasta")
a = readBStringSet("UP000005640_9606.fasta", format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
a = data.frame(a)
names = rownames(a)
names = sapply(strsplit(names,"GN="),function(x) x[2])
names = sapply(strsplit(names," "),function(x) x[1])
seqs = a$a
met_counts = str_count(seqs,pattern = "M")
cys_counts = str_count(seqs,pattern = "C")
result = data.frame(genes = names, met = met_counts, cys = cys_counts)
result$sum = result$met+result$cys
result = result[!is.na(result$genes),]
save(result,file = "met_cys_counts.RData")
#load(file = "met_cys_counts.RData")
#-------------------------------------------------------------------------------
# 第二个段落： 对每个具体的样本算一个统计量：每个表达的蛋白包含的cys和met个数，
#              乘以这个蛋白的表达量
#-------------------------------------------------------------------------------
load(file = "tcga33_gtex18_samples11348v5620_and_genes56189.RData")
common.genes = intersect(rownames(all.tpm.final),result$genes)
result = result[!duplicated(result$genes),]
rownames(result) = result$genes
result = result[,-1]
result = result[common.genes,]
expr = all.tpm.final[common.genes,]
identical(rownames(expr),rownames(result))
expr[expr<1] = 0
all_result = data.frame()
for(i in c(1:ncol(expr)))
{
  if(i %% 100 == 0){
    print(i)
  }
  single_expr = expr[,i]
  single_result = result*single_expr
  single_result = data.frame(colSums(single_result))
  if(dim(all_result)[1]==0)
  {
    all_result = single_result
  }else{
    all_result = cbind(all_result,single_result)
  }
}
colnames(all_result) = colnames(expr)
all_result = t(all_result)
# 这个结果是每个样本计算一个cys或met相关的统计量
save(all_result,file = "met_cys_stat.RData")


genelist = c("SLC7A11","SLC1A4","SLC1A5")
geneexpr = all.tpm.final[genelist,]
geneexpr = data.frame(t(geneexpr))
geneexpr$sum = rowSums(geneexpr)
save(geneexpr,file = "expr_SLC7A11.RData")
rm(all.tpm.final)
gc()

#-------------------------------------------------------------------------------
# 第三个段落： 数据统计和可视化
#              输入包括：all.colinfo,oxidative stress(17种可能),统计量，表达值
#-------------------------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))#画图辅助
load(file = "expr_SLC7A11.RData") #geneexpr
load(file = "met_cys_stat.RData") #all_result
load(file = "the_17_possible_oxidative_stress_level_of_all_samples.RData")#contrast_17 
contrast_17 = contrast_17[,2,drop = FALSE]# 使用ER1模块作为oxidative
load("tcga33_gtex18_samples11348v5620_and_genes56189_all.colinfo_only.RData")#all.colinfo


major_types = c("normal","Solid Tissue Normal","Primary Tumor","Metastatic")
cora  = cor(contrast_17[11349:16968,],all_result[11349:16968,],method = "spearman")
for(cancer_type in major_types)
{
  corb = cor(contrast_17[all.colinfo$Sample_type == cancer_type,],all_result[all.colinfo$Sample_type == cancer_type,],method = "spearman")
  cora = rbind(cora,corb)
}
cora = cora[-1,]
cora = data.frame(cora)
type = rep(major_types,times = 3)
molecule = rep(c('met','cys','sum'),each = 4)
value = c(cora$met,cora$cys,cora$sum)
cora = data.frame(type = factor(type,levels =major_types ), molecule = factor(molecule), value = value)
  # 画折线图
ggplot(data = cora, mapping = aes(x = type, y = value,group = molecule, colour = molecule)) + 
  geom_line()+ 
  geom_point()+
  geom_text(aes(label= round(value,2)),position=position_dodge(width=0.1),hjust=1.5)
ggsave("sample_type_CCC.png")

result_to_plot = data.frame(all_result)
result_to_plot$type = all.colinfo$Sample_type
result_to_plot = result_to_plot[all.colinfo$Sample_type %in% major_types,]
result_to_plot$type = factor(result_to_plot$type,levels =major_types)
Varns = names(result_to_plot)[1:ncol(all_result)]
i = 1
for(varn in Varns)
{
  p <- ggboxplot(result_to_plot, x = "type", y = varn,
                  color = "type", 
                  palette = "jco",
  )+ 
    ggtitle('sample type') +
    stat_compare_means(method = "anova")
  assign(paste0("p",i,sep = ''),p)
  i = i + 1
}
ggarrange(p1,p2,p3,ncol =3)
png(filename = paste('D:/R/Oxidative_Stress/',"met统计量随类型分布",'.png',sep = ''),width = 3000,height = 1800,res = 200)
ggarrange(p1,p2,p3,ncol =3)
dev.off()