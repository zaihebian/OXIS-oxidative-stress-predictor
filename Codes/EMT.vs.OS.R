#  计算EMT与OS的相关性
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library("FactoMineR"))
suppressPackageStartupMessages(library(psych))
rm(list=ls())
options(stringsAsFactors = F)
workingdir = "C:/Users/Dell/Documents/R/Oxidative_stress_new"
setwd(workingdir)
memory.limit(400000)
# 加载EMT的主成分
load("~/R/Oxidative_stress_new/paper.module.genes.and.pca.3geneSets.RData")
all_to_plot = read.csv("C:/Users/Dell/AnacondaDamao/OXIS/all_samples_results.csv") #读入最开始的数据
load("~/R/Oxidative_stress_new/tcga33_gtex18_samples11348v5620_and_genes56189_all.colinfo_only.RData")

rownames(all_to_plot) = all_to_plot$X
all_to_plot = all_to_plot[,-1]

OUT.pca1 = OUT.pca1[all.colinfo$Sample_type %in% c("Solid Tissue Normal","normal","Primary Tumor","Metastatic"),]
cora  = cor(OUT.pca1,all_to_plot$OS)

library(Hmisc)
matrix = cbind(OUT.pca1[,1],all_to_plot$OS)
rcorr(matrix)

