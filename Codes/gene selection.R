# get hub genes
suppressPackageStartupMessages(library(WGCNA))#共表达分析
rm(list=ls())
options(stringsAsFactors = F)
workingdir = "C:/Users/Dell/Documents/R/Oxidative_stress_new"
setwd(workingdir)
memory.limit(400000)

load("~/R/Oxidative_stress_new/all_GO_geneSets.RData")
F_genes = all_GO_geneSets[[grep("GOBP_RNA_CATABOLIC_PROCESS",names(all_GO_geneSets))]]
F_genes = all_GO_geneSets[[grep("GOBP_DNA_REPAIR",names(all_GO_geneSets))]]
load(file = "tcga33_gtex18_samples11348v5620_and_genes56189.RData")#all.rowinfo,all.colinfo,all.tpm.final
F_genes = F_genes[F_genes %in% rownames(all.tpm.final)]
all.tpm.final = all.tpm.final[F_genes,]
gc()

all_relavent_data = t(all.tpm.final)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(all_relavent_data, powerVector = powers, verbose = 5)
sft$powerEstimate
if(is.na(sft$powerEstimate))
{selected_power =5}else{selected_power = sft$powerEstimate}
cor <- WGCNA::cor
net = blockwiseModules(all_relavent_data, power = selected_power,
                       TOMType = "unsigned", minModuleSize = 10,#从30开始逐步递减
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)
cor<-stats::cor
MEs = as.matrix(net$MEs)
MEs = MEs[,order(colnames(MEs)),drop = FALSE]
color = as.data.frame(net$colors)
selected_genes = c()
for(j in c(0:(ncol(MEs)-1)))
{
  genes = rownames(color)[color$`net$colors`==j]
  cormodule = cor(MEs[,j],all_relavent_data[,genes])
  selected_genes = append(selected_genes,colnames(cormodule)[which.max(cormodule)])
}
selected_genes