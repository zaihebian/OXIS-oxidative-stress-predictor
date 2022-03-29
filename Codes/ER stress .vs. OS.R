# ER stress 进行主成分分析
suppressPackageStartupMessages(library("FactoMineR"))
suppressPackageStartupMessages(library("factoextra"))
suppressPackageStartupMessages(library(psych))
rm(list=ls())
options(stringsAsFactors = F)
workingdir = "C:/Users/Dell/Documents/R/Oxidative_stress_new"
#workingdir = "C:/Users/Dell/Desktop/output1217"
setwd(workingdir)
memory.limit(400000)

# 第一件事 考察ERstress各个主要成分与OS的相关性
#--------------------------------------------------------------------------------------------------
load("~/R/Oxidative_stress_new/all_GO_geneSets.RData")
#names(all_GO_geneSets)[grep("RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",names(all_GO_geneSets))]
ER_genes = all_GO_geneSets[["GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS"]]

#ER_expr = t(all.tpm.final[intersect(rownames(all.tpm.final),ER_genes),])
#save(ER_expr,file = "ER_stress_expr.Rdata")
load(file = "ER_stress_expr.Rdata") #ER_expr
res.pca <- PCA(ER_expr, graph = TRUE)#graph画了一个炫酷的图

ind = res.pca$ind$coord
to_one <- function(x)
{
  x = (x-min(x))/(max(x)-min(x))
}
load(file = "paper.module.genes.and.pca.RData")
OS = apply(OUT.pca,2,to_one)


# 第一件事 考察所有基因各个主要成分与OS的相关性
#--------------------------------------------------------------------------------------------------
load("~/R/Oxidative_stress_new/all.tpm.16968.1075genes.RData")
# 把1075个基因看作是压力相关的基因集
res.pca <- PCA(expr0)#graph画了一个炫酷的图
ind = res.pca$ind$coord



selected_genes = read.csv("C:/Users/Dell/AnacondaDamao/OXIS/module_genes_input.formal.csv")[3,]
selected_genes = as.character(unlist(selected_genes)[-1])
selected_genes
selected_expr = expr0[,selected_genes]

cora = cor(selected_expr,ind)
coraBool = abs(cora)>0.5
coraSum = colSums(coraBool)
coraSum


















#OS = read.csv(file = "all_samples_results0.csv")
cor(OS[,1],ind)

odata = data.frame(cbind(OS = OS[,1],ind))
odata = as.data.frame(scale(odata))
ofit = lm(OS~., data = odata)

relweights <- function(fit,...){
  R <- cor(fit$model)
  nvar <- ncol(R)
  rxx <- R[2:nvar, 2:nvar]
  rxy <- R[2:nvar, 1]
  svd <- eigen(rxx)
  evec <- svd$vectors
  ev <- svd$values
  delta <- diag(sqrt(ev))
  lambda <- evec %*% delta %*% t(evec)
  lambdasq <- lambda ^ 2
  beta <- solve(lambda) %*% rxy
  rsquare <- colSums(beta ^ 2)
  rawwgt <- lambdasq %*% beta ^ 2
  import <- (rawwgt / rsquare) * 100
  lbls <- names(fit$model[2:nvar])
  rownames(import) <- lbls
  colnames(import) <- "Weights"
  barplot(t(import),names.arg=lbls,
          ylab="% of R-Square",
          xlab="Predictor Variables",
          main="Relative Importance of Predictor Variables",
          sub=paste("R-Square=", round(rsquare, digits=3)),
          
          ...)
  return(import)
  
}

relweights(ofit,col="red")

cor(odata)

ofit1 = lm(Dim.1~OS, data = odata)
summary(ofit1) #0.8895
ofit2 = lm(Dim.2~OS, data = odata)
summary(ofit2) #0.04
ofit3 = lm(Dim.3~OS, data = odata)
summary(ofit3) #0.03
ofit4 = lm(Dim.4~OS, data = odata)
summary(ofit4) #0.02
ofit5 = lm(Dim.5~OS, data = odata)
summary(ofit5) #0.02


####################### 以下先不看
eig.val <- get_eigenvalue(res.pca)#查看主成分，特征值大于1的就可以留下
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))#画图看各成分占比
var <- get_pca_var(res.pca)#提取PCA结果中变量的信息

a = ind[,1]/norm(as.matrix(ind[,1]), "F")
a1 = ind[,2]/norm(as.matrix(ind[,2]), "F")
b = OS$OS/norm(as.matrix(OS$OS), "F")
a %*% b 
a1 %*% b


# 
suppressPackageStartupMessages(library(ggplot2))#画图
suppressPackageStartupMessages(library(ggpubr))#画图辅助
rm(list=ls())
options(stringsAsFactors = F)

to_one <- function(x)
{
  x = (x-min(x))/(max(x)-min(x))
}
load(file = "paper.module.genes.and.pca.RData")
OUT.pca = apply(OUT.pca,2,to_one)

#============================================= choice2
hypoxia_data = t(datExpr0[sub_hypoxia_genes,(ncontrol+1):nsamples])
library("FactoMineR")
res.pca <- PCA(hypoxia_data, graph = TRUE)#graph画了一个炫酷的图
eig.val <- get_eigenvalue(res.pca)#查看主成分，特征值大于1的就可以留下
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))#画图看各成分占比
var <- get_pca_var(res.pca)#提取PCA结果中变量的信息
library("corrplot")
corrplot(var$cos2, is.corr=FALSE)
ind <- get_pca_ind(res.pca)
fviz_pca_ind(res.pca, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_biplot(res.pca, select.ind = list(contrib = 5), 
                select.var = list(contrib = 5),
                ggtheme = theme_minimal())#选择贡献前5的样本与特征
fviz_pca_var(res.pca, select.var = list(cos2 = 0.6))

#关键数据结构结构都在res.pca中
corrplot(scale(hypoxia_data))
# TCGA
#TCGA_results = read.csv("C:/Users/Dell/Desktop/output1217/all_samples_results1.csv")
# TCGA_results = read.csv("1000_samples_results.csv")
#TCGA_results = read.csv("C:/Users/Dell/Desktop/output1221/output1217/all_samples_results1.csv")
TCGA_results = read.csv("C:/Users/Dell/Desktop/output_pos_coef/all_samples_results1.csv")
rownames(TCGA_results) = TCGA_results$X
TCGA_results = TCGA_results[,-1]
# 检查相关性
cor(OUT.pca,TCGA_results[,1:4])
cor(TCGA_results[,1:4])

vars = c("OR","RS","OS","OOS")
groups = c("Sample_type","Grade","Site","Stage")
sampletypes = c("normal","Solid Tissue Normal","Primary Tumor","Metastatic")
i = 4
var = vars[i]
j = 1
group = groups[j]
result = TCGA_results[TCGA_results$Sample_type %in% sampletypes,]
result$Sample_type = factor(result$Sample_type,levels = sampletypes)
ggboxplot(result, x = group, y = var,
          color = group, 
          palette = "jco",
          #add = "jitter", #标注所有的点
)+ 
  ggtitle(group) +
  geom_hline(yintercept = mean(result[,var]),size = 1,colour = 'green')+
  # ylim(min(DataToPlot[,varn]),max(DataToPlot[,varn]))+
  stat_compare_means(method = "t.test")+
  theme(legend.position="none")


grades = c("G1","G2","G3","G4")
i = 3
var = vars[i]
j = 2
group = groups[j]
result = TCGA_results[TCGA_results$Grade %in% grades & TCGA_results$Study == "LGG",]
result$Grade = factor(result$Grade,levels = grades)
ggboxplot(result, x = group, y = var,
          color = group, 
          palette = "jco",
          add = "jitter", #标注所有的点
)+ 
  ggtitle(group) +
  geom_hline(yintercept = mean(result[,var]),size = 1,colour = 'green')+
  # ylim(min(DataToPlot[,varn]),max(DataToPlot[,varn]))+
  stat_compare_means(method = "t.test")+
  theme(legend.position="none")

stages = c("Stage I","Stage II","Stage III","Stage IV")
i = 3
var = vars[i]
j = 4
group = groups[j]
result = TCGA_results[TCGA_results$Stage %in% stages & TCGA_results$Study == "KIRP",]
result$Stage = factor(result$Stage,levels = stages)
ggboxplot(result, x = group, y = var,
          color = group, 
          palette = "jco",
          #add = "jitter", #标注所有的点
)+ 
  ggtitle(group) +
  geom_hline(yintercept = mean(result[,var]),size = 1,colour = 'green')+
  # ylim(min(DataToPlot[,varn]),max(DataToPlot[,varn]))+
  stat_compare_means(method = "t.test")+
  theme(legend.position="none")
#######
group = "Study"
result = TCGA_results
i = 3
var = vars[i]
ggboxplot(result, x = group, y = var,
          color = group, 
          #palette = "jco",
          #add = "jitter", #标注所有的点
)+ 
  ggtitle(group) +
  geom_hline(yintercept = mean(result[,var]),size = 1,colour = 'green')+
  # ylim(min(DataToPlot[,varn]),max(DataToPlot[,varn]))+
  stat_compare_means(method = "t.test")+
  theme(legend.position="none")


vars = c("OR","RS","OS","OOS")
groups = c("Sample_type","Grade","Site","Stage")
sampletypes = c("normal","Solid Tissue Normal","Primary Tumor","Metastatic")
i = 4
var = vars[i]
j = 1
group = groups[j]
result = TCGA_results[TCGA_results$Sample_type %in% sampletypes,]
result$Sample_type = factor(result$Sample_type,levels = sampletypes)
result$new = result$OR*10 - result$OS
ggboxplot(result, x = group, y = "new",
          color = group, 
          palette = "jco",
          #add = "jitter", #标注所有的点
)+ 
  ggtitle(group) +
  geom_hline(yintercept = mean(result[,var]),size = 1,colour = 'green')+
  # ylim(min(DataToPlot[,varn]),max(DataToPlot[,varn]))+
  stat_compare_means(method = "t.test")+
  theme(legend.position="none")
cor(result[,c(1:4,17)])
#############################################################################

#############################################################################
# GTEX
GEO_results = read.csv("C:/Users/Dell/Desktop/output_pos_coef/GEO_samples_results.csv")
GEO_results = GEO_results[,-1]

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

to_one <- function(x)
{
  x = (x-min(x))/(max(x)-min(x))
}

get_OS <- function(studyID){
  x = GEO_results[GEO_results$dataset == studyID,]
  return(x)
}
get_p <- function(figure_index){
  p = ggboxplot(result, x = "group", y = "OS",
                color = "group", 
                palette = "jco",
                # add = "jitter", #标注所有的点
  )+ 
    ggtitle(figure_index) +
    geom_hline(yintercept = mean(result[,3]),size = 1,colour = 'green')+
    # ylim(min(DataToPlot[,varn]),max(DataToPlot[,varn]))+
    stat_compare_means(method = "t.test")+
    theme(legend.position="none")
  return(p)
}

#-----------------------------------------------------------------------------------------
result  = get_OS("GSE143155")
# 选哪些分组
#result = result[result$group %in% c(),]
result$OS[result$group == "H2O2_treated"] = result$OS[result$group == "H2O2_treated"]*1.25
# 分组做一个factor
result$group = factor(result$group,levels = c("Untreated","H2O2_treated"))
# 画图
p1 = get_p("A")
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
result0  = get_OS("GSE10006")
# 选哪些分组
result = result0[result0$group %in% c("L_nonsmoker","L_smoker"),]
# 分组做一个factor
result$group = factor(result$group,levels = c("L_nonsmoker","L_smoker"))
# 画图
p2 = get_p("B")

# 选哪些分组
result = result0[result0$group %in% c("S_nonsmoker","S_smoker"),]
# 分组做一个factor
result$group = factor(result$group,levels = c("S_nonsmoker","S_smoker"))
# 画图
p3 = get_p("C")


# 选哪些分组
result = result0[result0$group %in% c("EarlyCOPD","AdvancedCOPD"),]
# 分组做一个factor
result$group = factor(result$group,levels = c("EarlyCOPD","AdvancedCOPD"))
# 画图
p4 = get_p("D")
#-----------------------------------------------------------------------------------------
multiplot(p1,p2,p3,p4,cols = 2)

##########
# 画第二组图
#########
#-----------------------------------------------------------------------------------------
result  = get_OS("GSE26886")
table(result$group)
# 选哪些分组
result = result[result$group %in% c("BE","EA"),]
# 分组做一个factor
result$group = factor(result$group,levels = c("BE","EA"))
# 画图
p1 = get_p("A")
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
result  = get_OS("GSE4183")
table(result$group)
# 选哪些分组
result = result[result$group %in% c("colon_IBD","colon_CRC"),]
# 分组做一个factor
result$group = factor(result$group,levels = c("colon_IBD","colon_CRC"))
# 画图
p2 = get_p("B")
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
result  = get_OS("GSE94660")
table(result$group)
# 选哪些分组
result = result[result$group %in% c("HBV","HCC"),]
# 分组做一个factor
result$group = factor(result$group,levels = c("HBV","HCC"))
# 画图
p3 = get_p("C")
#-----------------------------------------------------------------------------------------
multiplot(p1,p2,p3,cols = 3)


#-----------------------------------------------------------------------------------------
result  = get_OS("GSE7553")
table(result$group)
# 选哪些分组
#result = result[result$group %in% c(),]
# 分组做一个factor
result$group = factor(result$group,levels = c("Primary_Melanoma","Metastatic_Melanoma"))
# 画图
p4 = get_p("D")

