suppressPackageStartupMessages(library(ggplot2))#画图
suppressPackageStartupMessages(library(ggpubr))#画图辅助
rm(list=ls())
options(stringsAsFactors = F)
workingdir = "C:/Users/Dell/Documents/R/Oxidative_stress_new"
setwd(workingdir)
memory.limit(400000)

to_one <- function(x)
{
  x = (x-min(x))/(max(x)-min(x))
}


# GEO
GEO_results = read.csv("C:/Users/Dell/AnacondaDamao/OXIS/results.csv")
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
    layout <- t(matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols)))
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
get_p <- function(figure_index,xlabels){
  p = ggboxplot(result, x = "group", y = "OS",
                color = "group", 
                palette = "jco",
                # add = "jitter", #标注所有的点
                size = 1.5    #设置线条粗细
  )+ 
    ggtitle(figure_index) +
    #geom_hline(yintercept = mean(result[,2]),size = 1,colour = 'green')+
    stat_compare_means(aes(label = paste0("p-value = ",..p.format..)),method = "t.test", size = 8,label.y = 1.8)+ 
    xlab(NULL)+     #不显示x label
    ylab("Oxidative Stress Level")+
    ylim(c(0,2))+
    theme_light(base_size = 18)+
    theme(legend.position="none") +
    theme(legend.title=element_blank())+
    theme(axis.text.x = element_text(angle = 0,color="black", size=18),axis.text.y = element_text(color="black", size=20))+
    scale_x_discrete(labels = xlabels)
  return(p)
}

#-----------------------------------------------------------------------------------------
result  = get_OS("GSE143155")
# 选哪些分组
#result = result[result$group %in% c(),]
result$OS[result$group == "H2O2_treated"] = result$OS[result$group == "H2O2_treated"]
# 分组做一个factor
result$group = factor(result$group,levels = c("Untreated","H2O2_treated"))
# 画图
p1 = get_p("A",c("Untreated",expression("H"[2]*"O"[2]~"t"*"r"*"e"*"a"*"t"*"e"*"d")))
#-----------------------------------------------------------------------------------------
p1