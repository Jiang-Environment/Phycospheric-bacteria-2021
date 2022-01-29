library(vegan)

##读入文件
#OTU 丰度表
otu <- read.csv("C:/Users/mumua/Desktop/16s paper/andos/ASV_table_no0.csv", header=TRUE, quote="", sep=",",row.names=1,stringsAsFactors = FALSE)
otu <- data.frame(t(otu))

#样本分组文件
group <- read.csv("C:/Users/mumua/Desktop/16s paper/andos/group.csv", header=TRUE, quote="", sep=",",row.names=1,stringsAsFactors = FALSE)

#（1）直接输入 OTU 丰度表，在参数中指定距离类型
#使用 Bray-Curtis 距离测度
adonis_result <- adonis(otu~lake, group, distance = 'bray', permutations = 1000)

#查看结果
adonis_result



#################################
#######两两配对对比

#install.packages('devtools')
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(cluster)
library(pairwiseAdonis)

# This is a wrapper function for multilevel pairwise comparison 
# using adonis() from package 'vegan'. 
# The function returns adjusted p-values using p.adjust().
dune.pairwise.adonis <- pairwise.adonis(x=otu, factors=group$lake, sim.function = "vegdist",
                                        sim.method = "bray",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 1000)

dune.pairwise.adonis


############correlation heatmap##########
##########################################
###########################################
library(pheatmap)
library("RColorBrewer")
# 读取比较列表
asvmat <- read.csv("C:/Users/mumua/Desktop/16s paper/andos/bray_nofree.csv", header=TRUE, quote="", row.names=1)

# 读取实验设计注释列分组
sampledata <- read.csv("C:/Users/mumua/Desktop/16s paper/andos/group_nofree.csv", header=TRUE, quote="", row.names=1)

# 准备行/列注释
anno_col = data.frame(Group = sampledata$lake, col.names = colnames(sampledata))
anno_row = data.frame(Group = sampledata$lake, row.names = rownames(sampledata))
# 绘图
pheatmap(asvmat,
         treeheight_col=32,treeheight_row=32,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdGy")[2:10]))(100),
         annotation_names_row= T,annotation_names_col=T,border_color = NA,
         annotation_col = anno_row,annotation_row = anno_row,
         filename = paste("p9.Bray-Curtis.tiff", sep=""),width=7, height=6,fontsize=10,display_numbers=F)

