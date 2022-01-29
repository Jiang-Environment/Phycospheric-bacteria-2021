library(vegan)

##读入文件
#OTU 丰度表
otu <- read.csv("ASV_table.csv", header=TRUE, quote="", sep=",",row.names=1,stringsAsFactors = FALSE)
otu <- data.frame(t(otu))

#样本分组文件
group <- read.csv("sample_data_PCoA.csv", header=TRUE, quote="", sep=",",row.names=1,stringsAsFactors = FALSE)

#（1）直接输入 ASV 丰度表，在参数中指定距离类型
#使用 Bray-Curtis 距离测度
adonis_result <- adonis(otu~Type, group, distance = 'bray', permutations = 1000)

#查看结果
adonis_result
