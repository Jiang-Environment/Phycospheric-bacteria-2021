library(vegan)

##读入文件
#OTU 丰度表
otu <- read.csv("C:/Users/mumua/Desktop/16s paper/PCoA/ASV_table_PCoA_BW.csv", header=TRUE, quote="", sep=",",row.names=1,stringsAsFactors = FALSE)
otu <- data.frame(t(otu))

#样本分组文件
group <- read.csv("C:/Users/mumua/Desktop/16s paper/PCoA/sample_data_PCoA _BW.csv", header=TRUE, quote="", sep=",",row.names=1,stringsAsFactors = FALSE)

#（1）直接输入 OTU 丰度表，在参数中指定距离类型
#使用 Bray-Curtis 距离测度
adonis_result <- adonis(otu~Type, group, distance = 'bray', permutations = 1000)

#查看结果
adonis_result
