#install.packages("reshape2")
library(reshape2)
#输入数据
asv <- read.csv("C:/Users/mumua/Desktop/16s paper/network/asv_network.csv", header=TRUE, quote="", row.names=1)
tax <- read.csv("C:/Users/mumua/Desktop/16s paper/network/tax.csv", header=TRUE, quote="", row.names=1)

##转换为边列表，即 asv 和样本的对应关系
#结果中，所有非 0 的值代表该 asv 在该样本中存在丰度
edge <- asv
edge$ASV <- rownames(edge)
edge <- reshape2::melt(edge, id = 'ASV')

#删除不存在的边
#即去除 0 丰度的值，该 OTU 不在该样本中存在
edge <- subset(edge, value != 0)

#修改列名称（以便后续软件识别），其中权重可表示为 OTU 在样本中的丰度
names(edge) <- c('source', 'target', 'weight')

#添加一列“shared name”，以“->”连接 source 和 target
#便于 cytoscape 读取识别，并防止读取后的名称错乱
edge$'shared name' <- paste(edge$source, edge$target, sep = '->')

#输出边列表，后续可导入至 cytoscape 用于构建网络
write.table(edge, 'edge.txt', sep = '\t', quote = FALSE, row.names = FALSE)


##获取节点属性列表
#获得 asv 在各组中的分布状态（Venn 分布）
asv[asv>0] <- 1
asv$venn <- apply(asv, 1, function(x) paste(x, collapse = '-'))
asv <- asv[unique(edge$source), ]

#asv 属性列表
node_asv <- data.frame(
  'shared name' = rownames(asv),  #asv 名称，以“shared name”命名列名称，便于 cytoscape 读取识别，并防止读取后的名称错乱
  group1 = tax[unique(edge$source),'Phylum'],  #用于后续按指定分组赋值颜色，这里定义为 OTU 所属的门分类水平
  group2 = 'asv',  #用于后续按指定分组定义形状，这里统一分组为“asv”
  group3 = asv$venn,  #用于后续按指定分组调整聚群，按 asv 在各组中的分布状态定义
  stringsAsFactors = FALSE, 
  check.names = FALSE
)

#样本属性列表
asv <- asv[-ncol(asv)]

node_sample <- data.frame(
  'shared name' = names(asv),  #样本名称，以“shared name”命名列名称，便于 cytoscape 读取识别，并防止读取后的名称错乱
  group1 = names(asv),  #用于后续按指定分组赋值颜色，指定样本分组
  group2 = names(asv),  #用于后续按指定分组定义形状，指定样本分组
  group3 = names(asv),  #用于后续按指定分组调整聚群，指定样本分组
  stringsAsFactors = FALSE, 
  check.names = FALSE
)

#二者合并构建节点属性列表后输出，后续可导入至 cytoscape 用于调整节点可视化
node <- rbind(node_asv, node_sample)
write.table(node, 'node.txt', sep = '\t', quote = FALSE, row.names = FALSE)
