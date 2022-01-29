library(vegan)
asvmat <- read.csv("ASV_table_PCoA.csv", header=TRUE, quote="", row.names=1)
asvmat <- data.frame(t(asvmat))
sampledata <- read.csv("sample_data_PCoA.csv", header=TRUE, quote="", row.names=1)

#distance calculation
bray_dist<-vegdist(asvmat,method = "bray")

#PCoA analysis
library(ape)
asvmat.pcoa<-pcoa(bray_dist,correction = "none")

#data for visualization
asvmat.pcoa$vectors

#Axis values
asvmat.pcoa$values
asv.plot<-data.frame(asvmat.pcoa$vectors)
head(asv.plot)

#value matching
asv.plot = cbind(asv.plot, sampledata[,1])
asv.plot[,35]
colnames(asv.plot)[35] = 'group'


#vasualizing using ggplot2
library(ggplot2)
library(ggsci)
x_label<-round(asvmat.pcoa$values$Relative_eig[1]*100,2)
y_label<-round(asvmat.pcoa$values$Relative_eig[2]*100,2)
x_label
y_label

fig <- ggplot(asv.plot, aes(x=Axis.1,y=Axis.2,fill=group))+ 
  # 选择X轴Y轴并映射颜色和形状
  geom_point(shape=21,aes(fill = group,alpha = 0.5),size=2)+ # 画散点图并设置大小
  geom_hline(yintercept = 0,linetype="dashed",size = 0.3) + # 添加横线
  geom_vline(xintercept = 0,linetype="dashed",size = 0.3) + # 添加竖线
  scale_color_igv()+ # 设置颜色,此处为Integrative Genomics Viewer配色
  theme(panel.border = element_rect(fill="transparent",color="black", size=1, linetype="solid"),
        panel.background = element_rect(fill="transparent")) + # 加上边框,主题设置
  stat_ellipse(type = "t",level = 0.95,linetype = 2, size = 0.4)+ # 添加置信椭圆
  # 自动提取主成分解释度进行绘图
  labs(x=paste0("PCoA1 ",x_label,"%"),  
       y=paste0("PCoA2 ",y_label,"%"))+
  scale_fill_manual(values=c("#00c4ff", "#ff9289", "#e199ff", "#96ca00"))+
  theme(legend.position = "none") # 设置图例位置，此处为相对位置
ggsave("PCoA_lake.tiff", width = 8, height = 8, units = "cm", dpi = 300)
