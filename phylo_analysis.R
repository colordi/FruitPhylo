library(tidyverse)
library(V.PhyloMaker2)
library(ggtree)
library(tidytree)
library(ggtreeExtra)
library(ggnewscale)

# 读取数据
data <- read_csv("data/YunnanFruit9370sppData.csv")

# 构建树
treesets <- data %>%
  mutate(Species = species_name) %>%
  select(Species, Genus, Family) %>%
  phylo.maker(tree = GBOTB.extended.TPL, 
              nodes = nodes.info.1.TPL, 
              scenarios = c("S1","S2","S3"))
tree <- treesets$scenario.1
# 保存系统树到本地
write.tree(tree,"result/YunnanFruit9370sppTree.tre")

# 读取树(原始附件中的树文件)
tree <- read.tree("data/phyloa_sp_YNfruit9370tree.tre")

# 简单绘制树的图形(非必需)
plot(tree, cex=.01, show.tip.label = F)
axisPhylo()

# 绘制系统树
## 创建颜色分类
my_colors <- c("dry" = "cornflowerblue", "fleshy" = "coral",
               "wood" = "blue", "herb" = "yellow",
               "Absence" = "white", "Tropic" = "red", 
               "Subtropic" = "orange", "Temperate" = "lightgreen")

## 添加分组信息（果实类型）
newtree <- groupOTU(tree, split(data$species_name, data$Fruittype))

(p <- ggtree(newtree,layout = "fan",aes(color=group)) +
  # 添加生活型映射
  geom_fruit(data = data,geom = geom_point,
             pch = 15,size=0.2,
             mapping = aes(y = species_name,color=Growthform)) +
  # 添加热带映射
  geom_fruit(data = data,geom = geom_point,
             pch = 15,size=0.21,
             mapping = aes(y = species_name,color=Tropic)) +
  # 添加亚热带映射
  geom_fruit(data = data,geom = geom_point,
             pch = 15,size=0.21,
             mapping = aes(y = species_name,color=Subtropic)) +
  # 添加温带映射
  geom_fruit(data = data,geom = geom_point,
             pch = 15,size=0.21,
             mapping = aes(y = species_name,color=Temperate)) +
  scale_color_manual(values = my_colors)
)
# save plot
ggsave("result/phylogeny_plot.png",plot=p,dpi=300)


