library(phylolm) # 注意，这一行需要先在ggtree前载入
library(geiger)
library(rr2)


library(tidyverse)
library(V.PhyloMaker2)

library(ggtree)
library(tidytree)
library(ggtreeExtra)
library(ggnewscale)

# 读取数据
merged_data <- read_csv("merged_data.csv")

treesets <- merged_data %>%
  select(Species, Genus, Family,`species.relative`) %>%
  phylo.maker(tree = GBOTB.extended.TPL, 
              nodes = nodes.info.1.TPL, 
              scenarios = c("S1","S2","S3"))

tree <- treesets$scenario.3 # 根据作者的推荐使用S3情形

# 保存系统树到本地
write.tree(tree,"result/682.tre")

# 读取树
tree <- read.tree("result/682.tre")

# 绘制系统树
## 创建颜色分类
my_colors <- c("0" = "cornflowerblue", "1" = "coral",
               "T" = "#82954B", "S" = "#EFD345",
               "short" = "#C8B6A6", "middle" = "#A4907C", "high" = "#8D7B68",
               "D" = "orange", "E" = "lightgreen",
               "TRR" = "#F5E8C7", "TRS" = "#ECCCB2", "TMR" = "#DEB6AB", 
               "TMS" = "#AC7088", "WDS" = "#DB13C6","DES" = "#A31ACB", 
               "BOR" = "#39B5E0")

## 添加分组信息（也就是绘制时线条颜色对应的分类）
newtree <- groupOTU(tree, split(merged_data$Species, merged_data$eff_safe))
newtree %>% as_tibble()

# 线条色为gorup即eff_safe
(p <- ggtree(newtree,layout = "fan",aes(color=group)) +
    # 添加Life.form映射
    geom_fruit(data = merged_data,geom = geom_point,
               pch = 15,size=0.2,
               mapping = aes(y = Species,color=Life.form)) +
    # 添加Height_category映射
    geom_fruit(data = merged_data,geom = geom_point,
               pch = 15,size=0.21,
               mapping = aes(y = Species,color=Height_category)) +
    # 添加Leaf.form映射
    geom_fruit(data = merged_data,geom = geom_point,
               pch = 15,size=0.21,
               mapping = aes(y = Species,color=Leaf.form)) +
    # 添加Biome映射
    geom_fruit(data = merged_data,geom = geom_point,
               pch = 15,size=0.21,
               mapping = aes(y = Species,color=Biome)) +
    scale_color_manual(values = my_colors)
)
# save plot
ggsave("result/phylogeny_plot.png",plot=p,dpi=300)
