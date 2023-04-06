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

# 系统发育信号计算
# 数据转换为0/1
data_01 <- data %>%
  # 选择对应的数据列
  select(Fruittype, Growthform, Tropic, Subtropic, Temperate) %>%
  # 将其转换为0/1格式
  mutate(Fruittype = if_else(Fruittype == "fleshy",
                             true = 1,
                             false = 0),
         Growthform = if_else(Growthform == "herb",
                              true = 1,
                              false = 0),
         Tropic = if_else(Tropic == "Tropic",
                          true = 1,
                          false = 0),
         Temperate = if_else(Temperate == "Temperate",
                             true = 1,
                             false = 0),
         Subtropic = if_else(Subtropic == "Subtropic",
                             true = 1,
                             false = 0))

## 对树进行标准化
st_tree <- tree
st_tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree))
## 我们的数据是tibble格式的，没有rowname
## phyloglm必需依赖rowname，因此首先要添加rowname，并命名检查
row.names(data_01) <- data$species_name
name.check(st_tree,data_01) # 命名检查
## 果实类型结果测试
phyloglm(Fruittype ~ 1, data=data_01, phy=st_tree)
R2_lik(phyloglm(Fruittype ~ 1, data=data_01, phy=tree))
## 自定义函数对5个性状进行计算
phy_sign <- function () {
  z <- list(NULL)
  z[[1]] <- phyloglm(Fruittype ~ 1, data=data_01, phy=st_tree)
  z[[2]] <- phyloglm(Growthform ~ 1, data=data_01, phy=st_tree)
  z[[3]] <- phyloglm(Tropic ~ 1, data=data_01, phy=st_tree)
  z[[4]] <- phyloglm(Subtropic ~ 1, data=data_01, phy=st_tree)
  z[[5]] <- phyloglm(Temperate ~ 1, data=data_01, phy=st_tree)
  df <- data.frame(traits=c("fruit", "growthform","Tropic", "Subtropic", "Temperate"))
  for (i in 1:5) {
    df$alpha[i] <- z[[i]]$alpha
    df$r2[i] <- R2_lik(z[[i]])
  }
  return(df)
}

phy_sign()
#       traits      alpha        r2
# 1      fruit  0.7786089 0.8634150
# 2 growthform  0.7467027 0.8681852
# 3     Tropic  7.6515278 0.2004165
# 4  Subtropic 54.5207242 0.0000000
# 5  Temperate 11.0576667 0.1218708