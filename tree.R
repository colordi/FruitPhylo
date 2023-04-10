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
               "short" = "#C8B6A6", "middle" = "#A4907C", "tall" = "#8D7B68",
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

# 系统发育信号计算
# 数据转换为0/1
data_01 <- merged_data %>%
  # 选择对应的数据列
  select(eff_safe,Residualall,Height_category, Biomeclass, Life.form, Leaf.form) %>% 
  transmute(eff_safe,
            Residualall,
            Short = if_else(Height_category == "short",
                            true = 1,
                            false = 0,),
            Middle = if_else(Height_category == "middle",
                             true = 1,
                             false = 0,),
            Tall = if_else(Height_category == "tall",
                           true = 1,
                           false = 0,),
            Biomeclass = case_when(Biomeclass == "seasonal" ~ 1,
                                   Biomeclass == "non_seasonal" ~ 0,
                                   is.na(Biomeclass) ~ NA),
            LifeForm = case_when(Life.form == "S" ~ 1,
                                 Life.form == "T" ~ 0,
                                 is.na(Life.form) ~ NA),
            LeafForm = case_when(Leaf.form == "D" ~ 1,
                                 Leaf.form == "E" ~ 0,
                                 is.na(Leaf.form) ~ NA)
            )

## 对树进行标准化
st_tree <- tree
st_tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree))
## 我们的数据是tibble格式的，没有rowname
## phyloglm必需依赖rowname，因此首先要添加rowname，并命名检查
row.names(data_01) <- merged_data$Species
name.check(st_tree,data_01) # 命名检查
## 测试
phyloglm(eff_safe ~ 1, data=data_01, phy=st_tree)
R2_lik(phyloglm(eff_safe ~ 1, data=data_01, phy=tree))
## 测试一般线性模型
phylolm(Residualall ~ 1, data=data_01, phy=st_tree,model="lambda")
phylolm(Residualall ~ 1, data=data_01, phy=st_tree,model="lambda") %>% R2_lik()

## 自定义函数对5个性状进行计算
phy_sign <- function () {
  z <- list(NULL)
  z[[1]] <- phyloglm(eff_safe ~ 1, data=data_01, phy=st_tree)
  z[[2]] <- phyloglm(Short ~ 1, data=data_01, phy=st_tree)
  z[[3]] <- phyloglm(Middle ~ 1, data=data_01, phy=st_tree)
  z[[4]] <- phyloglm(Tall ~ 1, data=data_01, phy=st_tree)
  z[[5]] <- phyloglm(Biomeclass ~ 1, data=data_01, phy=st_tree)
  z[[6]] <- phyloglm(LifeForm ~ 1, data=data_01, phy=st_tree)
  z[[7]] <- phyloglm(LeafForm ~ 1, data=data_01, phy=st_tree)
  df <- data.frame(traits=c("eff_safe", "Short","Middle", "Tall", "Biomeclass",
                            "LifeForm","LeafForm"))
  for (i in 1:7) {
    df$alpha[i] <- z[[i]]$alpha
    df$r2[i] <- R2_lik(z[[i]])
  }
  return(df)
}

phy_sign()
#       traits     alpha         r2
# 1   eff_safe 35.627253 0.05428210
# 2      Short 13.006934 0.34544882
# 3     Middle 14.385842 0.26451427
# 4       Tall 16.510276 0.62757073
# 5 Biomeclass 27.364871 0.06872171
# 6   LifeForm  7.227416 0.41407835
# 7   LeafForm  8.185456 0.36921030