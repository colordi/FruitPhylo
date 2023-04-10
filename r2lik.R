library(phylolm)
library(rr2)
library(geiger)

library(tidyverse)

# 读取数据
data <- read_csv("data/YunnanFruit9370sppData.csv") # 性状数据
tree <-  read.tree("data/phyloa_sp_YNfruit9370tree.tre") # 系统树数据

# 为了进行系统发育分析，需要把字符串类型的性状数据转换为0/1的数值型数据
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

## 我们的数据是tibble格式的，没有rowname
## phyloglm必需依赖rowname，因此首先要添加rowname，并命名检查
row.names(data_01) <- data$species_name
name.check(tree,data_01) # 命名检查

# 我们希望通过构建两个模型之间的比较，来衡量不同因子之间的相对贡献
# 模型的构建
# 全模型构建，1表示常数项
mod_GF_CR_PHY_1 <- phyloglm(Fruittype ~ Growthform + Tropic + Subtropic + Temperate + 1, data = data_01, phy = tree)
# CR_PHY_1
mod_CR_PHY_1 <- phyloglm(Fruittype ~ Tropic + Subtropic + Temperate + 1, data = data_01, phy = tree)
# GF_PHY_1
mod_GF_PHY_1 <- phyloglm(Fruittype ~ Growthform + 1, data = data_01, phy = tree)
# PHY_1
mod_PHY_1 <- phyloglm(Fruittype ~ 1, data = data_01, phy = tree)
# GF_CR_1
mod_GF_CR_1 <- glm(Fruittype ~ Growthform + Tropic + Subtropic + Temperate + 1, data = data_01, family = "binomial")
# GF_1
mod_GF_1 <- glm(Fruittype ~ Growthform + 1, data = data_01, family = "binomial")
# CR_1
mod_CR_1 <- glm(Fruittype ~ Tropic + Subtropic + Temperate + 1, data = data_01, family = "binomial")
# 1
mod_1 <- glm(Fruittype ~ 1, data = data_01, family = "binomial")

# 计算模型的差异
get_r2like <- function(mod1, mod2, df) {
  # - mod1: 模型
  # - mod2: 需要减去的模型
  # - df: 因子的自由度，这里需要自定判断两个模型之间因子的差异量，例如本文中CR虽然名字只有两个字母，但其实表示了3个气候区
  # 输出结果是一个包含模型名、减去模型名、对应偏方系数、逻辑似然差值和p值的列表
  w <- c(as.character(substitute(mod1)),
         as.character(substitute(mod2)),
         as.character(R2_lik(mod1, mod2)),
         as.character(logLik(mod1)[[1]] - logLik(mod2)[[1]]),
         as.character(pchisq(2 * (logLik(mod1)[[1]] - logLik(mod2)[[1]]),
                             df = df, lower.tail = F))
  )
  return(w)
}

# 输出结果
df_fullphy <- data.frame(get_r2like(mod_GF_CR_PHY_1, mod_1, df = 5),
                         get_r2like(mod_GF_CR_PHY_1, mod_GF_CR_1, df = 1),
                         get_r2like(mod_GF_CR_PHY_1, mod_GF_PHY_1, df = 3),
                         get_r2like(mod_GF_CR_PHY_1, mod_CR_PHY_1, df = 1),
                         get_r2like(mod_GF_PHY_1, mod_GF_1, df = 1),
                         get_r2like(mod_GF_PHY_1, mod_PHY_1, df = 1),
                         get_r2like(mod_CR_PHY_1, mod_CR_1, df = 1),
                         get_r2like(mod_CR_PHY_1, mod_PHY_1, df = 3),
                         get_r2like(mod_GF_CR_1, mod_GF_1, df = 3),
                         get_r2like(mod_GF_CR_1, mod_CR_1, df = 1),
                         get_r2like(mod_GF_1, mod_1, df = 1),
                         get_r2like(mod_CR_1, mod_1, df = 3),
                         get_r2like(mod_PHY_1, mod_1, df = 1),
                         row.names = c("full_mod","reduce_mod","R2_lik","delta_log_lik","p")) %>%
  t() %>% as_tibble() %>%
  mutate(R2_lik = parse_number(R2_lik),
         delta_log_lik = parse_number(delta_log_lik),
         p = parse_number(p)) %>%
  mutate_if(is.numeric, round, 3)
df_fullphy
#    full_mod        reduce_mod   R2_lik delta_log_lik     p
#    <chr>           <chr>         <dbl>         <dbl> <dbl>
#  1 mod_GF_CR_PHY_1 mod_1         0.866       4317.   0
#  2 mod_GF_CR_PHY_1 mod_GF_CR_1   0.796       3034.   0
#  3 mod_GF_CR_PHY_1 mod_GF_PHY_1  0.009          9.92 0
#  4 mod_GF_CR_PHY_1 mod_CR_PHY_1  0.018         19.7  0
#  5 mod_GF_PHY_1    mod_GF_1      0.799       3097.   0
#  6 mod_GF_PHY_1    mod_PHY_1     0.01          11.1  0
#  7 mod_CR_PHY_1    mod_CR_1      0.855       4096.   0
#  8 mod_CR_PHY_1    mod_PHY_1     0.001          1.38 0.431
#  9 mod_GF_CR_1     mod_GF_1      0.015         73.1  0
# 10 mod_GF_CR_1     mod_CR_1      0.206       1081.   0
# 11 mod_GF_1        mod_1         0.228       1210.   0
# 12 mod_CR_1        mod_1         0.042        202.   0
# 13 mod_PHY_1       mod_1         0.863       4296.   0