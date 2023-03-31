library(tidyverse)
library(scales)

# 读取数据
data <- read_csv("data/YunnanFruit9370sppData.csv")

# 气候区
# 数据转换为0/1，并汇总求和绘图
(p <- data %>%
  mutate(across(c(Tropic, Temperate, Subtropic),
                ~ifelse(. == "Absence", 0, 1))) %>%
  group_by(Fruittype) %>%
  summarize(across(c(Tropic, Temperate, Subtropic), sum)) %>%
  pivot_longer(cols = c(Tropic, Temperate, Subtropic),
               values_to = "nums",
               names_to = "CR") %>%
  ggplot(mapping = aes(x = factor(CR, levels = c("Tropic", "Temperate", "Subtropic")),
                       y = nums,
                       fill = Fruittype)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = nums), position = position_fill()) +
  scale_y_continuous(labels = percent_format())
)
# save plot
ggsave("result/climate_regions_barplot.png",plot=p,dpi=300)

# 生长型
(p <- data %>%
  group_by(Fruittype, Growthform) %>%
  summarize(nums = n()) %>%
  ggplot(mapping = aes(x = Growthform, y = nums, fill = Fruittype)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = nums), position = position_fill()) +
  scale_y_continuous(labels = percent_format())
)
# save plot
ggsave("result/growth_forms,_barplot.png",plot=p,dpi=300)

# spearman相关分析
## 相关性分析
data %>%
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
                          false = 0)
         ) %>% 
  cor(method = "spearman") %>% 
  corrplot::corrplot(method = "square",addCoef.col = "black",type = "lower")
