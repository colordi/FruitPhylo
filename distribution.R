library(tidyverse)
library(scales)
library(corrplot)
library(dplyr)
library(scales)


# 读取数据
treeheight <- read.csv("data/1843.csv",header = TRUE)
data_682 <- read.csv("data/682.csv",header = TRUE)

##匹配两个数据集

treeheight_subset <- treeheight %>%
  group_by(Species) %>%
  slice(which.max(Hmax))

treeheight_subset_682 <- treeheight_subset %>% 
  subset(Species %in% data_682$Species)

merged_data <- merge(data_682, treeheight_subset_682[, c("Species", "Hmax")], 
                     by = "Species", all.x = TRUE)

##对被子植物和裸子植物的残差列进行合并
merged_data$Residualall <- ifelse(is.na(merged_data$ResidualAng), 
                                  merged_data$ResidualGym, merged_data$ResidualAng)
##将Residualall改为判断式的，Residualall为负时效率和安全性较低，为0，为正时效率和安全性较高，为1
merged_data$eff_safe <- ifelse(merged_data$Residualall < 0, 0, 1)

##画相对贡献
merged_data %>%
  filter(Life.form %in% c("T", "S")) %>%
  mutate(eff_safe_binary = ifelse(eff_safe == 0, "Low Efficiency and Safety", "High Efficiency and Safety")) %>%
  group_by(Life.form, eff_safe_binary) %>%
  summarize(n = n()) %>%
  group_by(Life.form) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(mapping = aes(x = factor(Life.form, labels = c("Tree", "Shrub")), y = prop, fill = eff_safe_binary)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  labs(x = "Life Form", y = "Distribution", fill = "Efficiency and Safety") +
  theme_bw() +
  scale_fill_manual(values = c("#519259", "#F0BB62")) +
  scale_y_continuous(labels = percent_format()) +
  ggtitle("Efficiency and Safety Distribution by Life Form")

merged_data %>%
  filter(Leaf.form %in% c("E", "D")) %>%
  mutate(eff_safe_binary = ifelse(eff_safe == 0, "Low Efficiency and Safety", "High Efficiency and Safety")) %>%
  group_by(Leaf.form, eff_safe_binary) %>%
  summarize(n = n()) %>%
  group_by(Leaf.form) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(mapping = aes(x = factor(Leaf.form, labels = c("Evergreen", "Deciduous")), y = prop, fill = eff_safe_binary)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  labs(x = "Leaf Form", y = "Distribution", fill = "Efficiency and Safety") +
  theme_bw() +
  scale_fill_manual(values = c("#519259", "#F0BB62")) +
  scale_y_continuous(labels = percent_format()) +
  ggtitle("Efficiency and Safety Distribution by Leaf Form")


merged_data %>%
  filter(Biome %in% c("TRR", "TRS", "TMR", "TMS", "WDS", "DES", "BOR")) %>%
  mutate(eff_safe_binary = ifelse(eff_safe == 0, "Low Efficiency and Safety", "High Efficiency and Safety")) %>%
  group_by(Biome, eff_safe_binary) %>%
  summarize(n = n()) %>%
  group_by(Biome) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(mapping = aes(x = factor(Biome, labels = c("Tropical Rainforest", "Tropical Seasonal Forest", "Temperate Rainforest", "Temperate Seasonal Forest", "Woodland/Shrubland", "Deserts", "Boreal Forest")), y = prop, fill = eff_safe_binary)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  labs(x = "Biome", y = "Distribution", fill = "Efficiency and Safety") +
  theme_bw() +
  scale_fill_manual(values = c("#519259", "#F0BB62")) +
  scale_y_continuous(labels = percent_format()) +
  ggtitle("Efficiency and Safety Distribution by Biome") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

merged_data %>%
  filter(!is.na(Soiltype)) %>%
  mutate(eff_safe_binary = ifelse(eff_safe == 0, "Low Efficiency and Safety", "High Efficiency and Safety")) %>%
  group_by(Soiltype, eff_safe_binary) %>%
  summarize(n = n()) %>%
  group_by(Soiltype) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(mapping = aes(x = factor(Soiltype), y = prop, fill = eff_safe_binary)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  labs(x = "Soil Type", y = "Distribution", fill = "Efficiency and Safety") +
  theme_bw() +
  scale_fill_manual(values = c("#519259", "#F0BB62")) +
  scale_y_continuous(labels = percent_format()) +
  ggtitle("Efficiency and Safety Distribution by Soil Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##土壤性状分类后

merged_data %>%
  filter(!is.na(Soiltype)) %>%
  mutate(Soil_class = case_when(
    Soiltype %in% c("Gleysol", "Histosol") ~ "Hydric Soils",
    Soiltype %in% c("Podzol", "Cambisol") ~ "Forest Soils",
    Soiltype %in% c("Chernozem", "Phaeozem", "Kastanozem", "Luvisol") ~ "Grassland Soils",
    TRUE ~ "Arid Soils"
  )) %>%
  mutate(eff_safe_binary = ifelse(eff_safe == 0, "Low Efficiency and Safety", "High Efficiency and Safety")) %>%
  group_by(Soil_class, eff_safe_binary) %>%
  summarize(n = n()) %>%
  group_by(Soil_class) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(mapping = aes(x = factor(Soil_class), y = prop, fill = eff_safe_binary)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  labs(x = "Soil Type", y = "Distribution", fill = "Efficiency and Safety") +
  theme_bw() +
  scale_fill_manual(values = c("#519259", "#F0BB62")) +
  scale_y_continuous(labels = percent_format()) +
  ggtitle("Efficiency and Safety Distribution by Soil Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##树高性状分类后

merged_data$Height_category <- cut(merged_data$Hmax, breaks = c(0, 45, 100, Inf), labels = c("short", "middle", "tall"))

merged_data %>%
  filter(!is.na(Height_category)) %>%
  mutate(eff_safe_binary = ifelse(eff_safe == 0, "Low Efficiency and Safety", "High Efficiency and Safety")) %>%
  group_by(Height_category, eff_safe_binary) %>%
  summarize(n = n()) %>%
  group_by(Height_category) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(na.omit(merged_data), mapping = aes(x = factor(Height_category, labels = c("Short", "Middle", "Tall")), y = prop, fill = eff_safe_binary)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
  labs(x = "Height Category", y = "Distribution", fill = "Efficiency and Safety") +
  theme_bw() +
  scale_fill_manual(values = c("#519259", "#F0BB62")) +
  scale_y_continuous(labels = percent_format()) +
  ggtitle("Efficiency and Safety Distribution by Height Category")


# 因为数据集可能存在重复的物种，需要添加后缀
merged_data <- merged_data %>%
  # 把物种名中的空格替换为_
  mutate(Species = str_replace(Species, " ", "_")) %>% 
  group_by(Species) %>% # 以 Species 列为依据对数据框进行分组
  mutate(suffix = ifelse(row_number() == 1, "", row_number()),
         `species.relative` = Species,
         Species = paste0(Species, suffix)
         )
## 导出到本地
write.csv(merged_data, file = "merged_data.csv", row.names = FALSE)

