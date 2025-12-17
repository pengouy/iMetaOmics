
###############################################################################
### Figure 1H

setwd("Figure1/BodyWeight/")
rm(list = ls())
options(encoding = "utf-8")
library(openxlsx)  
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)

# 导入数据
data <- read.xlsx("BodyWeight_transform.xlsx", sheet = 1, startRow = 21, 
                  colNames = TRUE, rowNames = TRUE,)
head(data)
# 转成长数据
m_data <- melt(data, value.name = "weight")
head(m_data)
m_data$dpi <- rep(c(0, seq(1,15,2)), 5) %>% as.character()
m_data$dpi <- factor(m_data$dpi,
                     levels = c("0", "1", "3", "5", "7", "9", "11", "13", "15"),
                     labels = c("0", "1", "3", "5", "7", "9", "11", "13", "15"))
m_data$group <- rep(c("Mock", "PC" ,"Abx-PC", "Young-FMT", "Adult-FMT"), each = 180)
head(m_data)
str(m_data)

# 绘制折线图
my_colors <- c("#ac75a0", "#c1282c", "#f4a637", "#008869", "#79aedf")
ggplot(m_data, aes(x = dpi, y = weight)) +
  expand_limits(x = c(0.5, 9.5), y = c(85, 355)) + # 扩展坐标轴范围
  geom_point(aes(color = group), alpha = 0.8, size = 3.5) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Days post infection (dpi)", y = "Body weight") +
  theme_bw() +
  theme(
    panel.border = element_rect(linewidth = 0.75),
    panel.grid.major = element_line(linewidth = 0.5),
    axis.line = element_line(linewidth = 0.5),
    legend.position.inside = c(0.15, 0.8),
    legend.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")) +
  geom_smooth(method = NULL, se = FALSE, aes(group = group, color = group))
ggsave("bodyweight_pointplot.pdf", width = 7.5, height = 4)


###############################################################################
### Figure S4E−J

wd <- "FigureS4/barplot/"
setwd(wd)
library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)

bardata <- list.files(wd, pattern = ".csv")

for (i in 1:length(bardata)){
  df <- read.csv(bardata[i])
  df <- df %>% 
  level <- stringr::str_split(bardata[i], pattern = "_")[[1]][1]
  
  mdt <- melt(df, OTU.ID = c("X"))
  head(mdt)
  
  p <- ggplot(mdt, aes(variable, value)) +
    geom_bar(aes(fill = X), 
             stat = "identity", position = "fill") +
    labs(x= "Group", y = "Proportion", 
         fill = paste0(level)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw()+
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 8, color = "black"),
          axis.text.x = element_text(angle = 45, 
                                     hjust = 1.1, vjust = 1.2),
          strip.text = element_text(size = 8),
          legend.title = element_text(size = 8))
  ggsave(paste0(level, "_Sample.pdf"), p, width = 8, height = 6)
  
}

###############################################################################
### Figure S9

iphopdir <- "iphop/"
iphopres <- read.csv(paste0(iphopdir, "Host_prediction_to_genus_m90.csv"), header = T)
head(iphopres)
colnames(iphopres)[1] <- "contigID"
colnames(iphopres)

iphopres <- merge(iphopres, viruscontig0, by = "contigID") %>% 
  select(contigID, AAI.to.closest.RaFAH.reference, Host.genus, 
         Confidence.score, List.of.methods, family)
head(iphopres)

# 提取第一个 ; 前的字符作为 method
iphopres$method <- sapply(strsplit(as.character(iphopres$List.of.methods), ";"), `[`, 1)

# 提取第一个 ; 后的第一个数字作为 method.score
iphopres$method.score <- sapply(strsplit(as.character(iphopres$List.of.methods), ";"), function(x) {
  # 提取第一个分号后的部分，并使用正则表达式提取第一个数字
  score <- ifelse(length(x) > 1, as.numeric(sub(".*?(\\d+\\.?\\d*).*", "\\1", x[2])), NA)
  return(score)
})

# 提取 Host.genus 列中最后一个分号后的字符，组成新列 extracted.host.genus
iphopres$extracted.host.genus <- sapply(strsplit(as.character(iphopres$Host.genus), ";"), function(x) {
  return(tail(x, n=1))  # 返回最后一个分号后的内容
})

# 生成唯一标识符，结合 contigID 和 extracted.host.genus
iphopres$unique_pair <- interaction(iphopres$contigID, iphopres$extracted.host.genus, sep = "_")

# 查看处理后的数据框
head(iphopres)
write.csv(iphopres, paste0(iphopdir, "iphopres.csv"))

library(ggplot2)
library(ggforce)
library(ggalluvial)
library(cowplot)
library(scales)

# 绘制桑基图，显示 family 和 contigID 以及 Host.genus
p_sankey <- ggplot(iphopres, 
                   aes(axis1 = family, axis2 = contigID, axis3 = extracted.host.genus, 
                       y = Confidence.score)) +
  geom_alluvium(aes(fill = family), width = 1/12, alpha = 0.25) +  # 设置链接颜色
  geom_stratum(width = 1/12, aes(fill = factor(family))) +  # 为每个 axis 设置不同颜色
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +  # 文字旋转
  scale_x_discrete(limits = c("family", "contigID", "Host.genus"), expand = c(0.15, 0.05)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")  # 旋转 x 轴文字


# 在桑基图右边用折线图表示 confidence.score
p_line_confidence <- ggplot(iphopres, aes(x = unique_pair, y = Confidence.score, group = 1)) +
  geom_line(size = 1, color = "#f4a637") +
  geom_point(size = 2, color = "#f4a637") +
  ylim(75, 100) +
  theme_minimal() +
  labs(y = NULL, x = "Total confidence score") +
  theme(axis.text.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  coord_flip()

# 在桑基图右侧显示 method.score，并用 method 颜色填充
p_bar_method_score <- ggplot(iphopres, 
                             aes(x = unique_pair, y = method.score, fill = method)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("blast" = "#7aafdf", "CRISPR" = "#c1282c", "iPHoP-RF" = "#008969")) +
  scale_y_continuous(limits = c(75, 100), oob = rescale_none) +  # 要加载rescale包并加入后面的参数
  theme_minimal() +
  labs(y = NULL, x = "Max method score") +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank()) +
  coord_flip()

# 调整图形长宽比及相对位置，保证y轴方向对齐
p_sankey <- p_sankey + theme(plot.margin = margin(0, 0, 0, 0))
p_line_confidence <- p_line_confidence + theme(plot.margin = margin(25, 5, 20, -30))
p_bar_method_score <- p_bar_method_score + theme(plot.margin = margin(25, 0, 20, 5))

# 将三个图形水平排列
combined_plot <- plot_grid(
  p_sankey,
  p_line_confidence,
  p_bar_method_score,
  ncol = 3,  # 设置列数为3，即水平排列
  rel_widths = c(4.5, 0.5, 1.25) # 调整宽度比例
)

# 打印最终的组合图
print(combined_plot)

ggsave(paste0(iphopdir, "iphop_plot.pdf"), combined_plot, width = 12, height = 8)


###############################################################################
### Figure S10

###先在 http://cloudtutu.com.cn 中跑出结果，再用R画
# 加载必要的包
library(igraph)
library(dplyr)
library(purrr)
library(stringr)

data_source <- "Young"  # 数据来源

# 读取边文件和顶点文件
data <- read.csv(paste0(data_source, "/edge2.csv"), header = TRUE, stringsAsFactors = FALSE)
nodes <- read.csv(paste0(data_source, "/node2.csv"), header = TRUE, stringsAsFactors = FALSE)
groups <- read.csv(paste0(data_source, "/group.csv"), header = TRUE, stringsAsFactors = FALSE) %>% 
  mutate(node = str_replace_all(node, "[ /]", "_"))  # 替换空格和斜杠为下划线

# 确保每种 type 提取前 10 的节点
filtered_groups <- nodes %>%
  inner_join(groups, by = c("Id" = "node")) %>% # 合并 nodes 和 groups 数据
  group_by(type) %>%
  slice_max(order_by = Abundance, n = 10) %>% # 确保每组取前 10
  ungroup()


# 创建图对象（无向图）
g <- graph_from_data_frame(d = data, directed = FALSE, vertices = nodes)

# 检查图的基本信息
print(g)

# 设置顶点大小：根据 Abundance 调整
V(g)$size <- log(nodes$Abundance[match(V(g)$name, nodes$Id)] + 1, 10)

# 设置顶点颜色：根据 module 分组
module_colors <- rainbow(length(unique(nodes$module)))  # 自动生成颜色
names(module_colors) <- unique(nodes$module)
V(g)$color <- module_colors[as.character(nodes$module[match(V(g)$name, nodes$Id)])]

# 设置顶点形状：根据 type 信息
shape_mapping <- c("Microbiome" = "circle", "Virome" = "square")  # 定义形状映射
V(g)$shape <- shape_mapping[groups$type[match(V(g)$name, groups$node)]]  # 根据 type 设置形状
V(g)$shape[is.na(V(g)$shape)] <- "square"  # 缺失值设置为圆形

# 设置顶点标签：仅保留筛选后的节点
V(g)$label <- ifelse(V(g)$name %in% filtered_groups$Id, V(g)$name, NA)

# 设置边的颜色：根据权重正负性
E(g)$color <- ifelse(E(g)$Weight > 0, "#ff878c", "#5ea6c2")  # 红色表示正相关，蓝色表示负相关

# 设置边的宽度：根据权重绝对值
E(g)$width <- abs(E(g)$Weight) ^3 #让线粗细的差别更明显

# 创建网络图布局
layout_graphopt <- layout_with_kk(g)

# 保存为 PDF 文件
pdf(paste0(data_source, "_network.pdf"), width = 10, height = 8)

# 绘制网络图
plot.igraph(
  g,
  layout = layout_graphopt,          # 使用优化布局
  vertex.color = V(g)$color,         # 顶点颜色
  vertex.size = V(g)$size,           # 顶点大小
  vertex.frame.color = "white",      # 顶点框颜色
  vertex.frame.width = 0.1,          # 顶点框宽度
  vertex.label = V(g)$label,         # 顶点标签，仅筛选后的节点
  vertex.label.family = "sans",      # 顶点标签字体
  vertex.label.font = 1,             # 顶点标签字体样式
  vertex.label.cex = 0.4,            # 顶点标签大小
  vertex.label.dist = 0.1,           # 顶点标签距离
  vertex.label.degree = 0,           # 顶点标签角度
  vertex.label.color = "black",      # 顶点标签颜色
  vertex.shape = V(g)$shape,         # 顶点形状
  edge.width = E(g)$width            # 边宽度
)

# 添加模块图例
legend(
  x = 0.8, y = 1, 
  title = "Module", 
  legend = names(module_colors), 
  fill = module_colors, 
  border = "white", 
  bty = "n",
  cex = 0.5
)

# 添加节点类型图例
legend(
  x = 0.8, y = 0.5, 
  title = "Node Type", 
  legend = c("Microbiome", "Virome"), 
  pch = c(21, 22),  # 圆形和方形
  pt.bg = c("white", "black"),  # 填充颜色
  pt.cex = 1,  # 点大小
  border = "white", 
  bty = "n",
  cex = 0.5
)

# 添加边权重图例
legend(
  x = 0.8, y = 0.25, 
  title = "Weight", 
  legend = c("Positive", "Negative"), 
  col = c("#ff878c", "#5ea6c2"), 
  lty = 1, 
  lwd = 2, 
  bty = "n",
  cex = 0.5
)

# 添加节点丰度大小图例
legend(
  x = 0.8, y = 0.1, 
  title = "Node Size (Abundance)",
  legend = c("Low", "Medium", "High"),
  pt.cex = c(0.2, 0.4, 0.6),  # 节点大小表示丰度
  pch = 21,             # 圆形节点
  pt.bg = "white",        # 节点填充颜色
  border = "black",
  bty = "n",
  cex = 0.5
)

dev.off()


#############################################################################
### Figure S11D

rm(list=ls())
setwd("FigureS11/mantel/")
library(linkET)
library(ggplot2)
library(ggtext)
library(dplyr)
library(cols4all)

#组合网络热图绘制
#读入数据：
varespec <- read.csv("virome_metabolites.csv", header = T, row.names = 1)

varechem <- read.csv("16s_family.csv", header = T, row.names = 1)
#计算环境因子相关性系数：
cor2 <- correlate(varechem)
corr2 <- cor2 %>% as_md_tbl() %>%
  filter(is.na(r) == F)
write.csv(corr2, file = "pearson_correlate(env&env).csv", row.names = TRUE)

head(corr2)
#mantel test:
mantel <- mantel_test(varespec, varechem,
                      mantel_fun = 'mantel', #支持4种："mantel"使用vegan::mantel()；"mantel.randtest"使用ade4::mantel.randtest()；"mantel.rtest"使用ade4::mantel.rtest()；"mantel.partial"使用vegan::mantel.partial()
                      spec_select = list(Virome = 1:50,
                                         Metabolism = 51:3100)) #对应两个数据集的行索引
head(mantel)
write.csv(mantel, file = "mantel_result(bio&env).csv", row.names = TRUE)
#对mantel的r和P值重新赋值（设置绘图标签）：
mantel2 <- mantel %>%
  mutate(r = cut(r, breaks = c(-Inf, 0.25, 0.5, 0.75, Inf),
                 labels = c("<0.25", "0.25-0.5", "0.5-0.75", ">= 0.75")),
         p = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                 labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">= 0.05")))
head(mantel2)
#首先，绘制相关性热图(和上文相同):
p4 <- qcorrplot(cor2,
                grid_col = "grey50",
                grid_size = 0.2,
                type = "upper",
                diag = FALSE) +
  geom_square() +
  scale_fill_gradientn(colours = c4a('rd_bu',30),
                       limits = c(-1, 1))
p4

#添加显著性标签：
p5 <- p4 +
  geom_mark(size = 4,
            only_mark = T,
            sig_level = c(0.05, 0.01, 0.001),
            sig_thres = 0.05,
            colour = 'white')
p5
#在相关性热图上添加mantel连线：
p6 <- p5 +
  geom_couple(data = mantel2,
              aes(colour = p, size = r),
              curvature = nice_curvature())
p6
#继续美化连线：
p7 <- p6 +
  scale_size_manual(values = c(0.5, 1, 1.5, 2)) + #连线粗细
  scale_colour_manual(values = c("#57994d", "#69638d", "#aeb04c", "#cfcfce")) + #连线配色
  #修改图例：
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"),
                             order = 2),
         colour = guide_legend(title = "Mantel's p",
                               override.aes = list(size = 3),
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
p7

ggsave("mantel.pdf", p7, width = 10, height = 10)
