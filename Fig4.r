
# 载入所需的R包
library(ropls)
library(ggplot2)
library(openxlsx)
library(tidyverse)


setwd("01_correlation\\PLS-DA")
data <- read.table("input_raw.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
dim(data)
data_new <- data[c(2:8, 10), ]

group     <- factor(c(rep("Control", 5), rep("PS", 5)))
group_new <- factor(c(rep("Control", 4), rep("PS", 4)))

# ---- 1. 原始数据建模与绘图 ----
pdf(file = "PLS-DA_Output_Raw.pdf", height = 8, width = 8)
plsda <- opls(data, group, predI = 2)
dev.off()

pdf(file = "OPLS-DA_Output_Raw.pdf", height = 6, width = 6)
oplsda <- opls(data, group, predI = 1, orthoI = NA) # orthoI = NA 执行 OPLS
dev.off()

# 导出原始得分图
pdf(file = "PLS-DA-group_Raw.pdf", height = 6, width = 6)
plot(plsda, typeVc = "x-score", parAsColFcVn = group, parPaletteVc = c("green4", "orange"))
dev.off()

pdf(file = "OPLS-DA-group_Raw.pdf", height = 6, width = 6)
plot(oplsda, typeVc = "x-score", parAsColFcVn = group, parPaletteVc = c("green4", "orange"))
dev.off()

# 合并原始 VIP 并导出
df1 <- data.frame(t(data), opls_VIP = oplsda@vipVn, pls_VIP = plsda@vipVn)
write.csv(df1, file = "dataMatrix_VIP_raw.csv")


# ---- 2. 过滤后数据建模与绘图 ----
pdf(file = "PLS-DA_Output_Filtered.pdf", height = 8, width = 8)
plsda_new <- opls(data_new, group_new, predI = 2)
dev.off()

pdf(file = "OPLS-DA_Output_Filtered.pdf", height = 6, width = 6)
oplsda_new <- opls(data_new, group_new, predI = 1, orthoI = NA)
dev.off()

# 导出过滤后得分图
pdf(file = "PLS-DA-group-filtered.pdf", height = 6, width = 6)
plot(plsda_new, typeVc = "x-score", parAsColFcVn = group_new, parPaletteVc = c("green4", "orange"))
dev.off()

pdf(file = "OPLS-DA-group-filtered.pdf", height = 6, width = 6)
plot(oplsda_new, typeVc = "x-score", parAsColFcVn = group_new, parPaletteVc = c("green4", "orange"))
dev.off()

# 合并过滤后 VIP 并导出
df_new_1 <- data.frame(t(data_new), opls_VIP_new = oplsda_new@vipVn, pls_VIP_new = plsda_new@vipVn)
write.csv(df_new_1, file = "dataMatrix_VIP_filtered.csv")

# ---- 1. 原始数据差异分析 ----
df1$mean_control <- rowMeans(df1[, 1:5], na.rm = TRUE)
df1$mean_PS50    <- rowMeans(df1[, 6:10], na.rm = TRUE)

# 矩阵化执行 T 检验
df1$pvalue  <- apply(df1, 1, function(x) t.test(as.numeric(x[6:10]), as.numeric(x[1:5]))$p.value)
df1$pvalue1 <- apply(df1, 1, function(x) t.test(as.numeric(x[6:10]), as.numeric(x[1:5]), var.equal = TRUE)$p.value)

# 计算 Fold Change 并处理特殊值
df1$foldchange     <- df1$mean_PS50 / df1$mean_control
df1$log2foldchange <- log2(df1$foldchange)
df1$log2foldchange[is.na(df1$log2foldchange) | is.infinite(df1$log2foldchange)] <- 0
df1$pvalue[is.na(df1$pvalue) | is.infinite(df1$pvalue)]   <- 1
df1$pvalue1[is.na(df1$pvalue1) | is.infinite(df1$pvalue1)] <- 1

# 标注上下调状态
df1$up_down <- case_when(
  df1$log2foldchange > 1 ~ "up",
  df1$log2foldchange < -1 ~ "down",
  TRUE ~ "none"
)
df1$sig <- ifelse(df1$pvalue < 0.05, "yes", "no")

# 导出原始结果
write.xlsx(df1, file = "diff_metabolic_result_raw.xlsx", rowNames = TRUE)
write.table(df1, file = "diff_metabolic_result_raw.xls", sep = "\t", row.names = TRUE)


# ---- 2. 过滤后数据差异分析 ----
df_new_1$mean_con <- rowMeans(df_new_1[, 1:4], na.rm = TRUE)
df_new_1$mean_PS  <- rowMeans(df_new_1[, 5:8], na.rm = TRUE)

df_new_1$pvalue  <- apply(df_new_1, 1, function(x) t.test(as.numeric(x[5:8]), as.numeric(x[1:4]))$p.value)
df_new_1$pvalue1 <- apply(df_new_1, 1, function(x) t.test(as.numeric(x[5:8]), as.numeric(x[1:4]), var.equal = TRUE)$p.value)

df_new_1$foldchange     <- df_new_1$mean_PS / df_new_1$mean_con
df_new_1$log2foldchange <- log2(df_new_1$foldchange)
df_new_1$log2foldchange[is.na(df_new_1$log2foldchange) | is.infinite(df_new_1$log2foldchange)] <- 0
df_new_1$pvalue[is.na(df_new_1$pvalue) | is.infinite(df_new_1$pvalue)]   <- 1
df_new_1$pvalue1[is.na(df_new_1$pvalue1) | is.infinite(df_new_1$pvalue1)] <- 1

df_new_1$up_down <- case_when(
  df_new_1$log2foldchange > 1 ~ "up",
  df_new_1$log2foldchange < -1 ~ "down",
  TRUE ~ "none"
)
df_new_1$sig <- ifelse(df_new_1$pvalue < 0.05, "yes", "no")

# 导出过滤后结果
write.table(df_new_1, file = "diff_metabolic_result_filtered.xls", sep = "\t", row.names = TRUE)


setwd("02_diff_metabolics\\2_volcano_plot")
df_new_1$group <- case_when(
  df_new_1$pvalue1 < 0.05 & df_new_1$pls_VIP_new > 1 & df_new_1$log2foldchange > 0 ~ "up",
  df_new_1$pvalue1 < 0.05 & df_new_1$pls_VIP_new > 1 & df_new_1$log2foldchange <= 0 ~ "down",
  TRUE ~ "NS"
)
table(df_new_1$group)

# 2. 基础火山图绘制
p <- ggplot(df_new_1, aes(x = log2foldchange, y = -log10(pvalue1))) + 
  geom_point(aes(color = group), alpha = 0.8, size = 2) +
  scale_color_manual(values = c('down' = '#679AFD', 'NS' = 'grey', 'up' = '#EF98A1')) +
  geom_vline(xintercept = 0, lty = 3, color = 'black', lwd = 0.5) + 
  geom_hline(yintercept = -log10(0.05), lty = 3, color = 'black', lwd = 0.5) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  labs(x = 'log2 fold change', y = '-log10 pvalue')
p

# 3. 高亮并标记目标代谢物 (CAMP)
df_new_1$gene <- row.names(df_new_1)
df_new_1$group <- as.character(df_new_1$group)
df_new_1[df_new_1$gene == 'M328T159_2', 'group'] <- 'CAMP'
df_new_1$group <- as.factor(df_new_1$group)

p1 <- ggplot(df_new_1, aes(x = log2foldchange, y = -log10(pvalue1))) + 
  geom_point(aes(color = group), alpha = 0.8, size = 2) +
  scale_color_manual(values = c('CAMP' = 'red', 'down' = '#679AFD', 'NS' = 'grey', 'up' = '#EF98A1')) +
  geom_vline(xintercept = 0, lty = 3, color = 'black', lwd = 0.5) + 
  geom_hline(yintercept = -log10(0.05), lty = 3, color = 'black', lwd = 0.5) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  labs(x = 'log2 fold change', y = '-log10 pvalue')
ggsave("volplot_ps50_con_CAMP_1.pdf", p1, width = 5, height = 4.7)

# 4. 剔除离群极端点后的精细火山图
row_names_to_remove <- c("M225T424", "M272T229")
df_filtered_volcano <- subset(df_new_1, !gene %in% row_names_to_remove)  

p2 <- ggplot(df_filtered_volcano, aes(x = log2foldchange, y = -log10(pvalue1))) + 
  geom_point(aes(color = group), alpha = 0.9, size = 1.2) +
  scale_color_manual(values = c('CAMP' = 'red', 'down' = '#679AFD', 'NS' = 'grey', 'up' = '#EF98A1')) +
  geom_vline(xintercept = 0, lty = 3, color = 'black', lwd = 0.5) + 
  geom_hline(yintercept = -log10(0.05), lty = 3, color = 'black', lwd = 0.5) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  labs(x = 'log2 fold change', y = '-log10 pvalue')
ggsave("volplot_ps50_con_CAMP.pdf", p2, width = 5, height = 4)
