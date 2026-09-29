# 加载必要的包
library(vegan)
library(dplyr)

# 定义函数：读取距离矩阵和元数据，运行 PERMANOVA
run_permanova <- function(beta_file, group_file, source) {
  # 检查文件是否存在
  if (!file.exists(beta_file)) {
    stop("β多样性文件不存在：", beta_file)
  }
  if (!file.exists(group_file)) {
    stop("元数据文件不存在：", group_file)
  }
  
  # 读取距离矩阵
  beta_data <- read.table(beta_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  # 检查是否为方阵
  if (nrow(beta_data) != ncol(beta_data)) {
    stop("距离矩阵不是方阵：", beta_file)
  }
  
  # 转换为 dist 对象
  beta_dist <- as.dist(beta_data)
  
  # 读取元数据
  group_data <- read.table(group_file, header = TRUE, sep = "\t")
  
  # 检查必要列
  if (!all(c("Sample", "Group") %in% colnames(group_data))) {
    stop("group.txt 缺少 'Sample' 或 'Group' 列")
  }
  
  # 确保样本名匹配
  dist_samples <- labels(beta_dist)
  group_data <- group_data %>% filter(Sample %in% dist_samples)
  
  if (nrow(group_data) == 0) {
    stop("元数据中没有与距离矩阵样本匹配的样本：", beta_file)
  }
  
  # 检查样本是否完全匹配
  if (!all(dist_samples %in% group_data$Sample)) {
    warning("部分距离矩阵样本在元数据中缺少分组信息：", 
            paste(setdiff(dist_samples, group_data$Sample), collapse = ", "))
  }
  
  # 设置行名为样本名，确保顺序与距离矩阵一致
  rownames(group_data) <- group_data$Sample
  group_data <- group_data[dist_samples, ]
  
  # 运行 PERMANOVA
  permanova_result <- adonis2(beta_dist ~ Group, data = group_data, permutations = 999)
  
  # 打印结果并添加数据来源
  cat("\nPERMANOVA 结果 for", source, ":\n")
  print(permanova_result)
  
  # 返回结果
  return(permanova_result)
}

# 处理三个 β 多样性文件
beta_files <- c("beta_div_ff.txt", "beta_div_xy.txt", "beta_div_mt.txt")
sources <- c("ff", "xy", "mt")
group_file <- "group.txt"

# 存储所有结果
permanova_results <- list()

# 循环运行 PERMANOVA
for (i in seq_along(beta_files)) {
  cat("\n处理文件：", beta_files[i], "\n")
  permanova_results[[sources[i]]] <- run_permanova(beta_files[i], group_file, sources[i])
}

# 可选：保存所有结果到文件
for (source in sources) {
  sink(paste0("permanova_result_", source, ".txt"))
  cat("PERMANOVA 结果 for", source, ":\n")
  print(permanova_results[[source]])
  sink()
}

# 可选：提取关键统计信息汇总
summary_table <- do.call(rbind, lapply(names(permanova_results), function(source) {
  res <- permanova_results[[source]]
  data.frame(
    Source = source,
    Df_Group = res$Df[1],
    F_Value = res$F[1],
    R2 = res$R2[1],
    P_Value = res$`Pr(>F)`[1]
  )
}))

# 打印汇总表
cat("\nPERMANOVA 汇总表：\n")
print(summary_table)

# 保存汇总表
write.table(summary_table, "permanova_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)