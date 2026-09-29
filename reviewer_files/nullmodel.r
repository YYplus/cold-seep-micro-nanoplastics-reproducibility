# 基于相对丰度划分稀有、丰富、条件稀有或丰富微生物类群
# 2024.11.11
# 读取OTU 绝对丰度表
otu <- read.delim('bacterial ASV.txt', row.names = 1)

# 计算每个样本的总读数（总丰度）
total_counts <- colSums(otu, na.rm = TRUE)

# 将绝对丰度转换为相对丰度
otu_relative <- otu / total_counts[col(otu)]

#（i）稀有类群（rare taxa，RT），在所有样本中丰度均 ≤0.1% 的 OTU
otu_ART <- otu_relative[apply(otu_relative, 1, function(x) max(x, na.rm = TRUE) <= 0.001), ]

#（ii）丰富类群（abundant taxa，AT），在所有样本中丰度均 ≥1% 的 OTU
otu_AAT <- otu_relative[apply(otu_relative, 1, function(x) min(x, na.rm = TRUE) >= 0.01), ]

#（iii）中间类群（moderate taxa，MT），在所有样本中丰度均 >0.1% 且 <1% 的 OTU
otu_MT <- otu_relative[apply(otu_relative, 1, function(x) min(x, na.rm = TRUE) > 0.001 & max(x, na.rm = TRUE) < 0.01), ]

#（iv）条件稀有类群（conditionally rare taxa，CRT），在所有样本中丰度均 <1%，且仅在部分样本中丰度 <0.1% 的 OTU
otu_CRT <- otu_relative[apply(otu_relative, 1, function(x) min(x, na.rm = TRUE) < 0.001 & max(x, na.rm = TRUE) < 0.01), ]
otu_CRT <- otu_CRT[which(! rownames(otu_CRT) %in% rownames(otu_ART)), ]  # CRT 和 ART 是没有重叠的，在所有样本中丰度均 ≤0.1% 的 OTU 是不可取的

#（v）条件丰富类群（conditionally abundant taxa，CAT），在所有样本中丰度均 >0.1%，且仅在部分样本中丰度 >1% 的 OTU
otu_CAT <- otu_relative[apply(otu_relative, 1, function(x) min(x, na.rm = TRUE) > 0.001 & max(x, na.rm = TRUE) > 0.01), ]
otu_CAT <- otu_CAT[which(! rownames(otu_CAT) %in% rownames(otu_AAT)), ]  # CAT 和 AAT 是没有重叠的，在所有样本中丰度均 ≥1% 的 OTU 是不可取的

#（vi）条件稀有或丰富类群（conditionally rare or abundant taxa，CRAT），丰度跨越从稀有（最低丰度 ≤0.1%）到丰富（最高丰度 ≥1%）的 OTU
otu_CRAT <- otu_relative[apply(otu_relative, 1, function(x) min(x, na.rm = TRUE) <= 0.001 & max(x, na.rm = TRUE) >= 0.01), ]

# 备注：这 6 个类群没有重叠，总数即等于 OTU 表的总数
otu[which(rownames(otu) %in% rownames(otu_ART)), 'taxa'] <- 'ART'
otu[which(rownames(otu) %in% rownames(otu_AAT)), 'taxa'] <- 'AAT'
otu[which(rownames(otu) %in% rownames(otu_MT)), 'taxa'] <- 'MT'
otu[which(rownames(otu) %in% rownames(otu_CRT)), 'taxa'] <- 'CRT'
otu[which(rownames(otu) %in% rownames(otu_CAT)), 'taxa'] <- 'CAT'
otu[which(rownames(otu) %in% rownames(otu_CRAT)), 'taxa'] <- 'CRAT'

# 将相对丰度转换回绝对丰度
otu_absolute <- otu_relative * total_counts[col(otu_relative)]

# 将绝对丰度和分类信息合并
otu_absolute <- cbind(otu_absolute, taxa = otu$taxa)

write.table(otu_absolute, 'otu.taxa.txt', col.names = NA, sep = '\t', quote = FALSE)