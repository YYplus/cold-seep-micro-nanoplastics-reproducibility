library(Hmisc)
library(minpack.lm)
library(stats4)
library(extrafont)  # 用于使用自定义字体
loadfonts()  # 加载字体

# 定义函数分析单个文件
analyze_neutral_model <- function(input_file) {
  # 提取文件名（不含扩展名）作为输出文件名前缀
  file_prefix <- tools::file_path_sans_ext(basename(input_file))
  cat("\n开始处理文件:", input_file, "\n")
  
  # 读取数据
  spp_full <- read.csv(input_file, head=T, stringsAsFactors=F, row.names=1, sep="\t")
  
  # 去掉最后一列(taxa列)
  spp <- spp_full[, 1:(ncol(spp_full)-1)]
  
  # 转置矩阵
  spp <- t(spp)
  
  # 计算平均个体数
  N <- mean(apply(spp, 1, sum))
  p.m <- apply(spp, 2, mean)
  p.m <- p.m[p.m != 0]
  p <- p.m/N
  
  # 计算频率
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  # 合并数据
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
  
  # 检查是否有足够的数据点进行拟合
  if(nrow(C.0) < 5) {
    cat("警告: 数据点太少，无法进行可靠拟合\n")
    return(list(file = input_file, Nm = NA, Rsqr = NA, success = FALSE))
  }
  
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  d = 1/N
  
  # 尝试多个初始参数值
  initial_values <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5)
  m.fit <- NULL
  
  for(m_init in initial_values) {
    cat("尝试初始参数 m =", m_init, "\n")
    
    # 尝试拟合模型，捕获任何错误
    fit_attempt <- tryCatch({
      # 增加最大迭代次数和容差
      nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), 
            start=list(m=m_init),
            control=nls.lm.control(maxiter=500, ftol=1e-8, ptol=1e-8))
    }, error = function(e) {
      cat("  拟合失败:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if(!is.null(fit_attempt)) {
      m.fit <- fit_attempt
      cat("  拟合成功，参数 m =", coef(m.fit), "\n")
      break
    }
  }
  
  # 检查是否成功拟合模型
  if(is.null(m.fit)) {
    cat("错误: 所有拟合尝试均失败，跳过此文件\n")
    return(list(file = input_file, Nm = NA, Rsqr = NA, success = FALSE))
  }
  
  # 计算拟合质量
  m.ci <- tryCatch({
    confint(m.fit, 'm', level=0.95)
  }, error = function(e) {
    cat("  无法计算置信区间:", conditionMessage(e), "\n")
    return(c(NA, NA))
  })
  
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  
  # 如果预测值接近0或1，可能导致binconf函数出错，添加保护
  freq.pred.safe <- pmin(pmax(freq.pred, 0.001), 0.999)
  
  pred.ci <- tryCatch({
    binconf(freq.pred.safe*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  }, error = function(e) {
    cat("  无法计算预测置信区间:", conditionMessage(e), "\n")
    # 创建一个模拟的置信区间数据框
    data.frame(
      PointEst = freq.pred,
      Lower = freq.pred - 0.1,
      Upper = freq.pred + 0.1
    )
  })
  
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  cat("R² = ", Rsqr, "\n")
  
  # 创建bacnlsALL数据框，包含所有绘图所需的数据
  bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[,2:3])
  names(bacnlsALL)[4:5] <- c("Lower", "Upper")  # 重命名置信区间列名
  
  # 确保置信区间在[0,1]范围内
  bacnlsALL$Lower <- pmax(0, bacnlsALL$Lower)
  bacnlsALL$Upper <- pmin(1, bacnlsALL$Upper)
  
  # 设置点的颜色
  inter.col <- rep('black', nrow(bacnlsALL))  # 默认点颜色为黑色
  inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- '#A52A2A'  # 将低于预测下限的点设为棕色
  inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- '#29A6A6'  # 将高于预测上限的点设为青色
  
  # 设置PDF输出
  pdf_file <- paste0("ncm_", file_prefix, ".pdf")
  pdf(pdf_file, width=7, height=7)  # 创建PDF文件，设置宽高为7英寸
  
  # 绘图代码
  library(grid)  
  grid.newpage()  
  pushViewport(viewport(h=0.6, w=0.6))  
  pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(-0.02,1.02), extension=c(0.02,0)))  
  grid.rect() 
  grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch=20, gp=gpar(col=inter.col, cex=0.7))  
  grid.yaxis(at=c(0, 0.25, 0.50, 0.75, 1.00), label=c("0", "0.25", "0.50", "0.75", "1.00"), gp=gpar(fontsize=16))  
  grid.xaxis(gp=gpar(fontsize=16))  
  grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp=gpar(col='blue', lwd=2), default='native') 
  
  grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp=gpar(col='blue', lwd=2, lty=2), default='native')  
  grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp=gpar(col='blue', lwd=2, lty=2), default='native')  
  #grid.text(y=unit(0,'npc')-unit(2.5,'lines'), label='Mean Relative Abundance (log10)', gp=gpar(fontface=2, fontfamily="Times New Roman"))  
  #grid.text(x=unit(0,'npc')-unit(3,'lines'), label='Frequency of Occurance', gp=gpar(fontface=2, fontfamily="Times New Roman"), rot=90)  
  
  # 添加R²和Nm值到右下方位置
  # 分别显示R²和Nm值，使用expression来正确显示上标
  grid.text(bquote(R^2 == .(round(Rsqr, 3))), 
            x=unit(0.85, "npc"), y=unit(0.15, "npc"), just=c("centre", "bottom"), 
            gp=gpar(fontfamily="Times New Roman", fontsize=17))
  grid.text(paste("Nm=", round(coef(m.fit)*N, 3)), 
            x=unit(0.85, "npc"), y=unit(0.08, "npc"), just=c("centre", "bottom"), 
            gp=gpar(fontfamily="Times New Roman", fontsize=17))
  
  # 关闭PDF设备并保存文件
  dev.off()
  
  cat("已生成图形文件:", pdf_file, "\n")
  
  # 创建饼图显示三种过程
  neutral_count <- sum(inter.col == 'black')
  negative_selection_count <- sum(inter.col == '#A52A2A')
  positive_selection_count <- sum(inter.col == '#29A6A6')
  
  # 计算百分比
  total_count <- length(inter.col)
  neutral_pct <- round(neutral_count / total_count * 100, 1)
  negative_pct <- round(negative_selection_count / total_count * 100, 1)
  positive_pct <- round(positive_selection_count / total_count * 100, 1)
  
  # 创建饼图数据
  pie_data <- data.frame(
    process = c("中性过程", "负向选择", "正向选择"),
    count = c(neutral_count, negative_selection_count, positive_selection_count),
    color = c("black", "#A52A2A", "#29A6A6"),
    percentage = c(neutral_pct, negative_pct, positive_pct)
  )
  
  # 只绘制有数据的过程
  pie_data <- pie_data[pie_data$count > 0, ]
  
  if(nrow(pie_data) > 0) {
    # 设置饼图PDF输出
    pie_pdf_file <- paste0("pie_", file_prefix, ".pdf")
    pdf(pie_pdf_file, width=6, height=6)
    
    # 绘制饼图
    library(grid)
    grid.newpage()
    pushViewport(viewport(h=0.8, w=0.8))
    
    # 计算角度
    angles <- cumsum(pie_data$count) / sum(pie_data$count) * 2 * pi
    start_angles <- c(0, angles[-length(angles)])
    
    # 绘制每个扇形
    for(i in 1:nrow(pie_data)) {
      n_points <- 100
      theta <- seq(start_angles[i], angles[i], length.out = n_points)
      x_coords <- c(0.5, 0.5 + 0.4 * cos(theta), 0.5)
      y_coords <- c(0.5, 0.5 + 0.4 * sin(theta), 0.5)
      
      grid.polygon(x = unit(x_coords, "npc"), y = unit(y_coords, "npc"),
                   gp = gpar(fill = pie_data$color[i], col = "white", lwd = 2))
    }
    
    # 添加百分比标签
    for(i in 1:nrow(pie_data)) {
      # 计算标签位置（稍微向外移动以避免遮挡）
      mid_angle <- (start_angles[i] + angles[i]) / 2
      label_x <- 0.5 + 0.25 * cos(mid_angle)  # 从0.25增加到0.32
      label_y <- 0.5 + 0.25 * sin(mid_angle)  # 从0.25增加到0.32
      
      # 只绘制百分比
      grid.text(paste0(pie_data$percentage[i], "%"), 
                x = unit(label_x, "npc"), y = unit(label_y, "npc"),
                gp = gpar(fontfamily = "Times New Roman", fontsize = 22, col = "white", fontface = 2))
    }
    
    # 关闭饼图PDF
    dev.off()
    
    cat("已生成饼图文件:", pie_pdf_file, "\n")
  }
  
  # 返回结果
  return(list(
    file = input_file,
    Nm = coef(m.fit)*N,
    Rsqr = Rsqr,
    success = TRUE
  ))
}

# 处理两个输入文件
input_files <- c("ff.txt", "xy.txt")
results <- list()

# 循环处理每个文件
for (i in 1:length(input_files)) {
  # 捕获整个文件处理过程中可能出现的错误
  results[[i]] <- tryCatch({
    analyze_neutral_model(input_files[i])
  }, error = function(e) {
    cat("\n处理文件时出错:", input_files[i], "\n")
    cat("错误信息:", conditionMessage(e), "\n")
    return(list(file = input_files[i], Nm = NA, Rsqr = NA, success = FALSE))
  })
}

# 输出汇总结果
cat("\n分析结果汇总:\n")
cat("---------------------------\n")
for (i in 1:length(results)) {
  if(results[[i]]$success) {
    cat(sprintf("文件: %-10s  Nm: %.4f  R²: %.4f\n", 
                basename(results[[i]]$file), 
                results[[i]]$Nm, 
                results[[i]]$Rsqr))
  } else {
    cat(sprintf("文件: %-10s  [拟合失败]\n", basename(results[[i]]$file)))
  }
}
cat("---------------------------\n")
cat("\n所有文件处理完成！\n")

