#' @title 创建带散点的组合条形图
#' @description 该函数使用ggplot2和ggsignif创建一个组合图表，包括条形图（表示均值）、
#'              误差线（表示标准误）、散点（表示原始数据）和显著性标记。
#'
#' @param data 一个数据框（data.frame）。数据必须是“长格式”，即：
#'             - 包含一个表示分组的列（例如，'DMSO', 'Drug A'）。
#'             - 包含一个表示测量值的数值列（例如，150.5, 98.2）。
#' @param x_col 一个字符串，指定`data`中用作X轴分组的列名。
#' @param y_col 一个字符串，指定`data`中用作Y轴数值的列名。
#' @param group_order 一个字符向量，用于指定X轴上组的显示顺序。如果为NULL（默认），
#'                    将使用因子水平的默认顺序或字母顺序。
#' @param colors_vec 一个命名的字符向量，用于指定每个组的颜色。
#'                   例如: `c("Group1" = "black", "Group2" = "red")`。
#' @param fills_vec 一个命名的字符向量，用于指定每个组的条形图填充色。
#' @param shapes_vec 一个命名的数值向量，用于指定每个组的散点形状。
#' @param y_label 字符串，Y轴的标题。
#' @param comparisons_list 一个列表，其中每个元素都是一个包含两个组名的向量，用于`ggsignif`进行比较。
#'                         例如: `list(c("Group1", "Group2"), c("Group1", "Group3"))`。
#' @param annotations_vec 一个字符向量，包含要显示在显著性条上的文本（p值）。
#'                        其顺序必须与`comparisons_list`中的比较一一对应。
#' @param y_positions_vec 一个数值向量，指定每个显著性条的Y轴位置。
#'                        其顺序也必须与`comparisons_list`一一对应。
#' @param y_limits 一个数值向量 `c(min, max)`，用于设置Y轴的范围。
#'
#' @return 返回一个ggplot对象，可以直接打印显示。
#'
#' @examples
#' # plot <- create_bar_scatter_plot(...)
#' # print(plot)

create_bar_scatter_plot <- function(data, 
                                    x_col, 
                                    y_col,
                                    group_order = NULL,
                                    colors_vec,
                                    fills_vec,
                                    shapes_vec,
                                    y_label = "Value",
                                    comparisons_list = NULL,
                                    annotations_vec = NULL,
                                    y_positions_vec = NULL,
                                    y_limits = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请安装 'ggplot2' 包: install.packages('ggplot2')")
  }
  if (!requireNamespace("ggsignif", quietly = TRUE)) {
    stop("请安装 'ggsignif' 包: install.packages('ggsignif')")
  }
  if (!is.data.frame(data)) {
    stop("输入 'data' 必须是 data.frame 类型。")
  }
  if (!is.null(group_order)) {
    data[[x_col]] <- factor(data[[x_col]], levels = group_order)
  }
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
    
    ggplot2::stat_summary(
      fun = mean,
      geom = "bar",
      ggplot2::aes(fill = .data[[x_col]], color = .data[[x_col]]),
      width = 0.7,
      alpha = 0.5
    ) +
    
    ggplot2::stat_summary(
      fun.data = ggplot2::mean_se,
      geom = "errorbar",
      ggplot2::aes(color = .data[[x_col]]),
      width = 0.3,
      linewidth = 0.8
    ) +
    
    ggplot2::geom_jitter(
      ggplot2::aes(color = .data[[x_col]], shape = .data[[x_col]]),
      size = 3.5,
      width = 0.2
    ) +
    
    ggplot2::scale_color_manual(values = colors_vec) +
    ggplot2::scale_fill_manual(values = fills_vec) +
    ggplot2::scale_shape_manual(values = shapes_vec) +
    
    ggplot2::scale_y_continuous(
      limits = y_limits,
      expand = ggplot2::expansion(mult = c(0, 0.05)) # Y轴从0开始，顶部留5%空间
    ) +
    
    ggplot2::labs(
      y = y_label,
      x = NULL
    ) +
    
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 16, color = "black"),
      axis.text.y = ggplot2::element_text(size = 16, color = "black"),
      axis.title.y = ggplot2::element_text(size = 20, face = "plain"),
      axis.line = ggplot2::element_line(linewidth = 0.8),
      axis.ticks = ggplot2::element_line(linewidth = 0.8, color = "black"),
      axis.ticks.length = ggplot2::unit(0.2, "cm"),
      legend.position = "none"
    )
  
  if (!is.null(comparisons_list) && !is.null(annotations_vec) && !is.null(y_positions_vec)) {
    p <- p + ggsignif::geom_signif(
      comparisons = comparisons_list,
      annotations = annotations_vec,
      y_position = y_positions_vec,
      tip_length = 0.01,
      vjust = -0.2,
      textsize = 6,
      linewidth = 0.8
    )
  }
  
  return(p)
}



library(ggplot2)
library(ggsignif)

set.seed(123)
original_df <- data.frame(
  treatment = factor(
    rep(c("DMSO", "Infigratinib-24h", "Infigratinib-48h"), each = 8)
  ),
  lipid_droplets = c(
    rnorm(8, mean = 165, sd = 15),
    rnorm(8, mean = 110, sd = 25),
    rnorm(8, mean = 90, sd = 20)
  )
)

original_df <- data.frame(
  group = factor(
    rep(c("Parents", "F1"), each = 10)
  ),
  gene_expression = c(
    rnorm(10, mean = 10, sd = 3),
    rnorm(10, mean = 25, sd = 6)
  )
)

library(ggplot2)
p <- ggplot(data = original_df, aes(x = .data$group, y = .data$gene_expression)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    aes(fill = .data$group, colour = .data$group),
    width = 0.3,
    alpha = 0.5
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    mapping = aes(color = group),
    width = 0.1,
    linewidth = 0.4
  ) +
  geom_jitter(
    aes(color = group, shape = group),
    width = 0.1,
    size = 0.35
  ) +
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  scale_shape_manual(values = c(1, 2)) +
  scale_y_continuous(limits = c(NA, 30),
                     expand = expansion(mult = c(0, 0.05))) +
  
  labs(title = "The test figure",
       y = "Gene expression",
       x = "Group") +
  
  theme_classic()+ theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 16),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 20, face = "plain"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "none"
  )


p

library(ggpubr)


ggdensity(data = original_df$gene_expression)


my_order <- c("DMSO", "Infigratinib-24h", "Infigratinib-48h")

# 颜色、填充和形状的映射
my_colors <- c("DMSO" = "black", "Infigratinib-24h" = "#F8766D", "Infigratinib-48h" = "red")
my_fills <- c("DMSO" = "white", "Infigratinib-24h" = "#F8766D", "Infigratinib-48h" = "red")
my_shapes <- c("DMSO" = 16, "Infigratinib-24h" = 15, "Infigratinib-48h" = 17)

my_comparisons <- list(c("DMSO", "Infigratinib-24h"), c("DMSO", "Infigratinib-48h"))
my_annotations <- c("p=0.0017", "p<0.0001") #可以自己算p值然后自己填
my_y_positions <- c(195, 215) #y轴范围



plot1 <- create_bar_scatter_plot(
  data = original_df,
  x_col = "treatment",
  y_col = "lipid_droplets",
  group_order = my_order,
  colors_vec = my_colors,
  fills_vec = my_fills,
  shapes_vec = my_shapes,
  y_label = "Number of lipid droplets",
  comparisons_list = my_comparisons,
  annotations_vec = my_annotations,
  y_positions_vec = my_y_positions,
  y_limits = c(0, 230)
)


print(plot1)






