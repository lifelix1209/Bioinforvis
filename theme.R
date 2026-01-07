'''

'''

# ==============================================================================
# 1. import libraries
# ==============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(grid)      
library(patchwork) 
library(ggsci)


#---- Fit for UMAP and t-SNE ----#
theme_blank <- function(add_coord = TRUE, xlen_npc = 0.15, ylen_npc = 0.15, xlab = "", ylab = "", lab_size = 12, ...) {
  args1 <- list(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.background = element_blank(),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(10, "pt"),
    plot.margin = margin(lab_size + 2, lab_size + 2, lab_size + 2, lab_size + 2, unit = "points")
  )
  args2 <- as.list(match.call())[-1]
  call.envir <- parent.frame(1)
  args2 <- lapply(args2, function(arg) {
    if (is.symbol(arg)) eval(arg, envir = call.envir) else if (is.call(arg)) eval(arg, envir = call.envir) else arg
  })
  for (n in names(args2)) args1[[n]] <- args2[[n]]
  args <- args1[names(args1) %in% formalArgs(theme)]
  out <- do.call(what = theme, args = args)
  
  if (isTRUE(add_coord)) {
    g <- grobTree(gList(
      linesGrob(x = unit(c(0, xlen_npc), "npc"), y = unit(c(0, 0), "npc"), arrow = arrow(length = unit(0.02, "npc")), gp = gpar(lwd = 2)),
      textGrob(label = xlab, x = unit(0, "npc"), y = unit(0, "npc"), vjust = 4 / 3, hjust = 0, gp = gpar(fontsize = lab_size)),
      linesGrob(x = unit(c(0, 0), "npc"), y = unit(c(0, ylen_npc), "npc"), arrow = arrow(length = unit(0.02, "npc")), gp = gpar(lwd = 2)),
      textGrob(label = ylab, x = unit(0, "npc"), y = unit(0, "npc"), vjust = -2 / 3, hjust = 0, rot = 90, gp = gpar(fontsize = lab_size))
    ))
    return(list(list(annotation_custom(g)), list(theme_scp() + out), list(coord_cartesian(clip = "off"))))
  } else {
    return(list(list(theme_hz() + out)))
  }
}


theme_hz <- function(aspect.ratio = NULL, fig_scaler = 16, ...) {
  scaler <- fig_scaler / 16
  
  # ---- defaults ----
  args_panel <- list(
    aspect.ratio = aspect.ratio,
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )
  
  args_axis <- list(
    axis.title = element_text(size = 16 * scaler),
    axis.text  = element_text(size = 14 * scaler),
    axis.line  = element_blank()
  )
  
  args_legend <- list(
    legend.position   = "right",
    legend.title      = element_text(size = 12 * scaler),
    legend.text       = element_text(size = 11 * scaler),
    legend.background = element_blank()
  )
  
  args_strip <- list(
    strip.text            = element_text(size = 12 * scaler, margin = margin(3, 3, 3, 3)),
    strip.placement       = "outside",
    strip.background      = element_rect(fill = "transparent", linetype = 0),
    strip.switch.pad.grid = unit(-1, "pt"),
    strip.switch.pad.wrap = unit(-1, "pt")
  )
  
  args_plot <- list(
    text            = element_text(size = 12 * scaler, color = "black", family = "sans"),
    plot.title      = element_text(size = 14 * scaler, colour = "black", vjust = 1, face = "bold"),
    plot.subtitle   = element_text(size = 0 * scaler, hjust = 0, margin = margin(b = 3)),
    plot.background = element_rect(fill = "white", color = "white")
  )
  
  args1 <- c(args_panel, args_axis, args_legend, args_strip, args_plot)
  args2 <- as.list(match.call())[-1]
  call.envir <- parent.frame(1)
  args2 <- lapply(args2, function(arg) {
    if (is.symbol(arg)) eval(arg, envir = call.envir)
    else if (is.call(arg)) eval(arg, envir = call.envir)
    else arg
  })
  
  for (n in names(args2)) args1[[n]] <- args2[[n]]
  args <- args1[names(args1) %in% formalArgs(theme)]
  
  do.call(theme, args)
}


head(plot_data)
p1 <- ggplot(plot_data, aes(x = orig.ident, y = nFeature_RNA)) +
  geom_violin(linewidth = 1.1, width = 1, trim = FALSE, fill = "red", alpha = 0.5, color = "red") +
  geom_boxplot(width = 0.08, linewidth = 0.8, fill = "white",
               outlier.shape = NA, color = "red") +
  theme_hz()
p1


p2 <- ggplot(plot_data, aes(x = orig.ident, y = nCount_RNA)) +
  geom_violin(linewidth = 1.1, width = 1, trim = FALSE, fill = "blue", alpha = 0.5, color = "blue") +
  geom_boxplot(width = 0.08, linewidth = 0.8, fill = "white",
               outlier.shape = NA, color = "blue") +
  theme_hz()
p2

p3 <- ggplot(plot_data, aes(x = orig.ident, y = percent.mt)) +
  geom_violin(linewidth = 1.1, width = 1, trim = FALSE, fill = "green", alpha = 0.5, color = "green") +
  geom_boxplot(width = 0.08, linewidth = 0.8, fill = "white",
               outlier.shape = NA, color = "green") +
  theme_hz()
p3


p_final <- (p1 | p2 | p3)
p_final & theme(plot.margin = margin(6, 6, 6, 6))



