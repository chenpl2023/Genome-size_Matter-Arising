# 大图：2×2 布局，整体宽:高 = 3:1
#  (1,1) 热图1 （fengbase：行Sample EC降序，列Genome max→sum）
#  (1,2) 散点图2 （fengbase：((Genome_size*Eff_cov)/Total_bp)*100 按 Sample_file 求和）
#  (2,1) 热图3 （All_results：行Sample EC降序，列accession max→sum，Tax注释）
#  (2,2) 散点图4 （All_results：上述公式分组求和）
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(ggplot2)
  library(grid)
  library(data.table)          # fread 稳健解析含内嵌换行的 CSV
})

ht_opt(message = FALSE)   # 关闭 use_raster / magick 等例行提示

# ================= 公共工具：离散三档 =================
disc3 <- function(m) {
  matrix(ifelse(m == 0, "0", ifelse(m > 1, "high", "low")),
         nrow = nrow(m), dimnames = dimnames(m))
}

# ================= 热图1：Eff_cov（fengbase） =================
build_ht1 <- function() {
  hc <- read.csv("HEAT_MAP_col.csv", check.names = FALSE)
  gs <- read.csv("Genome_Size.csv", check.names = FALSE)
  fb <- read.csv("fengbase_result.csv", row.names = NULL, check.names = FALSE)
  genomes <- unique(hc$Genome_file)
  samples <- as.character(gs$Sample_ID)
  ec      <- setNames(gs[["EC (dS/m)"]], samples)

  fb$Sample_file <- as.character(fb$Sample_file); fb$Genome_file <- as.character(fb$Genome_file)
  fb$Eff_cov <- suppressWarnings(as.numeric(as.character(fb$Eff_cov)))

  m <- matrix(0, length(genomes), length(samples), dimnames = list(genomes, samples))
  r <- match(fb$Genome_file, genomes); c <- match(fb$Sample_file, samples)
  ok <- !is.na(r) & !is.na(c) & !is.na(fb$Eff_cov)
  m[cbind(r[ok], c[ok])] <- fb$Eff_cov[ok]

  m <- t(m)                                        # 行=样品，列=基因组
  # 行：EC 降序
  m <- m[order(-ec[rownames(m)]), , drop = FALSE]
  # 列：max → sum 降序
  m <- m[, order(-apply(m, 2, max), -apply(m, 2, sum)), drop = FALSE]

  k_vec <- setNames(hc[["#k"]], hc$Genome_file)
  k_arr <- k_vec[colnames(m)]
  ha <- HeatmapAnnotation(K = k_arr,
    col = list(K = c("Bacteria" = "#8B4513", "Archaea" = "#F5F5DC")),
    show_annotation_name = FALSE,
    annotation_legend_param = list(title = "Kingdom",
      title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)))

  Heatmap(disc3(m), name = "coverage",
    col = c("0" = "#F0F0F0", "low" = "#1F77B4", "high" = "#D62728"),
    cluster_rows = FALSE, row_order = rownames(m),
    cluster_columns = FALSE, column_order = colnames(m),
    show_row_names = FALSE, show_column_names = FALSE,
    column_title = "Available reference genome for taxa in Dong et al's Supplementary-Data2 (n = 356)",
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 15, fontface = "bold"),
    row_title = "Samples",
    row_title_gp = gpar(fontsize = 15, fontface = "bold"),
    top_annotation = ha,
    heatmap_legend_param = list(title = "Coverage",
      title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12),
      at = c("0", "low", "high"), labels = c("0", "0-1", ">1")))
}

# ================= 热图3：Eff_cov（All_results，accession列） =================
build_ht3 <- function() {
  d <- fread("All_results_tax_gtdb.csv", fill = TRUE)
  samp <- unique(as.character(d$Sample_file)); acc <- unique(as.character(d$accession))
  m <- matrix(0, length(samp), length(acc), dimnames = list(samp, acc))
  r <- match(as.character(d$Sample_file), samp); c <- match(as.character(d$accession), acc)
  ef <- suppressWarnings(as.numeric(as.character(d$Eff_cov)))
  ok <- !is.na(r) & !is.na(c) & !is.na(ef)
  m[cbind(r[ok], c[ok])] <- ef[ok]

  ecmap <- setNames(as.numeric(as.character(d$EC)), as.character(d$Sample_file))
  m <- m[order(-ecmap[rownames(m)]), , drop = FALSE]
  m <- m[, order(-apply(m, 2, max), -apply(m, 2, sum)), drop = FALSE]

  tax_map <- setNames(as.character(d$Tax), as.character(d$accession))
  tax_arr <- tax_map[colnames(m)]
  col_pal <- c(Archaea = "#F5F5DC", Bacteria = "#8B4513")
  ha <- HeatmapAnnotation(Taxon = tax_arr, col = list(Taxon = col_pal),
    show_annotation_name = FALSE,
    annotation_legend_param = list(title = "Kingdom",
      title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13)))

  Heatmap(disc3(m), name = "coverage",
    col = c("0" = "#F0F0F0", "low" = "#1F77B4", "high" = "#D62728"),
    cluster_rows = FALSE, row_order = rownames(m),
    cluster_columns = FALSE, column_order = colnames(m),
    show_row_names = FALSE, show_column_names = FALSE,
    raster_quality = 5,                                  # PDF 中光栅化更清晰
    column_title = "GTDB-derived reference genomes (n = 2486)",
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 15, fontface = "bold"),
    row_title = "Samples",
    row_title_gp = gpar(fontsize = 15, fontface = "bold"),
    top_annotation = ha,
    heatmap_legend_param = list(title = "Coverage",
      title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13),
      at = c("0", "low", "high"), labels = c("0", "0-1", ">1")))
}

# ================= 散点数据：((Genome_size*Eff_cov)/Total_bp)*100 按 Sample_file 求和 =================
scatter_data <- function(d, gs_col) {
  gs <- suppressWarnings(as.numeric(as.character(d[[gs_col]])))
  ef <- suppressWarnings(as.numeric(as.character(d$Eff_cov)))
  tb <- suppressWarnings(as.numeric(as.character(d$Total_bp)))
  ec <- suppressWarnings(as.numeric(as.character(d$EC)))
  sm <- as.character(d$Sample_file)
  val <- (gs * ef / tb) * 100
  dat <- data.frame(Sample_file = sm, val = val, EC = ec, stringsAsFactors = FALSE)
  dat <- dat[is.finite(dat$val), ]
  agg <- aggregate(val ~ Sample_file + EC, dat, sum)
  agg <- agg[order(agg$EC), ]                        # EC 升序
  agg$Sample_file <- factor(agg$Sample_file, levels = agg$Sample_file)
  agg
}

scatter_plot <- function(agg) {
  ggplot(agg, aes(x = Sample_file, y = val)) +
    geom_point(color = "#1F77B4", shape = 1, size = 2.4) +
    theme_bw(base_size = 12) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10),
                       name = "Metagenomic reads mapped (%)") +
    labs(x = "Samples") +
    theme(axis.ticks.x = element_blank(),
          axis.text.x  = element_blank(),            # 隐藏刻度标签
          axis.title.x = element_text(size = 15, face = "bold"),   # 保留标题
          axis.title.y = element_text(size = 15, face = "bold"),
          axis.text.y  = element_text(size = 13),
          plot.title   = element_blank())
}

# ================= 一并输出：2×2 网格（图1=图3、图2=图4 同尺寸） =================
ht1 <- build_ht1()
ht3 <- build_ht3()
agg2 <- scatter_data(fread("fengbase_result.csv", fill = TRUE), "Genome_size")
agg4 <- scatter_data(fread("All_results_tax_gtdb.csv", fill = TRUE), "genome_size")
p2 <- scatter_plot(agg2)
p4 <- scatter_plot(agg4)

b <- ggplotGrob(p2)
d <- ggplotGrob(p4)

pdf("combined_figure.pdf", width = 16, height = 8)
grid.newpage()
# 左列热图宽 0.6、右列散点宽 0.2；两行等行高 → 图1=图3、图2=图4 同宽同高
pushViewport(viewport(layout = grid.layout(2, 2,
  widths = unit(c(0.6, 0.2), "null"),
  heights = unit(c(1, 1), "null"))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht1, merge_legend = TRUE, newpage = FALSE)
grid.text("a", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(b)
grid.text("b", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht3, merge_legend = TRUE, newpage = FALSE)
grid.text("c", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(d)
grid.text("d", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 18, fontface = "bold"))
popViewport()

dev.off()
cat("输出：combined_figure.pdf (16 × 8 in)\n")