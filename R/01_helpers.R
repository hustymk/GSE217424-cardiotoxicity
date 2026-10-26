## ============================================================
## 01_helpers.R：三张主讲图共用的函数
## ============================================================
source("R/00_setup.R")

## ---------- 药物类别配色（对齐论文 Fig 1D 图例）----------
## 橙：小分子激酶抑制剂  红：抗激酶单抗  紫：蒽环类
## 蓝：心脏作用药        青：非心脏作用药
ANNO <- list(
  KI_mab = c("BEV", "TRA", "RTX", "TRS"),                       # 抗激酶单抗
  anthra = c("DOX", "DAU", "IDA", "EPI"),                       # 蒽环类
  KI_sm  = c("PAZ","DAB","AXI","BOS","CAB","CER","CRI","CTX","DAS","DIC",
             "ERL","GEF","IMA","LAP","NIL","PIO","PON","REG","RUX","SAX",
             "SOR","SUN","TOF","VAN","VEM","AFA","CAR","DOB","FLE","MIL","ROS"),
  cardiac   = c("AMI","AZA","DEC","DEL","EDN","ISO","OLM","PHP","PRE","TNF",
                "VRP","CYC","EST","IGF","BOR")
)
CLASS_COL <- c(KI_sm = "#ff7f0e", KI_mab = "#d62728", anthra = "#9467bd",
               cardiac = "#1f77b4", Other = "#17becf")
CLASS_LAB <- c(KI_sm = "小分子激酶抑制剂", KI_mab = "抗激酶单抗",
               anthra = "蒽环类", cardiac = "心脏作用药", Other = "非心脏作用药")

drug_class <- function(dr) {
  out <- rep("Other", length(dr))
  for (k in c("KI_mab", "anthra", "KI_sm", "cardiac"))
    out[dr %in% ANNO[[k]]] <- k
  out
}

## ---------- 相关 -> 聚类（论文距离：(1 - cor) / 2）----------
clust_order <- function(M, method = "complete") {
  cm <- cor(M, use = "pairwise.complete.obs")
  cm[is.na(cm)] <- 0
  hc <- stats::hclust(stats::as.dist((1 - cm) / 2), method = method)
  list(cor = cm, hc = hc, order = hc$order)
}

## ---------- F1 聚类效率（按作者源码移植）----------
## 来源：SVD_singleDrug_validations_calculate_and_add.R 第 24 至 79 行
##   clusters = partition_leaves(col_dend)  即枚举树上每一个内部节点
##   precision = 簇内同药样本数 / 簇大小
##   recall    = 簇内同药样本数 / 该药样本总数
##   F1        = 调和平均，对每个药取所有有效节点的最大值
## 有效条件：簇内至少 2 个样本、至少 2 种药；该药在簇内 >= 2 个样本（或该药总共只有 1 个样本）
partition_leaves_hc <- function(hc) {
  n <- length(hc$order)
  members <- vector("list", n - 1L)
  for (i in seq_len(n - 1L)) {
    a <- hc$merge[i, 1]; b <- hc$merge[i, 2]
    la <- if (a < 0) -a else members[[a]]
    lb <- if (b < 0) -b else members[[b]]
    members[[i]] <- c(la, lb)
  }
  members
}

f1_per_drug <- function(hc, labels, ks = NULL) {
  total <- table(labels)
  best  <- stats::setNames(rep(0, length(total)), names(total))
  for (m in partition_leaves_hc(hc)) {
    cs <- length(m)
    if (cs < 2) next                        # 簇内至少 2 个样本
    tab <- table(labels[m])
    if (length(tab) >= cs) next             # 簇内至少 2 种药
    for (dr in names(tab)) {
      cnt <- as.integer(tab[[dr]])
      if (!(cnt >= 2 || as.integer(total[[dr]]) == 1)) next
      p  <- cnt / cs
      r  <- cnt / as.integer(total[[dr]])
      f1 <- 2 * p * r / (p + r)
      if (f1 > best[[dr]]) best[[dr]] <- f1
    }
  }
  tibble::tibble(drug = names(best), f1 = as.numeric(best))
}

## ---------- 上下拼接（base grid，零依赖）----------
combine_v <- function(p_top, p_main, file, w, h, top_frac = 0.10) {
  png(file, width = w * FIG_DPI, height = h * FIG_DPI, res = FIG_DPI)
  on.exit(dev.off(), add = TRUE)
  grid::grid.newpage()
  lay <- grid::grid.layout(2, 1,
          heights = grid::unit(c(top_frac, 1 - top_frac), "null"))
  grid::pushViewport(grid::viewport(layout = lay))
  print(p_top,  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_main, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::popViewport()
}

## ---------- 左右拼接 ----------
combine_h <- function(p_left, p_right, file, w, h, left_frac = 0.55) {
  png(file, width = w * FIG_DPI, height = h * FIG_DPI, res = FIG_DPI)
  on.exit(dev.off(), add = TRUE)
  grid::grid.newpage()
  lay <- grid::grid.layout(1, 2,
          widths = grid::unit(c(left_frac, 1 - left_frac), "null"))
  grid::pushViewport(grid::viewport(layout = lay))
  print(p_left,  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::popViewport()
}

## ---------- 相关热图（Fig 1B 用）----------
cor_heatmap <- function(M, ord, meta_ord, title, subtitle) {
  hm <- tibble::tibble(
    r = rep(seq_along(ord), times = length(ord)),
    c = rep(seq_along(ord), each  = length(ord)),
    v = as.vector(M[ord, ord])
  )
  p_hm <- ggplot(hm, aes(c, r, fill = v)) +
    geom_raster() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, name = "Pearson r") +
    coord_fixed() +
    labs(x = NULL, y = NULL) +
    theme(axis.text = element_blank(), panel.grid = element_blank())

  strip <- tibble::tibble(i = seq_along(ord),
                          cls = drug_class(meta_ord$drug[ord]))
  p_strip <- ggplot(strip, aes(i, 1, fill = cls)) +
    geom_raster() +
    scale_fill_manual(values = CLASS_COL, labels = CLASS_LAB, name = NULL) +
    labs(title = title, subtitle = subtitle) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.key.width = unit(5, "mm"),
          legend.text = element_text(size = 7),
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8))
  list(strip = p_strip, heatmap = p_hm)
}
