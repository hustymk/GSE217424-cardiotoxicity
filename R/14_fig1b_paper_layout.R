## ============================================================
## 14_fig1b_paper_layout.R：按原文版式复现 Fig 1B
##
## 参照原图（ref/figs/41467_2024_52145_Fig1_HTML.png）的实际做法：
##   上：层次聚类树，**被压得很扁**（纵轴只占约 1/3 面板高）
##   下：一条厚色带，**每个叶子一根竖条，54 种药各一种颜色**，
##       药名竖排写在条内（只在空间够的地方写）
##   色带左侧标注 "266 条细胞系x药物组合" 与 "54 drugs (54 colors)"
##
## 与"第一版"（纵轴拉满、只标 5 类颜色）形成错误→改进的对比。
## ============================================================
source("R/01_helpers.R")

fm <- readRDS(file.path(DIR$cache, "fullmatrix.rds"))
M  <- fm$M
meta <- tibble::tibble(drug = fm$drug, cell_line = fm$cell)
cm <- cor(M); cm[is.na(cm)] <- 0
hc  <- stats::hclust(stats::as.dist((1 - cm) / 2), method = "average")   # 作者源码用平均连锁
ord <- hc$order
n <- length(ord)

## ---------- 药物颜色：照作者 SVD_colors.R 的手工配色 ----------
## 作者按机制成组指定颜色（蒽环类紫系、心血管药蓝系、TKI 绿橙系），
## 文件在 ref/github/SVD_colors.R；缺失的药退回灰色
drugs <- sort(unique(meta$drug))
pal <- stats::setNames(rep("grey60", length(drugs)), drugs)
frac <- "ref/github/SVD_colors.R"
if (file.exists(frac)) {
  src <- readLines(frac, warn = FALSE)
  pair <- paste0('drug_highlight_colors = c\\(drug_highlight_colors,"([^"]+)"\\);',
                 'drugs_for_highlight_colors = c\\(drugs_for_highlight_colors,"([^"]+)"\\);')
  hit <- unlist(regmatches(src, gregexpr(pair, src)))
  if (length(hit)) {
    pal[sub(pair, "\\2", hit)] <- sub(pair, "\\1", hit)
  }
}
cols <- unname(pal[meta$drug[ord]])
cols[is.na(cols)] <- "grey60"
if (any(cols == "grey60")) {
  warning("有药没在作者色表里，已用灰色：",
          paste(meta$drug[ord][cols == "grey60"], collapse = ", "))
}

## ---------- 画：上树 + 下色带 ----------
png(file.path(DIR$fig, "fig1b_paper_layout.png"),
    width = 16 * FIG_DPI, height = 7 * FIG_DPI, res = FIG_DPI)
layout(matrix(1:2, ncol = 1), heights = c(1, 1.30))

## 上：树（压扁：只让它占上面那块，不拉满整页）
par(mar = c(0.3, 9, 2.6, 0.5))
plot(hc, labels = FALSE, hang = -1, main = "", sub = "",
     xlab = "", ylab = "", axes = TRUE, yaxt = "n")
axis(2, las = 1, cex.axis = 0.75)
mtext("Fig 1B | 266 个「细胞系 × 药物」组合的层次聚类（hclust 平均连锁）",
      side = 3, adj = 0, line = 1.5, cex = 0.85, font = 2)
mtext("266\ndrug/cell line\ncombinations", side = 2, line = 5.5,
      cex = 0.62, las = 0)

## 下：色带（每叶一根竖条，54 色）
par(mar = c(0.3, 9, 0.3, 0.5))
plot.new(); plot.window(xlim = c(0.5, n + 0.5), ylim = c(0, 1))
rect(seq_len(n) - 0.5, 0, seq_len(n) + 0.5, 1,
     col = cols, border = NA)

## 药名竖排：只在同一药物的"连续段"长度 >= 3 处写，避免糊
r <- rle(meta$drug[ord])
ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
run_ids  <- ave(seq_along(r$values), r$values, FUN = seq_along)
n_runs   <- as.integer(table(r$values)[r$values])
labs_run <- ifelse(n_runs > 1, paste0(r$values, " ", run_ids), r$values)
lab_idx  <- which(r$lengths >= 3)
for (k in lab_idx) {
  xc <- (starts[k] + ends[k]) / 2
  text(xc, 0.5, labels = labs_run[k], srt = 90, cex = 0.62,
       col = "white", font = 2)
}
mtext("54 drugs\n(54 colors)", side = 2, line = 5.5, cex = 0.62, las = 0)
par(xpd = NA)
mtext(sprintf("共 %d 根竖条 = 266 个样本；同色相邻 = 同一种药的样本聚在一起",
              n), side = 1, line = -0.6, cex = 0.6)
dev.off()

## ---------- 定量：同色相邻比例（与图对应的一句话）----------
adj <- mean(meta$drug[ord][-1] == meta$drug[ord][-n])
message(sprintf("[1B-paper] 相邻竖条同色比例 = %.1f%%（随机期望 %.1f%%）",
                100 * adj, 100 / length(drugs)))
readr::write_tsv(tibble(metric = c("n_leaves", "n_drugs", "adjacent_same_drug",
                                   "expected_random"),
                        value = c(n, length(drugs), adj, 1 / length(drugs))),
                 file.path(DIR$tab, "fig1b_paper_layout.tsv"))
message("[1B-paper done] -> figures/fig1b_paper_layout.png")
