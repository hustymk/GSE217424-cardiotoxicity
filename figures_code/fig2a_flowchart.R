## ============================================================
## 16_fig2a_flowchart.R：按原文版式复刻 Fig 2A（流程图）
##
## 原文 Fig 2A 的元素（见 ref/figs/..._Fig2_HTML.png）：
##   橙色圆角框 + 深蓝箭头 + 深蓝粗体说明文字
##   左侧一列：upregulated 通路；右侧一列：downregulated 通路
##   两列汇聚 -> "Combine and rank by decreasing AUC" -> 底部结果框
##   右侧一条紫色回环箭头 + 竖排紫字 "Repeat for non-cardiotoxic drugs"
##
## 该 panel 是示意图，数据来自上游分析，不涉及绘图数据；
## 这里用 grid 精确复刻画布布局。
## ============================================================
source("R/00_setup.R")
library(grid)

BOX_FILL <- "#f8dcae"
BOX_LINE <- "#e08a2e"
BLUE     <- "#1f4e8c"
PURPLE   <- "#7b2d8e"

## 圆角框（R >= 4.1 有 roundrectGrob，否则退回矩形）
box_grob <- function(x, y, w, h, label, cex = 0.78) {
  g <- if (exists("roundrectGrob")) {
    roundrectGrob(x = unit(x, "native"), y = unit(y, "native"),
                  width = unit(w, "native"), height = unit(h, "native"),
                  r = unit(0.9, "mm"),
                  gp = gpar(fill = BOX_FILL, col = BOX_LINE, lwd = 1.2))
  } else {
    rectGrob(x = unit(x, "native"), y = unit(y, "native"),
             width = unit(w, "native"), height = unit(h, "native"),
             gp = gpar(fill = BOX_FILL, col = BOX_LINE, lwd = 1.2))
  }
  gTree(children = gList(
    g,
    textGrob(label, x = unit(x, "native"), y = unit(y, "native"),
             gp = gpar(cex = cex, col = "black"))
  ))
}

arrow_grob <- function(x0, y0, x1, y1, col = BLUE, lwd = 1.6) {
  segmentsGrob(x0 = unit(x0, "native"), y0 = unit(y0, "native"),
               x1 = unit(x1, "native"), y1 = unit(y1, "native"),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
               gp = gpar(col = col, lwd = lwd, fill = col))
}

png(file.path(DIR$fig, "fig2a_flowchart.png"),
    width = 11 * FIG_DPI, height = 5.2 * FIG_DPI, res = FIG_DPI)
grid.newpage()
pushViewport(viewport(xscale = c(0, 100), yscale = c(0, 100), name = "vp"))

## ---------- 框 ----------
b_tl <- box_grob(28, 91, 36, 9, "Ranked list of upregulated\nSCPs for each TKI-treated sample")
b_tr <- box_grob(72, 91, 36, 9, "Ranked list of downregulated\nSCPs for each TKI-treated sample")
b_ml <- box_grob(28, 63, 36, 9, "SCPs upregulated at higher\nranks by cardiotoxic TKIs")
b_mr <- box_grob(72, 63, 36, 9, "SCPs downregulated at higher\nranks by cardiotoxic TKIs")
b_bt <- box_grob(50, 21, 46, 9,
                 "Top 25 combined up- or downregulated\nlevel-3 SCPs for cardiotoxic TKIs")

## ---------- 箭头 ----------
a1 <- arrow_grob(28, 86.5, 28, 67.8)
a2 <- arrow_grob(72, 86.5, 72, 67.8)

## 两列汇聚
merge <- gTree(children = gList(
  segmentsGrob(unit(c(28, 28), "native"), unit(c(58.5, 50), "native"),
               unit(c(28, 50), "native"), unit(c(50, 50), "native"),
               gp = gpar(col = "black", lwd = 1.3)),
  segmentsGrob(unit(c(72, 72), "native"), unit(c(58.5, 50), "native"),
               unit(c(72, 50), "native"), unit(c(50, 50), "native"),
               gp = gpar(col = "black", lwd = 1.3)),
  arrow_grob(50, 50, 50, 25.8, col = "black", lwd = 1.3)
))

## 紫色回环：底部框右 -> 右外 -> 上方 -> 指回中右框
loop <- gTree(children = gList(
  segmentsGrob(unit(c(73, 95), "native"), unit(c(21, 21), "native"),
               unit(c(95, 95), "native"), unit(c(21, 63), "native"),
               gp = gpar(col = PURPLE, lwd = 1.5)),
  segmentsGrob(unit(95, "native"), unit(63, "native"),
               unit(90.4, "native"), unit(63, "native"),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
               gp = gpar(col = PURPLE, lwd = 1.5, fill = PURPLE))
))

## ---------- 文字 ----------
txt <- gTree(children = gList(
  textGrob("F1 score and AUC statistics", x = unit(50, "native"), y = unit(77, "native"),
           gp = gpar(cex = 0.82, col = BLUE, fontface = "bold")),
  textGrob("Combine and rank\nby decreasing AUC", x = unit(57, "native"), y = unit(45, "native"),
           gp = gpar(cex = 0.82, col = BLUE, fontface = "bold")),
  textGrob("Repeat for\nnon-cardiotoxic drugs", x = unit(98, "native"), y = unit(42, "native"),
           rot = 90, gp = gpar(cex = 0.78, col = PURPLE, fontface = "bold")),
  textGrob("A", x = unit(1, "native"), y = unit(99, "native"), just = c("left", "top"),
           gp = gpar(cex = 1.25, fontface = "bold"))
))

grid.draw(b_tl); grid.draw(b_tr); grid.draw(b_ml); grid.draw(b_mr); grid.draw(b_bt)
grid.draw(a1);   grid.draw(a2);   grid.draw(merge); grid.draw(loop); grid.draw(txt)

popViewport(); dev.off()
message("[2A done] -> figures/fig2a_flowchart.png")
