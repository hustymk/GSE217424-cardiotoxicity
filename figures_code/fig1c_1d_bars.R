## ============================================================
## 17_fig1cd_bars_paper_layout.R：按原文版式复现 Fig 1C 与 1D
##
## 原文 Fig 1C 与 1D 都是灰色柱状图，C 在上、D 在下叠放：
##   C：Complete dataset，54 根灰柱，x 轴无标签，只画刻度
##   D：After projection into final drug-selective subspaces，54 根灰柱，
##      x 轴药名竖排并按药物类别着色
##   两图都在柱子上叠一条深蓝色的中位数线与四分位区间
##
## 数据：Fig 1C 用作者 Data 04（未分解 DEG），Fig 1D 用 Data 09（药物选择性 DEG）
## ============================================================
source("R/01_helpers.R")

drug_tbl <- utils::read.delim(file.path(DIR$data, "suppl", "drug_table.tsv"),
                              check.names = FALSE, stringsAsFactors = FALSE)
names(drug_tbl)[1:3] <- c("abbr", "drug", "class")
drug_tbl$drug  <- tolower(drug_tbl$drug)
drug_tbl$class <- gsub("_", " ", drug_tbl$class)

CLASS5 <- c("Kinase inhibitor"        = "#ff7f0e",
            "Monoclonal antibody"     = "#d62728",
            "Anthracycline"           = "#9467bd",
            "Cardiovascular drug"     = "#1f77b4",
            "Non cardiovascular drug" = "#17becf")
CLASS_ORDER <- names(CLASS5)
GREY <- "grey62"
BLUE <- "#123c8c"

b  <- readRDS(file.path(DIR$cache, "hi_fig1b.rds"))    # 含 ord、f1_before
de <- readRDS(file.path(DIR$cache, "hi_fig1de.rds"))   # 含 cmp

## ---------- C：改用作者完整 16345 基因矩阵（与原文同源）----------
fm <- readRDS(file.path(DIR$cache, "fullmatrix.rds"))
cm_f <- cor(fm$M); cm_f[is.na(cm_f)] <- 0
hc_f <- stats::hclust(stats::as.dist((1 - cm_f) / 2), method = "average")
f1_full <- f1_per_drug(hc_f, fm$drug)
## 作者源码：Entity_factor 用 drugs_in_drugClass_specific_order，
## 即 C 与 D 共用同一套顺序（先按官方五类，组内按药名字母序）
c_df <- f1_full %>%
  left_join(drug_tbl %>% select(abbr, class), by = c("drug" = "abbr")) %>%
  filter(!is.na(f1)) %>%
  mutate(class = factor(class, levels = CLASS_ORDER),
         name  = as.character(drug)) %>%
  arrange(class, name)

## ---------- D：按药物类别分组，组内按药名排序 ----------
d_df <- de$cmp %>%
  mutate(class = factor(class, levels = CLASS_ORDER),
         abbr_chr = as.character(abbr)) %>%
  arrange(class, abbr_chr)

med_iqr <- function(v) c(median = median(v), q1 = quantile(v, .25), q3 = quantile(v, .75))

## ---------- 画：上下两个面板 ----------
png(file.path(DIR$fig, "fig1c_1d_f1_bars.png"),
    width = 13 * FIG_DPI, height = 8 * FIG_DPI, res = FIG_DPI)
layout(matrix(1:2, ncol = 1), heights = c(1, 1.15))

## 面板 C
par(mar = c(2.2, 5.5, 3, 1))
bp <- barplot(c_df$f1, col = GREY, border = NA, ylim = c(0, 1.05),
              names.arg = rep("", nrow(c_df)), space = 0.18,
              axes = FALSE, xpd = FALSE)
axis(2, las = 1, cex.axis = .8); box(col = "grey80")
mtext("F1 score", side = 2, line = 3.4, cex = .85)
title(main = "C   Complete dataset", adj = 0, cex.main = 1, font.main = 1, line = 1.6)
mtext(sprintf("F1 中位数 %.3f，n = %d 种药", median(c_df$f1), nrow(c_df)),
      side = 3, adj = 0, line = .4, cex = .68)
s <- med_iqr(c_df$f1)
rect(min(bp) - .4, s["q1"], max(bp) + .4, s["q3"], col = adjustcolor(BLUE, .25), border = NA)
segments(min(bp) - .4, s["median"], max(bp) + .4, s["median"], col = BLUE, lwd = 3)

## 面板 D
par(mar = c(6.5, 5.5, 3, 1))
bp2 <- barplot(d_df$after, col = GREY, border = NA, ylim = c(0, 1.05),
               names.arg = rep("", nrow(d_df)), axes = FALSE,
               cex.names = .62, xpd = FALSE)
axis(2, las = 1, cex.axis = .8); box(col = "grey80")
mtext("F1 score", side = 2, line = 3.4, cex = .85)
title(main = "D   After projection into final drug-selective subspaces",
      adj = 0, cex.main = 1, font.main = 1, line = 1.6)
mtext(sprintf("F1 中位数 %.3f，n = %d 种药", median(d_df$after), nrow(d_df)),
      side = 3, adj = 0, line = .4, cex = .68)
s2 <- med_iqr(d_df$after)
rect(min(bp2) - .4, s2["q1"], max(bp2) + .4, s2["q3"], col = adjustcolor(BLUE, .25), border = NA)
segments(min(bp2) - .4, s2["median"], max(bp2) + .4, s2["median"], col = BLUE, lwd = 3)

## x 轴药名按类别着色（原文做法）
cols <- unname(CLASS5[as.character(d_df$class)])
text(bp2, par("usr")[3] - .035, labels = d_df$abbr, srt = 90, adj = c(1, .5),
     xpd = NA, cex = .62, col = cols)

## 图例
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("bottom", horiz = TRUE, bty = "n", cex = .72,
       legend = c("Kinase inhibitor", "Monoclonal antibody", "Anthracycline",
                  "Cardiovascular drug", "Non cardiovascular drug"),
       fill = unname(CLASS5), border = NA)
dev.off()

message(sprintf("[1C/1D bars] 中位数 %.3f 到 %.3f，上升 %d / %d",
                median(c_df$f1), median(d_df$after),
                sum(d_df$after > d_df$before), nrow(d_df)))
readr::write_tsv(bind_rows(
  tibble(panel = "C", metric = "median_F1", value = median(c_df$f1)),
  tibble(panel = "D", metric = "median_F1", value = median(d_df$after)),
  tibble(panel = "D", metric = "n_increased", value = sum(d_df$after > d_df$before)),
  tibble(panel = "D", metric = "n_total", value = nrow(d_df))
), file.path(DIR$tab, "fig1cd_bars_key_numbers.tsv"))

## ---------- 报告用的汇总表 ----------
## 只放 C / D 两个面板的数字
readr::write_tsv(
  tibble(metric = c("median_F1_before", "median_F1_after",
                    "n_drugs_increased", "n_drugs_total"),
         value = c(median(c_df$f1), median(d_df$after),
                   sum(d_df$after > d_df$before), nrow(d_df))),
  file.path(DIR$tab, "fig1c_1d_key_numbers_v2.tsv"))
message("[1C/1D bars done] -> figures/fig1c_1d_f1_bars.png")
