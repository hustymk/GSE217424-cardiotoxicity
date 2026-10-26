## ============================================================
## 15_fig2b_paper_layout.R：按原文版式复现 Fig 2B
##
## 原文 Fig 2B 的版式A见 ref/figs/..._Fig2_HTML.png）：
##   竖排 25 个 level-3 SCP：右侧对齐的 SCP 名称
##   → 紧跟一个彩色小方块A红=上调 / 深蓝=下调），方块内白字写 AUC 排名
##   → 右侧用花括号把 25 行分成 7 组，组名横排ASarcomere dynamics …）
##
## 数据：Supplementary Data 14C (MBCOL3)
##   列：Scp / Association / Entry_type(Diffrna_up|down) / Selection_rank / Selected
## ============================================================
source("R/00_setup.R")

f <- file.path(DIR$data, "suppl", "MBCOL3_cardiotoxic_SCPs.txt")
if (!file.exists(f)) {
  ## 从 zip 里解出来
  z <- "ref/suppl/41467_2024_52145_MOESM4_ESM.zip"
  py <- file.path(DIR$cache, "unz14c.py")
  writeLines(sprintf(
    "import zipfile\nz=zipfile.ZipFile(r'%s')\nn='Supplementary_data/Supplementary Data 14C - Potential cardiotoxic MBCOL3 SCPs.txt'\nopen(r'%s','wb').write(z.read(n))\n",
    z, f), py)
  system2("D:/anaconda/anan/python.exe", shQuote(py))
}
d <- utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE,
                       quote = "", comment.char = "")
cat("Data14C 行数:", nrow(d), "\n")
cat("Association 取值:", paste(unique(d$Association), collapse = " | "), "\n")
cat("Entry_type 取值:", paste(unique(d$Entry_type), collapse = " | "), "\n")

## ---------- 原文 Fig 2B 的行顺序与分组A照抄原图）----------
layout <- tibble::tribble(
  ~scp, ~rank, ~group,
  "Thin myofilament organization", 21, "Sarcomere\ndynamics",
  "Myofibril formation",            6, "Sarcomere\ndynamics",
  "Potassium TM transport",         7, "Electric\ntransmission",
  "Gap junction organization",     24, "Electric\ntransmission",
  "Chloride TM transport",         17, "Transmembrane\ntransport",
  "Water TM transport",             9, "Transmembrane\ntransport",
  "Citric acid cycle",              1, "Metabolism",
  "Polyol pathway",                 4, "Metabolism",
  "Serine and glycine metabolism", 19, "Metabolism",
  "Transamination pathways",       25, "Metabolism",
  "Desaturation of fatty acids",   11, "Metabolism",
  "Cholesterol-sensitive control of SREBP activation", 15, "Metabolism",
  "Cholesterol synthesis",         22, "Metabolism",
  "Cellular cholesterol uptake and efflux", 23, "Metabolism",
  "Adrenergic receptor signaling",  8, "Signaling",
  "Hippo signaling",               13, "Signaling",
  "Thyroid hormone receptor signaling", 10, "Signaling",
  "PDGF receptor signaling",       14, "Signaling",
  "VEGF receptor signaling",       12, "Signaling",
  "HIF-1 receptor signaling pathway", 20, "Signaling",
  "Leptin receptor signaling",      5, "Signaling",
  "GCSF receptor signaling",       16, "Signaling",
  "Biglycan synthesis",            18, "Extracellular\nmatrix",
  "Lysosomal glycoprotein degradation", 2, "Cell cycle",
  "Restriction point",              3, "Cell cycle"
)

## ---------- 把原文的 rank 与数据里的方向对上 ----------
key <- d %>%
  filter(Association == "cardiotoxic TKIs", Selected %in% c(TRUE, "TRUE")) %>%
  select(scp = Scp, direction = Entry_type, rank_data = Selection_rank, auc = ValueForSelection)
cat("cardiotoxic 且 Selected 的 SCP 数:", nrow(key), "\n")

mm <- layout %>% left_join(key, by = "scp")
miss <- mm$scp[is.na(mm$direction)]
if (length(miss)) cat("!! 数据里没匹配到的名称:", paste(miss, collapse = "; "), "\n")
cat("方向分布:\n"); print(table(mm$direction, useNA = "ifany"))

## ---------- 画 ----------
UP   <- "#c0392b"   # 红：被 cardiotoxic TKIs 上调
DOWN <- "#1f3b73"   # 深蓝：下调
n <- nrow(mm)

png(file.path(DIR$fig, "fig2b_paper_layout.png"),
    width = 9.5 * FIG_DPI, height = 8.2 * FIG_DPI, res = FIG_DPI)
par(mar = c(1, 1, 2, 1))
plot.new(); plot.window(xlim = c(0, 1), ylim = c(0.5, n + 0.5))
y <- rev(seq_len(n))                     # 第一行在最上面
text(0.545, y, labels = mm$scp, adj = c(1, 0.5), cex = 0.72)

box_x <- c(0.555, 0.60)
cols <- ifelse(mm$direction == "Diffrna_up", UP, DOWN)
rect(box_x[1], y - 0.42, box_x[2], y + 0.42, col = cols, border = NA)
text(mean(box_x), y, labels = mm$rank, col = "white", cex = 0.62, font = 2)

## 花括号：每组的连续区间
grp_rle <- rle(mm$group)
ends <- cumsum(grp_rle$lengths); starts <- ends - grp_rle$lengths + 1
bx <- 0.645
for (k in seq_along(grp_rle$values)) {
  y1 <- y[ends[k]] - 0.42; y2 <- y[starts[k]] + 0.42
  segments(bx, y1, bx, y2, lwd = 1.1)                 # 竖线
  segments(bx, y1, bx + 0.012, y1, lwd = 1.1)         # 上下两钩
  segments(bx, y2, bx + 0.012, y2, lwd = 1.1)
  text(bx + 0.025, mean(c(y1, y2)), labels = grp_rle$values[k],
       adj = c(0, 0.5), cex = 0.68)
}
mtext("Fig 2B | 对 cardiotoxic TKI 有指示作用的 level-3 亚细胞过程AAUC 排名 1–25）",
      side = 3, adj = 0, line = 1.2, cex = 0.85, font = 2)
par(xpd = NA)
legend(0.02, 0.9, xjust = 0, yjust = 0, bty = "n", cex = 0.72,
       legend = c("Up by cardiotoxic TKIs", "Down by cardiotoxic TKIs"),
       fill = c(UP, DOWN), border = NA)
dev.off()

readr::write_tsv(mm, file.path(DIR$tab, "fig2b_paper_layout.tsv"))
message("[2B] 写出 ", n, " 行 -> figures/fig2b_paper_layout.png")
