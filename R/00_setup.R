## ============================================================
## 00_setup.R：全局设置（所有脚本开头 source 本文件）
## 项目：GSE217424 / Hansen et al. 2024 Nat Commun
##       "Multiscale mapping of transcriptomic signatures for cardiotoxic drugs"
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

## ---- 路径 ----
PROJ  <- "C:/Users/ymk/PycharmProjects/GSE217424_cardiotoxicity"
DIR   <- list(
  data  = file.path(PROJ, "data"),      # 下载的原始文件放这里
  cache = file.path(PROJ, "cache"),     # 中间结果（矩阵、RDS）
  fig   = file.path(PROJ, "figures"),   # 出图
  tab   = file.path(PROJ, "tables")     # 出表
)
invisible(lapply(DIR, dir.create, showWarnings = FALSE, recursive = TRUE))

## ---- 绘图主题 ----
theme_set(theme_minimal(base_size = 10))
FIG_DPI <- 300

## ---- 可选依赖检测 ----
## edgeR 装了就用真正的负二项检验（与论文一致）；
## 没装则自动退回 base R 的 Welch t 检验 + BH（见 03_de.R）
HAVE_EDGER <- requireNamespace("edgeR", quietly = TRUE)
message("[setup] edgeR 可用 = ", HAVE_EDGER)

## ---- 论文关键常数（写报告时引用）----
PAPER <- list(
  doi          = "10.1038/s41467-024-52145-4",
  journal      = "Nature Communications 2024; 15:7968",
  pmid         = "39261481",
  fdr_cutoff   = 0.10,     # 论文：FDR <= 10%
  n_genes      = 16345,    # 论文：每份 DEG 列表含 16345 个基因
  pzero_replace= 1e-64,    # 论文：p = 0 替换为 1e-64
  f1_median    = 0.116,    # 论文：未分解时药物聚类 F1 中位数
  dist_fun     = "(1 - cor) / 2"
)

## ---- 小工具 ----
gz <- function(f) gzfile(f, "rt")

## 从 GEO 样本名解析元信息
## 真实样本名形如："CTRL MSN09-04R-CM Rep 10 (Conv)"、"PAZ MSN08-13R-CM Rep 2 (DGE)"
parse_meta <- function(samples) {
  m <- stringr::str_match(
    samples,
    "^([^.]+)\\.([^.]+)\\.([0-9]+)\\.([0-9]+)\\.([0-9]+)\\.([A-Z]+[0-9]+)$")
  tibble(
    sample     = samples,
    drug       = m[, 2],
    cell_line  = m[, 3],
    experiment = m[, 4],
    dish       = m[, 5],
    plate      = m[, 6],
    well       = m[, 7]
  )
}
