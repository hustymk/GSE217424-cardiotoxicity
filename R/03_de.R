## ============================================================
## 03_de.R：每个「细胞系 x 药物」组合 vs 同细胞系 CTRL -> 差异表达
##             生成论文 Fig1B 的输入：signed -log10(p) 矩阵
##
## 论文做法：edgeR exactTest，FDR<=0.10，p=0 记作 1e-64，方向取符号。
## 本机：无 Bioconductor；且每个组合已合并为一条伪文库，两组各一条，
##       成组 t 检验不可用。改用**二项检验**：
##         t ~ Binomial(t + c, p0)，p0 由两条伪文库的文库大小定，
##       等价于 edgeR exactTest 在"单样本、不估离散度"时的极限。
## 若装了 edgeR，USE_EDGER 置 TRUE 即切回论文原方法。
## ============================================================
source("R/00_setup.R")

USE_EDGER <- HAVE_EDGER

qc_obj <- readRDS(file.path(DIR$cache, "01_qc.rds"))
counts <- qc_obj$counts
meta   <- qc_obj$meta

## ---------- 基因过滤（相当于 QC 的一部分，写进报告）----------
## 保留在 >= 20% 的样本里 CPM > 1 的基因，去掉几乎不表达的基因
cpm <- t(t(counts) / colSums(counts)) * 1e6
keep_gene <- rowMeans(cpm > 1) >= 0.2
counts <- counts[keep_gene, , drop = FALSE]
message("[de] 基因过滤: ", sum(keep_gene), " / ", length(keep_gene),
        " 个基因保留（CPM>1 且 >=20% 样本）")

lib <- colSums(counts)

## ---------- 二项检验（向量化）----------
binom_rows <- function(t_vec, c_vec, lib_t, lib_c) {
  n  <- t_vec + c_vec
  p0 <- lib_t / (lib_t + lib_c)
  lo <- stats::pbinom(t_vec, n, p0)
  hi <- stats::pbinom(t_vec - 1, n, p0, lower.tail = FALSE)
  p  <- pmin(1, 2 * pmin(lo, hi))
  p[n == 0] <- 1
  p
}

de_one <- function(s, ctrl) {
  if (USE_EDGER) {
    y  <- edgeR::DGEList(counts = counts[, c(s, ctrl), drop = FALSE])
    y  <- edgeR::calcNormFactors(y)
    et <- edgeR::exactTest(y, pair = c(2, 1))
    tb <- edgeR::topTags(et, n = Inf, sort.by = "none")$table
    return(list(p = tb$PValue, logfc = tb$logFC, gene = rownames(tb)))
  }
  t_vec <- counts[, s]; c_vec <- counts[, ctrl]
  p <- binom_rows(t_vec, c_vec, lib[s], lib[ctrl])
  lfc <- log2(((t_vec + 1) / (lib[s] / 1e6)) / ((c_vec + 1) / (lib[ctrl] / 1e6)))
  list(p = p, logfc = lfc, gene = rownames(counts))
}

## ---------- 主循环 ----------
treated <- which(meta$drug != "CTRL")
message("[de] 待计算组合数 = ", length(treated), "；方法 = ",
        if (USE_EDGER) "edgeR exactTest" else "二项检验（fallback）")

signed <- matrix(NA_real_, nrow(counts), length(treated),
                 dimnames = list(rownames(counts), meta$sample[treated]))
fdr <- signed; lfc <- signed

t0 <- Sys.time(); nskip <- 0
for (k in seq_along(treated)) {
  s  <- meta$sample[treated[k]]
  cl <- meta$cell_line[treated[k]]
  ctl <- meta$sample[meta$drug == "CTRL" & meta$cell_line == cl]
  if (length(ctl) < 1) { nskip <- nskip + 1; next }
  out <- de_one(s, ctl[1])
  g <- match(rownames(counts), out$gene)
  p <- out$p[g]; lf <- out$logfc[g]
  p[!is.finite(p)] <- 1; lf[!is.finite(lf)] <- 0
  signed[, k] <- sign(lf) * (-log10(pmax(p, PAPER$pzero_replace)))
  fdr[, k]    <- stats::p.adjust(p, method = "BH")
  lfc[, k]    <- lf
  if (k %% 50 == 0) message("  ..", k, "/", length(treated),
      " (", round(as.numeric(difftime(Sys.time(), t0, units = "min")), 1), " min)")
}
message("[de] 用时 ", round(as.numeric(difftime(Sys.time(), t0, units = "min")), 2),
        " 分钟；跳过 ", nskip, " 个（该细胞系无 CTRL）")

ok <- colSums(!is.na(signed)) > nrow(signed) / 2
signed <- signed[, ok, drop = FALSE]; fdr <- fdr[, ok, drop = FALSE]
lfc <- lfc[, ok, drop = FALSE]; meta_t <- meta[treated[ok], ]

sig <- colSums(fdr < PAPER$fdr_cutoff, na.rm = TRUE)
message("[de] 每个组合的显著基因数（FDR<0.1）：中位数 = ", median(sig),
        "，范围 ", min(sig), "-", max(sig))
readr::write_tsv(tibble(sample = colnames(signed), drug = meta_t$drug,
                        cell_line = meta_t$cell_line, n_sig_FDR0.1 = sig),
                 file.path(DIR$tab, "02_deg_counts.tsv"))

saveRDS(list(signed = signed, fdr = fdr, logfc = lfc, meta = meta_t),
        file.path(DIR$cache, "02_signed_log10p.rds"))
message("[de done] -> cache/02_signed_log10p.rds  (", nrow(signed), " x ", ncol(signed), ")")
