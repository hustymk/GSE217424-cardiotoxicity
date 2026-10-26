## ============================================================
## 02_qc.R：读入计数矩阵 -> 文库级质控 -> 按「细胞系 x 药物」合并重复
##
## 真实数据结构（已核实）：
##   文件 39.5 MB，38478 个基因 x 1171 个文库
##   列名 = State.Cell.Experiment.Dish.Plate.Well，如 PAZ.MSN08-13R-CM.116.2.0.A1
##   55 个 State（54 种药 + CTRL），6 个细胞系，共 305 个「细胞系 x 药物」组合
##
## 论文的输入是「每个组合一条 DEG 列表」（266 条）。原始数据里每个组合有
## 2-12 个重复文库，所以这里先把同一组合的重复文库**按计数求和合并**成一条
## 伪文库，再与同细胞系的 CTRL 比较，这样每个组合恰好得到一条 DEG 列表，
## 与论文的口径一致。
## ============================================================
source("R/00_setup.R")

## ---------- 读取 ----------
## 注意：该文件的表头**少了第一列（基因列）的名字**：
##   header 行只有 1171 个字段（全是样本名），而数据行有 1172 个字段
##   （基因 + 1171 个样本）。直接 read.delim(header=TRUE) 会整体错位，
##   实测只能读出 2522 个基因，所以这里手工修正表头。
read_counts <- function(f) {
  con <- gzfile(f, "rt"); hdr <- readLines(con, n = 1); close(con)
  hdr_fields <- strsplit(hdr, "\t", fixed = TRUE)[[1]]

  x <- utils::read.delim(gzfile(f), header = FALSE, skip = 1,
                         check.names = FALSE, stringsAsFactors = FALSE,
                         quote = "", comment.char = "")
  if (ncol(x) == length(hdr_fields) + 1) {
    names(x) <- c("gene_id", hdr_fields)          # 补上缺失的基因列名
  } else if (ncol(x) == length(hdr_fields)) {
    names(x) <- hdr_fields
    names(x)[1] <- "gene_id"
  } else {
    stop("表头与数据列数无法对齐: header=", length(hdr_fields),
         " data=", ncol(x))
  }
  gene <- as.character(x$gene_id)
  m <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- gene
  m[!duplicated(gene), , drop = FALSE]
}

counts_file <- file.path(DIR$data, "GSE217421_Conv-RNAseq-Read-Counts.tsv.gz")
if (!file.exists(counts_file)) stop("找不到计数矩阵：", counts_file,
                                    "\n请先运行 R/01_download.py")
message("[qc] 正在读取计数矩阵（39 MB，约需 30-60 秒）...")
counts <- read_counts(counts_file)
message("[qc] counts: ", nrow(counts), " genes x ", ncol(counts), " 文库")

## ---------- 元信息 ----------
meta <- parse_meta(colnames(counts))
if (any(is.na(meta$drug))) {
  message("[qc] 警告：", sum(is.na(meta$drug)), " 个列名无法解析")
}
meta <- meta %>% mutate(combo = paste(cell_line, drug, sep = "|"))
message("[qc] State 数 = ", dplyr::n_distinct(meta$drug),
        "（含 CTRL）；细胞系 = ", dplyr::n_distinct(meta$cell_line),
        "；组合数 = ", dplyr::n_distinct(meta$combo))

## ---------- 文库级 QC ----------
lib <- colSums(counts)
det <- colSums(counts > 0)
qc <- tibble(sample = colnames(counts),
             combo = meta$combo, drug = meta$drug, cell_line = meta$cell_line,
             lib_size = lib, n_detected = det) %>%
  group_by(combo) %>% mutate(combo_median = median(lib_size)) %>% ungroup()

## 文库大小明显偏离同组中位数（< 1/3 或 > 3 倍）的文库视为可疑
qc <- qc %>% mutate(flagged = lib_size < combo_median / 3 | lib_size > combo_median * 3)
message("[qc] 标记可疑文库: ", sum(qc$flagged), " / ", nrow(qc))
readr::write_tsv(qc, file.path(DIR$tab, "01_qc_metrics.tsv"))

## ---------- 合并重复文库 -> 每个组合一条伪文库 ----------
keep_samples <- qc$sample[!qc$flagged]
counts <- counts[, keep_samples, drop = FALSE]
meta <- meta[match(keep_samples, meta$sample), ]

combos <- unique(meta$combo)
agg <- matrix(0, nrow = nrow(counts), ncol = length(combos),
              dimnames = list(rownames(counts), combos))
for (cb in combos) {
  idx <- which(meta$combo == cb)
  agg[, cb] <- if (length(idx) == 1) counts[, idx] else rowSums(counts[, idx, drop = FALSE])
}
meta_combo <- meta %>%
  group_by(combo) %>%
  summarise(drug = dplyr::first(drug),
            cell_line = dplyr::first(cell_line),
            n_rep = dplyr::n(), .groups = "drop") %>%
  mutate(sample = combo)
message("[qc] 合并后: ", ncol(agg), " 个组合；每组合重复数中位数 = ",
        median(meta_combo$n_rep))

## ---------- QC 图 ----------
p1 <- ggplot(qc, aes(x = reorder(sample, lib_size), y = lib_size / 1e6)) +
  geom_col(aes(fill = flagged)) +
  scale_fill_manual(values = c(`FALSE` = "grey35", `TRUE` = "#d62728"), name = "可疑") +
  labs(title = "文库大小（每个文库）", x = NULL, y = "Library size (M reads)") +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
ggsave(file.path(DIR$fig, "qc_libsize.png"), p1, width = 9, height = 3.2, dpi = FIG_DPI)

p2 <- ggplot(qc, aes(lib_size / 1e6, n_detected, color = drug == "CTRL")) +
  geom_point(size = 1.2, alpha = .7) +
  scale_color_manual(values = c(`FALSE` = "grey45", `TRUE` = "#1f77b4"),
                     labels = c("药物处理", "对照 CTRL"), name = NULL) +
  labs(title = "文库大小 vs 检出基因数", x = "Library size (M)", y = "Detected genes")
ggsave(file.path(DIR$fig, "qc_libsize_vs_detected.png"), p2, width = 5.2, height = 4, dpi = FIG_DPI)

## ---------- 保存 ----------
saveRDS(list(counts = agg, meta = meta_combo, qc = qc),
        file.path(DIR$cache, "01_qc.rds"))
message("[qc done] -> cache/01_qc.rds （", nrow(agg), " genes x ", ncol(agg), " 组合）")
