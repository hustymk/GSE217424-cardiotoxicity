
source("R/01_helpers.R")

## ---------- 官方药物注释 ----------
drug_tbl <- utils::read.delim(file.path(DIR$data, "suppl", "drug_table.tsv"),
                              check.names = FALSE, stringsAsFactors = FALSE)
names(drug_tbl)[1:3] <- c("abbr", "drug", "class")
drug_tbl$drug <- tolower(drug_tbl$drug)
drug_tbl$class <- gsub("_", " ", drug_tbl$class)   # 官方表用下划线，DEG 长表用空格

CLASS5 <- c("Kinase inhibitor"        = "#ff7f0e",
            "Monoclonal antibody"     = "#d62728",
            "Anthracycline"           = "#9467bd",
            "Cardiovascular drug"     = "#1f77b4",
            "Non cardiovascular drug" = "#17becf")
CLASS5_LAB <- c("Kinase inhibitor"        = "激酶抑制剂",
                "Monoclonal antibody"     = "抗激酶单抗",
                "Anthracycline"           = "蒽环类",
                "Cardiovascular drug"     = "心血管药",
                "Non cardiovascular drug" = "非心血管药")

## ---------- 读作者 DEG 长表，转宽矩阵 ----------
read_long_matrix <- function(f) {
  x <- utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE,
                         quote = "", comment.char = "")
  drug <- tolower(x[[1]]); cell <- x[[4]]; gene <- x[[6]]; val <- x[[9]]
  key  <- paste(cell, drug, sep = "|")
  genes   <- unique(gene)
  samples <- sort(unique(key))
  M <- matrix(0, length(genes), length(samples), dimnames = list(genes, samples))
  M[cbind(match(gene, genes), match(key, samples))] <- val
  meta <- tibble(sample = samples) %>%
    mutate(cell_line = sub("\\|.*$", "", sample),
           drug      = sub("^.*\\|", "", sample)) %>%
    left_join(drug_tbl %>% select(drug, abbr, class), by = "drug")
  list(M = M, meta = meta)
}

d4 <- read_long_matrix(file.path(DIR$data, "suppl", "DEGs_each_cell_line.txt"))
M <- d4$M; meta <- d4$meta
message("[1B-hi] 矩阵: ", nrow(M), " genes x ", ncol(M), " 组合；药物 ",
        dplyr::n_distinct(meta$drug), " 种")

## ---------- 聚类：作者源码用 hclust 平均连锁 + 距离 (1-cor)/2 ----------
cm <- cor(M); cm[is.na(cm)] <- 0
dd <- stats::as.dist((1 - cm) / 2)
hc <- stats::hclust(dd, method = "average")
ord <- hc$order
message("[1B-hi] 树高 = ", sprintf("%.3f", max(hc$height)))

## ---------- 画树：叶子标签 = 药物缩写，颜色 = 官方类别 ----------
n <- length(ord)
cols <- unname(CLASS5[meta$class[ord]])
labs <- meta$abbr[ord]

png(file.path(DIR$fig, "appendix", "fig1b_dendrogram_full.png"),
    width = 13 * FIG_DPI, height = 6.5 * FIG_DPI, res = FIG_DPI)
op <- par(mar = c(6.5, 4.5, 4, 1))
plot(hc, labels = FALSE, hang = -1,
     main = "", sub = "", ylab = "距离 (1 - cor) / 2", xlab = "")
mtext("Fig 1B | 作者 DEG 数据的层次聚类树（hclust 平均连锁，266 组合）",
      side = 3, adj = 0, line = 2.2, cex = 0.9, font = 2)
mtext("叶端标签 = 药物缩写，颜色 = 官方药物类别；距离 = (1-cor)/2",
      side = 3, adj = 0, line = 0.9, cex = 0.6)
text(seq_len(n), par("usr")[3] * 0.02, labels = labs, srt = 90,
     adj = c(0, 0.5), cex = 0.32, col = cols, xpd = NA)
par(xpd = NA)
legend("bottom", inset = c(0, -0.42), horiz = TRUE, bty = "n",
       legend = CLASS5_LAB[names(CLASS5)], fill = CLASS5, border = NA, cex = 0.7)
par(op); dev.off()

## ---------- 定量 ----------
adj_cl   <- mean(meta$cell_line[ord][-1] == meta$cell_line[ord][-n])
adj_drug <- mean(meta$drug[ord][-1] == meta$drug[ord][-n])
message(sprintf("[1B-hi] 相邻叶子：同细胞系 %.1f%%，同药物 %.1f%%",
                100 * adj_cl, 100 * adj_drug))

## ---------- Fig 1C：每药 F1（未分解）----------
f1_before <- f1_per_drug(hc, meta$drug) %>%
  left_join(drug_tbl %>% select(drug, abbr, class), by = "drug")
med_before <- median(f1_before$f1)
message(sprintf("[1C-hi] F1 中位数 = %.3f（论文 0.116）", med_before))

pA <- f1_before %>% mutate(abbr = fct_reorder(abbr, f1)) %>%
  ggplot(aes(abbr, f1, fill = class)) +
  geom_col() +
  geom_hline(yintercept = med_before, linetype = 2, colour = "grey30") +
  annotate("text", x = 1, y = med_before, hjust = 0, vjust = -0.6, size = 2.8,
           label = sprintf("中位数 = %.3f", med_before)) +
  scale_fill_manual(values = CLASS5, labels = CLASS5_LAB, name = NULL) +
  labs(title = "Fig 1C | 各药物聚类效率 F1（未分解）",
       x = NULL, y = "F1 score") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size = 6),
        legend.position = "bottom", legend.text = element_text(size = 7),
        panel.grid.major.x = element_blank())
ggsave(file.path(DIR$fig, "appendix", "fig1c_f1_bar.png"), pA, width = 11, height = 4.2,
       dpi = FIG_DPI)

## ---------- 存 ----------
saveRDS(list(M = M, meta = meta, hc = hc, ord = ord, f1_before = f1_before,
             med_before = med_before),
        file.path(DIR$cache, "hi_fig1b.rds"))
readr::write_tsv(tibble(metric = c("n_genes", "n_samples", "n_drugs",
                                   "adjacent_same_cell_line", "adjacent_same_drug",
                                   "median_F1_before", "clustering"),
                        value = c(nrow(M), ncol(M), dplyr::n_distinct(meta$drug),
                                  adj_cl, adj_drug, med_before, "hclust_average")),
                 file.path(DIR$tab, "hi_fig1b_numbers.tsv"))
message("[1B-hi done] -> figures/fig1b_dendrogram_full.png, fig1c_f1_bar.png")
