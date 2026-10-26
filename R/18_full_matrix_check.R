## 用作者发布的完整 16345 基因矩阵重算 Fig 1B / 1C
source("R/01_helpers.R")

f <- file.path(DIR$data, "suppl", "DEG_matrix_fullgenes.txt")
x <- utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE,
                       quote = "", comment.char = "")
cat("矩阵:", nrow(x), "genes x", ncol(x) - 1, "samples\n")

M <- as.matrix(x[, -1]); rownames(M) <- x[[1]]

## 列名形如 Diffrna-H48-Anthracycline-Cell_line.MSN01_03R-DAU-Plate.0
cn <- colnames(M)
part <- strsplit(cn, "-", fixed = TRUE)
drug  <- vapply(part, function(v) v[[5]], "")
cls   <- vapply(part, function(v) v[[3]], "")
cell  <- vapply(part, function(v) sub("Cell_line\\.", "", v[[4]]), "")
cat("药物数:", length(unique(drug)), " 细胞系:", length(unique(cell)), "\n")

## 作者的做法：距离 (1-cor)/2 + hclust 平均连锁
cm <- cor(M); cm[is.na(cm)] <- 0
hc <- stats::hclust(stats::as.dist((1 - cm) / 2), method = "average")
ord <- hc$order

adj  <- mean(drug[ord][-1] == drug[ord][-length(ord)])
adjc <- mean(cell[ord][-1] == cell[ord][-length(ord)])
cat(sprintf("相邻叶子：同药物 %.1f%%，同细胞系 %.1f%%（随机 %.1f%%）\n",
            100 * adj, 100 * adjc, 100 / length(unique(drug))))

f1 <- f1_per_drug(hc, drug)
cat(sprintf("F1 中位数 = %.3f（论文 0.116）\n", median(f1$f1)))
cat(sprintf("F1 范围 = %.3f 到 %.3f\n", min(f1$f1), max(f1$f1)))

readr::write_tsv(tibble(metric = c("n_genes", "n_samples", "n_drugs",
                                   "adjacent_same_drug", "adjacent_same_cell_line",
                                   "median_F1"),
                        value = c(nrow(M), ncol(M), length(unique(drug)),
                                  adj, adjc, median(f1$f1))),
                 file.path(DIR$tab, "fig1b_fullmatrix_numbers.tsv"))
saveRDS(list(M = M, drug = drug, cell = cell, cls = cls, hc = hc,
             f1 = f1, cm = cm),
        file.path(DIR$cache, "fullmatrix.rds"))
cat("[done] -> tables/fig1b_fullmatrix_numbers.tsv\n")
