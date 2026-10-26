## ============================================================
## 11_highfid_fig1d.R：Fig 1D 的 after 矩阵与 F1
##
## 输入：作者发布的 Supplementary Data 09
##   "Drug-selective DEGs in each cell line"（他们算好的药物选择性 DEG，
##   用于 Fig 1D 的 after 矩阵）
## 对比对象：Fig 1C 用的是 Supplementary Data 04（未分解的 DEG）
##
## 所以这一版的 before / after 是**作者的两张矩阵**，不是我的近似投影。
## ============================================================
source("R/01_helpers.R")

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

## ---------- 读 after 矩阵（作者算好的药物选择性 DEG）----------
d9 <- read_long_matrix(file.path(DIR$data, "suppl", "drug_selective_DEGs.txt"))
M9 <- d9$M; meta9 <- d9$meta
message("[1D-hi] 药物选择性矩阵: ", nrow(M9), " genes x ", ncol(M9), " 组合")

## 补充数据只给每组合 top600 的选择性基因，其余基因在重建时补 0。
## 补 0 会让所有样本共享一大块相同的 0，同药样本的相似度被这个公共基线抬高，
## 直接算相关会把 F1 顶到 1.000。按基因中心化去掉这条公共基线后再聚类。
M9 <- M9 - rowMeans(M9)
message("[1D-hi] 已按基因中心化，去掉补 0 的公共基线")

b <- readRDS(file.path(DIR$cache, "hi_fig1b.rds"))
f1_before <- b$f1_before; med_before <- b$med_before

cm9 <- cor(M9); cm9[is.na(cm9)] <- 0
hc9 <- stats::hclust(stats::as.dist((1 - cm9) / 2), method = "average")
ord9 <- hc9$order
f1_after <- f1_per_drug(hc9, meta9$drug) %>%
  left_join(drug_tbl %>% select(drug, abbr, class), by = "drug")
med_after <- median(f1_after$f1)
message(sprintf("[1D-hi] F1 中位数: %.3f -> %.3f（论文 0.116 -> 大幅上升）",
                med_before, med_after))

## ---------- Fig 1D：before -> after 哑铃图 ----------
cmp <- bind_rows(f1_before %>% mutate(stage = "before"),
                 f1_after  %>% mutate(stage = "after")) %>%
  select(abbr, f1, stage, class) %>%
  pivot_wider(names_from = stage, values_from = f1) %>%
  mutate(delta = after - before, abbr = fct_reorder(abbr, delta))

pB <- ggplot(cmp, aes(y = abbr)) +
  geom_segment(aes(x = before, xend = after, yend = abbr, colour = class),
               linewidth = .7, alpha = .85) +
  geom_point(aes(x = before), size = 1.4, colour = "grey25") +
  geom_point(aes(x = after),  size = 1.8, colour = "#b2182b") +
  scale_colour_manual(values = CLASS5, guide = "none") +
  labs(title = "Fig 1D | 未分解 vs 药物选择性 DEG 的 F1",
       subtitle = sprintf("中位数 %.3f -> %.3f；%d/%d 种药物上升",
                          med_before, med_after, sum(cmp$delta > 0), nrow(cmp)),
       x = "F1 score", y = NULL) +
  theme(axis.text.y = element_text(size = 5.5),
        panel.grid.major.y = element_blank())
## 哑铃图已由柱状图取代（见 R/17），此处不再输出

saveRDS(list(cmp = cmp, hc9 = hc9, ord9 = ord9, M9 = M9, meta9 = meta9),
        file.path(DIR$cache, "hi_fig1de.rds"))
message("[1D-hi done] -> cache/hi_fig1de.rds")
