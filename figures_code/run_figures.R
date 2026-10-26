## ============================================================
## run_figures.R：一键生成 Fig 1B、1C、1D、2A、2B 五张复现图
##
## 用法（在仓库根目录执行）：
##   Rscript figures_code/run_figures.R
##
## 会先检查 cache，缺了就先跑 R/18、R/10、R/11 建 cache；
## cache 建好后依次跑 figures_code 下的四个出图脚本。
## ============================================================

need_data <- file.path("data", "suppl",
                       c("DEG_matrix_fullgenes.txt",
                         "DEGs_each_cell_line.txt",
                         "drug_selective_DEGs.txt"))
miss_data <- need_data[!file.exists(need_data)]
if (length(miss_data)) {
  stop("缺输入数据：\n  ", paste(miss_data, collapse = "\n  "),
       "\n见 README 第 2 节，或先跑 R/01_download.py")
}

need_cache <- file.path("cache",
                        c("fullmatrix.rds", "hi_fig1b.rds", "hi_fig1de.rds"))
if (!all(file.exists(need_cache))) {
  message("[run_figures] cache 不全，先建 cache（这一步较慢）")
  source("R/00_setup.R")
  source("R/18_full_matrix_check.R")   # 完整 16345 矩阵 -> fullmatrix.rds
  source("R/10_highfid_fig1b.R")       # 1B/1C 的 before -> hi_fig1b.rds
  source("R/11_highfid_fig1d.R")       # 1D 的 after -> hi_fig1de.rds
} else {
  message("[run_figures] cache 已就绪，直接出图")
}

source("figures_code/fig1b_paper_layout.R")   # Fig 1B
source("figures_code/fig1c_1d_bars.R")        # Fig 1C + 1D
source("figures_code/fig2a_flowchart.R")      # Fig 2A
source("figures_code/fig2b_paper_layout.R")   # Fig 2B

message("[run_figures] 完成，五张图在 figures/ 下")
