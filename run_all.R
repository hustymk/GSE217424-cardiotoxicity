## 一键运行（真实数据）
## 前提：已运行 R/01_download.py 下载 GEO 数据，并把作者数据包中的完整矩阵放到
##       data/suppl/DEG_matrix_fullgenes.txt
##
## 用法：在项目根目录执行 source("run_all.R")

source("R/02_qc.R")                 # 读矩阵 -> 文库质控 -> 合并重复
source("R/03_de.R")                 # 差异表达 -> signed -log10(p) 矩阵

source("R/10_highfid_fig1b.R")      # Fig 1B 与 1C 的基础数据
source("R/11_highfid_fig1d.R")      # Fig 1D 的 after 矩阵与 F1
source("R/18_full_matrix_check.R")  # 作者完整 16345 矩阵的聚类与 F1
source("R/17_fig1cd_bars_paper_layout.R")  # Fig 1C 与 1D 柱状图（原文版式）

source("R/14_fig1b_paper_layout.R")        # Fig 1B 作者版式
source("R/15_fig2b_paper_layout.R")        # Fig 2B
source("R/16_fig2a_flowchart.R")           # Fig 2A

message("\n全部完成。图在 figures/，关键数字在 tables/。")
