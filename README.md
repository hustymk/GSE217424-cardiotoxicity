# 心毒性药物转录组复现

复现 Hansen et al., *Multiscale mapping of transcriptomic signatures for cardiotoxic drugs*,
Nature Communications 2024; 15:7968（DOI 10.1038/s41467-024-52145-4）的五张主图：

| 图 | 输出 | 说明 |
|---|---|---|
| Fig 1B | `figures/fig1b_paper_layout.png` | 266 个组合的聚类树与 54 色药物色带，同药相邻 15.5%，随机 1.9% |
| Fig 1C + 1D | `figures/fig1c_1d_f1_bars.png` | F1 中位数 0.116 到 0.909，54 种药里 45 种上升 |
| Fig 2A | `figures/fig2a_flowchart.png` | 从转录组响应到毒性相关通路的筛选流程 |
| Fig 2B | `figures/fig2b_paper_layout.png` | 25 个 level-3 亚细胞过程，按 7 类功能分组 |

`figures/参考_原文Fig1.png` 是论文原图，用于版式对照。

## 环境

- R 4.x，`tidyverse`、`readr`、`stringr`
- 画图只用 `ggplot2` 与 base graphics，不需要 Bioconductor
- 下载脚本用 Python 3

## 运行

```bash
# 1) 取 GEO 数据（约 40 MB），放到 data/
python R/01_download.py

# 2) 把作者数据包里的完整矩阵放到 data/suppl/DEG_matrix_fullgenes.txt
#    GEO 计数矩阵：https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217421/suppl/
#    论文补充材料：https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-52145-4/MediaObjects/41467_2024_52145_MOESM4_ESM.zip
#    完整矩阵：https://iyengarlab.org/dtoxs/files/LINCS_DToxS_SVD.zip
#    路径 Results_to_run_R_SVD1to17_in_isolation/1_DEGs/SVD_DEGenes_iPSCdCMs_P0_Signed_minus_log10pvalue_topall.txt

# 3) 一键出图
Rscript -e 'source("run_all.R")'
```

只出五张主图可以走：

```bash
Rscript figures_code/run_figures.R
```

## 目录

```
R/               01_download.py、00_setup.R、01_helpers.R、02_qc.R、03_de.R、
                 10、11、14、15、16、17、18
figures_code/    Fig 1B/1C/1D/2A/2B 的出图脚本与 run_figures.R
figures/         交付的五张主图与论文原图
tables/          数据清单与两张关键数字表
ref/             论文原图、图注、作者 R 源码副本
run_all.R        一键运行
```

## 与作者源码的一致性

| 环节 | 作者源码 | 本项目 |
|---|---|---|
| 距离 | `as.dist((1 - cor(x)) / 2)` | 同 |
| 聚类 | `hclust(x, method = "average")` | 同 |
| F1 有效簇 | 簇内至少 2 个同药样本 | 同，枚举全部内部节点取最大 |

## 关键数字

- 1B：相邻同药 15.5%，随机 1.9%；相邻同细胞系 69.4%
- 1C：F1 中位数 0.116，与原文一致
- 1D：F1 中位数 0.909，45/54 上升
- 2B：25 个通路，AUC 排名 1 到 25，7 个功能分组

## 说明

`figures/` 与 `tables/` 放的是交付用到的子集，脚本运行后还会生成其他中间图与表。
