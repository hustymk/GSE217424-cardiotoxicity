# 生成复现图的 R 代码

这个文件夹只放**画图**的代码：Fig 1B、1C、1D、2A、2B 五张图。

四个脚本与 `R/` 下的流水线脚本内容完全相同（`R/` 里那份负责顺序执行），
这里按"论文图号"重命名，方便对照。数据准备与建 cache 的代码仍在 `R/`。

## 图、脚本、输出、输入

| 论文图 | 脚本 | 输出 | 输入 |
|---|---|---|---|
| **Fig 1B** 树 + 54 色带 | `fig1b_paper_layout.R` | `figures/fig1b_paper_layout.png` | `cache/fullmatrix.rds`（完整 16345 基因矩阵） |
| **Fig 1C** 完整数据集的 F1 | `fig1c_1d_bars.R` 上半个面板 | `figures/fig1c_1d_f1_bars.png` | `cache/fullmatrix.rds` |
| **Fig 1D** 投影后的 F1 | `fig1c_1d_bars.R` 下半个面板 | 同上（论文原图也是 C 上 D 下叠放，故共用一个脚本） | `cache/hi_fig1de.rds`（Data 09 药物选择性矩阵，按基因中心化） |
| **Fig 2A** 筛选流程 | `fig2a_flowchart.R` | `figures/fig2a_flowchart.png` | 无，纯图示复刻 |
| **Fig 2B** 25 个 level-3 通路 | `fig2b_paper_layout.R` | `figures/fig2b_paper_layout.png` | `data/suppl/MBCOL3_cardiotoxic_SCPs.txt`（Supplementary Data 14C） |

## 怎么跑

在仓库根目录执行（脚本里的相对路径都以仓库根目录为基准）：

```bash
Rscript figures_code/run_figures.R
```

或者单独跑某一张：

```bash
Rscript figures_code/fig2b_paper_layout.R
```

`run_figures.R` 会先检查 cache，缺了就先跑 `R/18`、`R/10`、`R/11` 把 cache 建起来。
建 cache 需要 `data/suppl/` 下的三张表：

- `DEG_matrix_fullgenes.txt`（完整 16345 基因矩阵，Fig 1B 与 1C）
- `DEGs_each_cell_line.txt`（Data 04，1B/1C 的对照口径）
- `drug_selective_DEGs.txt`（Data 09，Fig 1D）

这些文件都在总材料包的 `data/suppl/` 里；如果缺，先用 `R/01_download.py`
从 GEO 与作者门户取（见 README 第 2 节）。

## 与其他图的关系

`R/` 里另有辅助脚本，不属于交付的五张图：

| 脚本 | 输出 | 用途 |
|---|---|---|


## 一致性

这四个文件由 `tools/sync_packages.py` 每次从 `R/` 重新拷贝，两份内容始终一致，
改代码只改 `R/` 那一份即可。