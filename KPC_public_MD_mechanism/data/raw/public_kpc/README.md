# 原始公开数据说明

本目录由 `scripts/download_public_data.py` 填充。大型上游文件不重复存入 GitHub；这里保留来源、校验值和获取方法，使仓库保持轻量且可审计。

## Zenodo 输入

记录：[10.5281/zenodo.7114981](https://doi.org/10.5281/zenodo.7114981)，许可证：CC BY 4.0。

| 归档 | MD5 | 本分析用途 |
|---|---|---|
| `2.enes.tar.gz` | `2bdb00cb6bd463713359a911dec28a11` | 四体系各 200 条 QM/MM 去酰化能垒 |
| `3.datasets.tar.gz` | `76f646c858d2946135b00f152a9c129e` | 16 个活性位点距离特征及 SHAP 值 |

上游的 `4.code.tar.gz` 未被本分析调用，故下载脚本不获取它。上游约 GB 级反应路径坐标也不是本汇总统计的必需输入。

## RCSB PDB 输入

下载 5UL8、7TB7、7TBX、7TC1、4ZBE 和 8AKL 的 PDB 格式坐标，用于 Ω-loop 建模完整性、179 位残基身份和 R164–179 最短侧链距离的复核。

## 生成方式

在项目子目录根部执行：

```bash
python scripts/download_public_data.py
```

下载后会生成 `2.enes/`、`3.datasets/` 和六个 `*.pdb` 文件。它们被 `.gitignore` 排除；可提交的派生表位于 `outputs/tables/`。

