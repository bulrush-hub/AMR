# KPC 公开 MD/QM-MM 数据耐药机制分析

本目录把 KPC 耐药机制分析所需的代码、公开数据来源、派生数据和机制解释放在同一个可复现单元中。分析对象包括 KPC-2 WT/F72Y 与亚胺培南 Δ1/Δ2 互变异构体，以及 WT、D179N、D179Y 等公开晶体结构。

> 结论边界：Zenodo 中的是从 MD 轨迹抽取的反应物构象、活性位点特征及 QM/MM 最小能量路径结果；它们不是完整的经典 MD 坐标轨迹。晶体结构无序可支持构象异质性，但不能单独当作 MD 观测量。

## 目录说明

```text
KPC_public_MD_mechanism/
├── README.md                         # 本说明与复现入口
├── KPC_RESISTANCE_MECHANISM.md       # 中文机制报告及证据等级
├── LICENSE                           # 本项目代码许可证（MIT）
├── requirements.txt                  # Python 依赖
├── data/raw/public_kpc/README.md      # 原始数据来源、许可证与校验值
├── scripts/
│   ├── download_public_data.py       # 下载、校验并解压公开输入
│   └── analyze_kpc_public_data.py    # 统计、结构测量与绘图
└── outputs/
    ├── kpc_public_data_metadata.json # 参数、来源、局限性和校验值
    ├── figures/                      # 汇总图
    └── tables/                       # 可审计的派生 CSV 数据
```

## 主要结果

- WT 的平均去酰化能垒为 37.61 kcal/mol（Δ1）和 32.43 kcal/mol（Δ2）。
- F72Y 将平均能垒分别提高 11.10 和 5.28 kcal/mol，说明 72 位—Glu166 氢键重排会破坏 Glu166—去酰化水质子网络。
- 四个体系中，去酰化水 O 到亚胺培南 C7 的距离 `d9` 均为平均绝对 SHAP 值最高的特征，直接指向亲核水的近攻击构型。
- WT 5UL8 的 R164–D179 最短侧链距离为 2.76 Å，D179N 7TB7 为 5.70 Å；D179Y 7TBX 的 Ω-loop 164–178 在链 A 中未建模，支持 Ω-loop 网络被破坏。

完整解释、表型边界和文献链接见 [KPC_RESISTANCE_MECHANISM.md](KPC_RESISTANCE_MECHANISM.md)。

## 一键复现

要求 Python 3.10 或更高版本。在本目录执行：

```bash
python -m venv .venv
# Windows PowerShell: .venv\Scripts\Activate.ps1
# Linux/macOS: source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/download_public_data.py
python scripts/analyze_kpc_public_data.py
```

下载脚本只获取分析所需的紧凑数组（约 800 个 MD 派生构象）和 6 个 PDB 文件，不下载约 GB 级的 QM/MM 路径坐标归档。重复运行时，校验正确的已有文件会被保留。

## 派生数据字典

| 文件 | 内容 | 关键单位/说明 |
|---|---|---|
| `kpc_qmmm_energy_summary.csv` | 四体系去酰化能垒汇总 | kcal/mol；含均值、SD、分位数和 bootstrap 95% CI |
| `kpc_qmmm_contrasts.csv` | F72Y 效应与 Δ1–Δ2 差异 | kcal/mol；正值表示前一体系能垒更高 |
| `kpc_active_site_feature_summary.csv` | 16 个活性位点距离、SHAP 与能垒秩相关 | 距离为 Å；SHAP 为 kcal/mol |
| `kpc_structure_summary.csv` | Ω-loop 建模完整性及 R164–179 几何 | 距离为 Å；缺失残基以分号分隔 |

## 数据来源与许可

- MD 派生数组和 QM/MM 结果：[Zenodo 10.5281/zenodo.7114981](https://doi.org/10.5281/zenodo.7114981)，原记录标注 CC BY 4.0。
- 结构坐标：[RCSB PDB](https://www.rcsb.org/)，条目 5UL8、7TB7、7TBX、7TC1、4ZBE、8AKL。
- 本目录新增分析代码采用 MIT 许可证。转载原始或派生数据时仍需按各上游数据库的许可要求署名。

## 质量控制

- 下载时校验 Zenodo 归档 MD5；PDB 文件要求包含有效 `ATOM` 记录。
- 分析脚本检查每个体系恰有 200 条能垒和 `16 × 200` 个特征值。
- 固定随机种子 `20260820`，默认使用 20,000 次 bootstrap。
- 输出元数据明确记录方法限制，避免把势能垒误写为溶液自由能垒。

