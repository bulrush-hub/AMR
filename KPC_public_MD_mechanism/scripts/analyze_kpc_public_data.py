"""Analyze public KPC-2 MD/QM-MM features and resistance-related structures.

Inputs
------
Zenodo 10.5281/zenodo.7114981 (extracted ``2.enes`` and ``3.datasets``)
RCSB PDB entries 5UL8, 7TB7, 7TBX, 7TC1, 4ZBE, and 8AKL

The script deliberately analyzes the compact, deposited feature/energy arrays rather
than downloading the roughly 3.2 GB coordinate-path archives. The arrays comprise
200 MD-derived reactant conformations per system and their corresponding QM/MM
minimum-energy-path barriers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio.PDB import PDBParser

plt.switch_backend("Agg")


SYSTEMS = {
    "WT-Delta1": ("kpc_wt_d1_16_features.npy", "enes_all_paths_kpc_wt_d1.npz"),
    "WT-Delta2": ("kpc_wt_d2_16_features.npy", "enes_all_paths_kpc_wt_d2.npz"),
    "F72Y-Delta1": ("kpc_y72_d1_16_features.npy", "enes_all_paths_kpc_f72y_d1.npz"),
    "F72Y-Delta2": ("kpc_y72_d2_16_features.npy", "enes_all_paths_kpc_f72y_d2.npz"),
}

FEATURES = [
    ("d1", "Phe72 H-zeta (Tyr72 H-eta) - Glu166 O-epsilon2"),
    ("d2", "Lys73 H-zeta2 - Glu166 O-epsilon2"),
    ("d3", "Deacylating-water H1 - Glu166 O-epsilon2"),
    ("d4", "Deacylating-water H1 - Glu166 O-epsilon1"),
    ("d5", "Asn170 H-delta2 - Glu166 O-epsilon1"),
    ("d6", "Deacylating-water H2 - Asn170 O-delta"),
    ("d7", "Lys73 H-zeta1 - deacylating-water O"),
    ("d8", "Lys73 H-zeta2 - Asn132 O-delta"),
    ("d9", "Deacylating-water O - imipenem C7"),
    ("d10", "Lys73 H-zeta1 - Ser70 O-gamma"),
    ("d11", "Imipenem 6-alpha-OH - deacylating-water O"),
    ("d12", "Imipenem 6-alpha-OH - Asn132 O-delta"),
    ("d13", "Imipenem 6-alpha-OH - Glu166 O-epsilon1"),
    ("d14", "Imipenem 6-alpha-OH - Glu166 O-epsilon2"),
    ("d15", "Lys73 H-zeta1 - Ser130 O-gamma"),
    ("d16", "Ser130 H-gamma - imipenem N4"),
]

SHAP_FILES = {
    "WT-Delta1": "SHAP_values_KPC_WT_IPM_D1.npy",
    "WT-Delta2": "SHAP_values_KPC_WT_IPM_D2.npy",
    "F72Y-Delta1": "SHAP_values_KPC_F72Y_IPM_D1.npy",
    "F72Y-Delta2": "SHAP_values_KPC_F72Y_IPM_D2.npy",
}


@dataclass(frozen=True)
class Paths:
    raw: Path
    output: Path

    @property
    def energies(self) -> Path:
        return self.raw / "2.enes"

    @property
    def features(self) -> Path:
        return self.raw / "3.datasets" / "0.rawdata"

    @property
    def shap(self) -> Path:
        return self.raw / "3.datasets" / "1.shap_values"

    @property
    def tables(self) -> Path:
        return self.output / "tables"

    @property
    def figures(self) -> Path:
        return self.output / "figures"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-dir", type=Path, default=Path("data/raw/public_kpc"), help="Extracted inputs"
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--bootstrap", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=20260820)
    return parser.parse_args()


def require_inputs(paths: Paths) -> None:
    expected = []
    for feature_file, energy_file in SYSTEMS.values():
        expected.extend([paths.features / feature_file, paths.energies / energy_file])
    expected.extend(paths.shap / file_name for file_name in SHAP_FILES.values())
    expected.extend(paths.raw / f"{pdb_id}.pdb" for pdb_id in PDB_STRUCTURES)
    missing = [path for path in expected if not path.exists()]
    if missing:
        joined = "\n  ".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing public-data inputs:\n  {joined}")


def bootstrap_mean_ci(values: np.ndarray, n: int, rng: np.random.Generator) -> tuple[float, float]:
    indices = rng.integers(0, len(values), size=(n, len(values)))
    means = values[indices].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return float(low), float(high)


def bootstrap_difference_ci(
    first: np.ndarray, second: np.ndarray, n: int, rng: np.random.Generator
) -> tuple[float, float]:
    first_indices = rng.integers(0, len(first), size=(n, len(first)))
    second_indices = rng.integers(0, len(second), size=(n, len(second)))
    differences = first[first_indices].mean(axis=1) - second[second_indices].mean(axis=1)
    low, high = np.quantile(differences, [0.025, 0.975])
    return float(low), float(high)


def load_public_arrays(paths: Paths) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    features: dict[str, np.ndarray] = {}
    barriers: dict[str, np.ndarray] = {}
    for system, (feature_file, energy_file) in SYSTEMS.items():
        feature_values = np.load(paths.features / feature_file)
        energy_archive = np.load(paths.energies / energy_file)
        barrier_values = energy_archive["ene_barrier"]
        if feature_values.shape != (len(FEATURES), 200):
            raise ValueError(f"Unexpected feature shape for {system}: {feature_values.shape}")
        if barrier_values.shape != (200,):
            raise ValueError(f"Unexpected barrier shape for {system}: {barrier_values.shape}")
        features[system] = feature_values.astype(float)
        barriers[system] = barrier_values.astype(float)
    return features, barriers


def summarize_energies(
    barriers: dict[str, np.ndarray], bootstrap: int, rng: np.random.Generator
) -> pd.DataFrame:
    rows = []
    for system, values in barriers.items():
        low, high = bootstrap_mean_ci(values, bootstrap, rng)
        rows.append(
            {
                "system": system,
                "n_paths": len(values),
                "mean_kcal_mol": values.mean(),
                "mean_ci95_low": low,
                "mean_ci95_high": high,
                "sd_kcal_mol": values.std(ddof=1),
                "median_kcal_mol": np.median(values),
                "q25_kcal_mol": np.quantile(values, 0.25),
                "q75_kcal_mol": np.quantile(values, 0.75),
                "min_kcal_mol": values.min(),
                "max_kcal_mol": values.max(),
            }
        )
    return pd.DataFrame(rows)


def summarize_contrasts(
    barriers: dict[str, np.ndarray], bootstrap: int, rng: np.random.Generator
) -> pd.DataFrame:
    comparisons = [
        ("F72Y effect in Delta1", "F72Y-Delta1", "WT-Delta1"),
        ("F72Y effect in Delta2", "F72Y-Delta2", "WT-Delta2"),
        ("Delta1 minus Delta2 in WT", "WT-Delta1", "WT-Delta2"),
        ("Delta1 minus Delta2 in F72Y", "F72Y-Delta1", "F72Y-Delta2"),
    ]
    rows = []
    for label, first_name, second_name in comparisons:
        first = barriers[first_name]
        second = barriers[second_name]
        low, high = bootstrap_difference_ci(first, second, bootstrap, rng)
        rows.append(
            {
                "contrast": label,
                "first_system": first_name,
                "second_system": second_name,
                "mean_difference_kcal_mol": first.mean() - second.mean(),
                "difference_ci95_low": low,
                "difference_ci95_high": high,
            }
        )
    return pd.DataFrame(rows)


def summarize_features(
    paths: Paths, features: dict[str, np.ndarray], barriers: dict[str, np.ndarray]
) -> pd.DataFrame:
    rows = []
    for system in SYSTEMS:
        shap_values = np.load(paths.shap / SHAP_FILES[system]).astype(float)
        if shap_values.shape != (200, len(FEATURES)):
            raise ValueError(f"Unexpected SHAP shape for {system}: {shap_values.shape}")
        for index, (feature, description) in enumerate(FEATURES):
            values = features[system][index]
            barrier_ranks = pd.Series(barriers[system]).rank().to_numpy()
            value_ranks = pd.Series(values).rank().to_numpy()
            rows.append(
                {
                    "system": system,
                    "feature": feature,
                    "description": description,
                    "mean_distance_A": values.mean(),
                    "sd_distance_A": values.std(ddof=1),
                    "mean_shap_kcal_mol": shap_values[:, index].mean(),
                    "mean_abs_shap_kcal_mol": np.abs(shap_values[:, index]).mean(),
                    "spearman_r_with_barrier": np.corrcoef(value_ranks, barrier_ranks)[0, 1],
                }
            )
    result = pd.DataFrame(rows)
    result["importance_rank_within_system"] = (
        result.groupby("system")["mean_abs_shap_kcal_mol"]
        .rank(method="min", ascending=False)
        .astype(int)
    )
    return result


PDB_STRUCTURES = {
    # The last value maps an Ambler position to the PDB author's residue number.
    "5UL8": ("WT apo", "A", 0),
    "7TB7": ("D179N apo", "A", 0),
    "7TBX": ("D179Y apo", "A", 0),
    "7TC1": ("D179N + vaborbactam", "A", 0),
    "4ZBE": ("WT + avibactam", "A", -1),
    "8AKL": ("E166Q + meropenem", "A", 0),
}


def residue(chain, number: int):
    key = (" ", number, " ")
    return chain[key] if key in chain else None


def minimum_atom_distance(first, first_names: set[str], second, second_names: set[str]) -> float:
    distances = [
        atom_a - atom_b
        for atom_a in first
        for atom_b in second
        if atom_a.name in first_names and atom_b.name in second_names
    ]
    return float(min(distances)) if distances else float("nan")


def summarize_structures(paths: Paths) -> pd.DataFrame:
    parser = PDBParser(QUIET=True)
    rows = []
    for pdb_id, (description, chain_id, author_delta) in PDB_STRUCTURES.items():
        pdb_path = paths.raw / f"{pdb_id}.pdb"
        structure = parser.get_structure(pdb_id, pdb_path)
        chain = structure[0][chain_id]
        modeled = {
            ambler_number
            for ambler_number in range(164, 180)
            if residue(chain, ambler_number + author_delta) is not None
        }
        missing = sorted(set(range(164, 180)) - modeled)
        res164 = residue(chain, 164 + author_delta)
        res179 = residue(chain, 179 + author_delta)
        salt_bridge = float("nan")
        if res164 is not None and res179 is not None:
            salt_bridge = minimum_atom_distance(
                res164,
                {"NE", "CZ", "NH1", "NH2"},
                res179,
                {"CG", "OD1", "OD2", "ND2", "CZ", "OH"},
            )
        rows.append(
            {
                "pdb_id": pdb_id,
                "description": description,
                "chain": chain_id,
                "ambler_to_author_number_delta": author_delta,
                "omega_loop_modeled_residues_164_179": len(modeled),
                "omega_loop_missing_residues": ";".join(map(str, missing)),
                "residue_179_identity": res179.resname if res179 is not None else "missing",
                "R164_to_179_sidechain_min_distance_A": salt_bridge,
                "sha256": hashlib.sha256(pdb_path.read_bytes()).hexdigest(),
            }
        )
    return pd.DataFrame(rows)


def plot_summary(
    barriers: dict[str, np.ndarray],
    energy_summary: pd.DataFrame,
    feature_summary: pd.DataFrame,
    output_path: Path,
) -> None:
    long_energy = pd.DataFrame(
        [(system, value) for system, values in barriers.items() for value in values],
        columns=["system", "barrier"],
    )
    heatmap = feature_summary.pivot(
        index="feature", columns="system", values="mean_abs_shap_kcal_mol"
    )
    heatmap = heatmap.loc[[name for name, _ in FEATURES], list(SYSTEMS)]

    sns.set_theme(context="paper", style="whitegrid", font_scale=0.9)
    figure = plt.figure(figsize=(12, 8.2), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=[1.0, 1.15])
    axis_box = figure.add_subplot(grid[0, 0])
    axis_mean = figure.add_subplot(grid[0, 1])
    axis_heat = figure.add_subplot(grid[1, :])

    palette = ["#4C78A8", "#72B7B2", "#E45756", "#F2CF5B"]
    sns.boxplot(
        data=long_energy,
        x="system",
        y="barrier",
        order=list(SYSTEMS),
        hue="system",
        palette=palette,
        legend=False,
        showfliers=False,
        ax=axis_box,
    )
    sns.stripplot(
        data=long_energy,
        x="system",
        y="barrier",
        order=list(SYSTEMS),
        color="black",
        size=1.6,
        alpha=0.35,
        jitter=0.25,
        ax=axis_box,
    )
    axis_box.set_title("A  Deposited QM/MM barrier distributions")
    axis_box.set_xlabel("")
    axis_box.set_ylabel("Deacylation barrier (kcal mol$^{-1}$)")
    axis_box.tick_params(axis="x", rotation=20)

    summary = energy_summary.set_index("system").loc[list(SYSTEMS)]
    means = summary["mean_kcal_mol"].to_numpy()
    lower = means - summary["mean_ci95_low"].to_numpy()
    upper = summary["mean_ci95_high"].to_numpy() - means
    positions = np.arange(len(summary))
    axis_mean.bar(positions, means, color=palette, width=0.72)
    axis_mean.errorbar(
        positions,
        means,
        yerr=np.vstack([lower, upper]),
        color="black",
        capsize=4,
        fmt="none",
        linewidth=1.2,
    )
    axis_mean.set_xticks(positions, summary.index, rotation=20)
    axis_mean.set_ylabel("Mean barrier (kcal mol$^{-1}$)")
    axis_mean.set_title("B  Mean barriers with bootstrap 95% CI")
    axis_mean.set_ylim(0, max(upper + means) * 1.18)
    for x, mean in zip(positions, means, strict=True):
        axis_mean.text(x, mean + 1.3, f"{mean:.1f}", ha="center", va="bottom")

    sns.heatmap(
        heatmap,
        cmap="mako",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Mean |SHAP| (kcal mol$^{-1}$)"},
        ax=axis_heat,
    )
    axis_heat.set_title("C  Active-site features driving the calculated barriers")
    axis_heat.set_xlabel("")
    axis_heat.set_ylabel("Feature (definitions in output table)")
    axis_heat.tick_params(axis="x", rotation=0)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def write_metadata(paths: Paths, args: argparse.Namespace) -> None:
    archive_checksums = {}
    for archive in sorted(paths.raw.glob("*.tar.gz")):
        archive_checksums[archive.name] = hashlib.md5(archive.read_bytes()).hexdigest()  # noqa: S324
    metadata = {
        "source": {
            "zenodo_doi": "10.5281/zenodo.7114981",
            "zenodo_license": "CC BY 4.0",
            "rcsb_pdb_ids": list(PDB_STRUCTURES),
        },
        "method": {
            "bootstrap_resamples": args.bootstrap,
            "random_seed": args.seed,
            "barrier_field": "ene_barrier",
            "feature_units": "angstrom",
            "barrier_units": "kcal/mol",
        },
        "archive_md5": archive_checksums,
        "limitations": [
            (
                "The deposited compact arrays are MD-derived snapshots, not full "
                "classical-MD trajectories."
            ),
            (
                "The coordinate archives contain QM/MM reaction paths and were not "
                "required for these summaries."
            ),
            (
                "Crystal disorder supports conformational heterogeneity but is not "
                "itself an MD observable."
            ),
            (
                "The reported barriers are single-point QM/MM potential-energy "
                "differences, not solution free-energy barriers."
            ),
            (
                "The simple bootstrap treats snapshots as independent and may "
                "underestimate uncertainty if the MD-derived samples are correlated."
            ),
        ],
    }
    output = paths.output / "kpc_public_data_metadata.json"
    output.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def main() -> None:
    args = parse_args()
    paths = Paths(raw=args.raw_dir, output=args.output_dir)
    paths.tables.mkdir(parents=True, exist_ok=True)
    paths.figures.mkdir(parents=True, exist_ok=True)
    require_inputs(paths)

    rng = np.random.default_rng(args.seed)
    features, barriers = load_public_arrays(paths)
    energy_summary = summarize_energies(barriers, args.bootstrap, rng)
    contrasts = summarize_contrasts(barriers, args.bootstrap, rng)
    feature_summary = summarize_features(paths, features, barriers)
    structure_summary = summarize_structures(paths)

    energy_summary.to_csv(paths.tables / "kpc_qmmm_energy_summary.csv", index=False)
    contrasts.to_csv(paths.tables / "kpc_qmmm_contrasts.csv", index=False)
    feature_summary.to_csv(paths.tables / "kpc_active_site_feature_summary.csv", index=False)
    structure_summary.to_csv(paths.tables / "kpc_structure_summary.csv", index=False)
    plot_summary(
        barriers,
        energy_summary,
        feature_summary,
        paths.figures / "kpc_public_data_mechanism.png",
    )
    write_metadata(paths, args)

    print(energy_summary.round(3).to_string(index=False))
    print("\nContrasts (first minus second):")
    print(contrasts.round(3).to_string(index=False))
    top_features = feature_summary.query("importance_rank_within_system <= 3").sort_values(
        ["system", "importance_rank_within_system"]
    )
    print("\nTop three SHAP features by system:")
    print(
        top_features[
            ["system", "feature", "description", "mean_abs_shap_kcal_mol"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
