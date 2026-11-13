#!/usr/bin/env python
"""Plot expression and strict validation status of the locked 12-gene signature."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = PROJECT_ROOT / "results" / "tables" / "validation" / "validation_gene_results.csv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "figures" / "validation"

REQUIRED_COLUMNS = {
    "cohort",
    "gene",
    "expected_direction",
    "observed_direction",
    "logFC",
    "adj_p_value",
    "AUC",
}

STATUS_ORDER = [
    "Validated",
    "Not significant",
    "Significant; AUC < 0.70",
    "Direction mismatch",
    "Missing",
]
STATUS_COLORS = {
    "Validated": "#238b57",
    "Not significant": "#e69f00",
    "Significant; AUC < 0.70": "#56b4e9",
    "Direction mismatch": "#cc79a7",
    "Missing": "#bdbdbd",
}
STATUS_LABELS = {
    "Validated": "VALIDATED",
    "Not significant": "NOT SIGNIFICANT",
    "Significant; AUC < 0.70": "SIG.; AUC < 0.70",
    "Direction mismatch": "DIRECTION MISMATCH",
    "Missing": "MISSING",
}
COHORTS = [
    {
        "label": "GSE53845 (Microarray)",
        "expression_file": "GSE53845_expr_gene_level.csv",
        "metadata_file": "GSE53845_metadata.csv",
        "y_label": "Normalized microarray expression",
    },
    {
        "label": "GSE150910 (RNA-seq)",
        "expression_file": "GSE150910_expr_gene_level.csv",
        "metadata_file": "GSE150910_metadata.csv",
        "y_label": "log2(CPM + 1)",
    },
]
GROUP_COLORS = {"IPF": "#e84a5f", "Control": "#168ac2"}


def validation_status(row: pd.Series) -> str:
    if pd.isna(row["logFC"]):
        return "Missing"
    if row["observed_direction"] != row["expected_direction"]:
        return "Direction mismatch"
    if pd.isna(row["adj_p_value"]) or row["adj_p_value"] >= 0.05:
        return "Not significant"
    if pd.isna(row["AUC"]) or row["AUC"] < 0.70:
        return "Significant; AUC < 0.70"
    return "Validated"


def significance_label(q_value: float) -> str:
    if pd.isna(q_value):
        return "ND"
    if q_value < 0.001:
        return "***"
    if q_value < 0.01:
        return "**"
    if q_value < 0.05:
        return "*"
    return "ns"


def create_expression_plot(data: pd.DataFrame, genes: list[str]) -> None:
    rng = np.random.default_rng(20261102)
    fig, axes = plt.subplots(1, len(COHORTS), figsize=(16, 7), squeeze=False)

    for ax, cohort_info in zip(axes[0], COHORTS):
        expr_path = PROJECT_ROOT / "data" / "processed" / cohort_info["expression_file"]
        meta_path = PROJECT_ROOT / "data" / "processed" / cohort_info["metadata_file"]
        if not expr_path.is_file() or not meta_path.is_file():
            raise FileNotFoundError(
                f"Validation expression inputs not found: {expr_path} or {meta_path}. "
                "Run scripts/13_external_validation.R first."
            )

        expression = pd.read_csv(expr_path).set_index("gene")
        metadata = pd.read_csv(meta_path)
        if "sample_id" not in metadata or "group" not in metadata:
            raise ValueError(f"Metadata must contain sample_id and group columns: {meta_path}")
        if metadata["sample_id"].duplicated().any():
            raise ValueError(f"Metadata contains duplicate sample IDs: {meta_path}")
        sample_ids = metadata["sample_id"].tolist()
        missing_samples = set(sample_ids).difference(expression.columns)
        if missing_samples:
            raise ValueError(
                f"Expression matrix {expr_path} is missing {len(missing_samples)} metadata samples."
            )

        cohort_results = data[data["cohort"] == cohort_info["label"]].set_index("gene")
        if cohort_results.index.duplicated().any():
            raise ValueError(f"Validation results contain duplicate genes for {cohort_info['label']}.")

        group_counts = metadata["group"].value_counts()
        present_count = 0
        x_positions = list(range(len(genes)))
        offsets = {"Control": -0.19, "IPF": 0.19}
        expression_values = []

        for gene_index, gene in enumerate(genes):
            if gene not in expression.index:
                ax.text(
                    gene_index,
                    0.03,
                    "ND",
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="bottom",
                    color="#8295ad",
                    fontweight="bold",
                )
                continue

            present_count += 1
            for group in ("IPF", "Control"):
                group_samples = metadata.loc[metadata["group"] == group, "sample_id"].tolist()
                values = pd.to_numeric(expression.loc[gene, group_samples], errors="coerce").dropna().to_numpy()
                if values.size == 0:
                    continue
                expression_values.extend(values)
                position = gene_index + offsets[group]
                box = ax.boxplot(
                    values,
                    positions=[position],
                    widths=0.32,
                    patch_artist=True,
                    showfliers=False,
                    manage_ticks=False,
                    medianprops={"color": GROUP_COLORS[group], "linewidth": 1.8},
                    whiskerprops={"color": GROUP_COLORS[group], "linewidth": 1.1},
                    capprops={"color": GROUP_COLORS[group], "linewidth": 1.1},
                    boxprops={"edgecolor": GROUP_COLORS[group], "linewidth": 1.3},
                )
                box["boxes"][0].set_facecolor(GROUP_COLORS[group])
                box["boxes"][0].set_alpha(0.18)
                jitter = rng.uniform(-0.055, 0.055, size=values.size)
                ax.scatter(
                    position + jitter,
                    values,
                    s=18,
                    color=GROUP_COLORS[group],
                    alpha=0.65,
                    edgecolors="white",
                    linewidths=0.35,
                    zorder=3,
                )

            result = cohort_results.loc[gene] if gene in cohort_results.index else None
            q_value = result["adj_p_value"] if result is not None else float("nan")
            label = significance_label(q_value)
            gene_values = pd.to_numeric(expression.loc[gene, sample_ids], errors="coerce").dropna()
            ax.annotate(
                label,
                (gene_index, float(gene_values.max())),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
                color="#e31a1c" if label != "ns" and label != "ND" else "#718096",
                fontweight="bold",
            )

        ax.set_title(
            f"{cohort_info['label']}\n{int(group_counts.get('Control', 0))} Control, "
            f"{int(group_counts.get('IPF', 0))} IPF | {present_count}/{len(genes)} detected",
            fontsize=11,
            fontweight="bold",
        )
        ax.set_xticks(x_positions, genes, rotation=35, ha="right", rotation_mode="anchor")
        ax.set_ylabel(cohort_info["y_label"], fontweight="bold")
        ax.set_xlim(-0.65, len(genes) - 0.35)
        ax.yaxis.grid(True, linestyle=":", color="#d9d9d9", linewidth=0.9)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        all_values = np.asarray(expression_values)
        y_min = float(all_values.min())
        y_max = float(all_values.max())
        y_span = max(y_max - y_min, 0.25)
        ax.set_ylim(y_min - 0.04 * y_span, y_max + 0.18 * y_span)
        legend_handles = [
            Patch(
                facecolor=GROUP_COLORS[group],
                edgecolor=GROUP_COLORS[group],
                alpha=0.2,
                label=f"{group} (n={int(group_counts.get(group, 0))})",
            )
            for group in ("IPF", "Control")
        ]
        ax.legend(handles=legend_handles, loc="upper right", frameon=True, fontsize=8)
    fig.suptitle(
        "Validation Expression of Discovery Signature Genes",
        x=0.04,
        y=0.99,
        ha="left",
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.04,
        0.945,
        "IPF vs Control; stars show FDR-adjusted limma significance (* q<0.05, ** q<0.01, *** q<0.001; ns=not significant; ND=not detected).",
        ha="left",
        fontsize=9,
        color="#555555",
    )
    fig.tight_layout(rect=(0.03, 0.04, 0.98, 0.90))
    fig.savefig(OUTPUT_DIR / "validation_expression_boxplots.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "validation_expression_boxplots.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not INPUT_FILE.is_file():
        raise FileNotFoundError(
            f"Validation results not found: {INPUT_FILE}. "
            "Run scripts/13_external_validation.R first."
        )

    data = pd.read_csv(INPUT_FILE)
    missing_columns = REQUIRED_COLUMNS.difference(data.columns)
    if missing_columns:
        raise ValueError(
            "Validation results are missing required columns: "
            + ", ".join(sorted(missing_columns))
        )
    if data.empty:
        raise ValueError(f"Validation results contain no rows: {INPUT_FILE}")

    data["status"] = data.apply(validation_status, axis=1)
    genes = data["gene"].drop_duplicates().tolist()
    cohorts = data["cohort"].drop_duplicates().tolist()
    if data.duplicated(["gene", "cohort"]).any():
        raise ValueError("Validation results contain duplicate gene/cohort rows.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    create_expression_plot(data, genes)
    fig, ax = plt.subplots(figsize=(10, max(6.5, len(genes) * 0.48)))
    ax.set_xlim(0, len(cohorts))
    ax.set_ylim(len(genes), 0)

    for row_index, gene in enumerate(genes):
        for column_index, cohort in enumerate(cohorts):
            matches = data[(data["gene"] == gene) & (data["cohort"] == cohort)]
            if matches.empty:
                status = "Missing"
                q_value = None
            else:
                row = matches.iloc[0]
                status = row["status"]
                q_value = row["adj_p_value"]

            ax.add_patch(
                plt.Rectangle(
                    (column_index, row_index),
                    1,
                    1,
                    facecolor=STATUS_COLORS[status],
                    edgecolor="white",
                    linewidth=2,
                )
            )
            label = STATUS_LABELS[status]
            if q_value is not None and pd.notna(q_value):
                label += f"\nq={q_value:.1e}"
            text_color = "white" if status == "Validated" else "#202020"
            ax.text(
                column_index + 0.5,
                row_index + 0.5,
                label,
                ha="center",
                va="center",
                fontsize=8.5,
                color=text_color,
                fontweight="bold" if status == "Validated" else "normal",
            )

    ax.set_xticks([index + 0.5 for index in range(len(cohorts))], cohorts)
    ax.set_yticks([index + 0.5 for index in range(len(genes))], genes)
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, length=0, pad=10)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.xaxis.set_label_position("top")
    fig.suptitle(
        "12-Gene Signature: Independent Cohort Validation",
        fontsize=15,
        fontweight="bold",
        x=0.02,
        y=0.985,
        ha="left",
    )
    fig.text(
        0.02,
        0.94,
        "Strict validation requires concordant direction, FDR-adjusted p < 0.05, and AUC ≥ 0.70.",
        ha="left",
        va="bottom",
        fontsize=9.5,
        color="#444444",
    )
    for spine in ax.spines.values():
        spine.set_visible(False)

    legend_handles = [
        Patch(facecolor=STATUS_COLORS[status], label=STATUS_LABELS[status])
        for status in STATUS_ORDER
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        ncol=3,
        frameon=False,
        fontsize=8.5,
    )
    ax.set_position((0.12, 0.13, 0.86, 0.70))
    fig.savefig(OUTPUT_DIR / "gene_validation_status.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "gene_validation_status.pdf", bbox_inches="tight")
    plt.close(fig)

    counts = data.groupby("cohort")["status"].apply(lambda values: (values == "Validated").sum())
    print("Saved validation_expression_boxplots.png and validation_expression_boxplots.pdf")
    print("Saved gene_validation_status.png and gene_validation_status.pdf")
    print("Strictly validated genes by cohort:")
    for cohort in cohorts:
        print(f"  {cohort}: {counts[cohort]}/{len(genes)}")


if __name__ == "__main__":
    main()
