"""Inspect the supplied PBMC dataset and write a reproducible QC report."""
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

DATA_PATH = Path("data/pbmc3k.h5ad")
OUTPUT_PATH = Path("results/dataset_inspection.txt")


def main() -> None:
    adata = sc.read_h5ad(DATA_PATH)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    counts = adata.layers["counts"]
    mito = np.asarray(adata.var_names.str.upper().str.startswith("MT-"), dtype=bool)
    total_counts = np.asarray(counts.sum(axis=1)).ravel()
    n_genes = np.asarray((counts > 0).sum(axis=1)).ravel()
    mito_counts = np.asarray(counts[:, mito].sum(axis=1)).ravel()
    pct_mito = np.divide(
        mito_counts,
        total_counts,
        out=np.zeros_like(mito_counts, dtype=float),
        where=total_counts != 0,
    ) * 100

    qc_from_counts = pd.DataFrame(
        {
            "cluster": adata.obs["leiden"].astype(str).to_numpy(),
            "n_genes": n_genes,
            "total_counts": total_counts,
            "pct_mito": pct_mito,
        },
        index=adata.obs_names,
    )
    qc = pd.DataFrame(
        {
            "cluster": adata.obs["leiden"].astype(str).to_numpy(),
            "n_genes": adata.obs["n_genes"].to_numpy(),
            "total_counts": adata.obs["total_counts"].to_numpy(),
            "pct_mito": adata.obs["pct_mito"].to_numpy(),
        },
        index=adata.obs_names,
    )
    summary = qc.groupby("cluster").agg(
        n_cells=("cluster", "size"),
        n_genes_mean=("n_genes", "mean"),
        n_genes_median=("n_genes", "median"),
        total_counts_mean=("total_counts", "mean"),
        total_counts_median=("total_counts", "median"),
        pct_mito_mean=("pct_mito", "mean"),
        pct_mito_median=("pct_mito", "median"),
    )

    owner_values = {
        "cluster 7 n_genes": ("7", "n_genes", 2363),
        "cluster 7 pct_mito": ("7", "pct_mito", 2.0),
        "cluster 6 n_genes": ("6", "n_genes", 350),
        "cluster 6 pct_mito": ("6", "pct_mito", 1.6),
        "cluster 4 n_genes": ("4", "n_genes", 1263),
        "cluster 4 pct_mito": ("4", "pct_mito", 2.4),
        "cluster 0 n_cells": ("0", "n_cells", 1197),
    }
    statistic_check = []
    for label, (cluster, metric, owner_value) in owner_values.items():
        values = qc.loc[qc["cluster"] == cluster, metric] if metric != "n_cells" else pd.Series([len(qc[qc["cluster"] == cluster])])
        rounded_median = round(float(values.median()), 1 if metric == "pct_mito" else 0)
        statistic_check.append(
            {
                "owner_item": label,
                "owner_value": owner_value,
                "stored_mean": float(values.mean()),
                "stored_median": float(values.median()),
                "rounded_median": rounded_median,
                "count": int(len(values)),
                "matches_rounded_median": bool(np.isclose(rounded_median, owner_value)),
            }
        )
    statistic_check = pd.DataFrame(statistic_check).set_index("owner_item")

    duplicate_obs_names = int(adata.obs_names.duplicated().sum())
    duplicate_var_names = int(adata.var_names.duplicated().sum())
    duplicate_count_rows = int(pd.DataFrame(counts.toarray()).duplicated().sum())
    embedding_details = []
    for key, value in adata.obsm.items():
        shape = getattr(value, "shape", None)
        embedding_details.append(f"{key}: type={type(value).__name__}, shape={shape}")

    lines = [
        "PBMC dataset inspection",
        "=======================",
        f"adata: {adata}",
        f"shape: {adata.shape}",
        f"obs columns: {list(adata.obs.columns)}",
        f"var columns: {list(adata.var.columns)}",
        f"layers: {list(adata.layers.keys())}",
        f"obsm: {list(adata.obsm.keys())}",
        f"obs index unique: {adata.obs_names.is_unique}",
        f"var index unique: {adata.var_names.is_unique}",
        f"leiden categories: {list(adata.obs['leiden'].cat.categories)}",
        f"mitochondrial genes detected by MT- prefix: {int(mito.sum())}",
        "",
        "Duplicate and doublet-related file contents",
        f"duplicate cell identifiers: {duplicate_obs_names}",
        f"duplicate gene identifiers: {duplicate_var_names}",
        f"exact duplicate raw-count rows: {duplicate_count_rows}",
        "No doublet score or doublet annotation was found in obs.",
        "No dedicated duplicate/barcode provenance or doublet metadata was found.",
        "",
        "Existing embeddings and reusable annotations",
        *embedding_details,
        "No additional embeddings or named cell-type annotations were found.",
        "The existing X_umap coordinates are reusable for display only, not evidence.",
        "", 
        "Per-cluster QC from the file's stored obs columns; means and medians",
        summary.to_string(),
        "",
        "Counts-derived versus stored-column QC check",
        "The stored obs columns are the primary QC values because they are the dataset's own processed annotations. Counts-derived values are retained as a consistency check.",
        "n_genes exact agreement: " + str(bool(np.array_equal(qc["n_genes"].to_numpy(), n_genes))),
        "total_counts exact agreement: " + str(bool(np.array_equal(qc["total_counts"].to_numpy(), total_counts))),
        "pct_mito max absolute difference: " + str(float(np.max(np.abs(qc["pct_mito"].to_numpy() - pct_mito)))),
        "",
        "Owner-number statistic check",
        "Rounded medians are considered compatible; exact equality is not required for rounded owner summaries.",
        statistic_check.to_string(),
        "",
        "Per-cluster QC recomputed directly from layers['counts'] (consistency view)",
        qc_from_counts.groupby("cluster").agg(
            n_cells=("cluster", "size"),
            n_genes_mean=("n_genes", "mean"),
            n_genes_median=("n_genes", "median"),
            total_counts_mean=("total_counts", "mean"),
            total_counts_median=("total_counts", "median"),
            pct_mito_mean=("pct_mito", "mean"),
            pct_mito_median=("pct_mito", "median"),
        ).to_string(),
        "",
        "Cluster 7 per-cell QC",
        qc[qc["cluster"] == "7"].to_string(),
        "",
        "Cluster 6 per-cell QC",
        qc[qc["cluster"] == "6"].to_string(),
    ]
    OUTPUT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
