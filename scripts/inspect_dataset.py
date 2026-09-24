"""Inspect the supplied PBMC dataset and write a plain-text inspection report."""
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

    qc = pd.DataFrame(
        {
            "cluster": adata.obs["leiden"].astype(str).to_numpy(),
            "n_genes": n_genes,
            "total_counts": total_counts,
            "pct_mito": pct_mito,
        },
        index=adata.obs_names,
    )
    summary = qc.groupby("cluster").agg(
        n_cells=("cluster", "size"),
        n_genes_mean=("n_genes", "mean"),
        total_counts_mean=("total_counts", "mean"),
        pct_mito_mean=("pct_mito", "mean"),
    )

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
        "Recomputed per-cluster QC (from layers['counts'])",
        summary.to_string(),
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
