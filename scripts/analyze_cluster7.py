"""Reproduce cluster-7 marker, lineage co-expression, and QC headline results."""
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

DATA = Path("data/pbmc3k.h5ad")
OUT = Path("results/cluster7_analysis")
OUT.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(DATA)
labels = adata.obs["leiden"].astype(str).to_numpy()
cluster = labels == "7"
counts = adata.layers["counts"]

# Markers use the app/specification's log-normalised ad.X.
sc.tl.rank_genes_groups(adata, "leiden", groups=["7"], reference="rest", method="wilcoxon", use_raw=False)
r = adata.uns["rank_genes_groups"]
ranked = pd.DataFrame({
    "gene": r["names"]["7"],
    "score": r["scores"]["7"],
    "logfoldchange": r["logfoldchanges"]["7"],
    "pvals_adj": r["pvals_adj"]["7"],
}).head(50)

# Use raw counts for detection prevalence and lineage co-expression.
genes = adata.var_names.to_numpy()
gene_index = {g: i for i, g in enumerate(genes)}
def present(marker_set):
    idx = [gene_index[g] for g in marker_set if g in gene_index]
    if not idx:
        return np.zeros(adata.n_obs, dtype=bool)
    return np.asarray((counts[:, idx] > 0).sum(axis=1)).ravel() >= 2

lineages = {
    "T_cell": {"CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R", "LTB"},
    "B_cell": {"MS4A1", "CD79A", "CD37", "HLA-DRA", "CD74"},
    "myeloid": {"LYZ", "S100A8", "S100A9", "FCN1", "CTSS"},
    "NK_cell": {"NKG7", "GNLY", "KLRD1", "TRBC1", "CD3D"},
    "platelet": {"PPBP", "PF4", "GNG11", "RGS18"},
}
flags = {name: present(markers) for name, markers in lineages.items()}
raw = pd.DataFrame({"cell": adata.obs_names.astype(str), "cluster": labels})
for name, values in flags.items():
    raw[f"{name}_positive"] = values
raw["cluster7"] = cluster
raw["cluster7_lineage_count"] = sum(values for values in flags.values())
raw["cluster7_double_lineage"] = raw["cluster7_lineage_count"] >= 2

qc = pd.DataFrame({
    "cell": adata.obs_names.astype(str), "cluster": labels,
    "n_genes": adata.obs["n_genes"].to_numpy(),
    "total_counts": adata.obs["total_counts"].to_numpy(),
    "pct_mito": adata.obs["pct_mito"].to_numpy(),
})

ranked.to_csv(OUT / "cluster7_markers.csv", index=False)
raw[raw["cluster7"]].to_csv(OUT / "cluster7_lineage_flags.csv", index=False)
qc[qc["cluster"] == "7"].to_csv(OUT / "cluster7_qc.csv", index=False)
pd.DataFrame([{
    "cluster7_cells": int(cluster.sum()),
    "marker_rows": len(ranked),
    "double_lineage_cells": int(raw.loc[cluster, "cluster7_double_lineage"].sum()),
    "double_lineage_fraction": float(raw.loc[cluster, "cluster7_double_lineage"].mean()),
    **{f"{name}_positive_cells": int(raw.loc[cluster, f"{name}_positive"].sum()) for name in lineages},
}]).to_csv(OUT / "cluster7_headline.csv", index=False)
print(f"Cluster 7 cells: {cluster.sum()}")
print(ranked.head(10).to_string(index=False))
print(raw[raw["cluster7"]][["cell", *[f"{n}_positive" for n in lineages], "cluster7_double_lineage"]].to_string(index=False))
