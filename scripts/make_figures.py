"""Generate the two owner-facing cluster-7 figures."""
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc

DATA = Path("data/pbmc3k.h5ad")
OUT = Path("results")
adata = sc.read_h5ad(DATA)
labels = adata.obs["leiden"].astype(str)

fig, ax = plt.subplots(figsize=(8, 6))
for cl in sorted(labels.unique()):
    mask = labels == cl
    ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1], s=12 if cl != "7" else 60, alpha=0.35 if cl != "7" else 1, label=f"Cluster {cl}", edgecolors="black" if cl == "7" else "none", linewidths=0.8)
ax.set_title("PBMC UMAP — cluster 7 highlighted (display only)")
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.legend(ncol=2, fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "cluster7_umap.png", dpi=180); plt.close(fig)

mask = (labels == "7").to_numpy()
cells = adata.obs_names[mask].astype(str)
gene = "MKI67"
idx = list(adata.var_names).index(gene)
expr = adata.X[mask, idx].toarray().ravel() if hasattr(adata.X[mask, idx], "toarray") else adata.X[mask, idx].ravel()
q = adata.obs.loc[mask]
fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True)
series = [(f"{gene} expression", expr), ("Genes detected", q["n_genes"]), ("Total counts", q["total_counts"]), ("Mitochondrial %", q["pct_mito"])]
for ax, (title, values) in zip(axes, series):
    ax.scatter(range(len(cells)), values, s=35)
    ax.set_ylabel(title); ax.grid(alpha=0.2)
axes[-1].set_xticks(range(len(cells))); axes[-1].set_xticklabels(cells, rotation=60, ha="right", fontsize=7)
fig.suptitle("Cluster 7 per-cell proliferation expression and quality", y=0.995)
fig.tight_layout(); fig.savefig(OUT / "cluster7_per_cell_expression_quality.png", dpi=180); plt.close(fig)
