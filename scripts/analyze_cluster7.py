"""Reproduce cluster-7 marker, calibrated doublet, proliferation, and QC results."""
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

DATA = Path("data/pbmc3k.h5ad")
OUT = Path("results/cluster7_analysis")
OUT.mkdir(parents=True, exist_ok=True)
adata = sc.read_h5ad(DATA)
labels = adata.obs["leiden"].astype(str).to_numpy()
cluster7 = labels == "7"
counts = adata.layers["counts"]
X = adata.X

def dense_col(matrix, idx):
    v = matrix[:, idx]
    return np.asarray(v.toarray() if hasattr(v, "toarray") else v).ravel()

genes = adata.var_names.to_numpy()
idx = {g: i for i, g in enumerate(genes)}

def raw_positive(marker_set, min_genes=2):
    cols = [idx[g] for g in marker_set if g in idx]
    return np.asarray((counts[:, cols] > 0).sum(axis=1)).ravel() >= min_genes if cols else np.zeros(adata.n_obs, bool)

def log_score(marker_set):
    cols = [idx[g] for g in marker_set if g in idx]
    return np.asarray(X[:, cols].mean(axis=1)).ravel() if cols else np.zeros(adata.n_obs)

# Marker ranking uses log-normalised ad.X, as used by the app.
sc.tl.rank_genes_groups(adata, "leiden", groups=["7"], reference="rest", method="wilcoxon", use_raw=False)
r = adata.uns["rank_genes_groups"]
markers = pd.DataFrame({"gene": r["names"]["7"], "score": r["scores"]["7"], "logfoldchange": r["logfoldchanges"]["7"], "pvals_adj": r["pvals_adj"]["7"]})
cluster_values = X[cluster7, :]
rest_values = X[~cluster7, :]
markers["cluster_detection_pct"] = [float((dense_col(cluster_values, idx[g]) > 0).mean() * 100) for g in markers.gene]
markers["rest_detection_pct"] = [float((dense_col(rest_values, idx[g]) > 0).mean() * 100) for g in markers.gene]
markers["cluster_median"] = [float(np.median(dense_col(cluster_values, idx[g]))) for g in markers.gene]
markers["rest_median"] = [float(np.median(dense_col(rest_values, idx[g]))) for g in markers.gene]
# Owner screening: meaningful effect and substantial prevalence, not ubiquitous housekeeping.
markers["passes_owner_screen"] = (markers.logfoldchange.abs() >= 1) & ((markers.cluster_detection_pct - markers.rest_detection_pct).abs() >= 25) & (markers.pvals_adj < .05)
markers.to_csv(OUT / "cluster7_markers.csv", index=False)
markers[markers.passes_owner_screen].head(50).to_csv(OUT / "cluster7_markers_owner_screen.csv", index=False)

# Mutually specific lineage sets: no marker appears in another set.
lineages = {
    "T_cell": {"TRAC", "CD3D", "CD3E", "TRBC1", "TRBC2"},
    "B_cell": {"MS4A1", "CD79A", "CD37", "CD74", "HLA-DRA"},
    "myeloid": {"LYZ", "FCN1", "S100A8", "S100A9", "CTSS"},
    "NK_cell": {"NKG7", "GNLY", "KLRD1", "FCER1G", "TYROBP"},
    "platelet": {"PPBP", "PF4", "RGS18", "GNG11"},
}
lineage_flags = {name: raw_positive(markers_set) for name, markers_set in lineages.items()}
lineage_scores = {name: log_score(markers_set) for name, markers_set in lineages.items()}
flag_table = pd.DataFrame({"cell": adata.obs_names.astype(str), "cluster": labels})
for name in lineages:
    flag_table[f"{name}_positive"] = lineage_flags[name]
    flag_table[f"{name}_score"] = lineage_scores[name]
flag_table["lineage_count"] = sum(lineage_flags.values())
flag_table["multi_lineage"] = flag_table.lineage_count >= 2
flag_table["cluster7"] = cluster7
flag_table.to_csv(OUT / "lineage_flags_all_cells.csv", index=False)

# Calibrate cluster 7 against cluster-level and QC-matched backgrounds.
qc = pd.DataFrame({"cell": adata.obs_names.astype(str), "cluster": labels, "n_genes": adata.obs.n_genes.to_numpy(), "total_counts": adata.obs.total_counts.to_numpy(), "pct_mito": adata.obs.pct_mito.to_numpy()})
flag_table = flag_table.join(qc.set_index("cell"), on="cell", rsuffix="_qc")
cluster_calibration = flag_table.groupby("cluster").agg(n_cells=("cell", "size"), multi_lineage_pct=("multi_lineage", "mean"), median_lineage_count=("lineage_count", "median"), median_n_genes=("n_genes", "median"), median_total_counts=("total_counts", "median"), median_pct_mito=("pct_mito", "median"))
cluster_calibration.to_csv(OUT / "lineage_calibration_by_cluster.csv")
# Match each cluster-7 cell to up to 20 non-7 cells by standardized QC distance.
features = qc[["n_genes", "total_counts", "pct_mito"]].to_numpy(float)
scale = features.std(axis=0); scale[scale == 0] = 1
z = features / scale
matched = []
for i in np.flatnonzero(cluster7):
    candidates = np.flatnonzero(~cluster7)
    distances = np.sqrt(((z[candidates] - z[i]) ** 2).sum(axis=1))
    take = candidates[np.argsort(distances)[:20]]
    matched.append(flag_table.iloc[take])
matched = pd.concat(matched, ignore_index=True)
match_summary = pd.DataFrame([{"cluster7_multi_lineage_pct": float(flag_table.loc[cluster7, "multi_lineage"].mean() * 100), "qc_matched_multi_lineage_pct": float(matched.multi_lineage.mean() * 100), "qc_matched_cells": len(matched)}])
match_summary.to_csv(OUT / "qc_matched_lineage_calibration.csv", index=False)

# Scrublet on raw counts.
scr = sc.AnnData(X=counts.copy(), obs=adata.obs.copy(), var=adata.var.copy())
try:
    sc.pp.scrublet(scr, random_state=0, threshold=0.25)
    scrub = pd.DataFrame({"cell": adata.obs_names.astype(str), "cluster": labels, "doublet_score": scr.obs["doublet_score"].to_numpy(), "predicted_doublet": scr.obs["predicted_doublet"].to_numpy()})
    scrub_status = "completed"
except Exception as exc:
    scrub = pd.DataFrame({"cell": adata.obs_names.astype(str), "cluster": labels, "error": [str(exc)] * adata.n_obs})
    scrub_status = f"failed: {exc}"
scrub.to_csv(OUT / "scrublet_all_cells.csv", index=False)

# Proliferation and per-cell lineage assignment.
cell_cycle = {"MKI67", "TOP2A", "TYMS", "PCNA", "STMN1", "TUBA1B", "HMGB2", "UBE2C", "BIRC5"}
per_cell = flag_table.loc[cluster7, ["cell", "cluster", "n_genes", "total_counts", "pct_mito", "lineage_count", "multi_lineage", *[f"{n}_positive" for n in lineages], *[f"{n}_score" for n in lineages]]].copy()
per_cell["cell_cycle_score"] = log_score(cell_cycle)[cluster7]
per_cell["cell_cycle_positive"] = raw_positive(cell_cycle, min_genes=2)[cluster7]
per_cell["dominant_lineage"] = per_cell[[f"{n}_score" for n in lineages]].idxmax(axis=1).str.replace("_score", "", regex=False)
per_cell.to_csv(OUT / "cluster7_per_cell_lineage_cycle_qc.csv", index=False)

headline = {"cluster7_cells": int(cluster7.sum()), "owner_screen_markers": int(markers.passes_owner_screen.sum()), "cluster7_multi_lineage_pct": float(per_cell.multi_lineage.mean()*100), "all_cluster_multi_lineage_pct_min": float(cluster_calibration.multi_lineage_pct.min()*100), "all_cluster_multi_lineage_pct_max": float(cluster_calibration.multi_lineage_pct.max()*100), "cluster7_cycle_positive": int(per_cell.cell_cycle_positive.sum()), "scrublet_status": scrub_status}
if "predicted_doublet" in scrub:
    headline.update({"cluster7_scrublet_doublets": int(scrub.loc[cluster7, "predicted_doublet"].sum()), "cluster7_scrublet_doublet_pct": float(scrub.loc[cluster7, "predicted_doublet"].mean()*100)})
pd.DataFrame([headline]).to_csv(OUT / "cluster7_headline.csv", index=False)
print(pd.Series(headline).to_string())
print("Owner-screen markers:")
print(markers[markers.passes_owner_screen].head(20).to_string(index=False))
