"""Cluster-7 identity, proliferation, doublet, controls, and stability analysis."""
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
genes = adata.var_names.to_numpy()
idx = {g: i for i, g in enumerate(genes)}

def vals(matrix, cols):
    a = matrix[:, cols]
    return np.asarray(a.toarray() if hasattr(a, "toarray") else a)

def raw_positive(marker_set, min_genes=2):
    cols = [idx[g] for g in marker_set if g in idx]
    return np.asarray((vals(counts, cols) > 0).sum(axis=1)).ravel() >= min_genes if cols else np.zeros(adata.n_obs, bool)

def log_score(marker_set):
    cols = [idx[g] for g in marker_set if g in idx]
    return vals(X, cols).mean(axis=1) if cols else np.zeros(adata.n_obs)

# Single shared marker ranking on log-normalised X.
sc.tl.rank_genes_groups(adata, "leiden", groups=["7"], reference="rest", method="wilcoxon", use_raw=False)
r = adata.uns["rank_genes_groups"]
markers = pd.DataFrame({"gene": r["names"]["7"], "score": r["scores"]["7"], "logfoldchange": r["logfoldchanges"]["7"], "pvals_adj": r["pvals_adj"]["7"]})
cluster_values, rest_values = X[cluster7, :], X[~cluster7, :]
markers["cluster_detection_pct"] = [float((vals(cluster_values, [idx[g]])[:, 0] > 0).mean()*100) for g in markers.gene]
markers["rest_detection_pct"] = [float((vals(rest_values, [idx[g]])[:, 0] > 0).mean()*100) for g in markers.gene]
markers["cluster_median"] = [float(np.median(vals(cluster_values, [idx[g]])[:, 0])) for g in markers.gene]
markers["rest_median"] = [float(np.median(vals(rest_values, [idx[g]])[:, 0])) for g in markers.gene]
# Common/depth-sensitive genes are not treated as identity markers.
markers["rest_detection_high"] = markers.rest_detection_pct >= 50
markers["identity_candidate"] = (markers.logfoldchange.abs() >= 1) & ((markers.cluster_detection_pct - markers.rest_detection_pct).abs() >= 25) & (markers.pvals_adj < .05) & ~markers.rest_detection_high
markers.to_csv(OUT / "cluster7_markers.csv", index=False)
markers[markers.identity_candidate].to_csv(OUT / "cluster7_identity_candidates.csv", index=False)

# Strict lineage sets: only canonical markers assigned to one lineage here.
lineages = {
    "T_cell": {"TRAC", "TRBC1", "TRBC2", "CD3D", "CD3E"},
    "B_cell": {"MS4A1", "CD79A", "CD37"},
    "myeloid": {"LYZ", "FCN1", "S100A8", "S100A9"},
    "NK_cell": {"NKG7", "GNLY", "KLRD1"},
    "platelet": {"PPBP", "PF4", "RGS18"},
}
flags = {n: raw_positive(s) for n, s in lineages.items()}
scores = {n: log_score(s) for n, s in lineages.items()}
lineage = pd.DataFrame({"cell": adata.obs_names.astype(str), "cluster": labels})
for n in lineages:
    lineage[f"{n}_positive"] = flags[n]
    lineage[f"{n}_score"] = scores[n]
lineage["lineage_count"] = sum(flags.values())
lineage["multi_lineage"] = lineage.lineage_count >= 2
lineage["assigned_lineage"] = lineage[[f"{n}_score" for n in lineages]].idxmax(axis=1).str.replace("_score", "", regex=False)

qc = pd.DataFrame({"cell": adata.obs_names.astype(str), "cluster": labels, "n_genes": adata.obs.n_genes.to_numpy(), "total_counts": adata.obs.total_counts.to_numpy(), "pct_mito": adata.obs.pct_mito.to_numpy()})
lineage = lineage.join(qc.set_index("cell"), on="cell", rsuffix="_qc")
lineage.to_csv(OUT / "lineage_flags_all_cells.csv", index=False)
lineage[lineage.cluster == "7"].to_csv(OUT / "cluster7_per_cell_lineage.csv", index=False)

# DNA replication/cell-cycle genes, prevalence across clusters and QC-matched cells.
replication = {"MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "TYMS", "RRM1", "RRM2", "STMN1", "TK1", "PCLAF", "UBE2C", "TOP2A", "MKI67", "BIRC5"}
rep = pd.DataFrame({"gene": sorted(replication & set(idx))})
for cl in sorted(set(labels)):
    mask = labels == cl
    rep[f"cluster_{cl}_detection_pct"] = [float((vals(counts[mask, :], [idx[g]])[:, 0] > 0).mean()*100) for g in rep.gene]
# QC-matched background: 20 nearest non-7 cells per cluster-7 cell.
f = qc[["n_genes", "total_counts", "pct_mito"]].to_numpy(float); scale = f.std(axis=0); scale[scale == 0] = 1
z = f / scale; non7 = np.flatnonzero(~cluster7); matched_idx = []
for i in np.flatnonzero(cluster7):
    d = np.sqrt(((z[non7] - z[i])**2).sum(axis=1)); matched_idx.extend(non7[np.argsort(d)[:20]])
matched_idx = np.asarray(matched_idx)
rep["qc_matched_detection_pct"] = [float((vals(counts[matched_idx, :], [idx[g]])[:, 0] > 0).mean()*100) for g in rep.gene]
rep.to_csv(OUT / "replication_marker_prevalence.csv", index=False)

# Scrublet on raw counts.
scr = sc.AnnData(X=counts.copy(), obs=adata.obs.copy(), var=adata.var.copy())
try:
    sc.pp.scrublet(scr, random_state=0, threshold=0.25)
    scrub = pd.DataFrame({"cell": adata.obs_names.astype(str), "cluster": labels, "doublet_score": scr.obs.doublet_score.to_numpy(), "predicted_doublet": scr.obs.predicted_doublet.to_numpy()})
    scrub_status = "completed"
except Exception as exc:
    scrub = pd.DataFrame({"cell": adata.obs_names.astype(str), "cluster": labels, "error": [str(exc)]*adata.n_obs}); scrub_status = f"failed: {exc}"
scrub.to_csv(OUT / "scrublet_all_cells.csv", index=False)

# Stability: same identity proxy = >=8 cells positive for the same dominant strict lineage and >=8 cells positive for >=2 replication genes.
def stability(keep):
    m = cluster7.copy(); m[np.flatnonzero(cluster7)[~keep]] = False
    sub = lineage[m]
    dominant = sub.assigned_lineage.value_counts().index[0] if len(sub) else "none"
    dominant_n = int((sub.assigned_lineage == dominant).sum()) if len(sub) else 0
    rep_positive = np.asarray((vals(counts[m, :], [idx[g] for g in replication if g in idx]) > 0).sum(axis=1)).ravel() >= 2
    return dominant, dominant_n, int(rep_positive.sum()), int(m.sum())
base = np.ones(cluster7.sum(), dtype=bool)
rows = [{"removed": "none", "removed_n": 0, "dominant_lineage": stability(base)[0], "dominant_n": stability(base)[1], "replication_positive_n": stability(base)[2], "retained_n": stability(base)[3]}]
for n in (1, 2):
    from itertools import combinations
    for rem in combinations(range(cluster7.sum()), n):
        keep = base.copy(); keep[list(rem)] = False; s = stability(keep)
        rows.append({"removed": ",".join(map(str, rem)), "removed_n": n, "dominant_lineage": s[0], "dominant_n": s[1], "replication_positive_n": s[2], "retained_n": s[3]})
pd.DataFrame(rows).to_csv(OUT / "cluster7_leave_out_stability.csv", index=False)

# Controls specified in spec: known-type cluster 1 and heterogeneous random groups.
rng = np.random.default_rng(0); controls = []
for draw in range(100):
    sample = rng.choice(np.flatnonzero(labels == "1"), size=max(1, (labels == "1").sum()), replace=True)
    t = lineage.iloc[sample]; controls.append({"control":"known_cluster_1", "draw":draw, "coherent": bool(t.T_cell_positive.mean() >= .8 and t.lineage_count.mean() < 2)})
for draw in range(100):
    sample = rng.choice(np.arange(len(labels)), size=10, replace=False)
    t = lineage.iloc[sample]; controls.append({"control":"random_10", "draw":draw, "coherent": bool(t.lineage_count.mean() >= 2)})
pd.DataFrame(controls).to_csv(OUT / "control_draws.csv", index=False)

headline = {
    "cluster7_cells": 10, "identity_candidates": int(markers.identity_candidate.sum()), "cluster7_multi_lineage_pct": float(lineage.loc[cluster7, "multi_lineage"].mean()*100), "scrublet_status": scrub_status,
    "cluster7_scrublet_doublets": int(scrub.loc[cluster7, "predicted_doublet"].sum()) if "predicted_doublet" in scrub else -1,
    "known_control_pass_rate": float(pd.DataFrame(controls).query("control=='known_cluster_1'").coherent.mean()),
    "negative_control_no_identity_rate": float((~pd.DataFrame(controls).query("control=='random_10'").coherent).mean()),
}
pd.DataFrame([headline]).to_csv(OUT / "cluster7_headline.csv", index=False)
print(pd.Series(headline).to_string())
