from contextlib import asynccontextmanager
from pathlib import Path

import numpy as np
import scanpy as sc
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles

DATA_PATH = Path("data/pbmc3k.h5ad")
adata = None


def _cluster_mask(cluster: str):
    if cluster not in set(adata.obs["leiden"].astype(str)):
        raise HTTPException(status_code=404, detail=f"Unknown cluster: {cluster}")
    return adata.obs["leiden"].astype(str).to_numpy() == cluster


def _expression_vector(gene: str):
    matches = np.flatnonzero(adata.var_names.to_numpy() == gene)
    if len(matches) == 0:
        raise HTTPException(status_code=404, detail=f"Unknown gene: {gene}")
    values = adata[:, matches[0]].X
    return np.asarray(values.toarray() if hasattr(values, "toarray") else values).ravel()


@asynccontextmanager
async def lifespan(app: FastAPI):
    global adata
    adata = sc.read_h5ad(DATA_PATH)
    yield
    adata = None


app = FastAPI(title="PBMC cluster explorer", lifespan=lifespan)


@app.get("/", response_class=HTMLResponse)
def index():
    return HTMLResponse((Path(__file__).parent / "index.html").read_text())


@app.get("/api/meta")
def meta():
    return {"clusters": sorted(adata.obs["leiden"].astype(str).unique()), "genes": adata.var_names.tolist()}


@app.get("/api/umap")
def umap(color: str = Query("cluster")):
    coords = np.asarray(adata.obsm["X_umap"])
    labels = adata.obs["leiden"].astype(str).to_numpy()
    if color == "cluster":
        values = labels.tolist()
        numeric = False
    elif color in {"n_genes", "total_counts", "pct_mito"}:
        values = adata.obs[color].astype(float).tolist()
        numeric = True
    else:
        values = _expression_vector(color).tolist()
        numeric = True
    return {"color": color, "numeric": numeric, "points": [{"cell": str(cell), "x": float(x), "y": float(y), "cluster": label, "value": value} for cell, (x, y), label, value in zip(adata.obs_names, coords, labels, values)]}


@app.get("/api/expression/{gene}")
def expression(gene: str):
    values = _expression_vector(gene)
    labels = adata.obs["leiden"].astype(str).to_numpy()
    return {"gene": gene, "values": [{"cell": str(cell), "cluster": cluster, "value": float(value)} for cell, cluster, value in zip(adata.obs_names, labels, values)]}


@app.get("/api/cluster/{cluster}")
def cluster(cluster: str, gene: str | None = None):
    mask = _cluster_mask(cluster)
    idx = np.flatnonzero(mask)
    result = {
        "cluster": cluster,
        "quality": [{"cell": str(adata.obs_names[i]), "n_genes": float(adata.obs.iloc[i]["n_genes"]), "total_counts": float(adata.obs.iloc[i]["total_counts"]), "pct_mito": float(adata.obs.iloc[i]["pct_mito"])} for i in idx],
        "cells": [str(adata.obs_names[i]) for i in idx],
    }
    if gene:
        values = _expression_vector(gene)
        result["expression"] = [{"cell": str(adata.obs_names[i]), "value": float(values[i])} for i in idx]
    if "rank_genes_groups" not in adata.uns:
        sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=False, layer="counts")
    ranked = adata.uns["rank_genes_groups"]["names"][cluster][:20]
    scores = adata.uns["rank_genes_groups"].get("scores")
    logfoldchanges = adata.uns["rank_genes_groups"].get("logfoldchanges")
    pvals_adj = adata.uns["rank_genes_groups"].get("pvals_adj")
    marker_rows = []
    cluster_values = adata[mask, :].X
    rest_values = adata[~mask, :].X
    for i, gene in enumerate(ranked):
        gene = str(gene)
        gene_idx = int(np.flatnonzero(adata.var_names.to_numpy() == gene)[0])
        in_cluster = np.asarray(cluster_values[:, gene_idx].toarray() if hasattr(cluster_values[:, gene_idx], "toarray") else cluster_values[:, gene_idx]).ravel()
        in_rest = np.asarray(rest_values[:, gene_idx].toarray() if hasattr(rest_values[:, gene_idx], "toarray") else rest_values[:, gene_idx]).ravel()
        marker_rows.append({
            "gene": gene,
            "score": float(scores[gene][i]) if scores is not None else None,
            "logfoldchange": float(logfoldchanges[gene][i]) if logfoldchanges is not None else None,
            "pvals_adj": float(pvals_adj[gene][i]) if pvals_adj is not None else None,
            "cluster_detection_pct": float((in_cluster > 0).mean() * 100),
            "rest_detection_pct": float((in_rest > 0).mean() * 100),
            "cluster_median": float(np.median(in_cluster)),
            "rest_median": float(np.median(in_rest)),
            "pvals_adj": float(pvals_adj[gene][i]) if pvals_adj is not None else None,
        })
    result["top_markers"] = marker_rows
    return result
