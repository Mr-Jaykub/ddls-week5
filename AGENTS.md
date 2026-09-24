# AGENTS.md

## GOAL
Make a defensible **sequence-or-set-aside** call on cluster 7, naming its cell type and the marker genes supporting that identity. If the markers do not support a coherent identity, say so explicitly. Never claim novelty.

## DATA
Use `data/pbmc3k.h5ad`, loaded once at startup with `sc.read_h5ad`. Key meanings are in `spec.md` and `data/ABOUT_THIS_FILE.txt`. Every number from the owner is unverified: print `adata` and recompute per-cluster QC before using any of them. Use `.venv/bin/python`; never create an environment.

## MUST NOT
- Do not use the UMAP as evidence; show it only.
- Do not choose comparison clusters by map distance; compare cluster 7 against all 7 other clusters.
- Do not claim novelty.
- Do not report a cluster average without the corresponding per-cell view.
- Do not rely on a p-value alone.
- Do not check co-expression on `ad.X`; use `ad.layers["counts"]`.
- Do not reload the `.h5ad` inside an endpoint.
- Do not put a number in the app or README that is not in a `results/` file.
- Do not commit `.env` or `data/`.
- Do not overwrite `obs["leiden"]`.
- Do not claim to have seen the app render.
- Do not open or inspect the `.h5ad` outside the startup load, or bypass the documented loading approach.
- Do not omit per-cell consistency, quality measures, uncertainty, limitations, or direct comparisons.
- Do not ignore unequal cluster sizes, multiple testing, small clusters, duplicates, doublets, mixed profiles, or QC-driven separation.
- Do not present owner-supplied counts or QC summaries as verified before recomputing them.
- Do not infer unavailable donor, batch, treatment, collection, provenance, filtering, preprocessing, replication, cost, or follow-up cell-number facts.
- Do not automatically delete near-identical cells; investigate duplicates and document exclusions.
- Do not recommend a prospective surface marker without a coherent program and verified surface-protein counterpart.
- Do not deviate from the analysis stack specified in `spec.md`.
- Do not put scripts anywhere except `scripts/`.
- Do not put generated outputs anywhere except `results/`.

## VERSION CONTROL
Commit after every change that works. Stage only the files changed, by name; never use `git add -A` or `git add .`. Write a one-line commit message saying what changed and why, then show `git log --oneline -3`. Never commit `.env`, `data/`, or `.venv/`. Never push, and never rewrite history: no amend, rebase, or reset.

Everything detailed (criteria, methods, build contract, done-condition) is in spec.md. Read it before planning.
