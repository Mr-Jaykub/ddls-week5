# Specification: Cluster 7 follow-up decision

## Objective

Analyze the supplied processed human PBMC single-cell RNA-seq dataset and produce a defensible recommendation for Alex Rios and the PI: **sequence cluster 7 in an independent follow-up run, or set it aside**. The work should determine whether cluster 7 is a coherent, biologically meaningful signal rather than a technical or sampling artefact. It must not present cluster 7 as proven novel.

## Dataset contract

Use `data/pbmc3k.h5ad`, loaded once at startup with `sc.read_h5ad`, using `.venv/bin/python`. At startup, print `adata` and recompute per-cluster QC before using any owner-supplied number. Do not create another environment. Use the existing `pbmc3k.h5ad` file only through its documented contents:

- 2,700 human PBMC cells in clusters `0`–`7`.
- Rows conceptually represent cells and columns genes.
- `ad.X`: log-normalised expression for 13,714 genes.
- `ad.layers["counts"]`: raw UMI counts.
- `ad.var_names`: gene symbols.
- `ad.obs["leiden"]`: numeric cluster labels (`"0"` through `"7"`).
- `ad.obs["n_genes"]`, `ad.obs["total_counts"]`, `ad.obs["pct_mito"]`: per-cell QC.
- `ad.obsm["X_umap"]`: existing two-dimensional coordinates.
- Do not recompute the existing clusters or map.

The summary does not provide donor, batch, treatment, collection-time, provenance, filtering, preprocessing, duplicate, or biological-replicate information. These are limitations, not facts to infer.

## Primary scientific questions

1. Is cluster 7’s expression program coherent across its 10 cells?
2. Is it distinct from clusters 0–6, particularly its nearest UMAP neighbours and QC-matched cells?
3. Are differences explained by `n_genes`, `total_counts`, or `pct_mito` rather than biology?
4. Does cluster 7 resemble a known PBMC population, or is no reliable identity supported?
5. Is there evidence strong enough to justify one follow-up sequencing run?

Cluster 4 is an optional comparator: it has 163 cells, 1,263 genes/cell, and 2.4% mitochondrial RNA. The original cluster-4 novelty question is secondary to the cluster-7 decision.

## Analysis stack and build contract

Use the existing Python virtual environment and the documented Scanpy/H5AD workflow; do not substitute another analysis stack. Load the dataset once at startup, print `adata`, recompute per-cluster QC, and preserve a reproducible record of software/settings, exact cluster-7 cell identifiers, exclusions and reasons, and analysis choices. Put scripts only in `scripts/` and generated outputs only in `results/`.

## Analysis requirements

### Expression and differential comparisons

- Rank top genes for every cluster to describe expression programs.
- Compare cluster 7 with every other cluster, prioritizing UMAP-nearest groups.
- Use normalisation/library-size adjustment appropriate to the expression data and a method that handles unequal cluster sizes.
- Report marker detection prevalence, median expression, log2 fold-change, uncertainty/confidence intervals where possible, and multiple-testing-adjusted results.
- Treat absolute log2FC < 0.25 as small, approximately 0.5 as moderate, and approximately 1.0 as sizeable.
- Use adjusted p < 0.05 only as an exploratory screening criterion, not proof. For a stronger exploratory signal, look for approximately log2FC 1 and a 25–30 percentage-point detection-prevalence difference.
- Compare against QC-matched cells or otherwise assess whether library size and quality explain apparent separation.
- Include a whole-dataset baseline across all 2,700 cells and, where possible, a well-annotated same-species/comparable-assay PBMC reference. Reference matching is a naming aid, not proof of identity.
- If using a null comparison, pre-specify it (for example, label shuffling or matched-cell comparisons).

### Coherence and robustness

- Evaluate individual cells rather than cluster averages alone.
- As a working criterion, leading markers should be present above a defined threshold in at least 8 of 10 cluster-7 cells, with no one cell driving the result. Report the actual threshold and results.
- Perform leave-one-cell-out analysis. If removal of one or two cells eliminates the signal, classify it as unstable.
- Inspect all cluster-7 cells for mixed profiles, doublet-like behaviour, unusually high RNA, and extreme QC.
- Check exact barcodes/records and repeated expression profiles for duplicates if the source data permit it. Investigate before removal; document any removal and rerun clustering if records are removed. Do not automatically remove merely similar cells.
- Conduct reasonable sensitivity checks for filtering, dimensions, and clustering settings where supported by the available data.

### Quality assessment

Compare cluster 7 with every relevant cluster and the full dataset for:

- cell count;
- `n_genes`;
- `total_counts`;
- `pct_mito`;
- distributions and individual-cell extremes, not merely means.

Cluster 7’s known summary values are 10 cells, 2,363 genes/cell, and 2.0% mitochondrial RNA. Treat its high gene count as potentially biological, high capture, or doublet-related until checked. Clusters 6 (13 cells; 350 genes/cell; 1.6% mitochondrial RNA) and 7 are especially uncertain; cluster 5 is also exploratory. Cluster 0 has 1,197 cells, illustrating the uneven sizes.

## Decision rule

Recommend **sequence** only if all of the following are reasonably supported:

- a coherent marker program is present across most cluster-7 cells;
- the distinction from relevant neighbouring/QC-matched clusters is meaningful, not merely statistically significant;
- quality metrics, high RNA capture, duplicates, doublets, or mixed profiles do not plausibly explain it;
- the conclusion is not overturned by leave-one-cell-out or reasonable sensitivity checks;
- the limitations of 10 cells, absent biological replicates, and incomplete provenance are clearly acknowledged.

Otherwise recommend **set aside**. This means the cluster is not sufficiently defensible for the one run, not that it is disproven. If it matches a known cell type, name it only with appropriate evidence; sequence only if it remains unusually coherent or biologically important, not merely because it is known.

Select one alternative only if it beats cluster 7 under the same criteria: coherent markers, meaningful separation, adequate cell number, non-explanatory QC, and fewer unresolved technical concerns. Cluster 4 is a plausible candidate, not an automatic alternative. If none qualifies, state “no stronger alternative identified.”

## Done-condition

The work is done only when the memo contains a defensible sequence-or-set-aside decision for cluster 7, names its supported cell type and marker genes (or explicitly says no coherent identity is supported), compares it with the other clusters and relevant baselines, reports per-cell coherence and QC/technical checks, states uncertainty and missing provenance, and includes the reproducibility record and required figures. It must distinguish “worth testing” from “proven to be a new population.”

## Deliverable format

Produce a one-page budget-meeting memo with:

1. **Decision:** bold “SEQUENCE” or “SET ASIDE”.
2. **Evidence table:** cluster 7 versus clusters 0–6, including cell count, marker effects, adjusted results, `n_genes`, `total_counts`, and `pct_mito`.
3. **Biological evidence:** top marker programs, nearest-cluster comparisons, known-reference comparison if available, and observed versus interpreted findings clearly separated.
4. **Consistency/quality:** per-cell marker prevalence, QC distributions, duplicate/doublet checks, and leave-one-cell-out results.
5. **Uncertainty and limitations:** effect-size uncertainty, small-cell-count caveat, missing metadata/provenance, and lack of replication.
6. **Alternative:** one stronger alternative with rationale, or the explicit no-alternative statement.
7. **Figures:** an existing-map figure highlighting cluster 7, plus a per-cell expression/QC figure showing whether the signal is consistent across all 10 cells.
8. **Reproducibility appendix/record:** software and settings, exact cluster-7 cell identifiers, exclusions and their reasons, and analysis choices.

The conclusion must say whether cluster 7 is **worth testing**, not claim that it is **proven to be a new population**. Do not specify a prospective surface marker or lab-isolation strategy without first establishing a coherent program and confirming usable surface-protein counterparts.
