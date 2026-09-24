# Specification: cluster 7 follow-up decision

## Objective

Produce a defensible **SEQUENCE** or **SET ASIDE** recommendation for cluster 7. The memo must name a supported cell type and its marker genes, or explicitly say **“no coherent identity”**; it must never claim novelty. A coherent identity with only a small effect is **SET ASIDE**, not SEQUENCE.

## Data and build contract

Use `data/pbmc3k.h5ad`, loaded with `sc.read_h5ad` via `.venv/bin/python`. Analysis scripts and the later app may load the file as needed; do not reload it inside an endpoint request. Print `adata` and recompute per-cluster QC before using any owner-supplied number; record verified values in `results/`. Never create an environment, overwrite `obs["leiden"]`, or recompute the supplied labels/UMAP.

The file contains expression in `ad.X`, raw UMI counts in `ad.layers["counts"]`, gene names in `ad.var_names`, cluster labels in `ad.obs["leiden"]`, QC columns in `ad.obs`, and the existing display coordinates in `ad.obsm["X_umap"]`. Full meanings are documented in `data/ABOUT_THIS_FILE.txt`.

The owner’s summary is not evidence until checked against the file. Preserve these owner-reported values as **[UNVERIFIED]** until recomputed: cluster 7 reportedly has 10 cells, 2,363 genes per cell, and 2.0% mitochondrial RNA; cluster 6 reportedly has 13 cells, 350 genes per cell, and 1.6% mitochondrial RNA; cluster 4 reportedly has 163 cells, 1,263 genes per cell, and 2.4% mitochondrial RNA; and cluster 0 reportedly has 1,197 cells. Report any mismatch explicitly. The owner selected cluster 4 as the named comparator, but it must still be compared with all clusters under the same rules. Unknowns include donor, batch, treatment, collection time, provenance, filtering, preprocessing, duplicate status, sampling history, and biological replication. The dataset may also contain ambient RNA, doublets, duplicated records, quality-driven separation, and cluster boundaries dependent on analysis settings. “One donor, no replicates” means observations are not independent biological replicates: cell-level p-values cannot establish population-level reproducibility, and uncertainty/generalisation must be stated conservatively.

Use the documented Scanpy/H5AD stack only. `README.md` is the final one-page memo. Put scripts in `scripts/`; put every table, verified number, statistical result, and figure in `results/`. The two required figures must be committed image files in `results/`. An app will be built later under separate app instructions; this specification defines the analysis and memo, not the app implementation.

## Questions and comparisons

1. Is cluster 7 coherent across its cells?
2. Is it distinct from **all seven** other clusters, with cluster 4 as the owner’s named comparator and without selecting or excluding comparisons by UMAP distance?
3. Are differences explained by `n_genes`, `total_counts`, `pct_mito`, capture, duplicates, ambient RNA, or doublets?
4. Does it match a known PBMC type, or is there no coherent identity?
5. Is it worth the one follow-up run?

## Pre-specified analysis

### Markers and multiplicity

Rank genes for all clusters, then use Scanpy `sc.tl.rank_genes_groups` as the single marker-ranking method for both the offline analysis and the app’s live marker table, with one pre-specified cluster-7-versus-rest analysis and identical settings. Use the documented rank-test method and library-size handling supported by the stack; do not fit a fragile negative-binomial model to ten cells. Correct that one gene-testing family with Benjamini–Hochberg FDR. Compare cluster 7 with each other cluster descriptively and with pre-specified effect/QC checks; do not create uncontrolled independent test families for every pair or plot.

For every reported marker, provide gene, comparison, detection prevalence in both groups, median expression, log2 fold-change, raw and BH-adjusted p-value, and a bootstrap uncertainty interval where supported; also report the number of genes tested and correction family. The app and analysis must use the same `rank_genes_groups` outputs/settings so a displayed marker means the same thing in both places, avoiding contradictory results, duplicated logic, and irreproducible review. Treat absolute log2FC below 0.25 as small, about 0.5 as moderate, and about 1 as sizeable. Adjusted p < 0.05 is screening evidence only; a stronger exploratory marker has about log2FC 1 and a 25–30 percentage-point prevalence difference.

### Coherence and identity

Define “expressed” before inspecting results as at least one raw UMI in `ad.layers["counts"]`. Report the raw-count distribution and marker prevalence to check that this threshold is supported rather than driven by transformed noise. Leading markers must be expressed in at least 8/10 cluster-7 cells, with no single cell driving the program.

Perform leave-one-cell-out and every relevant leave-two-cell-out check. The result **survives** if the same named identity and leading marker program remain supported in at least 80% of retained cells, leading effects do not reverse, and the decision is unchanged for every removal. If one or two removals eliminate the signal, call it unstable and SET ASIDE.

Identity controls are written before results: (a) use cluster 1 as the known-type control, draw 100 bootstrap resamples of its cells with replacement, and require the expected canonical PBMC type and marker program in at least 95/100 draws with no contradictory lineage; (b) use 100 random groups of 10 cells sampled from the full dataset as the negative control, and require “no coherent identity” in at least 95/100 draws. Failure of either control invalidates the pipeline and prevents naming cluster 7.

### Required traps and technical checks

- **Doublet:** define incompatible lineage marker sets in advance. A lineage is positive when at least two of its markers have ≥1 raw UMI in `ad.layers["counts"]`; flag cells positive for both lineages, especially with extreme `total_counts` or `n_genes`. Raw counts are required because co-expression is a molecule-detection question; `ad.X` is log-normalised. Report sets, thresholds, flagged cells, and interpretation. A standard doublet score is secondary only.
- **Cluster 6’s low gene count and mitochondrial context:** the owner reported **[UNVERIFIED]** 350 genes per cell and 1.6% mitochondrial RNA for cluster 6. Low detected-gene counts make its expression estimates fragile even when mitochondrial percentage is low; recompute both, show every cell, and keep cluster 6 exploratory. A few cells can dominate its estimates, so it cannot provide stable population evidence.
- **Separated does not mean novel:** UMAP is display only. Require markers, per-cell coherence, all-cluster comparison, QC/doublet checks, effect sizes, and robustness. Name a known type or say “no coherent identity”; never claim novelty.
- Check exact barcodes and repeated profiles for duplicates; investigate before removal, document exclusions, and rerun only the sensitivity analysis if records are removed. Do not delete merely similar cells.

### Nulls and distrust

Rank and pre-specify nulls: (1) an existing Scanpy doublet detector, `sc.pp.scrublet`, run on the raw-count layer as required by its API, with its simulated-doublet score and flagged cells; (2) a QC-matched non-cluster-7 comparison matched on `total_counts`, `n_genes`, and `pct_mito`; (3) a label-shuffle sensitivity check preserving the observed cluster sizes. `scrublet` is the doublet-mixture null: it simulates artificial doublets from the observed count profiles and compares them with observed cells, so do not invent a second bespoke synthetic-mixture null unless the function is unavailable. The QC-matched comparison tests quality explanation. Label shuffling is not a valid novelty null here because the observed clusters were produced from expression-derived structure: shuffling labels destroys the clustering relationship and asks only whether arbitrary labels produce a signal, not whether the existing cluster-7 assignment is biologically coherent. Use it only as a negative pipeline check, never as evidence for or against identity.

Distrust the result if identity changes after one or two removals, markers occur in a minority, one cell/QC extreme drives it, either primary null explains it, marker directions reverse, intervals include negligible effects, results depend on arbitrary filtering, duplicate/barcode checks fail, ambient RNA or mixed profiles are plausible, controls fail, or missing provenance prevents defensible interpretation.

## Sensitivity and labels

First inspect the file for unique cell identifiers, `obs["leiden"]` integrity, raw-count layer availability, QC columns, duplicate profiles/barcodes, and the existing UMAP. Analysis scripts may load the file directly with `sc.read_h5ad`; the app should load it once at startup, not per request. Run the core analysis first. Prefer a neighbour-purity check using the existing representation and labels: report, for each cluster-7 cell, the fraction of its existing-space neighbours carrying each label. Do not select comparisons by map distance. Only if that check is unavailable or inadequate may a separate-copy re-clustering sensitivity be attempted after the core analysis; never alter the original labels or UMAP. If neither is possible, state it as a limitation.

## Decision and alternative rules

**SEQUENCE** requires a coherent marker program across most cells, meaningful effects against all clusters, no persuasive QC/doublet/duplicate explanation, survival under removals, and explicit limitations. If identity is coherent but effects are small, unstable, technically explained, or not sufficiently reproducible, choose **SET ASIDE**. SET ASIDE does not disprove the identity.

Before seeing results, define “stronger alternative” as a better lead for the one follow-up run, not the largest, cleanest, or best-behaved cluster. Do not invent a numeric weighted score. First require the candidate to pass the same coherence bar as cluster 7: at least 80% of cells express the leading marker program, it remains stable under one- and two-cell removal, and it has no contradictory doublet signal. Then compare candidates using this owner-facing evidence ladder, in order: (1) coherent identity and marker support; (2) meaningful distinction from all clusters, including cluster 7; (3) technical explanations ruled out; (4) robustness and uncertainty; and only then (5) enough cells to make the proposed run informative. A large cluster cannot compensate for failing identity or coherence, and a small but coherent lead is not automatically disqualified. Choose at most one candidate only when it is better supported as a follow-up experiment under this ladder and has a stronger SEQUENCE/SET ASIDE case; otherwise report **“no stronger alternative identified.”**

## Deliverable and done-condition

`README.md` must contain:

1. Bold **SEQUENCE** or **SET ASIDE** and the named identity plus markers, or “no coherent identity”.
2. A table comparing cluster 7 with clusters 0–6: verified cell count, marker effects/adjusted results, `n_genes`, `total_counts`, and `pct_mito`.
3. Per-cell marker/QC views, doublet and duplicate checks, leave-one/two-out results, null results, controls, uncertainty, one-donor/no-replicate limitation, and distrust triggers.
4. At most one stronger alternative justified by the pre-specified evidence ladder, or the no-alternative statement.
5. Two committed images in `results/`: the UMAP as display only, and a per-cell expression/QC figure.
6. Reproducibility details: software/settings, exact cluster-7 identifiers, exclusions/reasons, verified-number tables, and analysis choices. Every README number must correspond to a `results/` file.

The work is done only when the memo, tables, figures, controls, sensitivity results, and reproducibility record support the decision without claiming novelty.
