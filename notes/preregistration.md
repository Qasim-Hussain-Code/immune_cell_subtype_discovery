# Pre-registration

This file restates, in plain prose, the rules published on Days 62 to 66 of the series and fixed in `config.yaml`. It was committed before any data were downloaded. Where this file and a published post disagree, the post and `config.yaml` govern.

## 1. Question

Given only gene expression, can k-means, never shown a label, recover the immune populations defined by antibody-based bead enrichment? The labels are compared with the groupings only after every grouping has been frozen and committed.

## 2. Data

Zheng GXY et al. 2017. Massively parallel digital transcriptional profiling of single cells. Nature Communications 8:14049, doi `10.1038/ncomms14049`. Ten bead-enriched populations of peripheral blood mononuclear cells from one donor, released by 10x Genomics (Cell Ranger `1.1.0`, reference `hg19`, GemCode 3' v1 chemistry, licence CC BY 4.0). For each population the filtered gene-barcode matrix is used, that is, only the barcodes that Cell Ranger called as cells. No further cells are removed, no doublets are removed, no batch correction is applied and no reads are downsampled.

The populations were enriched with antibody-coated beads. Flow cytometry was used afterwards only to measure purity; the cells were not sorted by flow cytometry.

The ten populations and the six lineages they map to:

- CD19+ B Cells: B cells.
- CD14+ Monocytes: Monocytes.
- CD34+ Cells: CD34+ progenitors.
- CD4+ Helper T Cells, CD4+/CD25+ Regulatory T Cells, CD4+/CD45RA+/CD25- Naive T cells and CD4+/CD45RO+ Memory T Cells: CD4 T cells.
- CD56+ Natural Killer Cells: NK cells.
- CD8+ Cytotoxic T cells and CD8+/CD45RA+ Naive Cytotoxic T Cells: CD8 T cells.

Properties of the labels that are reported and not fixed:

- The helper population was enriched on CD4 alone, so it contains naive, memory and regulatory CD4 T cells.
- The naive cytotoxic population (CD8+/CD45RA+) is nested inside the cytotoxic population, and CD45RA+ CD8 cells include terminally differentiated effector memory cells that re-express CD45RA.
- CD45RA and CD45RO are splice isoforms of one gene, PTPRC, differing in exons 4 to 6 near the 5' end of the transcript. This 3' chemistry barely sees them. No isoform quantification is attempted.
- The CD34+ population was 45% pure by the dataset's own flow check.
- Each population was its own capture and sequencing run, so label and batch coincide. For that reason no batch correction is applied.
- The original authors downsampled reads and kept marker-matching clusters when building their references. Those label-using steps are not replicated.
- A 2020 re-analysis of these cells by the Stephens lab (single-cell-topics, fastTopics) informed prediction P4, as disclosed on Day 63.

## 3. Definitions

The adjusted Rand index (ARI) is `sklearn.metrics.adjusted_rand_score`.

For a label L and a cluster c, n(L, c) is the number of cells with label L in cluster c. The best cluster for a label is the cluster holding most of its cells, with ties going to the lowest cluster ID. Recovery of a label is the share of its cells that sit in its best cluster. Purity of a label is the share of its best cluster that carries the label.

The overlap of two labels A and B is the sum over clusters of the smaller of the two shares, n(A, c)/n(A) and n(B, c)/n(B). It is 0 when the two labels never share a cluster and 1 when they are spread across clusters identically.

The sub-ARI of two labels is the ARI between the partition and the labels, computed on only the cells carrying either label.

The truth of a cluster is its most common lineage if that lineage makes up at least half of the cluster, and "unresolved" otherwise.

Confidence intervals come from 1,000 resamples of cells with replacement (random stream 3). The same resamples are used for every arm and both keys, so every difference is paired. Intervals are the 2.5th and 97.5th percentiles. The groupings are held fixed while resampling, so an interval covers which cells were scored, not variation in the clustering itself; variation between random starts is reported separately.

The permutation null shuffles the labels 1,000 times (random stream 4). Reported for each arm and key: the null mean, its 2.5th and 97.5th percentiles, its maximum, and p = (1 + number of null values at or above the observed value) / 1001.

Rule-based naming uses log-normalised expression x (each cell scaled to 10,000 counts, then log1p, over all genes). For a gene g and a cluster c, d(g, c) is the mean of x inside c minus the mean of x outside c, divided by the standard deviation of x over all cells, and is set to 0 when that standard deviation is 0. A cluster is called T when the mean d over CD3D, CD3E and CD3G is above 0; it is then named CD8 T cells if its mean x of CD8A and CD8B (averaged) exceeds its mean x of CD4, and CD4 T cells otherwise. Any other cluster is scored for each lineage as the mean d over that lineage's markers (B cells: CD79A, CD79B, MS4A1, CD19; Monocytes: CD14, LYZ, S100A8, S100A9, FCN1; NK cells: GNLY, NKG7, KLRF1, NCAM1; CD34+ progenitors: CD34) and named by the highest score if that score is above 0, otherwise "unresolved". If a marker symbol is missing from the gene table, the run stops.

A name counts as correct when it equals the cluster's truth.

## 4. Predictions (Day 63)

The predictions concern ten labels, so they are scored on the run in which k-means was told to make ten groups (the label-informed k = 10 run). The same statistics for the primary k-means run and for Leiden are reported as description only.

- P1. B cells come back largely intact: recovery of the B cell label is at least 0.80.
- P2. Monocytes come back largely intact: recovery of the monocyte label is at least 0.80.
- P3. NK cells partly blur into cytotoxic T cells: the overlap of the NK and cytotoxic T labels is at least 0.10.
- P4. The CD34+ label mostly holds together: its recovery is at least 0.80.
- P5. The six T cell labels do not separate cleanly: none of them has a best cluster with both recovery and purity of at least 0.80.
- P6. The helper label cannot separate from the labels nested inside it: the sub-ARI of helper against naive CD4 T cells is no higher than the sub-ARI of naive against memory CD4 T cells.

## 5. Failure-mode checks (Day 66)

The published Day 66 post lists four of these checks, numbered 1, 2, 4 and 5; its point 2 is FM3 below. FM2 was declared on Day 63 (scoring against six lineages as well as ten labels) and is fixed in `config.yaml`, so all five checks run.

FM1, the run is the label. Whether gene expression beats the technical baseline (k-means on log10 total counts, log10 genes detected and percent mitochondrial counts, z-scored) under the paired bootstrap rule; the two sub-ARIs from P6 on the k = 10 run, which set a flag saying the runs are doing the sorting if P6 fails, with the primary and Leiden values as description; technical medians per population; and for CD34, its recovery under both technical baselines and the share of cells with at least one CD34 count inside and outside its best k = 10 cluster, over all cells and over CD34-labelled cells only.

FM2, the key's fault. Ten-label and six-lineage ARI side by side; the ten most overlapping label pairs for every arm; and, for each population, the share of cells in which each label gene (CD4, CD8A, CD8B, CD14, CD19, IL2RA, CD34, PTPRC, NCAM1) is detected.

FM3, counting versus grouping. The chosen k against 10; the ARI of the k = 10 run minus the primary run, with its interval; and the silhouette at the chosen k and at k = 10.

FM4, every cell gets a home. The share of cells with a negative silhouette and the silhouette of every cluster; stability over 20 reruns with seeds 101 to 120 for both k-means and Leiden; and the agreement between k-means and Leiden.

FM5, an unmatched group is not a discovery. Every cluster whose truth is "unresolved" in any arm, with its size, most common lineage and share, and top 10 marker genes; and every cluster the rule named "unresolved". No new cell types are claimed.

## 6. Order of work

Every grouping, the rule-based names and the second set of names are committed before the answer key is opened. The key is opened only by `scripts/common.py`, function `open_answer_key()`, which refuses unless it is called from script 07 or later, the freeze commit is on the history of HEAD, the frozen files still match `provenance/freeze_manifest.json`, the second set of names is complete and committed after the freeze, and the label file still matches its recorded hash.
