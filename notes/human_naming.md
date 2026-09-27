# Human naming: primary k-means partition

Filled by the owner from marker genes before the answer key opens. Cluster IDs follow column `kmeans_primary` of `results/assignments/assignments_main.csv.gz`, largest cluster first. Marker tables: `results/tables/06_markers_kmeans_primary.csv`. Figures: `figures/06_dotplot_kmeans_primary.png`, `figures/06_umap_kmeans_primary.png`.

Lineage call, exactly one of: `B cells`, `Monocytes`, `NK cells`, `CD34+ progenitors`, `CD4 T cells`, `CD8 T cells`, `unresolved`.

Finer call, exactly one of the ten dataset page names or `none`: `CD19+ B Cells`, `CD14+ Monocytes`, `CD34+ Cells`, `CD4+ Helper T Cells`, `CD4+/CD25+ Regulatory T Cells`, `CD4+/CD45RA+/CD25- Naive T cells`, `CD4+/CD45RO+ Memory T Cells`, `CD56+ Natural Killer Cells`, `CD8+ Cytotoxic T cells`, `CD8+/CD45RA+ Naive Cytotoxic T Cells`.

Note is free text and may be left empty.

| cluster | cells | rule-based name | top 10 markers | lineage call | finer call | note |
|---|---|---|---|---|---|---|
| 0 | 91,134 | CD8 T cells | RPL3, RPS27, RPL5, RPS3, RPL15, RPS6, LTB, RPSA, MALAT1, RPL10 |  |  |  |
| 1 | 3,521 | Monocytes | CST3, S100A9, S100A8, TYROBP, FTH1, FTL, LYZ, S100A4, S100A6, AIF1 |  |  |  |
