# Blind report

Written by script 06 before the answer key opens. No population label has been read.

Cells: 94,655. Highly variable genes: 999. Variance held by 50 components: 0.1947.

## Number of groups

Chosen k: 2 (mean silhouette 0.6336).

| k | mean silhouette | inertia (not used) |
|---|---|---|
| 2 | 0.6336 | 15850200.0 |
| 3 | 0.5116 | 13480990.0 |
| 4 | 0.4284 | 11480187.0 |
| 5 | 0.4500 | 9847018.0 |
| 6 | 0.4583 | 9324436.0 |
| 7 | 0.2647 | 8881436.0 |
| 8 | 0.2643 | 8638304.0 |
| 9 | 0.2095 | 8415227.0 |
| 10 | 0.2129 | 8204635.0 |
| 11 | 0.2714 | 8017030.0 |
| 12 | 0.2172 | 7793884.5 |
| 13 | 0.2286 | 7638136.5 |
| 14 | 0.2303 | 7488932.5 |
| 15 | 0.2195 | 7370553.0 |
| 16 | 0.2319 | 7217894.0 |
| 17 | 0.2325 | 7068544.5 |
| 18 | 0.2154 | 6982489.5 |
| 19 | 0.2165 | 6862709.0 |
| 20 | 0.2171 | 6805211.0 |

Chosen Leiden resolution: 0.8, giving 23 clusters (mean silhouette 0.1970).

| resolution | clusters | mean silhouette |
|---|---|---|
| 0.1 | 10 | 0.1890 |
| 0.2 | 12 | 0.1877 |
| 0.3 | 16 | 0.1861 |
| 0.4 | 18 | 0.1894 |
| 0.5 | 19 | 0.1847 |
| 0.6 | 21 | 0.1920 |
| 0.7 | 21 | 0.1873 |
| 0.8 | 23 | 0.1970 |
| 0.9 | 26 | 0.1910 |
| 1.0 | 27 | 0.1895 |
| 1.1 | 28 | 0.1849 |
| 1.2 | 34 | 0.1317 |
| 1.3 | 32 | 0.1343 |
| 1.4 | 37 | 0.1300 |
| 1.5 | 37 | 0.1186 |
| 1.6 | 43 | 0.0940 |
| 1.7 | 43 | 0.1046 |
| 1.8 | 47 | 0.1086 |
| 1.9 | 44 | 0.0940 |
| 2.0 | 50 | 0.0953 |

## Cluster sizes and rule-based names

### k-means, chosen k

| cluster | cells | rule-based name | top 10 markers |
|---|---|---|---|
| 0 | 91,134 | CD8 T cells | RPL3, RPS27, RPL5, RPS3, RPL15, RPS6, LTB, RPSA, MALAT1, RPL10 |
| 1 | 3,521 | Monocytes | CST3, S100A9, S100A8, TYROBP, FTH1, FTL, LYZ, S100A4, S100A6, AIF1 |

Marker table: `results/tables/06_markers_kmeans_primary.csv`. Rule scores: `results/tables/06_rule_names_kmeans_primary.csv`.

### k-means, k = 10 (label-informed)

| cluster | cells | rule-based name | top 10 markers |
|---|---|---|---|
| 0 | 24,476 | CD4 T cells | RPL3, EEF1A1, RPL10, RPS6, RPS2, RPL4, RPL13, RPS4X, RPS3, LDHB |
| 1 | 23,104 | CD8 T cells | RPL32, RPS12, RPS14, RPL31, RPLP2, RPS27, CD8B, RPL21, RPS25, RPS18 |
| 2 | 16,601 | CD4 T cells | IL32, S100A4, B2M, SH3BGRL3, HLA-A, HLA-B, PFN1, S100A10, HLA-C, TMSB4X |
| 3 | 10,472 | B cells | CD74, HLA-DRA, CD79A, HLA-DRB1, HLA-DPB1, HLA-DPA1, CD79B, CD37, HLA-DRB5, HLA-DQA1 |
| 4 | 8,579 | NK cells | GNLY, NKG7, GZMA, CST7, GZMB, TYROBP, FCER1G, MALAT1, CTSW, CLIC3 |
| 5 | 3,140 | Monocytes | FTL, S100A9, S100A8, TYROBP, CST3, FTH1, S100A4, LYZ, S100A6, AIF1 |
| 6 | 2,771 | CD34+ progenitors | SPINK2, ZFAS1, AIF1, LST1, RPS24, KIAA0125, GSTP1, SOX4, GPX1, CD74 |
| 7 | 2,731 | CD34+ progenitors | RPS24, ZFAS1, PRSS57, SPINK2, RP11-620J15.3, HINT1, SERPINB1, AIF1, RPS23, RPL37A |
| 8 | 2,337 | CD34+ progenitors | RP11-620J15.3, PRSS57, RPS24, RPL37A, ZFAS1, SERPINB1, CYTL1, SOX4, FCER1A, HSP90AB1 |
| 9 | 444 | Monocytes | CST3, HLA-DPB1, HLA-DRA, HLA-DPA1, LYZ, HLA-DRB1, GPX1, HLA-DQA1, HLA-DMA, CD74 |

Marker table: `results/tables/06_markers_kmeans_oracle_k10.csv`. Rule scores: `results/tables/06_rule_names_kmeans_oracle_k10.csv`.

### Leiden, chosen resolution

| cluster | cells | rule-based name | top 10 markers |
|---|---|---|---|
| 0 | 17,045 | CD8 T cells | CD8B, RPL32, RPS14, RPS12, RPL31, RPL21, RPLP2, RPS25, RPS27, RPS19 |
| 1 | 11,272 | CD4 T cells | EEF1A1, B2M, IL32, JUNB, RPL3, RPL10, LTB, RPS2, ACTB, VIM |
| 2 | 10,930 | CD4 T cells | RPL3, RPS3A, EEF1A1, RPS6, RPL10, RPL7, RPL4, RPS2, RPL13, RPS4X |
| 3 | 10,431 | B cells | CD74, HLA-DRA, CD79A, HLA-DRB1, HLA-DPB1, HLA-DPA1, CD79B, CD37, HLA-DRB5, HLA-DQA1 |
| 4 | 9,445 | CD4 T cells | S100A4, IL32, CD52, RPL41, RPS27, RPL36, RPS18, LTB, S100A6, S100A10 |
| 5 | 8,361 | NK cells | GNLY, NKG7, TYROBP, GZMA, CST7, FCER1G, GZMB, MALAT1, CLIC3, KLRB1 |
| 6 | 5,754 | CD4 T cells | RPL3, RPL10, RPS2, RPL13, EEF1A1, RPS6, RPL4, RPL5, EEF1D, RPS3 |
| 7 | 5,435 | CD8 T cells | CCL5, NKG7, GZMK, IL32, GZMA, LYAR, HLA-B, CST7, MALAT1, B2M |
| 8 | 4,677 | CD34+ progenitors | RPS24, PRSS57, RP11-620J15.3, ZFAS1, SERPINB1, RPL37A, AIF1, HSP90AB1, RPS23, HINT1 |
| 9 | 3,772 | CD4 T cells | IL32, HLA-A, B2M, ACTB, PFN1, SH3BGRL3, S100A4, HLA-B, ARHGDIB, TMSB4X |
| 10 | 3,425 | Monocytes | S100A9, FTL, S100A8, CST3, TYROBP, FTH1, S100A4, LYZ, S100A6, AIF1 |
| 11 | 2,781 | CD34+ progenitors | SPINK2, ZFAS1, AIF1, LST1, RPS24, GSTP1, KIAA0125, SOX4, GPX1, CD74 |
| 12 | 364 | Monocytes | CST3, GPX1, HLA-DPA1, HLA-DRA, HLA-DPB1, HLA-DRB1, CD74, HLA-DQA1, VIM, LYZ |
| 13 | 330 | CD34+ progenitors | LMO4, RPL37A, PRSS57, RP11-620J15.3, TPSAB1, RPS24, S100A6, RPS20, CD63, RPLP1 |
| 14 | 129 | CD8 T cells | SIRPB1, RPS3, NLRP3, JUNB, RPSA, LDHB, CD3E, RPL3, CD3D, RPL13 |
| 15 | 115 | CD8 T cells | CD101, RPS6, RPL3, RPL13, RPS5, RPS2, RPL32, RPS3, RPL21, EEF1A1 |
| 16 | 107 | CD8 T cells | USP35, CD3D, IL32, CD3E, JUNB, ACTB, UBC, CORO1A, EEF1A1, AES |
| 17 | 77 | CD8 T cells | LINC00937, IL32, BATF3, CD3D, CD3E, JUNB, B2M, TMSB4X, S100A6, SEPT9 |
| 18 | 75 | CD4 T cells | HOMER3, FXYD5, S100A4, IL32, RPL10, EEF1A1, VIM, LTB, B2M, HSPA8 |
| 19 | 39 | CD4 T cells | HES4, RPL4, RPS2, EEF1A1, RPL3, RPS6, JUNB, RPL5, TMEM66, RPL10 |
| 20 | 38 | CD4 T cells | NOL3, LTB, RPS3, IL32, TRADD, IFITM1, PSME1, RPL10, EEF1A1, EEF1D |
| 21 | 34 | unresolved | RPS27, RPS14, RPS5, RPS18, RPS15A, RPS12, RPL31, RPS19, RPS6, RPL18A |
| 22 | 19 | CD8 T cells | ADM, IL7R, CD3D, IER2, RPL29, LCK, RPSAP58, FXYD5, RPL13, RPS3 |

Marker table: `results/tables/06_markers_leiden_primary.csv`. Rule scores: `results/tables/06_rule_names_leiden_primary.csv`.

### Camera-settings baseline

| partition | cluster sizes |
|---|---|
| technical_k_chosen | 76,588, 18,067 |
| technical_k10 | 19,145, 17,801, 15,485, 10,569, 8,989, 8,616, 7,537, 5,541, 884, 88 |

## Stability, agreement and silhouette

| measure | min | median | max |
|---|---|---|---|
| kmeans pairwise ari | 0.4026 | 0.4026 | 1.0000 |
| kmeans ari to primary | 0.4026 | 0.4026 | 1.0000 |
| leiden pairwise ari | 0.7973 | 0.8845 | 0.9503 |
| leiden ari to primary | 0.8342 | 0.8980 | 0.9289 |

Distinct k-means solutions at the chosen k, among the primary run (seed 7) and the 20 reruns, tightest first:

| cluster sizes | inertia | silhouette (subsample) | seeds |
|---|---|---|---|
| 83,404, 11,251 | 15693770.8 | 0.5366 | 101, 104, 107, 108, 110, 111, 114, 115, 116, 117, 119, 120 |
| 91,134, 3,521 | 15850196.1 | 0.6336 | 7, 102, 103, 105, 106, 109, 112, 113, 118 |

| label-free agreement | ARI |
|---|---|
| oracle vs leiden | 0.5407 |
| primary vs leiden | 0.0174 |
| primary vs oracle | 0.0334 |

| partition | mean silhouette | share of cells with negative silhouette |
|---|---|---|
| kmeans_primary | 0.6336 | 0.0098 |
| kmeans_oracle_k10 | 0.2129 | 0.1178 |
| leiden_primary | 0.1970 | 0.1399 |

Refits match script 04: primary True, Leiden True. Oracle identical to primary: False.

## Runtimes (seconds)

Preprocessing 146; k and resolution selection 1139; clustering 295; naming 460.
