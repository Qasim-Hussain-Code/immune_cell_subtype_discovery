# Decisions log

Judgement calls made during the analysis that do not change the pre-registered plan. Departures from the plan are recorded separately in `deviations.md`.

## 2026-09-26

- The ten download links built from the pattern in `config.yaml` had not been verified in advance. All ten answered a header-only request with HTTP 200, and script 01 logs each download with its size and hash.
- The purity, cell counts and reads per cell recorded by script 00 were read by hand from the 10x Genomics dataset pages. I could not re-check them automatically, because the pages load these values with JavaScript and the server refused repeated requests (HTTP 429).

## 2026-09-28

- The analysis runs in Python `3.11.16` inside a conda prefix environment at `.venv`. All analysis packages were installed with pip and are pinned in `requirements.txt`; the full list is in `provenance/environment_freeze.txt`.
- Every script after 00 stops if less than 1 GB is free on the drive, and the intermediate `.h5ad` files are gzip-compressed. Neither choice affects a result.
- The quarantine test forbids the config key `labels_file` in scripts 03 to 06, in addition to the label file name, `data/raw`, `open_answer_key` and the population keys. This is stricter than the plan and changes no result.
- Standard deviations for the technical z-scores and for the naming statistic d(g, c) use the population form (ddof 0).
- The neighbour graph is computed once in script 04 and reused in script 05 for the primary Leiden run, the Leiden stability reruns and the UMAP picture, so all three rest on the same graph.
- Assignment files are gzipped with a fixed timestamp, so identical content always produces an identical hash.
- If two lineage scores tie for the highest value in rule-based naming, the first in config order wins and the tie is flagged in the rule table. If two lineages tie for plurality within a cluster, the first in config order is reported and the tie is flagged; such a cluster can reach the 0.5 share only when each holds exactly half. No ties occurred.
- Script 06 ranks Wilcoxon markers only for the three named partitions. For the FM5 check, script 08 ranks markers with the same settings for a technical-baseline partition, and only when that partition contains a cluster whose truth is unresolved.
- The finer calls in the second naming are compared with each cluster's most common population as description only. They are not part of any pre-registered score.
- Figures with more than 20 clusters cycle through the tab20 palette followed by tab20b.
- The posts in `posts/` are stored exactly as published on LinkedIn, including LinkedIn's "hashtag#" copy artefacts, trailing spaces, the numbering of the Day 66 list (1, 2, 4, 5) and the closing hashtag of Day 65.
- The published Day 62 describes the populations as "sorted". Day 63 corrects this: they were enriched with antibody-coated beads, and flow cytometry was used only to measure purity. The repository follows Day 63.
- `scanpy.pp.recipe_zheng17` was run with `n_top_genes=1000`, as pre-registered, and returned 999 genes. The selection is made inside the recipe by `filter_genes_dispersion(flavor="cell_ranger")`, whose source is saved in `provenance/recipe_zheng17_source.py`. I used the recipe's output as returned rather than altering the recipe to force 1,000 genes.
- The k-means stability reruns at the chosen k did not all reproduce the primary partition. Script 06 therefore also records, without reference to any label, each distinct solution among the primary run (seed 7) and the 20 reruns: its cluster sizes, its inertia recomputed in float64 on the 50 components, its silhouette on the same 10,000-cell subsample, and the seeds that reached it (`06_blind_naming.json`, key `kmeans_rerun_solutions`). This is description only; the primary partition remains the seed 7 run, as pre-registered.
- The answer key was opened only after the second set of cluster names had been committed.
