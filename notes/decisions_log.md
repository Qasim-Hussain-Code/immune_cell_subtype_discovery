# Decisions log

Judgement calls that do not change the pre-registered plan.

## 2026-09-26

- Git identity for all commits is the owner's global identity, Qasim Hussain, chosen by the owner and set in this repository's local config. No co-author or tool attribution lines are added (rule 0.6).
- The master prompt's opening line says the rules were published on Days 62 to 65, while config.yaml lists days 62 to 66 and rule 0.11 checks days 61 to 66. config.yaml and rule 0.11 are followed.
- All ten download URLs built from the unverified pattern in config.yaml returned HTTP 200 to a header-only request. Script 01 still records each download.
- The 10x dataset pages could not be cross-checked: their static HTML does not contain purity, cell counts or reads per cell (the values load with JavaScript), and a second fetch was refused with HTTP 429. The section 2 facts are used as read from the pages by the owner.

## 2026-09-28

- No Python 3.11 was installed, so `.venv` is a conda prefix environment with Python `3.11.16` from conda-forge; all analysis packages come from pip and are pinned in `requirements.txt`.
- The drive had little free space, so every script after 00 stops if less than 1 GB is free, and the intermediate `.h5ad` files are written with gzip compression. Neither changes a result.
- The quarantine test also forbids the config key `labels_file` in scripts 03 to 06, in addition to the strings listed in the prompt. This is stricter than the plan and changes no result.
- Standard deviations for the technical z-scores and for the naming statistic d(g, c) use the population form (ddof 0).
- The neighbour graph is computed once in script 04 and reused in script 05 for the primary Leiden run, the Leiden stability reruns and the UMAP picture.
- Assignment files are gzipped with a fixed timestamp so identical content always has an identical hash.
- If two lineage scores tie for the highest value in rule-based naming, the first in config order wins and the tie is flagged in the rule table. If two lineages tie for plurality in a cluster, the first in config order is reported and the tie is flagged; such a cluster can only reach the 0.5 share when both hold exactly half.
- Script 06 ranks Wilcoxon markers only for the three named partitions. For FM5, script 08 ranks markers for a technical-baseline partition with the same settings only if that partition has a cluster whose truth is unresolved.
- The owner's finer calls are compared with each cluster's most common population as description only; they are not part of any pre-registered score.
- Figures with more than 20 clusters cycle through tab20 followed by tab20b.
- Rule 0.11 says the owner pastes the posts. At the owner's request, the published text of Days 61 to 66 was saved verbatim into `posts/` from text the owner supplied in the conversation. The owner confirmed that Days 63 to 65 went out exactly as drafted and that Day 66 went out exactly as pasted. LinkedIn's "hashtag#" copy artefacts, trailing spaces, the Day 66 numbering (1, 2, 4, 5) and the Day 65 closing hashtag are kept as published.
- The published Day 62 says the populations were "sorted"; Day 63 corrects this to bead enrichment with flow cytometry used only to measure purity. The repository follows Day 63.
