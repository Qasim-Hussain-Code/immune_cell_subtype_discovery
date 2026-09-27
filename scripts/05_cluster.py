"""Script 05: final partitions, stability reruns and the technical baseline.

Writes the assignment files that the freeze commit locks. UMAP coordinates are
computed here for figures only; nothing downstream uses them for a choice,
metric or name.
"""
import itertools
import time

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_samples

from common import ROOT, check_free_space, check_posts, gzip_csv, load_config, read_metrics, relabel_by_size, write_metrics

INTERIM = ROOT / "data" / "interim"


def kmeans_labels(x, km_cfg, k, seed):
    model = KMeans(
        n_clusters=k,
        init=km_cfg["init"],
        n_init=km_cfg["n_init"],
        max_iter=km_cfg["max_iter"],
        tol=km_cfg["tol"],
        algorithm=km_cfg["algorithm"],
        random_state=seed,
    )
    return model.fit(x).labels_


def leiden_labels(graph, ld_cfg, resolution, seed):
    sc.tl.leiden(graph, resolution=resolution, flavor=ld_cfg["flavor"], n_iterations=ld_cfg["n_iterations"], directed=ld_cfg["directed"], random_state=seed, key_added="leiden_run")
    labels = graph.obs["leiden_run"].astype(int).to_numpy()
    del graph.obs["leiden_run"]
    return labels


def summary(values):
    values = np.asarray(values, dtype=float)
    return {"min": float(values.min()), "median": float(np.median(values)), "max": float(values.max()), "n": int(len(values))}


def stability(runs, primary):
    pairwise = [adjusted_rand_score(a, b) for a, b in itertools.combinations(runs, 2)]
    to_primary = [adjusted_rand_score(primary, r) for r in runs]
    return {"pairwise_ari": summary(pairwise), "ari_to_primary": summary(to_primary), "pairwise_values": pairwise, "to_primary_values": to_primary}


def silhouette_summary(x_sub, labels_sub):
    values = silhouette_samples(x_sub, labels_sub, metric="euclidean")
    per_cluster = [
        {"cluster": int(c), "n_subsample": int((labels_sub == c).sum()), "mean_silhouette": float(values[labels_sub == c].mean()), "negative_share": float((values[labels_sub == c] < 0).mean())}
        for c in np.unique(labels_sub)
    ]
    return {"mean_silhouette": float(values.mean()), "negative_share": float((values < 0).mean()), "per_cluster": per_cluster}


def sizes(labels):
    values, counts = np.unique(labels, return_counts=True)
    return [{"cluster": int(v), "n_cells": int(c)} for v, c in zip(values, counts)]


def main():
    check_posts()
    check_free_space()
    sc.settings.verbosity = 0
    cfg = load_config()
    km_cfg, ld_cfg = cfg["kmeans"], cfg["leiden"]
    started = time.time()
    timings = {}

    stored = np.load(INTERIM / "pca.npz", allow_pickle=True)
    x_pca, ids = stored["X_pca"], stored["cell_id"]
    selection = read_metrics("04_select_k")
    k, res = selection["chosen_k"], selection["chosen_resolution"]
    km04 = np.load(INTERIM / "04_kmeans_labels.npz")
    ld04 = np.load(INTERIM / "04_leiden_labels.npz")

    t0 = time.time()
    primary = kmeans_labels(x_pca, km_cfg, k, km_cfg["random_state"])
    primary_matches_04 = bool(np.array_equal(primary, km04[f"k{k}"]))
    oracle = kmeans_labels(x_pca, km_cfg, km_cfg["oracle_k"], km_cfg["random_state"])
    oracle_matches_04 = bool(np.array_equal(oracle, km04[f"k{km_cfg['oracle_k']}"]))
    oracle_identical_to_primary = bool(k == km_cfg["oracle_k"] and adjusted_rand_score(primary, oracle) == 1.0)
    timings["kmeans_primary_and_oracle"] = time.time() - t0

    graph = sc.read_h5ad(INTERIM / "04_neighbours.h5ad")
    t0 = time.time()
    leiden = leiden_labels(graph, ld_cfg, res, ld_cfg["random_state"])
    leiden_matches_04 = bool(np.array_equal(leiden, ld04[f"r{res:.1f}"]))
    timings["leiden_primary"] = time.time() - t0

    t0 = time.time()
    km_runs = []
    for seed in km_cfg["stability_seeds"]:
        km_runs.append(kmeans_labels(x_pca, km_cfg, k, seed))
        print(f"k-means rerun seed {seed} done", flush=True)
    timings["kmeans_stability"] = time.time() - t0
    t0 = time.time()
    ld_runs = [leiden_labels(graph, ld_cfg, res, seed) for seed in ld_cfg["stability_seeds"]]
    timings["leiden_stability"] = time.time() - t0

    technical = pd.read_csv(INTERIM / "technical.csv")
    if not np.array_equal(technical["cell_id"].to_numpy(), ids):
        raise SystemExit("Stop: technical.csv order differs from the PCA order")
    features = technical[cfg["technical_baseline"]["features"]].to_numpy(dtype=float)
    z = (features - features.mean(axis=0)) / features.std(axis=0)
    t0 = time.time()
    tech_chosen = kmeans_labels(z, km_cfg, k, km_cfg["random_state"])
    tech_10 = kmeans_labels(z, km_cfg, km_cfg["oracle_k"], km_cfg["random_state"])
    timings["technical_baseline"] = time.time() - t0

    main_parts = {
        "kmeans_primary": relabel_by_size(primary),
        "kmeans_oracle_k10": relabel_by_size(oracle),
        "leiden_primary": relabel_by_size(leiden),
        "technical_k_chosen": relabel_by_size(tech_chosen),
        "technical_k10": relabel_by_size(tech_10),
    }
    stab_parts = {}
    for seed, labels in zip(km_cfg["stability_seeds"], km_runs):
        stab_parts[f"kmeans_seed{seed}"] = relabel_by_size(labels)
    for seed, labels in zip(ld_cfg["stability_seeds"], ld_runs):
        stab_parts[f"leiden_seed{seed}"] = relabel_by_size(labels)

    subsample_ids = pd.read_csv(ROOT / "results" / "tables" / "04_silhouette_subsample.csv")["cell_id"].to_numpy()
    position = pd.Series(np.arange(len(ids)), index=ids)
    sub = position.loc[subsample_ids].to_numpy()
    silhouettes = {name: silhouette_summary(x_pca[sub], main_parts[name][sub]) for name in ["kmeans_primary", "kmeans_oracle_k10", "leiden_primary"]}

    agreement = {
        "primary_vs_leiden": adjusted_rand_score(main_parts["kmeans_primary"], main_parts["leiden_primary"]),
        "primary_vs_oracle": adjusted_rand_score(main_parts["kmeans_primary"], main_parts["kmeans_oracle_k10"]),
        "oracle_vs_leiden": adjusted_rand_score(main_parts["kmeans_oracle_k10"], main_parts["leiden_primary"]),
    }

    t0 = time.time()
    sc.tl.umap(graph, random_state=cfg["project"]["master_seed"])
    np.savez_compressed(INTERIM / "05_umap.npz", X_umap=np.asarray(graph.obsm["X_umap"]), cell_id=ids)
    timings["umap"] = time.time() - t0

    gzip_csv(pd.DataFrame({"cell_id": ids, **main_parts}), "results/assignments/assignments_main.csv.gz")
    gzip_csv(pd.DataFrame({"cell_id": ids, **stab_parts}), "results/assignments/assignments_stability.csv.gz")

    write_metrics(
        "05_cluster",
        {
            "chosen_k": k,
            "chosen_resolution": res,
            "oracle_k": km_cfg["oracle_k"],
            "primary_identical_to_04": primary_matches_04,
            "oracle_identical_to_04": oracle_matches_04,
            "leiden_identical_to_04": leiden_matches_04,
            "oracle_identical_to_primary": oracle_identical_to_primary,
            "n_clusters": {name: int(len(np.unique(v))) for name, v in main_parts.items()},
            "cluster_sizes": {name: sizes(v) for name, v in main_parts.items()},
            "stability": {"kmeans": stability([stab_parts[f"kmeans_seed{s}"] for s in km_cfg["stability_seeds"]], main_parts["kmeans_primary"]),
                          "leiden": stability([stab_parts[f"leiden_seed{s}"] for s in ld_cfg["stability_seeds"]], main_parts["leiden_primary"])},
            "stability_n_clusters": {name: int(len(np.unique(v))) for name, v in stab_parts.items()},
            "agreement_ari": agreement,
            "silhouette": silhouettes,
            "technical_zscore_ddof": 0,
            "timings_seconds": timings,
            "runtime_seconds": time.time() - started,
        },
    )
    print(f"partitions written; primary k = {k}, Leiden clusters = {len(np.unique(leiden))}")
    if not (primary_matches_04 and leiden_matches_04):
        print("note: a refit differs from script 04; see 05_cluster.json")


if __name__ == "__main__":
    main()
