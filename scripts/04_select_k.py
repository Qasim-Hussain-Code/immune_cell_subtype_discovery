"""Choose the number of k-means groups and the Leiden resolution.

Both choices follow the rule published on Day 65: the highest mean silhouette,
computed on the same fixed 10,000-cell subsample in the 50-component PCA
space, with ties going to the smaller k or the lower resolution. Inertia is
recorded for description only and plays no part in the choice.
"""
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from common import ROOT, apply_style, check_free_space, check_posts, load_config, rng, write_metrics

INTERIM = ROOT / "data" / "interim"


def kmeans_model(km_cfg, k, seed):
    return KMeans(
        n_clusters=k,
        init=km_cfg["init"],
        n_init=km_cfg["n_init"],
        max_iter=km_cfg["max_iter"],
        tol=km_cfg["tol"],
        algorithm=km_cfg["algorithm"],
        random_state=seed,
    )


def pick(values, scores):
    """Highest score after rounding to 4 decimal places; ties go to the smallest k or lowest resolution."""
    rounded = {v: round(s, 4) for v, s in zip(values, scores) if s is not None}
    best = max(rounded.values())
    return min(v for v, s in rounded.items() if s == best)


def main():
    check_posts()
    check_free_space()
    sc.settings.verbosity = 0
    cfg = load_config()
    ks_cfg, km_cfg, ld_cfg = cfg["k_selection"], cfg["kmeans"], cfg["leiden"]
    started = time.time()

    stored = np.load(INTERIM / "pca.npz", allow_pickle=True)
    x_pca, ids = stored["X_pca"], stored["cell_id"]
    n = x_pca.shape[0]

    sub = np.sort(rng(cfg["rng_streams"]["silhouette_subsample"]).choice(n, size=ks_cfg["subsample_size"], replace=False))
    pd.DataFrame({"cell_id": ids[sub]}).to_csv(ROOT / "results" / "tables" / "04_silhouette_subsample.csv", index=False, lineterminator="\n")
    x_sub = x_pca[sub]

    kmeans_curve, kmeans_labels = [], {}
    for k in ks_cfg["k_values"]:
        t0 = time.time()
        model = kmeans_model(km_cfg, k, km_cfg["random_state"]).fit(x_pca)
        labels = model.labels_
        score = float(silhouette_score(x_sub, labels[sub], metric=ks_cfg["metric"]))
        kmeans_labels[f"k{k}"] = labels
        kmeans_curve.append({"k": k, "silhouette": score, "inertia": float(model.inertia_), "n_iter": int(model.n_iter_), "seconds": time.time() - t0})
        print(f"k={k}: silhouette {score:.4f} ({time.time() - t0:.0f} s)", flush=True)
    np.savez_compressed(INTERIM / "04_kmeans_labels.npz", **kmeans_labels)

    graph = ad.AnnData(obs=pd.DataFrame(index=pd.Index(ids)))
    graph.obsm["X_pca"] = x_pca
    t0 = time.time()
    sc.pp.neighbors(graph, n_neighbors=ld_cfg["n_neighbors"], n_pcs=ld_cfg["n_pcs"], use_rep="X_pca", random_state=cfg["project"]["master_seed"])
    neighbours_seconds = time.time() - t0

    leiden_curve, leiden_labels = [], {}
    for res in ld_cfg["resolutions"]:
        t0 = time.time()
        sc.tl.leiden(graph, resolution=res, flavor=ld_cfg["flavor"], n_iterations=ld_cfg["n_iterations"], directed=ld_cfg["directed"], random_state=ld_cfg["random_state"], key_added="leiden_scan")
        labels = graph.obs["leiden_scan"].astype(int).to_numpy()
        n_clusters = int(len(np.unique(labels)))
        score = None
        if n_clusters > 1 and len(np.unique(labels[sub])) > 1:
            score = float(silhouette_score(x_sub, labels[sub], metric=ks_cfg["metric"]))
        leiden_labels[f"r{res:.1f}"] = labels
        leiden_curve.append({"resolution": res, "n_clusters": n_clusters, "silhouette": score, "seconds": time.time() - t0})
        print(f"resolution {res:.1f}: {n_clusters} clusters, silhouette {score}", flush=True)
    del graph.obs["leiden_scan"]
    np.savez_compressed(INTERIM / "04_leiden_labels.npz", **leiden_labels)
    graph.write_h5ad(INTERIM / "04_neighbours.h5ad", compression="gzip")

    chosen_k = pick([r["k"] for r in kmeans_curve], [r["silhouette"] for r in kmeans_curve])
    chosen_res = pick([r["resolution"] for r in leiden_curve], [r["silhouette"] for r in leiden_curve])
    chosen_k_row = next(r for r in kmeans_curve if r["k"] == chosen_k)
    chosen_res_row = next(r for r in leiden_curve if r["resolution"] == chosen_res)

    plt = apply_style()
    ks = [r["k"] for r in kmeans_curve]
    fig, ax = plt.subplots(figsize=(6, 3.2))
    ax.plot(ks, [r["silhouette"] for r in kmeans_curve], marker="o", color="#17456B")
    ax.axvline(chosen_k, color="#F07D22", linestyle="--", linewidth=1)
    ax.set_xticks(ks)
    ax.set_xlabel("k")
    ax.set_ylabel("mean silhouette (10,000-cell subsample)")
    ax.set_title(f"k-means: silhouette by k, chosen k = {chosen_k}")
    fig.savefig(ROOT / "figures" / "04_silhouette_kmeans.png")
    plt.close(fig)

    resolutions = [r["resolution"] for r in leiden_curve]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6, 5), sharex=True)
    a1.plot(resolutions, [np.nan if r["silhouette"] is None else r["silhouette"] for r in leiden_curve], marker="o", color="#17456B")
    a1.axvline(chosen_res, color="#F07D22", linestyle="--", linewidth=1)
    a1.set_ylabel("mean silhouette")
    a1.set_title(f"Leiden: silhouette by resolution, chosen resolution = {chosen_res:.1f}")
    a2.plot(resolutions, [r["n_clusters"] for r in leiden_curve], marker="o", color="#21A6B8")
    a2.axvline(chosen_res, color="#F07D22", linestyle="--", linewidth=1)
    a2.set_ylabel("number of clusters")
    a2.yaxis.get_major_locator().set_params(integer=True)
    a2.set_xlabel("resolution")
    fig.savefig(ROOT / "figures" / "04_silhouette_leiden.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 3.2))
    ax.plot(ks, [r["inertia"] for r in kmeans_curve], marker="o", color="#7A8899")
    ax.set_xticks(ks)
    ax.set_xlabel("k")
    ax.set_ylabel("inertia")
    ax.set_title("k-means inertia by k: not used for the choice")
    fig.savefig(ROOT / "figures" / "04_kmeans_inertia.png")
    plt.close(fig)

    write_metrics(
        "04_select_k",
        {
            "n_cells": int(n),
            "subsample_size": int(len(sub)),
            "kmeans_curve": kmeans_curve,
            "leiden_curve": leiden_curve,
            "chosen_k": int(chosen_k),
            "chosen_k_silhouette": chosen_k_row["silhouette"],
            "chosen_resolution": chosen_res,
            "chosen_resolution_n_clusters": chosen_res_row["n_clusters"],
            "chosen_resolution_silhouette": chosen_res_row["silhouette"],
            "neighbours_seconds": neighbours_seconds,
            "runtime_seconds": time.time() - started,
        },
    )
    print(f"chosen k = {chosen_k}; chosen resolution = {chosen_res:.1f} ({chosen_res_row['n_clusters']} clusters)")


if __name__ == "__main__":
    main()
