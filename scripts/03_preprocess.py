"""Script 03: technical covariates, marker matrix, recipe_zheng17 and PCA.

The only input is the anonymised pooled matrix.
"""
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from common import ROOT, apply_style, check_free_space, check_posts, load_config, write_metrics

INTERIM = ROOT / "data" / "interim"


def main():
    check_posts()
    check_free_space()
    sc.settings.verbosity = 0
    cfg = load_config()
    pp = cfg["preprocess"]
    started = time.time()

    pooled = sc.read_h5ad(ROOT / cfg["quarantine"]["pooled_file"])
    ids = pooled.obs_names.to_numpy()
    raw = pooled.X.tocsr()
    symbols = pooled.var["gene_symbols"]

    # Technical covariates, one row per cell.
    qc = ad.AnnData(X=raw.copy(), obs=pd.DataFrame(index=pooled.obs_names), var=pooled.var[["gene_symbols"]].copy())
    qc.var["mito"] = symbols.str.startswith(pp["mito_prefix"]).to_numpy()
    sc.pp.calculate_qc_metrics(qc, qc_vars=["mito"], percent_top=None, log1p=False, inplace=True)
    technical = pd.DataFrame(
        {
            "cell_id": ids,
            "total_counts": qc.obs["total_counts"].to_numpy(),
            "n_genes": qc.obs["n_genes_by_counts"].to_numpy(),
            "pct_mito": qc.obs["pct_counts_mito"].to_numpy(),
        }
    )
    technical["log10_total_counts"] = np.log10(technical["total_counts"])
    technical["log10_n_genes"] = np.log10(technical["n_genes"])
    technical.to_csv(INTERIM / "technical.csv", index=False, lineterminator="\n")
    n_mito = int(qc.var["mito"].sum())
    del qc

    # Log-normalised marker matrix over all genes.
    lognorm = ad.AnnData(X=raw.copy(), obs=pd.DataFrame(index=pooled.obs_names), var=pooled.var[["gene_symbols"]].copy())
    sc.pp.normalize_total(lognorm, target_sum=pp["marker_target_sum"])
    sc.pp.log1p(lognorm)
    lognorm.write_h5ad(INTERIM / "lognorm.h5ad", compression="gzip")
    del lognorm

    # Preprocessing recipe from the 2017 paper, as shipped by Scanpy.
    n_genes_kept = int((np.asarray(raw.sum(axis=0)).ravel() >= 1).sum())
    recipe = ad.AnnData(X=raw.astype(np.float32), obs=pd.DataFrame(index=pooled.obs_names), var=pooled.var[["gene_symbols"]].copy())
    del pooled, raw
    sc.pp.recipe_zheng17(recipe, n_top_genes=pp["n_top_genes"], log=pp["log"], plot=False)
    hvg = pd.DataFrame({"gene_id": recipe.var_names.to_numpy(), "gene_symbol": recipe.var["gene_symbols"].to_numpy()})
    hvg.to_csv(ROOT / "results" / "tables" / "03_hvg_genes.csv", index=False, lineterminator="\n")

    sc.pp.pca(recipe, n_comps=pp["n_pcs"], svd_solver=pp["pca_svd_solver"], random_state=cfg["project"]["master_seed"])
    x_pca = np.asarray(recipe.obsm["X_pca"])
    ratio = np.asarray(recipe.uns["pca"]["variance_ratio"], dtype=float)
    np.savez_compressed(INTERIM / "pca.npz", X_pca=x_pca, cell_id=np.asarray(recipe.obs_names))
    cumulative = np.cumsum(ratio)

    plt = apply_style()
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    pcs = np.arange(1, len(ratio) + 1)
    ax.bar(pcs, ratio, color="#2E6DA4", width=0.8, label="per component")
    ax.set_xlabel("principal component")
    ax.set_ylabel("share of variance")
    twin = ax.twinx()
    twin.plot(pcs, cumulative, color="#F07D22", marker="o", markersize=2, label="cumulative")
    twin.set_ylabel("cumulative share")
    twin.set_ylim(0, 1)
    twin.spines["right"].set_visible(True)
    ax.set_title(f"Variance captured by {len(ratio)} components of {hvg.shape[0]} scaled genes")
    fig.savefig(ROOT / "figures" / "03_pca_variance.png")
    plt.close(fig)

    write_metrics(
        "03_preprocess",
        {
            "n_cells": int(len(ids)),
            "n_genes_input": int(len(symbols)),
            "n_genes_kept": n_genes_kept,
            "n_mito_genes": n_mito,
            "n_hvg": int(hvg.shape[0]),
            "n_pcs": int(len(ratio)),
            "variance_ratio": ratio,
            "cumulative_variance": cumulative,
            "cumulative_variance_all_pcs": float(cumulative[-1]),
            "median_total_counts_all_cells": float(technical["total_counts"].median()),
            "median_n_genes_all_cells": float(technical["n_genes"].median()),
            "median_pct_mito_all_cells": float(technical["pct_mito"].median()),
            "runtime_seconds": time.time() - started,
        },
    )
    print(f"PCA done: {len(ratio)} components hold {cumulative[-1]:.4f} of the variance")


if __name__ == "__main__":
    main()
