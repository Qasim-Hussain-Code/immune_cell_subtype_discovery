"""Pool all cells, shuffle them, replace their names and set the labels aside.

This is the only script before the freeze that sees which population each
cell came from. It writes an anonymised pooled matrix, in which each cell
carries a random-order ID and nothing else, and a separate label file that no
clustering script may read. Both files are hashed, and nothing is printed per
population.
"""
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from common import ROOT, check_free_space, check_posts, load_config, rng, sha256_file, write_metrics, write_text


def main():
    check_posts()
    check_free_space()
    sc.settings.verbosity = 0
    cfg = load_config()
    qcfg = cfg["quarantine"]

    blocks, keys, barcodes = [], [], []
    var = None
    for key in cfg["data"]["populations"]:
        matrices = sorted((ROOT / "data" / "raw" / key).rglob("matrix.mtx"))
        if len(matrices) != 1:
            raise SystemExit(f"Stop: expected one matrix.mtx for {key}")
        part = sc.read_10x_mtx(matrices[0].parent, var_names="gene_ids", cache=False)
        if var is None:
            var = part.var[["gene_symbols"]].copy()
        elif not (part.var_names.equals(var.index) and (part.var["gene_symbols"].to_numpy() == var["gene_symbols"].to_numpy()).all()):
            raise SystemExit(f"Stop: gene table for {key} differs from the first matrix")
        blocks.append(sp.csr_matrix(part.X))
        keys.extend([key] * part.n_obs)
        barcodes.extend(part.obs_names.tolist())

    counts = sp.vstack(blocks, format="csr")
    order = rng(cfg["rng_streams"]["shuffle"]).permutation(counts.shape[0])
    counts = counts[order]
    keys = np.asarray(keys)[order]
    barcodes = np.asarray(barcodes)[order]
    ids = [f"c{i:06d}" for i in range(1, counts.shape[0] + 1)]

    obs = pd.DataFrame(index=pd.Index(ids))
    var_out = pd.DataFrame({"gene_symbols": var["gene_symbols"].to_numpy()}, index=pd.Index(var.index.to_numpy()))
    pooled = ad.AnnData(X=counts, obs=obs, var=var_out)
    pooled.uns.clear()
    for column in list(pooled.obs.columns):
        del pooled.obs[column]
    pooled_path = ROOT / qcfg["pooled_file"]
    pooled_path.parent.mkdir(parents=True, exist_ok=True)
    pooled.write_h5ad(pooled_path, compression="gzip")

    labels_path = ROOT / qcfg["labels_file"]
    pd.DataFrame({"cell_id": ids, "population_key": keys, "original_barcode": barcodes}).to_csv(labels_path, index=False, lineterminator="\n")

    write_text("provenance/quarantine_sha256.txt", f"{sha256_file(labels_path)}  {qcfg['labels_file']}\n")
    write_text("provenance/pooled_sha256.txt", f"{sha256_file(pooled_path)}  {qcfg['pooled_file']}\n")

    adjacent = float(np.mean(keys[1:] == keys[:-1]))
    passed = adjacent < qcfg["max_adjacent_same_label"]
    write_metrics(
        "02_build",
        {
            "n_cells": int(counts.shape[0]),
            "n_genes": int(counts.shape[1]),
            "shuffle_check_passed": passed,
            "adjacent_same_label_fraction": adjacent,
        },
    )
    if not passed:
        raise SystemExit("Stop: shuffle check failed")
    print(f"pooled {counts.shape[0]} cells x {counts.shape[1]} genes; shuffle check passed")


if __name__ == "__main__":
    main()
