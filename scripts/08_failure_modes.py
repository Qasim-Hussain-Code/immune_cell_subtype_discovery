"""Score the Day 63 predictions and the Day 66 failure-mode checks.

The predictions concern ten labels, so they are scored on the partition that
was told to make ten groups. The same statistics for the primary k-means and
Leiden partitions are reported as description only.
"""
import itertools
import re
import time

import numpy as np
import pandas as pd
import scanpy as sc

from common import (
    ROOT,
    check_free_space,
    check_posts,
    contingency,
    load_config,
    open_answer_key,
    overlap,
    read_metrics,
    sub_ari,
    write_metrics,
)

TABLES = ROOT / "results" / "tables"
INTERIM = ROOT / "data" / "interim"
ARMS = ["kmeans_primary", "kmeans_oracle_k10", "leiden_primary", "technical_k_chosen", "technical_k10"]
DESCRIPTIVE = ["kmeans_primary", "leiden_primary"]


def threshold(rule):
    match = re.match(r"\s*(>=|<=)\s*([0-9.]+)\s*$", rule)
    return match.group(1), float(match.group(2))


def compare(value, op, limit):
    return value >= limit if op == ">=" else value <= limit


def prediction_stats(cfg, arm, ten, clusters, rp):
    pcfg = cfg["predictions"]
    table = contingency(ten, clusters, list(cfg["data"]["populations"]))
    rec = {r["label"]: r for r in rp[arm]["ten_labels"]}
    out = {}
    for pid in ["P1", "P2", "P4"]:
        op, limit = threshold(pcfg[pid]["rule"])
        value = rec[pcfg[pid]["label"]]["recovery"]
        out[pid] = {"label": pcfg[pid]["label"], "stat": "recovery", "value": value, "operator": op, "threshold": limit, "passed": compare(value, op, limit), "best_cluster": rec[pcfg[pid]["label"]]["best_cluster"]}
    op, limit = threshold(pcfg["P3"]["rule"])
    a, b = pcfg["P3"]["labels"]
    value = overlap(table, a, b)
    out["P3"] = {"labels": [a, b], "stat": "overlap", "value": value, "operator": op, "threshold": limit, "passed": compare(value, op, limit)}
    rec_limit, pur_limit = (float(x) for x in re.findall(r">=\s*([0-9.]+)", pcfg["P5"]["rule"]))
    per_label = []
    for label in pcfg["P5"]["labels"]:
        r = rec[label]
        clean = r["recovery"] >= rec_limit and r["purity"] >= pur_limit
        per_label.append({"label": label, "best_cluster": r["best_cluster"], "recovery": r["recovery"], "purity": r["purity"], "clean": clean})
    out["P5"] = {"labels": pcfg["P5"]["labels"], "recovery_threshold": rec_limit, "purity_threshold": pur_limit, "per_label": per_label, "n_clean": int(sum(p["clean"] for p in per_label)), "passed": not any(p["clean"] for p in per_label)}
    nested = cfg["diagnostics"]["nested_pair"]
    distinct = cfg["diagnostics"]["distinct_pair"]
    s_nested = sub_ari(ten, clusters, *nested)
    s_distinct = sub_ari(ten, clusters, *distinct)
    out["P6"] = {"nested_pair": nested, "distinct_pair": distinct, "sub_ari_nested": s_nested, "sub_ari_distinct": s_distinct, "operator": "<=", "passed": s_nested <= s_distinct}
    return out


def main():
    check_posts()
    check_free_space()
    sc.settings.verbosity = 0
    cfg = load_config()
    started = time.time()
    ev = read_metrics("07_evaluate")
    selection = read_metrics("04_select_k")
    clustering = read_metrics("05_cluster")
    blind = read_metrics("06_blind_naming")
    rp = ev["recovery_purity"]

    key = open_answer_key()
    assignments = pd.read_csv(ROOT / "results" / "assignments" / "assignments_main.csv.gz")
    key = key.set_index("cell_id").loc[assignments["cell_id"]]
    populations = list(cfg["data"]["populations"])
    ten = key["population_key"].to_numpy()
    parts = {arm: assignments[arm].to_numpy() for arm in ARMS}

    # Predictions.
    scored_on = cfg["predictions"]["scored_on"]
    scored = prediction_stats(cfg, scored_on, ten, parts[scored_on], rp)
    descriptive = {arm: prediction_stats(cfg, arm, ten, parts[arm], rp) for arm in DESCRIPTIVE}
    write_metrics("08_predictions", {"scored_on": scored_on, "predictions": scored, "descriptive": descriptive, "n_passed": int(sum(v["passed"] for v in scored.values())), "n_predictions": len(scored)})

    # Raw counts for the label genes, from the anonymised pooled matrix.
    pooled = sc.read_h5ad(ROOT / cfg["quarantine"]["pooled_file"])
    if not np.array_equal(pooled.obs_names.to_numpy(), assignments["cell_id"].to_numpy()):
        raise SystemExit("Stop: pooled order differs from the assignments")
    symbols = pooled.var["gene_symbols"].to_numpy()
    label_genes = cfg["diagnostics"]["label_genes"]
    columns = {}
    for g in label_genes:
        hits = np.flatnonzero(symbols == g)
        if len(hits) != 1:
            raise SystemExit(f"Stop: label gene {g} found {len(hits)} times")
        columns[g] = int(hits[0])
    counts = pooled.X[:, [columns[g] for g in label_genes]].toarray()
    del pooled
    detected = pd.DataFrame(counts > 0, columns=label_genes)
    detected["population"] = ten
    detection = detected.groupby("population", sort=False)[label_genes].mean().reindex(populations)
    detection.to_csv(TABLES / "08_label_gene_detection.csv", lineterminator="\n")

    # FM1: the run is the label.
    cd34_key = cfg["predictions"]["P4"]["label"]
    oracle = parts[scored_on]
    best = next(r for r in rp[scored_on]["ten_labels"] if r["label"] == cd34_key)["best_cluster"]
    cd34_detected = counts[:, label_genes.index("CD34")] > 0
    inside = oracle == best
    is_cd34 = ten == cd34_key
    fm1 = {
        "beats_verdicts": {name: {k: v["verdict"] for k, v in ev["paired_differences"][name].items()} for name in ["primary_minus_technical_k_chosen", "oracle_minus_technical_k10"]},
        "sub_ari": {arm: {"nested": s["P6"]["sub_ari_nested"], "distinct": s["P6"]["sub_ari_distinct"]} for arm, s in [(scored_on, scored)] + list(descriptive.items())},
        "run_driven": not scored["P6"]["passed"],
        "technical_medians": ev["populations"],
        "cd34": {
            "recovery_technical_k_chosen": next(r for r in rp["technical_k_chosen"]["ten_labels"] if r["label"] == cd34_key)["recovery"],
            "recovery_technical_k10": next(r for r in rp["technical_k10"]["ten_labels"] if r["label"] == cd34_key)["recovery"],
            "oracle_best_cluster": best,
            "cd34_detected_share_inside_best_cluster_all_cells": float(cd34_detected[inside].mean()),
            "cd34_detected_share_outside_best_cluster_all_cells": float(cd34_detected[~inside].mean()),
            "cd34_detected_share_inside_best_cluster_cd34_labelled": float(cd34_detected[inside & is_cd34].mean()) if (inside & is_cd34).any() else None,
            "cd34_detected_share_outside_best_cluster_cd34_labelled": float(cd34_detected[~inside & is_cd34].mean()) if (~inside & is_cd34).any() else None,
            "n_cd34_labelled_inside": int((inside & is_cd34).sum()),
            "n_cd34_labelled_outside": int((~inside & is_cd34).sum()),
        },
    }

    # FM2: the key's fault.
    pair_tables = {}
    for arm in ARMS:
        table = contingency(ten, parts[arm], populations)
        pairs = [{"label_a": a, "label_b": b, "overlap": overlap(table, a, b)} for a, b in itertools.combinations(populations, 2)]
        frame = pd.DataFrame(pairs).sort_values("overlap", ascending=False).reset_index(drop=True)
        frame.to_csv(TABLES / f"08_label_pair_overlap_{arm}.csv", index=False, lineterminator="\n")
        pair_tables[arm] = frame.head(10).to_dict(orient="records")
    fm2 = {
        "ari_ten_vs_six": {arm: {"ten_labels": ev["ari"][arm]["ten_labels"]["observed"], "six_lineages": ev["ari"][arm]["six_lineages"]["observed"]} for arm in ARMS},
        "top10_overlapping_pairs": pair_tables,
        "label_gene_detection_share": {p: {g: float(detection.loc[p, g]) for g in label_genes} for p in populations},
    }

    # FM3: counting versus grouping.
    curve = {r["k"]: r["silhouette"] for r in selection["kmeans_curve"]}
    fm3 = {
        "chosen_k": selection["chosen_k"],
        "oracle_k": cfg["kmeans"]["oracle_k"],
        "oracle_minus_primary": ev["paired_differences"]["oracle_minus_primary"],
        "silhouette_chosen_k": curve[selection["chosen_k"]],
        "silhouette_k10": curve[cfg["kmeans"]["oracle_k"]],
    }

    # FM4: every cell gets a home.
    fm4 = {
        "silhouette": clustering["silhouette"],
        "stability_label_free": {m: {s: clustering["stability"][m][s] for s in ["pairwise_ari", "ari_to_primary"]} for m in ["kmeans", "leiden"]},
        "stability_vs_labels": ev["stability_vs_labels"],
        "agreement_ari": clustering["agreement_ari"],
    }

    # FM5: unmatched is not a discovery.
    unresolved = []
    need_markers = {}
    for arm in ARMS:
        for r in ev["cluster_truth"][arm]:
            if r["truth"] == "unresolved":
                unresolved.append(dict(r, arm=arm))
                if arm not in blind["top10_markers"]:
                    need_markers.setdefault(arm, []).append(r["cluster"])
    extra_markers = {}
    if need_markers:
        lognorm = sc.read_h5ad(INTERIM / "lognorm.h5ad")
        for arm in need_markers:
            labels = parts[arm]
            groups = [str(c) for c in np.unique(labels)]
            lognorm.obs["group"] = pd.Categorical([str(v) for v in labels], categories=groups)
            sc.tl.rank_genes_groups(lognorm, groupby="group", method="wilcoxon", use_raw=False, pts=True, n_genes=cfg["naming"]["wilcoxon_top_n"], key_added="rgg")
            table = sc.get.rank_genes_groups_df(lognorm, group=None, key="rgg")
            table["gene_symbol"] = lognorm.var["gene_symbols"].loc[table["names"]].to_numpy()
            table = table.rename(columns={"group": "cluster", "names": "gene_id"})
            table["cluster"] = table["cluster"].astype(int)
            table.to_csv(TABLES / f"08_markers_{arm}.csv", index=False, lineterminator="\n")
            extra_markers[arm] = {int(c): g["gene_symbol"].head(10).tolist() for c, g in table.groupby("cluster", sort=True)}
            del lognorm.uns["rgg"]
    for r in unresolved:
        source = blind["top10_markers"].get(r["arm"])
        r["top10_markers"] = source[str(r["cluster"])] if source is not None else extra_markers[r["arm"]][r["cluster"]]
    pd.DataFrame(unresolved).to_csv(TABLES / "08_unresolved_clusters.csv", index=False, lineterminator="\n")
    rule_unresolved = []
    for part, rows in blind["rule_names"].items():
        truth = {t["cluster"]: t for t in ev["cluster_truth"][part]}
        for r in rows:
            if r["name"] == "unresolved":
                t = truth[r["cluster"]]
                rule_unresolved.append({"partition": part, "cluster": r["cluster"], "n_cells": r["n_cells"], "truth": t["truth"], "plurality_lineage": t["plurality_lineage"], "plurality_share": t["plurality_share"], "top10_markers": blind["top10_markers"][part][str(r["cluster"])]})
    fm5 = {"unresolved_truth_clusters": unresolved, "n_unresolved_truth_clusters": len(unresolved), "rule_unresolved_names": rule_unresolved, "n_rule_unresolved_names": len(rule_unresolved), "new_cell_types_claimed": False}

    write_metrics("08_failure_modes", {"FM1": fm1, "FM2": fm2, "FM3": fm3, "FM4": fm4, "FM5": fm5, "runtime_seconds": time.time() - started})
    print(f"predictions passed: {sum(v['passed'] for v in scored.values())} of {len(scored)}")


if __name__ == "__main__":
    main()
