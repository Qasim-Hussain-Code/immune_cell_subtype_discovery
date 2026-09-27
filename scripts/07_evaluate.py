"""Script 07: open the answer key and score every frozen partition.

The first script allowed to read the labels. ARI with paired bootstrap
confidence intervals and a permutation null for every arm and both keys,
paired differences with verdicts, recovery and purity, truth per cluster,
naming scores and stability against the labels.
"""
import time

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from common import (
    NAMED_PARTITIONS,
    ROOT,
    TENX_NAMES,
    best_cluster,
    check_free_space,
    check_posts,
    ci_bounds,
    cluster_truth,
    contingency,
    load_config,
    open_answer_key,
    parse_human_naming,
    rng,
    write_metrics,
)

ARMS = ["kmeans_primary", "kmeans_oracle_k10", "leiden_primary", "technical_k_chosen", "technical_k10"]
TABLES = ROOT / "results" / "tables"


def verdict(kind, low, high):
    if kind == "beats":
        return "beats" if low > 0 else "does not beat"
    if kind == "leiden_vs_kmeans":
        return "Leiden better" if low > 0 else ("k-means better" if high < 0 else "no clear difference")
    return "oracle higher" if low > 0 else ("primary higher" if high < 0 else "no clear difference")


PAIRS = {
    "primary_minus_technical_k_chosen": ("kmeans_primary", "technical_k_chosen", "beats"),
    "oracle_minus_technical_k10": ("kmeans_oracle_k10", "technical_k10", "beats"),
    "leiden_minus_primary": ("leiden_primary", "kmeans_primary", "leiden_vs_kmeans"),
    "oracle_minus_primary": ("kmeans_oracle_k10", "kmeans_primary", "oracle_vs_primary"),
}


def naming_score(names, truth_rows):
    truth = {r["cluster"]: r for r in truth_rows}
    detail = []
    for cluster, name in names.items():
        t = truth[cluster]
        detail.append({"cluster": cluster, "n_cells": t["n_cells"], "name": name, "truth": t["truth"], "correct": name == t["truth"]})
    frame = pd.DataFrame(detail).sort_values("cluster")
    n_cells = frame["n_cells"].sum()
    score = {
        "n_clusters": int(len(frame)),
        "n_correct": int(frame["correct"].sum()),
        "cluster_accuracy": float(frame["correct"].mean()),
        "cell_weighted_accuracy": float(frame.loc[frame["correct"], "n_cells"].sum() / n_cells),
        "n_named_unresolved": int((frame["name"] == "unresolved").sum()),
    }
    return score, frame


def main():
    check_posts()
    check_free_space()
    cfg = load_config()
    ecfg = cfg["evaluation"]
    started = time.time()

    key = open_answer_key()
    assignments = pd.read_csv(ROOT / "results" / "assignments" / "assignments_main.csv.gz")
    stability = pd.read_csv(ROOT / "results" / "assignments" / "assignments_stability.csv.gz")
    key = key.set_index("cell_id").loc[assignments["cell_id"]]
    populations = list(cfg["data"]["populations"])
    lineage_of = cfg["data"]["populations"]
    lineages = list(dict.fromkeys(lineage_of.values()))
    ten = key["population_key"].to_numpy()
    six = np.array([lineage_of[p] for p in ten])
    label_sets = {"ten_labels": (ten, populations), "six_lineages": (six, lineages)}
    parts = {arm: assignments[arm].to_numpy() for arm in ARMS}
    n = len(ten)

    # Per-population technical medians.
    technical = pd.read_csv(ROOT / "data" / "interim" / "technical.csv").set_index("cell_id").loc[assignments["cell_id"]]
    technical["population"] = ten
    pop_rows = []
    for p in populations:
        sub = technical[technical["population"] == p]
        pop_rows.append({"population": p, "lineage": lineage_of[p], "n_cells": len(sub), "median_total_counts": sub["total_counts"].median(), "median_n_genes": sub["n_genes"].median(), "median_pct_mito": sub["pct_mito"].median()})
    pd.DataFrame(pop_rows).to_csv(TABLES / "07_population_technical.csv", index=False, lineterminator="\n")

    # Observed ARI, paired bootstrap and permutation null.
    observed = {arm: {k: float(adjusted_rand_score(labels, parts[arm])) for k, (labels, _) in label_sets.items()} for arm in ARMS}
    n_boot, n_perm = ecfg["n_bootstrap"], ecfg["n_permutations"]
    boot = {(arm, k): np.empty(n_boot) for arm in ARMS for k in label_sets}
    gen = rng(cfg["rng_streams"]["bootstrap"])
    t0 = time.time()
    for b in range(n_boot):
        ix = gen.integers(0, n, size=n)
        for k, (labels, _) in label_sets.items():
            resampled = labels[ix]
            for arm in ARMS:
                boot[(arm, k)][b] = adjusted_rand_score(resampled, parts[arm][ix])
        if (b + 1) % 100 == 0:
            print(f"bootstrap {b + 1}/{n_boot} ({time.time() - t0:.0f} s)", flush=True)
    null = {(arm, k): np.empty(n_perm) for arm in ARMS for k in label_sets}
    gen = rng(cfg["rng_streams"]["permutation"])
    t0 = time.time()
    for i in range(n_perm):
        order = gen.permutation(n)
        for k, (labels, _) in label_sets.items():
            shuffled = labels[order]
            for arm in ARMS:
                null[(arm, k)][i] = adjusted_rand_score(shuffled, parts[arm])
        if (i + 1) % 100 == 0:
            print(f"permutation {i + 1}/{n_perm} ({time.time() - t0:.0f} s)", flush=True)
    np.savez_compressed(ROOT / "data" / "interim" / "07_resamples.npz", **{f"boot__{a}__{k}": v for (a, k), v in boot.items()}, **{f"null__{a}__{k}": v for (a, k), v in null.items()})

    ari = {}
    for arm in ARMS:
        ari[arm] = {}
        for k in label_sets:
            low, high = ci_bounds(boot[(arm, k)], ecfg["ci_level"])
            nulls = null[(arm, k)]
            n_low, n_high = ci_bounds(nulls, ecfg["ci_level"])
            ari[arm][k] = {
                "observed": observed[arm][k],
                "ci_low": low,
                "ci_high": high,
                "null_mean": float(nulls.mean()),
                "null_p2_5": n_low,
                "null_p97_5": n_high,
                "null_max": float(nulls.max()),
                "p_value": float((1 + (nulls >= observed[arm][k]).sum()) / (n_perm + 1)),
            }

    differences = {}
    for name, (a, b, kind) in PAIRS.items():
        differences[name] = {}
        for k in label_sets:
            diff = boot[(a, k)] - boot[(b, k)]
            low, high = ci_bounds(diff, ecfg["ci_level"])
            differences[name][k] = {"arm": a, "minus": b, "observed": observed[a][k] - observed[b][k], "ci_low": low, "ci_high": high, "verdict": verdict(kind, low, high)}

    # Contingency tables, recovery, purity and truth per cluster.
    recovery_purity, truth = {}, {}
    for arm in ARMS:
        recovery_purity[arm] = {}
        rp_rows = []
        for k, (labels, order) in label_sets.items():
            table = contingency(labels, parts[arm], order)
            table.to_csv(TABLES / f"07_contingency_{arm}_{k}.csv", lineterminator="\n")
            stats = [best_cluster(table, label) for label in order]
            recovery_purity[arm][k] = stats
            rp_rows += [dict(s, key=k) for s in stats]
        pd.DataFrame(rp_rows).to_csv(TABLES / f"07_recovery_purity_{arm}.csv", index=False, lineterminator="\n")
        lineage_table = contingency(six, parts[arm], lineages)
        pop_table = contingency(ten, parts[arm], populations)
        rows = cluster_truth(lineage_table, lineages, ecfg["unresolved_share"])
        for r in rows:
            counts = pop_table[r["cluster"]]
            r["plurality_population"] = str(counts.idxmax())
            r["plurality_population_share"] = float(counts.max() / counts.sum())
        truth[arm] = rows
        pd.DataFrame(rows).to_csv(TABLES / f"07_cluster_truth_{arm}.csv", index=False, lineterminator="\n")

    # Naming scores: rule-based names for three partitions, human names for the primary partition.
    naming = {}
    for part in NAMED_PARTITIONS:
        rule = pd.read_csv(TABLES / f"06_rule_names_{part}.csv")
        score, frame = naming_score(dict(zip(rule["cluster"].astype(int), rule["name"])), truth[part])
        naming[f"rule_{part}"] = score
        frame.to_csv(TABLES / f"07_naming_rule_{part}.csv", index=False, lineterminator="\n")
    human_rows = parse_human_naming()
    score, frame = naming_score({int(r["cluster"]): r["lineage_call"] for r in human_rows}, truth["kmeans_primary"])
    finer = {int(r["cluster"]): r["finer_call"] for r in human_rows}
    plurality_population = {r["cluster"]: r["plurality_population"] for r in truth["kmeans_primary"]}
    frame["finer_call"] = frame["cluster"].map(finer)
    frame["plurality_population_name"] = frame["cluster"].map(lambda c: TENX_NAMES[plurality_population[c]])
    frame["finer_matches_plurality_population"] = frame["finer_call"] == frame["plurality_population_name"]
    score["descriptive_n_finer_named"] = int((frame["finer_call"] != "none").sum())
    score["descriptive_n_finer_match_plurality_population"] = int(frame["finer_matches_plurality_population"].sum())
    naming["human_kmeans_primary"] = score
    frame.to_csv(TABLES / "07_naming_human_kmeans_primary.csv", index=False, lineterminator="\n")

    # Stability reruns against the labels.
    stab = {}
    for method in ["kmeans", "leiden"]:
        columns = [c for c in stability.columns if c.startswith(f"{method}_seed")]
        stab[method] = {}
        for k, (labels, _) in label_sets.items():
            values = np.array([adjusted_rand_score(labels, stability[c].to_numpy()) for c in columns])
            stab[method][k] = {"min": float(values.min()), "median": float(np.median(values)), "max": float(values.max()), "n": int(len(values))}

    write_metrics(
        "07_evaluate",
        {
            "n_cells": n,
            "populations": {r["population"]: {k: v for k, v in r.items() if k != "population"} for r in pop_rows},
            "ari": ari,
            "paired_differences": differences,
            "recovery_purity": recovery_purity,
            "cluster_truth": truth,
            "naming_scores": naming,
            "stability_vs_labels": stab,
            "runtime_seconds": time.time() - started,
        },
    )
    print("evaluation written")


if __name__ == "__main__":
    main()
