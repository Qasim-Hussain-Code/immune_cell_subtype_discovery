"""Script 06: name every cluster from marker genes, before the answer key opens.

Rule-based names follow the fixed rule in config.yaml (T score first, then CD8
against CD4, otherwise the best lineage score above zero). Wilcoxon marker
tables, the blind figures, the human naming template and the blind report are
written here.
"""
import time

import numpy as np
import pandas as pd
import scanpy as sc

from common import (
    NAMED_PARTITIONS,
    NULL_COLOUR,
    ROOT,
    apply_style,
    check_free_space,
    check_posts,
    cluster_colours,
    load_config,
    read_metrics,
    sequential_cmap,
    write_metrics,
    write_text,
)

INTERIM = ROOT / "data" / "interim"

# Allowed finer calls for the human naming template: the ten dataset page names, or "none".
FINER_CALLS = [
    "CD19+ B Cells",
    "CD14+ Monocytes",
    "CD34+ Cells",
    "CD4+ Helper T Cells",
    "CD4+/CD25+ Regulatory T Cells",
    "CD4+/CD45RA+/CD25- Naive T cells",
    "CD4+/CD45RO+ Memory T Cells",
    "CD56+ Natural Killer Cells",
    "CD8+ Cytotoxic T cells",
    "CD8+/CD45RA+ Naive Cytotoxic T Cells",
]
PARTITION_TITLES = {"kmeans_primary": "k-means, chosen k", "kmeans_oracle_k10": "k-means, k = 10 (label-informed)", "leiden_primary": "Leiden, chosen resolution"}


def naming_genes(ncfg):
    genes = list(ncfg["t_markers"])
    for markers in ncfg["lineage_markers"].values():
        genes += [g for g in markers if g not in genes]
    genes += [g for g in ncfg["cd8_markers"] + [ncfg["cd4_marker"]] if g not in genes]
    return genes


def gene_columns(var, symbols):
    column_symbols = var["gene_symbols"].to_numpy()
    found, missing, repeated = {}, [], []
    for symbol in symbols:
        hits = np.flatnonzero(column_symbols == symbol)
        if len(hits) == 0:
            missing.append(symbol)
        elif len(hits) > 1:
            repeated.append(symbol)
        else:
            found[symbol] = int(hits[0])
    if missing or repeated:
        raise SystemExit(f"Stop: marker symbols missing {missing} or repeated {repeated} in genes.tsv; log this and ask the owner.")
    return found


def rule_names(x, labels, genes, ncfg):
    col = {g: i for i, g in enumerate(genes)}
    sd = x.std(axis=0)
    rows = []
    for cluster in np.unique(labels):
        inside = labels == cluster
        mean_in = x[inside].mean(axis=0)
        mean_out = x[~inside].mean(axis=0)
        d = np.divide(mean_in - mean_out, sd, out=np.zeros_like(sd), where=sd > 0)
        t_score = float(np.mean([d[col[g]] for g in ncfg["t_markers"]]))
        lineage = {name: float(np.mean([d[col[g]] for g in markers])) for name, markers in ncfg["lineage_markers"].items()}
        cd8 = float(np.mean([mean_in[col[g]] for g in ncfg["cd8_markers"]]))
        cd4 = float(mean_in[col[ncfg["cd4_marker"]]])
        best = max(lineage, key=lineage.get)
        if t_score > 0:
            branch = "T"
            name = "CD8 T cells" if cd8 > cd4 else "CD4 T cells"
        else:
            branch = "lineage"
            name = best if lineage[best] > 0 else "unresolved"
        row = {"cluster": int(cluster), "n_cells": int(inside.sum()), "branch": branch, "name": name, "t_score": t_score}
        row.update({f"score[{k}]": v for k, v in lineage.items()})
        row["lineage_tie"] = int(sum(v == lineage[best] for v in lineage.values()) > 1)
        row["cd8_mean_x"] = cd8
        row["cd4_mean_x"] = cd4
        row.update({f"d[{g}]": float(d[col[g]]) for g in genes})
        row.update({f"mean_x[{g}]": float(mean_in[col[g]]) for g in genes})
        rows.append(row)
    return pd.DataFrame(rows)


def wilcoxon_markers(lognorm, labels, top_n):
    groups = [str(c) for c in np.unique(labels)]
    lognorm.obs["group"] = pd.Categorical([str(v) for v in labels], categories=groups)
    sc.tl.rank_genes_groups(lognorm, groupby="group", method="wilcoxon", use_raw=False, pts=True, n_genes=top_n, key_added="rgg")
    table = sc.get.rank_genes_groups_df(lognorm, group=None, key="rgg")
    symbol = lognorm.var["gene_symbols"]
    table = table.rename(columns={"group": "cluster", "names": "gene_id"})
    table["cluster"] = table["cluster"].astype(int)
    table["rank"] = table.groupby("cluster").cumcount() + 1
    table["gene_symbol"] = symbol.loc[table["gene_id"]].to_numpy()
    del lognorm.obs["group"]
    del lognorm.uns["rgg"]
    columns = ["cluster", "rank", "gene_id", "gene_symbol", "scores", "logfoldchanges", "pvals", "pvals_adj", "pct_nz_group", "pct_nz_reference"]
    return table[[c for c in columns if c in table.columns]].sort_values(["cluster", "rank"]).reset_index(drop=True)


def scatter_clusters(plt, umap, labels, names, title, path):
    clusters = np.unique(labels)
    colours = cluster_colours(len(clusters))
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for colour, cluster in zip(colours, clusters):
        mask = labels == cluster
        ax.scatter(umap[mask, 0], umap[mask, 1], s=0.3, color=colour, linewidths=0, rasterized=True, label=f"{cluster}: {names[cluster]}")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1 (picture only)")
    ax.set_ylabel("UMAP 2 (picture only)")
    ax.set_title(title)
    ax.legend(markerscale=12, fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.savefig(path)
    plt.close(fig)


def dot_plot(plt, x, labels, genes, names, title, path):
    clusters = np.unique(labels)
    frac = np.array([(x[labels == c] > 0).mean(axis=0) for c in clusters])
    mean = np.array([x[labels == c].mean(axis=0) for c in clusters])
    fig, ax = plt.subplots(figsize=(0.38 * len(genes) + 2.5, 0.32 * len(clusters) + 1.4))
    yy, xx = np.meshgrid(np.arange(len(clusters)), np.arange(len(genes)), indexing="ij")
    points = ax.scatter(xx.ravel(), yy.ravel(), s=frac.ravel() * 120, c=mean.ravel(), cmap=sequential_cmap(), edgecolors=NULL_COLOUR, linewidths=0.3)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90)
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels([f"{c}: {names[c]}" for c in clusters])
    ax.invert_yaxis()
    ax.set_title(title)
    fig.colorbar(points, ax=ax, label="mean log-normalised expression", shrink=0.6)
    fig.savefig(path)
    plt.close(fig)


def fmt_int(value):
    return f"{int(value):,}"


def kmeans_rerun_solutions(cfg, primary, stability):
    """Distinct k-means solutions among the primary run and the 20 reruns, with inertia and silhouette (label-free)."""
    from sklearn.metrics import adjusted_rand_score, silhouette_score

    stored = np.load(INTERIM / "pca.npz", allow_pickle=True)
    x = stored["X_pca"].astype(np.float64)
    position = pd.Series(np.arange(len(stored["cell_id"])), index=stored["cell_id"])
    sub = position.loc[pd.read_csv(ROOT / "results" / "tables" / "04_silhouette_subsample.csv")["cell_id"]].to_numpy()
    runs = [(cfg["kmeans"]["random_state"], primary)] + [(s, stability[f"kmeans_seed{s}"].to_numpy()) for s in cfg["kmeans"]["stability_seeds"]]
    groups = []
    for seed, labels in runs:
        for group in groups:
            if adjusted_rand_score(group["labels"], labels) == 1.0:
                group["seeds"].append(seed)
                break
        else:
            groups.append({"labels": labels, "seeds": [seed]})
    out = []
    for group in groups:
        labels = group["labels"]
        inertia = float(sum(((x[labels == c] - x[labels == c].mean(axis=0)) ** 2).sum() for c in np.unique(labels)))
        out.append(
            {
                "seeds": group["seeds"],
                "n_reruns": int(sum(s != cfg["kmeans"]["random_state"] for s in group["seeds"])),
                "includes_primary": cfg["kmeans"]["random_state"] in group["seeds"],
                "sizes": sorted((int(v) for v in np.bincount(labels)), reverse=True),
                "inertia": inertia,
                "silhouette_subsample": float(silhouette_score(x[sub], labels[sub], metric="euclidean")),
            }
        )
    return sorted(out, key=lambda g: g["inertia"])


def main():
    check_posts()
    check_free_space()
    sc.settings.verbosity = 0
    cfg = load_config()
    ncfg = cfg["naming"]
    started = time.time()

    lognorm = sc.read_h5ad(INTERIM / "lognorm.h5ad")
    assignments = pd.read_csv(ROOT / "results" / "assignments" / "assignments_main.csv.gz")
    if not np.array_equal(assignments["cell_id"].to_numpy(), lognorm.obs_names.to_numpy()):
        raise SystemExit("Stop: assignment order differs from the log-normalised matrix")
    genes = naming_genes(ncfg)
    cols = gene_columns(lognorm.var, genes)
    x = lognorm.X[:, [cols[g] for g in genes]].toarray().astype(np.float64)

    names, markers, top10 = {}, {}, {}
    marker_seconds = {}
    for part in NAMED_PARTITIONS:
        labels = assignments[part].to_numpy()
        table = rule_names(x, labels, genes, ncfg)
        table.to_csv(ROOT / "results" / "tables" / f"06_rule_names_{part}.csv", index=False, lineterminator="\n")
        names[part] = table
        t0 = time.time()
        marker_table = wilcoxon_markers(lognorm, labels, ncfg["wilcoxon_top_n"])
        marker_seconds[part] = time.time() - t0
        marker_table.to_csv(ROOT / "results" / "tables" / f"06_markers_{part}.csv", index=False, lineterminator="\n")
        markers[part] = marker_table
        top10[part] = {int(c): g.sort_values("rank")["gene_symbol"].head(10).tolist() for c, g in marker_table.groupby("cluster")}
        print(f"named and ranked markers for {part} ({marker_seconds[part]:.0f} s)", flush=True)

    # Blind figures.
    plt = apply_style()
    umap = np.load(INTERIM / "05_umap.npz", allow_pickle=True)["X_umap"]
    for part in NAMED_PARTITIONS:
        lookup = dict(zip(names[part]["cluster"], names[part]["name"]))
        labels = assignments[part].to_numpy()
        scatter_clusters(plt, umap, labels, lookup, f"UMAP by cluster: {PARTITION_TITLES[part]}", ROOT / "figures" / f"06_umap_{part}.png")
        dot_plot(plt, x, labels, genes, lookup, f"Naming genes by cluster: {PARTITION_TITLES[part]}", ROOT / "figures" / f"06_dotplot_{part}.png")

    technical = pd.read_csv(INTERIM / "technical.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, column, label in zip(axes, ["log10_total_counts", "pct_mito"], ["log10 total counts", "percent mitochondrial counts"]):
        points = ax.scatter(umap[:, 0], umap[:, 1], c=technical[column].to_numpy(), s=0.3, cmap=sequential_cmap(), linewidths=0, rasterized=True)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"UMAP coloured by {label}")
        fig.colorbar(points, ax=ax, shrink=0.7)
    fig.savefig(ROOT / "figures" / "06_umap_technical.png")
    plt.close(fig)

    clustering = read_metrics("05_cluster")
    stability_table = pd.read_csv(ROOT / "results" / "assignments" / "assignments_stability.csv.gz")
    solutions = kmeans_rerun_solutions(cfg, assignments["kmeans_primary"].to_numpy(), stability_table)
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    series = [
        ("k-means\npairwise", clustering["stability"]["kmeans"]["pairwise_values"]),
        ("k-means\nto primary", clustering["stability"]["kmeans"]["to_primary_values"]),
        ("Leiden\npairwise", clustering["stability"]["leiden"]["pairwise_values"]),
        ("Leiden\nto primary", clustering["stability"]["leiden"]["to_primary_values"]),
    ]
    jitter = np.random.default_rng(0)
    for i, (label, values) in enumerate(series):
        ax.scatter(i + jitter.uniform(-0.18, 0.18, len(values)), values, s=8, color=["#17456B", "#2E6DA4", "#21A6B8", "#4CAF7D"][i], alpha=0.7, linewidths=0)
    ax.set_xticks(range(len(series)))
    ax.set_xticklabels([s[0] for s in series])
    ax.set_ylabel("adjusted Rand index")
    ax.set_title("Stability over 20 reruns (seeds 101 to 120)")
    fig.savefig(ROOT / "figures" / "06_stability.png")
    plt.close(fig)

    # Human naming template for the primary partition.
    primary = names["kmeans_primary"]
    lines = [
        "# Human naming: primary k-means partition",
        "",
        "Filled by the owner from marker genes before the answer key opens. Cluster IDs follow column `kmeans_primary` of `results/assignments/assignments_main.csv.gz`, largest cluster first. Marker tables: `results/tables/06_markers_kmeans_primary.csv`. Figures: `figures/06_dotplot_kmeans_primary.png`, `figures/06_umap_kmeans_primary.png`.",
        "",
        "Lineage call, exactly one of: " + ", ".join(f"`{v}`" for v in ncfg["human_vocabulary"]) + ".",
        "",
        "Finer call, exactly one of the ten dataset page names or `none`: " + ", ".join(f"`{v}`" for v in FINER_CALLS) + ".",
        "",
        "Note is free text and may be left empty.",
        "",
        "| cluster | cells | rule-based name | top 10 markers | lineage call | finer call | note |",
        "|---|---|---|---|---|---|---|",
    ]
    for _, row in primary.iterrows():
        genes_text = ", ".join(top10["kmeans_primary"][int(row["cluster"])])
        lines.append(f"| {int(row['cluster'])} | {fmt_int(row['n_cells'])} | {row['name']} | {genes_text} |  |  |  |")
    write_text("notes/human_naming.md", "\n".join(lines) + "\n")

    # Blind report, from the metrics files only.
    selection = read_metrics("04_select_k")
    pre = read_metrics("03_preprocess")
    report = [
        "# Blind report",
        "",
        "Written by script 06 before the answer key opens. No population label has been read.",
        "",
        f"Cells: {fmt_int(selection['n_cells'])}. Highly variable genes: {pre['n_hvg']}. Variance held by {pre['n_pcs']} components: {pre['cumulative_variance_all_pcs']:.4f}.",
        "",
        "## Number of groups",
        "",
        f"Chosen k: {selection['chosen_k']} (mean silhouette {selection['chosen_k_silhouette']:.4f}).",
        "",
        "| k | mean silhouette | inertia (not used) |",
        "|---|---|---|",
    ]
    report += [f"| {r['k']} | {r['silhouette']:.4f} | {r['inertia']:.1f} |" for r in selection["kmeans_curve"]]
    report += [
        "",
        f"Chosen Leiden resolution: {selection['chosen_resolution']:.1f}, giving {selection['chosen_resolution_n_clusters']} clusters (mean silhouette {selection['chosen_resolution_silhouette']:.4f}).",
        "",
        "| resolution | clusters | mean silhouette |",
        "|---|---|---|",
    ]
    report += [f"| {r['resolution']:.1f} | {r['n_clusters']} | {'n/a' if r['silhouette'] is None else format(r['silhouette'], '.4f')} |" for r in selection["leiden_curve"]]
    report += ["", "## Cluster sizes and rule-based names", ""]
    for part in NAMED_PARTITIONS:
        report += [f"### {PARTITION_TITLES[part]}", "", "| cluster | cells | rule-based name | top 10 markers |", "|---|---|---|---|"]
        for _, row in names[part].iterrows():
            report.append(f"| {int(row['cluster'])} | {fmt_int(row['n_cells'])} | {row['name']} | {', '.join(top10[part][int(row['cluster'])])} |")
        report += ["", f"Marker table: `results/tables/06_markers_{part}.csv`. Rule scores: `results/tables/06_rule_names_{part}.csv`.", ""]
    report += ["### Camera-settings baseline", "", "| partition | cluster sizes |", "|---|---|"]
    for part in ["technical_k_chosen", "technical_k10"]:
        report.append(f"| {part} | {', '.join(fmt_int(s['n_cells']) for s in clustering['cluster_sizes'][part])} |")
    stab = clustering["stability"]
    report += [
        "",
        "## Stability, agreement and silhouette",
        "",
        "| measure | min | median | max |",
        "|---|---|---|---|",
    ]
    for method in ["kmeans", "leiden"]:
        for measure in ["pairwise_ari", "ari_to_primary"]:
            s = stab[method][measure]
            report.append(f"| {method} {measure.replace('_', ' ')} | {s['min']:.4f} | {s['median']:.4f} | {s['max']:.4f} |")
    report += ["", "Distinct k-means solutions at the chosen k, among the primary run (seed 7) and the 20 reruns, tightest first:", "", "| cluster sizes | inertia | silhouette (subsample) | seeds |", "|---|---|---|---|"]
    report += [f"| {', '.join(fmt_int(v) for v in s['sizes'])} | {s['inertia']:.1f} | {s['silhouette_subsample']:.4f} | {', '.join(str(v) for v in s['seeds'])} |" for s in solutions]
    report += ["", "| label-free agreement | ARI |", "|---|---|"]
    report += [f"| {k.replace('_', ' ')} | {v:.4f} |" for k, v in clustering["agreement_ari"].items()]
    report += ["", "| partition | mean silhouette | share of cells with negative silhouette |", "|---|---|---|"]
    report += [f"| {p} | {clustering['silhouette'][p]['mean_silhouette']:.4f} | {clustering['silhouette'][p]['negative_share']:.4f} |" for p in NAMED_PARTITIONS]
    report += [
        "",
        f"Refits match script 04: primary {clustering['primary_identical_to_04']}, Leiden {clustering['leiden_identical_to_04']}. Oracle identical to primary: {clustering['oracle_identical_to_primary']}.",
        "",
        "## Runtimes (seconds)",
        "",
        f"Preprocessing {pre['runtime_seconds']:.0f}; k and resolution selection {selection['runtime_seconds']:.0f}; clustering {clustering['runtime_seconds']:.0f}; naming {time.time() - started:.0f}.",
        "",
    ]
    write_text("notes/blind_report.md", "\n".join(report))

    write_metrics(
        "06_blind_naming",
        {
            "naming_genes": genes,
            "rule_names": {part: [{"cluster": int(r["cluster"]), "n_cells": int(r["n_cells"]), "name": r["name"], "branch": r["branch"]} for _, r in names[part].iterrows()] for part in NAMED_PARTITIONS},
            "rule_name_counts": {part: {k: int(v) for k, v in names[part]["name"].value_counts().items()} for part in NAMED_PARTITIONS},
            "n_rule_unresolved": {part: int((names[part]["name"] == "unresolved").sum()) for part in NAMED_PARTITIONS},
            "n_lineage_ties": {part: int(names[part]["lineage_tie"].sum()) for part in NAMED_PARTITIONS},
            "top10_markers": {part: {str(c): g for c, g in top10[part].items()} for part in NAMED_PARTITIONS},
            "kmeans_rerun_solutions": solutions,
            "n_kmeans_rerun_solutions": len(solutions),
            "marker_seconds": marker_seconds,
            "runtime_seconds": time.time() - started,
        },
    )
    print("blind naming written")


if __name__ == "__main__":
    main()
