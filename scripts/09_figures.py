"""Script 09: figures drawn after the answer key opens.

UMAP coordinates appear only as pictures. Every figure is listed with its
SHA-256 in 09_figures.json.
"""
import numpy as np
import pandas as pd

from common import (
    NULL_COLOUR,
    PALETTE,
    ROOT,
    apply_style,
    check_free_space,
    check_posts,
    cluster_colours,
    load_config,
    open_answer_key,
    read_metrics,
    sequential_cmap,
    sha256_file,
    write_metrics,
)

FIG = ROOT / "figures"
TABLES = ROOT / "results" / "tables"
ARMS = ["kmeans_primary", "kmeans_oracle_k10", "leiden_primary", "technical_k_chosen", "technical_k10"]
ARM_TITLES = {
    "kmeans_primary": "k-means, chosen k",
    "kmeans_oracle_k10": "k-means, k = 10 (label-informed)",
    "leiden_primary": "Leiden, chosen resolution",
    "technical_k_chosen": "technical baseline, chosen k",
    "technical_k10": "technical baseline, k = 10",
}
SHORT = {
    "kmeans_primary": "k-means\nchosen k",
    "kmeans_oracle_k10": "k-means\nk = 10",
    "leiden_primary": "Leiden\nchosen res.",
    "technical_k_chosen": "technical\nchosen k",
    "technical_k10": "technical\nk = 10",
}
KEY_TITLES = {"ten_labels": "ten labels", "six_lineages": "six lineages"}


def main():
    check_posts()
    check_free_space()
    cfg = load_config()
    plt = apply_style()
    ev = read_metrics("07_evaluate")
    fm = read_metrics("08_failure_modes")
    populations = list(cfg["data"]["populations"])
    written = []

    def save(fig, name):
        path = FIG / name
        fig.savefig(path)
        plt.close(fig)
        written.append(path)

    key = open_answer_key()
    assignments = pd.read_csv(ROOT / "results" / "assignments" / "assignments_main.csv.gz")
    ten = key.set_index("cell_id").loc[assignments["cell_id"], "population_key"].to_numpy()

    # UMAP by population.
    umap = np.load(ROOT / "data" / "interim" / "05_umap.npz", allow_pickle=True)["X_umap"]
    colours = cluster_colours(len(populations))
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for colour, p in zip(colours, populations):
        mask = ten == p
        ax.scatter(umap[mask, 0], umap[mask, 1], s=0.3, color=colour, linewidths=0, rasterized=True, label=p)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1 (picture only)")
    ax.set_ylabel("UMAP 2 (picture only)")
    ax.set_title("UMAP coloured by enriched population (answer key)")
    ax.legend(markerscale=12, fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    save(fig, "09_umap_by_population.png")

    # Row-normalised contingency heatmaps.
    for arm in ARMS:
        for k in ["ten_labels", "six_lineages"]:
            table = pd.read_csv(TABLES / f"07_contingency_{arm}_{k}.csv", index_col=0)
            shares = table.div(table.sum(axis=1), axis=0)
            fig, ax = plt.subplots(figsize=(0.45 * shares.shape[1] + 3, 0.4 * shares.shape[0] + 1.5))
            image = ax.imshow(shares.to_numpy(), cmap=sequential_cmap(), vmin=0, vmax=1, aspect="auto")
            ax.set_xticks(range(shares.shape[1]))
            ax.set_xticklabels(shares.columns)
            ax.set_yticks(range(shares.shape[0]))
            ax.set_yticklabels(shares.index)
            ax.set_xlabel("cluster")
            ax.set_title(f"Share of each label per cluster: {ARM_TITLES[arm]}, {KEY_TITLES[k]}")
            fig.colorbar(image, ax=ax, label="share of the label's cells", shrink=0.8)
            save(fig, f"09_contingency_{arm}_{k}.png")

    # ARI summary with confidence intervals and the null band.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, k in zip(axes, ["ten_labels", "six_lineages"]):
        for i, arm in enumerate(ARMS):
            a = ev["ari"][arm][k]
            ax.fill_between([i - 0.3, i + 0.3], a["null_p2_5"], a["null_p97_5"], color=NULL_COLOUR, alpha=0.35, linewidth=0)
            ax.errorbar(i, a["observed"], yerr=[[a["observed"] - a["ci_low"]], [a["ci_high"] - a["observed"]]], fmt="o", color=PALETTE[i], capsize=3)
        ax.axhline(0, color=NULL_COLOUR, linewidth=0.8)
        ax.set_xticks(range(len(ARMS)))
        ax.set_xticklabels([SHORT[a] for a in ARMS], fontsize=7)
        ax.set_title(f"Adjusted Rand index, {KEY_TITLES[k]}")
    axes[0].set_ylabel("ARI (95% bootstrap CI; grey band: permutation null)")
    save(fig, "09_ari_summary.png")

    # Recovery and purity by population.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    width = 0.26
    for ax, stat in zip(axes, ["recovery", "purity"]):
        for j, arm in enumerate(["kmeans_primary", "kmeans_oracle_k10", "leiden_primary"]):
            values = {r["label"]: r[stat] for r in ev["recovery_purity"][arm]["ten_labels"]}
            ax.bar(np.arange(len(populations)) + (j - 1) * width, [values[p] for p in populations], width=width, color=PALETTE[j * 2], label=ARM_TITLES[arm])
        ax.set_xticks(range(len(populations)))
        ax.set_xticklabels(populations, rotation=60, ha="right")
        ax.set_title(f"{stat.capitalize()} of each label in its best cluster")
    axes[0].axhline(0.8, color=NULL_COLOUR, linestyle="--", linewidth=0.8)
    axes[0].legend(frameon=False, fontsize=7)
    save(fig, "09_recovery_purity.png")

    # Technical metrics by population.
    technical = pd.read_csv(ROOT / "data" / "interim" / "technical.csv").set_index("cell_id").loc[assignments["cell_id"]]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, column, label in zip(axes, ["log10_total_counts", "log10_n_genes", "pct_mito"], ["log10 total counts", "log10 genes detected", "percent mitochondrial counts"]):
        data = [technical[column].to_numpy()[ten == p] for p in populations]
        parts = ax.boxplot(data, showfliers=False, patch_artist=True)
        for patch, colour in zip(parts["boxes"], colours):
            patch.set_facecolor(colour)
        ax.set_xticks(range(1, len(populations) + 1))
        ax.set_xticklabels(populations, rotation=60, ha="right")
        ax.set_title(label)
    save(fig, "09_technical_by_population.png")

    # Nested-label test.
    sub = fm["FM1"]["sub_ari"]
    arms = list(sub)
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    x = np.arange(len(arms))
    nested, distinct = cfg["diagnostics"]["nested_pair"], cfg["diagnostics"]["distinct_pair"]
    ax.bar(x - 0.18, [sub[a]["nested"] for a in arms], width=0.36, color=PALETTE[5], label=f"nested: {nested[0]} vs {nested[1]}")
    ax.bar(x + 0.18, [sub[a]["distinct"] for a in arms], width=0.36, color=PALETTE[0], label=f"distinct: {distinct[0]} vs {distinct[1]}")
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT[a] for a in arms], fontsize=7)
    ax.set_ylabel("ARI within the pair")
    ax.set_title("Nested-label test: nested should separate less than distinct")
    ax.legend(frameon=False, fontsize=7)
    save(fig, "09_nested_label_test.png")

    # Label-gene detection by population.
    detection = pd.read_csv(TABLES / "08_label_gene_detection.csv", index_col=0)
    fig, ax = plt.subplots(figsize=(7, 5))
    image = ax.imshow(detection.to_numpy(), cmap=sequential_cmap(), vmin=0, vmax=1, aspect="auto")
    for i in range(detection.shape[0]):
        for j in range(detection.shape[1]):
            value = detection.iat[i, j]
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=6, color="white" if value > 0.6 else "#16222E")
    ax.set_xticks(range(detection.shape[1]))
    ax.set_xticklabels(detection.columns, rotation=45, ha="right")
    ax.set_yticks(range(detection.shape[0]))
    ax.set_yticklabels(detection.index)
    ax.set_title("Share of cells with at least one count of each label gene")
    fig.colorbar(image, ax=ax, shrink=0.8)
    save(fig, "09_label_gene_detection.png")

    write_metrics("09_figures", {"files": [{"path": p.relative_to(ROOT).as_posix(), "sha256": sha256_file(p)} for p in written], "n_files": len(written)})
    print(f"{len(written)} figures written")


if __name__ == "__main__":
    main()
