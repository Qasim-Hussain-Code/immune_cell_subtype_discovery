"""Shared helpers for the pipeline.

Configuration, random number streams, hashing, the metrics writer, figure
style and git checks live here, together with open_answer_key(), the single
guarded route to the population labels. The clustering scripts import only
helpers that cannot reach the labels, which tests/test_label_quarantine.py
verifies.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
POST_DAYS = range(61, 67)
MIN_FREE_GB = 1.0

# Partitions scored by the rule-based names, in a fixed order.
NAMED_PARTITIONS = ["kmeans_primary", "kmeans_oracle_k10", "leiden_primary"]

# Series palette and text colours for every figure.
PALETTE = ["#17456B", "#2E6DA4", "#21A6B8", "#4CAF7D", "#F5B841", "#F07D22", "#E85D2F", "#7A8899"]
TEXT_COLOUR = "#16222E"
NULL_COLOUR = "#7A8899"

FREEZE_SUBJECT = "freeze_cluster_assignments"
FREEZE_MANIFEST = "provenance/freeze_manifest.json"
HUMAN_NAMING = "notes/human_naming.md"

# The ten names used on the 10x Genomics dataset pages, keyed by population.
TENX_NAMES = {
    "b_cells": "CD19+ B Cells",
    "cd14_monocytes": "CD14+ Monocytes",
    "cd34": "CD34+ Cells",
    "cd4_t_helper": "CD4+ Helper T Cells",
    "regulatory_t": "CD4+/CD25+ Regulatory T Cells",
    "naive_t": "CD4+/CD45RA+/CD25- Naive T cells",
    "memory_t": "CD4+/CD45RO+ Memory T Cells",
    "cd56_nk": "CD56+ Natural Killer Cells",
    "cytotoxic_t": "CD8+ Cytotoxic T cells",
    "naive_cytotoxic": "CD8+/CD45RA+ Naive Cytotoxic T Cells",
}


# ---------------------------------------------------------------- config, RNG


def load_config():
    with open(ROOT / "config.yaml", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def rng(stream):
    """NumPy generator for one declared stream, seeded as [master_seed, stream]."""
    seed = load_config()["project"]["master_seed"]
    return np.random.default_rng([seed, stream])


# ---------------------------------------------------------------- files


def sha256_file(path, chunk=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def _plain(obj):
    """Convert NumPy scalars and containers to plain JSON types; NaN and inf become null."""
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_plain(v) for v in obj.tolist()]
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        value = float(obj)
        return value if math.isfinite(value) else None
    return obj


def write_metrics(name, payload):
    """Write results/metrics/<name>.json with sorted keys, indent 2 and full precision."""
    path = ROOT / "results" / "metrics" / f"{name}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(_plain(payload), fh, sort_keys=True, indent=2)
        fh.write("\n")
    return path


def read_metrics(name):
    with open(ROOT / "results" / "metrics" / f"{name}.json", encoding="utf-8") as fh:
        return json.load(fh)


def write_text(rel_path, text):
    path = ROOT / rel_path
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    return path


def gzip_csv(frame, rel_path):
    """Write a gzipped CSV with a fixed gzip timestamp, so equal content gives an equal hash."""
    path = ROOT / rel_path
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})
    return path


def relabel_by_size(labels):
    """Renumber clusters so the largest is 0; equal sizes keep their original order."""
    labels = np.asarray(labels)
    values, counts = np.unique(labels, return_counts=True)
    order = sorted(range(len(values)), key=lambda i: (-counts[i], values[i]))
    lookup = {values[i]: new for new, i in enumerate(order)}
    return np.array([lookup[v] for v in labels], dtype=np.int64)


# ---------------------------------------------------------------- guards


def check_posts():
    """Stop unless the six published posts, Days 61 to 66, are present and non-empty."""
    missing = []
    for day in POST_DAYS:
        path = ROOT / "posts" / f"day_{day}.md"
        if not path.is_file() or path.stat().st_size == 0:
            missing.append(path.relative_to(ROOT).as_posix())
    if missing:
        raise SystemExit("Stop: published posts missing or empty: " + ", ".join(missing) + ". They are the pre-registration and must be in place first.")


def check_free_space(min_gb=MIN_FREE_GB):
    free_gb = shutil.disk_usage(ROOT).free / 1e9
    if free_gb < min_gb:
        raise SystemExit(f"Stop: only {free_gb:.2f} GB free on the project drive; at least {min_gb} GB is required.")
    return free_gb


# ---------------------------------------------------------------- figures


def apply_style():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "text.color": TEXT_COLOUR,
            "axes.labelcolor": TEXT_COLOUR,
            "axes.edgecolor": TEXT_COLOUR,
            "axes.titlecolor": TEXT_COLOUR,
            "xtick.color": TEXT_COLOUR,
            "ytick.color": TEXT_COLOUR,
            "axes.prop_cycle": matplotlib.cycler(color=PALETTE),
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    return plt


def cluster_colours(n):
    """tab10 for up to ten groups, tab20 above that, extended with tab20b past twenty."""
    import matplotlib

    if n <= 10:
        base = list(matplotlib.colormaps["tab10"].colors)
    else:
        base = list(matplotlib.colormaps["tab20"].colors) + list(matplotlib.colormaps["tab20b"].colors)
    return [base[i % len(base)] for i in range(n)]


def sequential_cmap():
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("series", ["#F4F6F8", PALETTE[2], PALETTE[0]])


# ---------------------------------------------------------------- git


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)


def find_commit(subject):
    """Hash of the most recent commit on HEAD's history whose subject equals `subject`."""
    result = git("log", "--format=%H%x09%s", "HEAD")
    if result.returncode != 0:
        return None
    for line in result.stdout.splitlines():
        sha, _, subj = line.partition("\t")
        if subj.strip() == subject:
            return sha
    return None


def is_ancestor(ancestor, descendant="HEAD"):
    return git("merge-base", "--is-ancestor", ancestor, descendant).returncode == 0


def freeze_file_list():
    files = [
        "results/assignments/assignments_main.csv.gz",
        "results/assignments/assignments_stability.csv.gz",
        "results/metrics/04_select_k.json",
        "results/metrics/05_cluster.json",
        "results/metrics/06_blind_naming.json",
    ]
    files += [f"results/tables/06_rule_names_{p}.csv" for p in NAMED_PARTITIONS]
    files += ["notes/blind_report.md", "config.yaml"]
    return files


def write_freeze_manifest():
    parent = git("rev-parse", "HEAD").stdout.strip()
    manifest = {
        "parent_commit": parent,
        "files": {rel: sha256_file(ROOT / rel) for rel in freeze_file_list()},
    }
    path = ROOT / FREEZE_MANIFEST
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(manifest, fh, sort_keys=True, indent=2)
        fh.write("\n")
    return path


def verify_freeze_manifest():
    """Return a list of files whose current hash differs from the manifest."""
    with open(ROOT / FREEZE_MANIFEST, encoding="utf-8") as fh:
        manifest = json.load(fh)
    bad = []
    for rel, digest in manifest["files"].items():
        path = ROOT / rel
        if not path.is_file() or sha256_file(path) != digest:
            bad.append(rel)
    return bad


# ---------------------------------------------------------------- second naming


def parse_human_naming(path=None):
    """Rows of the naming table as dicts keyed by the template's column names."""
    path = Path(path) if path else ROOT / HUMAN_NAMING
    columns = ["cluster", "cells", "rule_name", "top_markers", "lineage_call", "finer_call", "note"]
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if not cells or not re.fullmatch(r"\d+", cells[0]):
            continue
        if len(cells) > len(columns):
            cells = cells[: len(columns) - 1] + [" | ".join(cells[len(columns) - 1 :])]
        cells += [""] * (len(columns) - len(cells))
        rows.append(dict(zip(columns, cells)))
    return rows


def validate_human_naming(path=None):
    """Problems with the naming table; an empty list means it is complete and valid."""
    cfg = load_config()
    vocabulary = set(cfg["naming"]["human_vocabulary"])
    finer_allowed = set(TENX_NAMES.values()) | {"none"}
    problems = []
    rows = parse_human_naming(path)
    if not rows:
        return ["no cluster rows found"]
    blind = read_metrics("06_blind_naming")
    expected = sorted(int(r["cluster"]) for r in blind["rule_names"]["kmeans_primary"])
    found = sorted(int(r["cluster"]) for r in rows)
    if found != expected:
        problems.append(f"cluster rows {found} do not match the primary partition clusters {expected}")
    for r in rows:
        if r["lineage_call"] not in vocabulary:
            problems.append(f"cluster {r['cluster']}: lineage call '{r['lineage_call']}' is not in the vocabulary")
        if r["finer_call"] not in finer_allowed:
            problems.append(f"cluster {r['cluster']}: finer call '{r['finer_call']}' is not a 10x name or 'none'")
    return problems


def human_naming_committed_after(freeze_sha):
    last = git("log", "-1", "--format=%H", "--", HUMAN_NAMING).stdout.strip()
    if not last or last == freeze_sha or not is_ancestor(freeze_sha, last):
        return False
    return git("status", "--porcelain", "--", HUMAN_NAMING).stdout.strip() == ""


# ---------------------------------------------------------------- scoring (scripts 07 and later)


def contingency(labels, clusters, label_order):
    """Counts n(L, c): rows are labels in `label_order`, columns are cluster IDs in ascending order."""
    import pandas as pd

    table = pd.crosstab(pd.Categorical(labels, categories=label_order), np.asarray(clusters), dropna=False)
    table.index = list(label_order)
    table.index.name = "label"
    table.columns = [int(c) for c in table.columns]
    table.columns.name = "cluster"
    return table[sorted(table.columns)]


def best_cluster(table, label):
    """c*(L), recovery(L) and purity(L); ties in n(L, c) go to the lowest cluster ID."""
    row = table.loc[label].to_numpy()
    idx = int(np.argmax(row))
    cluster = int(table.columns[idx])
    n_label = row.sum()
    n_cluster = table[cluster].sum()
    return {
        "label": label,
        "n_label": int(n_label),
        "best_cluster": cluster,
        "n_label_in_best": int(row[idx]),
        "recovery": float(row[idx] / n_label) if n_label else None,
        "purity": float(row[idx] / n_cluster) if n_cluster else None,
    }


def overlap(table, a, b):
    """Sum over clusters of min(n(A, c)/n(A), n(B, c)/n(B)); 0 is fully separated, 1 an identical spread."""
    ra = table.loc[a].to_numpy(dtype=float)
    rb = table.loc[b].to_numpy(dtype=float)
    return float(np.minimum(ra / ra.sum(), rb / rb.sum()).sum())


def sub_ari(labels, clusters, a, b):
    """ARI between partition and labels using only cells labelled A or B."""
    from sklearn.metrics import adjusted_rand_score

    labels = np.asarray(labels)
    mask = (labels == a) | (labels == b)
    return float(adjusted_rand_score(labels[mask], np.asarray(clusters)[mask]))


def cluster_truth(lineage_table, lineage_order, share_needed):
    """Per cluster: plurality lineage, its share, and truth(c) (plurality if share >= threshold, else unresolved)."""
    rows = []
    for cluster in lineage_table.columns:
        counts = lineage_table[cluster].reindex(lineage_order).to_numpy(dtype=float)
        total = counts.sum()
        idx = int(np.argmax(counts))
        share = counts[idx] / total
        rows.append(
            {
                "cluster": int(cluster),
                "n_cells": int(total),
                "plurality_lineage": lineage_order[idx],
                "plurality_share": float(share),
                "plurality_tie": bool((counts == counts[idx]).sum() > 1),
                "truth": lineage_order[idx] if share >= share_needed else "unresolved",
            }
        )
    return rows


def ci_bounds(values, level):
    tail = (1 - level) / 2 * 100
    low, high = np.percentile(values, [tail, 100 - tail])
    return float(low), float(high)


# ---------------------------------------------------------------- answer key


def calling_script_number():
    main = sys.modules.get("__main__")
    name = Path(getattr(main, "__file__", "") or "").name
    match = re.match(r"^(\d{2})_", name)
    return int(match.group(1)) if match else None


def open_answer_key():
    """The only reader of the quarantined labels. Raises unless every freeze condition holds."""
    import pandas as pd

    cfg = load_config()
    qcfg = cfg["quarantine"]
    number = calling_script_number()
    if number is None or number < qcfg["first_script_allowed"]:
        raise PermissionError(f"answer key refused: calling script number {number} is below {qcfg['first_script_allowed']}")
    freeze = find_commit(FREEZE_SUBJECT)
    if freeze is None or not is_ancestor(freeze):
        raise PermissionError("answer key refused: no freeze commit on the history of HEAD")
    if not (ROOT / FREEZE_MANIFEST).is_file():
        raise PermissionError("answer key refused: freeze manifest missing")
    changed = verify_freeze_manifest()
    if changed:
        raise PermissionError("answer key refused: frozen files changed: " + ", ".join(changed))
    problems = validate_human_naming()
    if problems:
        raise PermissionError("answer key refused: second naming incomplete: " + "; ".join(problems))
    if not human_naming_committed_after(freeze):
        raise PermissionError("answer key refused: second naming not committed after the freeze")
    path = ROOT / qcfg["labels_file"]
    expected = (ROOT / "provenance" / "quarantine_sha256.txt").read_text(encoding="utf-8").split()[0]
    if sha256_file(path) != expected:
        raise PermissionError("answer key refused: quarantine file hash does not match")
    return pd.read_csv(path, dtype=str)
