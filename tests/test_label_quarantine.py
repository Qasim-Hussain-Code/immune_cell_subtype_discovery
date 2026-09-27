"""Quarantine checks: blind scripts never reach the labels, the pooled file is anonymous,
and no evaluation metrics were committed before the freeze."""
import ast
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import common  # noqa: E402

BLIND_SCRIPTS = ["03_preprocess.py", "04_select_k.py", "05_cluster.py", "06_blind_naming.py"]


def forbidden_strings():
    cfg = common.load_config()
    return ["labels_quarantine", "labels_file", "data/raw", "open_answer_key", *cfg["data"]["populations"]]


def imported_common_names(tree):
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == "common":
            names.update(alias.name for alias in node.names)
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == "common":
            names.add(node.attr)
    return names


def common_sources(names):
    """Source of each imported object in common.py, following references to other objects in common."""
    source = (ROOT / "scripts" / "common.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    defs = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            defs[node.name] = node
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    defs[target.id] = node
    seen, todo = set(), list(names)
    while todo:
        name = todo.pop()
        if name in seen or name not in defs:
            continue
        seen.add(name)
        for sub in ast.walk(defs[name]):
            if isinstance(sub, ast.Name) and sub.id in defs:
                todo.append(sub.id)
    return {name: ast.get_source_segment(source, defs[name]) for name in seen}


@pytest.mark.parametrize("script", BLIND_SCRIPTS)
def test_blind_scripts_do_not_reach_labels(script):
    path = ROOT / "scripts" / script
    if not path.is_file():
        pytest.skip(f"{script} not written yet")
    text = path.read_text(encoding="utf-8")
    bad = forbidden_strings()
    assert [s for s in bad if s in text] == [], f"{script} mentions a forbidden string"
    for name, source in common_sources(imported_common_names(ast.parse(text))).items():
        assert [s for s in bad if s in source] == [], f"{script} imports common.{name}, which mentions a forbidden string"


def test_pooled_file_is_anonymous():
    path = ROOT / common.load_config()["quarantine"]["pooled_file"]
    if not path.is_file():
        pytest.skip("pooled file not built yet")
    import anndata as ad

    pooled = ad.read_h5ad(path, backed="r")
    try:
        assert list(pooled.obs.columns) == []
        assert all(re.fullmatch(r"c\d{6}", i) for i in pooled.obs_names)
        assert len(pooled.uns) == 0
    finally:
        pooled.file.close()


def test_no_evaluation_metrics_before_freeze():
    if common.git("rev-parse", "HEAD").returncode != 0:
        pytest.skip("no commits yet")
    freeze = common.find_commit(common.FREEZE_SUBJECT)
    revision = freeze if freeze else "HEAD"
    log = subprocess.run(["git", "log", "--diff-filter=A", "--name-only", "--format=%H", revision], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    added = [line.strip() for line in log.splitlines() if line.startswith("results/metrics/")]
    late = [f for f in added if re.match(r"results/metrics/(\d+)_", f) and int(re.match(r"results/metrics/(\d+)_", f).group(1)) >= 7]
    assert late == [], f"evaluation metrics committed before the freeze: {late}"
