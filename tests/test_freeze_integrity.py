"""Check that the frozen files still match their recorded hashes and that the
freeze commit lies on the history of HEAD.
"""
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import common  # noqa: E402


def manifest():
    path = ROOT / common.FREEZE_MANIFEST
    if not path.is_file():
        pytest.skip("no freeze manifest yet")
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def test_frozen_hashes_match_manifest():
    manifest()
    assert common.verify_freeze_manifest() == []


def test_freeze_commit_is_ancestor_of_head():
    data = manifest()
    freeze = common.find_commit(common.FREEZE_SUBJECT)
    assert freeze is not None, "manifest exists but no freeze commit was found"
    assert common.is_ancestor(freeze)
    parent = common.git("rev-parse", f"{freeze}^").stdout.strip()
    assert parent == data["parent_commit"]
