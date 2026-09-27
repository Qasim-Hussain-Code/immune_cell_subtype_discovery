"""Check that every number in README.md can be traced to a metrics file.

Fenced blocks, code spans, URLs and link targets are ignored. Every remaining
number must equal a value in results/metrics/*.json, written either to four
decimal places or in full. The check also fails on em dashes, en dashes,
emojis and American spellings, which keeps the prose in the house style.
"""
import json
import re
import sys

from common import ROOT, check_posts, write_metrics

AMERICAN = ["analyze", "normalize", "normalized", "color", "behavior", "center", "favor", "labeled", "modeling", "summarize", "visualize", "optimize"]
EMOJI = re.compile("[\U0001F000-\U0001FAFF\u2600-\u27bf\ufe0f\u200d\u2b50\u2b55\u2190-\u21ff\u2300-\u23ff\u2b00-\u2bff]")
NUMBER = re.compile(r"(?<![A-Za-z0-9])-?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?![A-Za-z0-9])")


def strip_markdown(text):
    text = re.sub(r"^\s*(```|~~~).*?^\s*\1[^\n]*$", " ", text, flags=re.S | re.M)
    text = re.sub(r"<!--.*?-->", " ", text, flags=re.S)
    text = re.sub(r"`[^`\n]*`", " ", text)
    text = re.sub(r"\]\([^)]*\)", "] ", text)
    text = re.sub(r"https?://\S+", " ", text)
    text = re.sub(r"^(\s*)(?:[-*+]|\d+[.)])\s+", r"\1", text, flags=re.M)
    return text


def leaves(obj, out):
    if isinstance(obj, dict):
        for value in obj.values():
            leaves(value, out)
    elif isinstance(obj, list):
        for value in obj:
            leaves(value, out)
    elif isinstance(obj, bool) or obj is None:
        return
    elif isinstance(obj, int):
        out.add(str(obj))
    elif isinstance(obj, float):
        out.add(f"{obj:.4f}")
        out.add(repr(obj))


def main():
    check_posts()
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    allowed = set()
    for path in sorted((ROOT / "results" / "metrics").glob("*.json")):
        if path.name == "10_readme_audit.json":
            continue
        with open(path, encoding="utf-8") as fh:
            leaves(json.load(fh), allowed)

    stripped = strip_markdown(readme)
    found = [m.group(0) for m in NUMBER.finditer(stripped)]
    unmatched = sorted({n for n in found if n.replace(",", "") not in allowed})
    dashes = [c for c in readme if c in "\u2014\u2013"]
    emojis = EMOJI.findall(readme)
    spellings = sorted({w.lower() for w in AMERICAN for _ in re.finditer(rf"\b{w}\b", stripped, flags=re.I)})
    passed = not unmatched and not dashes and not emojis and not spellings
    write_metrics(
        "10_readme_audit",
        {
            "n_numbers_checked": len(found),
            "n_unmatched": len(unmatched),
            "unmatched": unmatched,
            "n_em_or_en_dashes": len(dashes),
            "n_emojis": len(emojis),
            "american_spellings": spellings,
            "passed": passed,
        },
    )
    print(f"checked {len(found)} numbers; unmatched {len(unmatched)}; dashes {len(dashes)}; emojis {len(emojis)}; spellings {spellings}")
    if unmatched:
        print("unmatched: " + ", ".join(unmatched))
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
