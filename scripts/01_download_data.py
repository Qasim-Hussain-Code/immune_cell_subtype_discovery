"""Download and unpack the ten filtered gene-barcode matrices.

Each archive is fetched from the 10x Genomics server, hashed and logged in
provenance/downloads.tsv. The script stops at the first failed request,
because a changed link should be checked against the dataset page rather than
guessed. No per-population cell counts are written: which population each
cell came from stays hidden until every grouping is frozen.
"""
import datetime
import tarfile

import pandas as pd
import requests

from common import ROOT, check_free_space, check_posts, load_config, sha256_file, write_metrics


def locate(folder, name):
    hits = sorted(folder.rglob(name))
    if len(hits) != 1:
        raise SystemExit(f"Stop: expected one {name} under {folder}, found {len(hits)}")
    return hits[0]


def main():
    check_posts()
    cfg = load_config()
    rows = []
    gene_hashes = set()
    n_genes = None
    for key in cfg["data"]["populations"]:
        check_free_space()
        url = cfg["data"]["url_pattern"].format(key=key)
        folder = ROOT / "data" / "raw" / key
        folder.mkdir(parents=True, exist_ok=True)
        archive = folder / url.rsplit("/", 1)[-1]
        stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        response = requests.get(url, stream=True, timeout=120)
        if response.status_code != 200:
            raise SystemExit(f"Stop: download of {key} returned HTTP {response.status_code} for {url}. Check the link in the Output and supplemental files section of the dataset page.")
        with open(archive, "wb") as fh:
            for block in response.iter_content(chunk_size=1 << 20):
                fh.write(block)
        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(folder, filter="data")
        matrix = locate(folder, "matrix.mtx")
        genes = locate(folder, "genes.tsv")
        locate(folder, "barcodes.tsv")
        if genes.parent != matrix.parent:
            raise SystemExit(f"Stop: matrix.mtx and genes.tsv are in different folders for {key}")
        gene_hashes.add(sha256_file(genes))
        if n_genes is None:
            with open(genes, encoding="utf-8") as fh:
                n_genes = sum(1 for line in fh if line.strip())
        rows.append(
            {
                "key": key,
                "url": url,
                "status": response.status_code,
                "bytes": archive.stat().st_size,
                "sha256": sha256_file(archive),
                "utc_time": stamp,
                "path": archive.relative_to(ROOT).as_posix(),
            }
        )
        print(f"downloaded and unpacked {key}")
    pd.DataFrame(rows).to_csv(ROOT / "provenance" / "downloads.tsv", sep="\t", index=False, lineterminator="\n")
    genes_identical = len(gene_hashes) == 1
    if not genes_identical:
        raise SystemExit("Stop: genes.tsv differs between the ten matrices")
    write_metrics(
        "01_download",
        {
            "n_populations": len(rows),
            "total_bytes": int(sum(r["bytes"] for r in rows)),
            "genes_identical": genes_identical,
            "n_genes": n_genes,
        },
    )


if __name__ == "__main__":
    main()
