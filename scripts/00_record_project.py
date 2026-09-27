"""Record the fixed facts of the project and echo the configuration.

Some numbers in the README come from the published posts or from the 10x
Genomics dataset pages rather than from the analysis: the purity of each
population, the cells and reads reported on its page, and the size of the
blood sample in the 2017 paper. They are written to a metrics file here so
that the README audit can trace every number it finds to a file.
"""
from common import TENX_NAMES, load_config, write_metrics

# Facts read from the 10x Genomics dataset pages before any code existed.
# "approx" marks values shown on the page with a tilde.
PAGE_FACTS = {
    "b_cells": {"purity_percent": 100, "purity_approx": True, "cells_on_page": 10085, "cells_approx": False, "reads_per_cell": 25000},
    "cd14_monocytes": {"purity_percent": 98, "purity_approx": False, "cells_on_page": 2612, "cells_approx": False, "reads_per_cell": 100000},
    "cd34": {"purity_percent": 45, "purity_approx": False, "cells_on_page": 9000, "cells_approx": True, "reads_per_cell": 24700},
    "cd4_t_helper": {"purity_percent": 99, "purity_approx": False, "cells_on_page": 11000, "cells_approx": True, "reads_per_cell": 21000},
    "regulatory_t": {"purity_percent": 95, "purity_approx": False, "cells_on_page": 10000, "cells_approx": True, "reads_per_cell": 27000},
    "naive_t": {"purity_percent": 98, "purity_approx": False, "cells_on_page": 10479, "cells_approx": False, "reads_per_cell": 19000},
    "memory_t": {"purity_percent": 98, "purity_approx": False, "cells_on_page": 10224, "cells_approx": False, "reads_per_cell": 24000},
    "cd56_nk": {"purity_percent": 92, "purity_approx": False, "cells_on_page": 8000, "cells_approx": True, "reads_per_cell": 29000},
    "cytotoxic_t": {"purity_percent": 98, "purity_approx": False, "cells_on_page": 10000, "cells_approx": True, "reads_per_cell": 28600},
    "naive_cytotoxic": {"purity_percent": 99, "purity_approx": False, "cells_on_page": 11953, "cells_approx": False, "reads_per_cell": 20000},
}


def main():
    cfg = load_config()
    table = {}
    for key, facts in PAGE_FACTS.items():
        table[key] = dict(facts, tenx_name=TENX_NAMES[key], lineage=cfg["data"]["populations"][key], reads_per_cell_approx=True)
    payload = {
        "chapter": cfg["project"]["chapter"],
        "days_referenced": [2, 48, 61, 62, 63, 64, 65, 66],
        "citation": {"year": 2017, "volume": 8, "article_number": 14049, "doi": cfg["data"]["doi"]},
        "paper_pbmc_sample": {"n_cells": 68000, "k_used": 10},
        "dataset_pages": table,
        "ci_level_percent": 95,
        "config": cfg,
    }
    path = write_metrics("00_project", payload)
    print(f"wrote {path.relative_to(path.parents[2]).as_posix()}")


if __name__ == "__main__":
    main()
