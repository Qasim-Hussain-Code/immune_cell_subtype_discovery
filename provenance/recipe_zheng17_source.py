# inspect.getsource(scanpy.pp.recipe_zheng17), scanpy 1.11.5
@old_positionals("n_top_genes", "log", "plot", "copy")
def recipe_zheng17(
    adata: AnnData,
    *,
    n_top_genes: int = 1000,
    log: bool = True,
    plot: bool = False,
    copy: bool = False,
) -> AnnData | None:
    """Normalize and filter as of :cite:t:`Zheng2017`.

    Reproduces the preprocessing of :cite:t:`Zheng2017` – the Cell Ranger R Kit of 10x
    Genomics.

    Expects non-logarithmized data.
    If using logarithmized data, pass `log=False`.

    The recipe runs the following steps

    .. code:: python

        sc.pp.filter_genes(adata, min_counts=1)         # only consider genes with more than 1 count
        sc.pp.normalize_per_cell(                       # normalize with total UMI count per cell
             adata, key_n_counts='n_counts_all'
        )
        filter_result = sc.pp.filter_genes_dispersion(  # select highly-variable genes
            adata.X, flavor='cell_ranger', n_top_genes=n_top_genes, log=False
        )
        adata = adata[:, filter_result.gene_subset]     # subset the genes
        sc.pp.normalize_per_cell(adata)                 # renormalize after filtering
        if log: sc.pp.log1p(adata)                      # log transform: adata.X = log(adata.X + 1)
        sc.pp.scale(adata)                              # scale to unit variance and shift to zero mean


    Parameters
    ----------
    adata
        Annotated data matrix.
    n_top_genes
        Number of genes to keep.
    log
        Take logarithm.
    plot
        Show a plot of the gene dispersion vs. mean relation.
    copy
        Return a copy of `adata` instead of updating it.

    Returns
    -------
    Returns or updates `adata` depending on `copy`.

    """
    start = logg.info("running recipe zheng17")
    if copy:
        adata = adata.copy()
    # only consider genes with more than 1 count
    pp.filter_genes(adata, min_counts=1)
    # normalize with total UMI count per cell
    normalize_total(adata, key_added="n_counts_all")
    filter_result = filter_genes_dispersion(
        adata.X, flavor="cell_ranger", n_top_genes=n_top_genes, log=False
    )
    if plot:  # should not import at the top of the file
        from ..plotting import _preprocessing as ppp

        ppp.filter_genes_dispersion(filter_result, log=True)
    # actually filter the genes, the following is the inplace version of
    #     adata = adata[:, filter_result.gene_subset]
    adata._inplace_subset_var(filter_result.gene_subset)  # filter genes
    normalize_total(adata)  # renormalize after filtering
    if log:
        pp.log1p(adata)  # log transform: X = log(X + 1)
    pp.scale(adata)
    logg.info("    finished", time=start)
    return adata if copy else None
