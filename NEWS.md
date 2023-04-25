---
editor_options: 
  markdown: 
    wrap: 72
---

# 230424: Gene accessing, LISI

-   `get_PC_genes`: accesses BioMart database for genes annotated as
    protein coding
-   `filter_PC_genes`: filter dataframe to specified list of genes, with
    possible `Seurat::UpdateSymbolList` powered gene symbol correcting.
    Technically this isn't specific to protein-coding filtering only and
    words with list input.
-   `get_gene_desc`: accesses BioMart database for gene description w/
    list of genes input
    -   Technically not specific to getting description fields given
        list of BioMart attrs
    -   Should probably add a error/warning if gene names don't match,
        possible Seurat correction?
-   `run_LISI`: for Seurat obj, assess integration quantitatively by
    using LISI metric
    -   outputs visualization to help interpret resulting values

------------------------------------------------------------------------

# Rubrary 0.7.0

230329 function improvements

# 230323: new GSEA, Seurat DimPlot functions & improvements

-   plot_GSEA_pathway
-   plot_dim_set
-   hypergeo_coexp

# 230303: New scRNA hypergeometric p-value and plot function

-   genecoexp_scatter_hyper
-   phyper_df
-   pltAB
-   run_RRHO (WIP)

# 230203: new functions

-   plot_waterfall
-   plot_density
-   check_normal
-   calc_sumZscore + various fixes

# New DESeq functions + extras

-   output_DESeq
-   filter_DESeq_PC
-   plot_DESeq_volcano
-   plot_scatter_compare
-   plot_screeplot
-   plot_scatter_mtx (WIP)

# Rubrary 0.2.0

-   plot_distribution: Density plot with histo

# Rubrary 0.1.5

-   Save option for "plot_scatter"

# Rubrary 0.1.4

-   Optional KS pval for plot_waterfall_hl

# Rubrary 0.1.3

-   plot_waterfall_hl label outside

# Rubrary 0.1.2

-   Fixed export of plot_waterfall_hl.R to namespace

# Rubrary 0.1.1

# Rubrary 0.1.0

-   Added a `NEWS.md` file to track changes to the package.
