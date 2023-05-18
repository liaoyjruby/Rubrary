General

- PCA improvements + associated vignettes
- scRNA integration assessment metrics
- GSEA barplot + associated helper functions

New
- `plot_PCA_biplot`
- `plot_GSEA_barplot`: function-ified version of horizontal barplot originally made for Jack
  - `format_GSEA_name`: cleans up MSigDB underscore + all caps names a bit
  - `split_line`: splits long text to multiple lines based on number of char or lines
- `rotate_varimax`: Varimax rotation for prcomp output
- `assess_integration`: `run_LISI` but more general, incorporating `Seurat::MixingMetric` + `CellMixS::cms`
  - `run_LISI`: wrapper for `assess_integration` specific to `LISI`
  - `run_MixingMetric`: wrapper for `assess_integration` specific to `MixingMetric`
  - `run_CellMixS`: wrapper for `assess_integration` specific to `cms`

Improved

- `run_PCA`: fix eigenvector / loadings distinction
-   `plot_volcano`: rewritten terms + args to be more general and align better w/ `EnhancedVolcano` documentation
-   `genes` functions: added `Mart` arg to functions that try to query BioMart
    -   `plot_PC_genes` -\> `plot_genes`, still defaults to filtering by PC genes but technically any gene list can be passed in to filter dataframe by
-   `plot_waterfall`: a little smarter on where to place high/low value labels but still quite dumb
- `plot_scatter`: color by group, more text options

Vignettes

- `PCA_Walkthrough`: based on Lindsey Smith's PCA tutorial, a more step by step toy dataset for PCA including some of the math behind it
  - TODO: Varimax section... once I get visualization working for plotting varimax
- `PCA_Quickstart`: very basic Palmer Penguin demonstration of plotting scripts
-   `DESeq_GSEA`: fgsea section + loads of properly cited citations
-   possibly split vignette into DE genes and GSEA separately? can do comparative DE methods if so and RRHO it out

# Rubrary 0.8.0

## General

-   Gene accessing, LISI
-   PCA functions improvements
-   Functional vignettes

## New

-   `get_PC_genes`: accesses BioMart database for genes annotated as protein coding
-   `filter_PC_genes`: filter dataframe to specified list of genes, with possible `Seurat::UpdateSymbolList` powered gene symbol correcting. Technically this isn't specific to protein-coding filtering only and words with list input.
-   `get_gene_desc`: accesses BioMart database for gene description w/ list of genes input
    -   Technically not specific to getting description fields given list of BioMart attrs
    -   Should probably add a error/warning if gene names don't match, possible Seurat correction?
-   `convert_genes`: access BioMart database to convert gene name format
-   `run_LISI`: for Seurat obj, assess integration quantitatively by using LISI metric
    -   outputs visualization to help interpret resulting values

## Improved

-   `run_PCA`: additional documentation recommending standardization
-   `plot_PCA`: accommodation for passing in `prcomp` object, custom colors

## Vignettes

-   `PCA`: changed to `palmerpenguins` dataset and have it actually run and output
-   `DESeq_GSEA`: working DESeq example with `airway` dataset

# Rubrary 0.7.0

230329 function improvements

# Rubrary 0.6.0

230323: new GSEA and Seurat DimPlot functions & improvements

-   plot_GSEA_pathway
-   plot_dim_set
-   hypergeo_coexp

# Rubrary 0.5.0

230303: New scRNA hypergeometric p-value and plot function

-   genecoexp_scatter_hyper
-   phyper_df
-   pltAB
-   run_RRHO (WIP)

# Rubrary 0.4.0

230203: new functions

-   plot_waterfall
-   plot_density
-   check_normal
-   calc_sumZscore
-   various fixes

# Rubrary 0.3.0

New DESeq functions + extras

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
