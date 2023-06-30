
- removed a LOT of old functions that weren't generalized very well so were probably not used anyways
- added more examples to various functions (DE genes / GSEA)
- reorganization of function reference/index

# Rubrary 0.10.3

- theme adjustments

# Rubrary 0.10.2

- Removed old DESeq functions, reorganized, and fixed various bugs

# Rubrary 0.10.1

- `GSEAsq`uared functions bug fixes

# Rubrary 0.10.0

- `GSEAsq`uared functions
  - `run_GSEA_squared`: main function, replicating [glab.library::gsea_squared](https://github.com/graeberlab-ucla/glab.library/blob/master/R/gsea_squared.R) with `keyword_plot_method == 1` setting
  - `get_GSEAsq_terms`: terms/keywords
  - and `plot_GSEAsq_density` moved to be with this set of functions
- `get_kspval` signed enrichment score implementation and visualization based on [franapoli/signed-ks-test](https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R)
- Various bug fixes
  

# Rubrary 0.9.6

## New

- `rwrite` & `rread`: fast `data.table` based read-in and write-out functions

## Improved

- `get_log_coding_geneexp` improvements - can specify log base now

# Rubrary 0.9.5

## New

- complete `plot_DESeq_boxplot` rewrite + rename to `plot_DEgene_boxplot`
  - finally !! able to remove `gridExtra` dependency

# Rubrary 0.9.4

## Improved

- rename `project_PCA` to the more descriptive `predict_PCA`

# Rubrary 0.9.3

## New

- `rread`: wrapper for `data.table::fread()` that can pull in rownames
- `project_PCA`: long awaited *massive* PCA projection script incorporating varimax

## Improved

- Switch to `Seurat::DiscretePalette` for `pals`'s `alphabet` color palette

# Rubrary 0.9.2

## Improved

- `plot_GSEA_barplot`: `sig_cutoff` format by p value significance parameter

# Rubrary 0.9.1

## Improved

- GSEA functions: improved name formatting for enrichment plot
- `plot_PCA`: fixed subtitle for non-annotated plots (thanks Duncan :) )
- `plot_GSEA_barplot`: formatting flexibility improvements

# Rubrary 0.9.0

## General

- PCA improvements + associated vignettes
- scRNA integration assessment metrics
- GSEA barplot, GSEAsq density plot comparison, and associated helper functions

## New

- `plot_GSEAsq_density`: plot categorical percentile rank comparison between two GSEA squared signatures
- `plot_PCA_biplot`: conventional biplot with loadings and standardized scores
- `plot_GSEA_barplot`: function-ified version of horizontal barplot originally made for Jack
  - `format_GSEA_name`: cleans up MSigDB underscore + all caps names a bit
  - `split_line`: splits long text to multiple lines based on number of char or lines
- `rotate_varimax`: Varimax rotation for prcomp output
- `assess_integration`: `run_LISI` but more general, incorporating `Seurat::MixingMetric` + `CellMixS::cms`
  - `run_LISI`: wrapper for `assess_integration` specific to `LISI`
  - `run_MixingMetric`: wrapper for `assess_integration` specific to `MixingMetric`
  - `run_CellMixS`: wrapper for `assess_integration` specific to `cms`

## Improved

- `run_PCA`: fix eigenvector / loadings distinction
-   `plot_volcano`: rewritten terms + args to be more general and align better w/ `EnhancedVolcano` documentation
-   `genes` functions: added `Mart` arg to functions that try to query BioMart
    -   `plot_PC_genes` -\> `plot_genes`, still defaults to filtering by PC genes but technically any gene list can be passed in to filter dataframe by
-   `plot_waterfall`: a little smarter on where to place high/low value labels but still quite dumb
- `plot_scatter`: color by group, more text options

## Vignettes

- `PCA_Walkthrough`: based on Lindsey Smith's PCA tutorial, a more step by step toy dataset for PCA including some of the math behind it
- `PCA_Quickstart`: very basic Palmer Penguin demonstration of plotting scripts + replication of Graeber Beltran PCA
-   `DE_Genes` & `GSEA`: previous combined tutorial nows plit into two

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
