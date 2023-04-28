## code to prepare `airway_deseq_res` dataset goes here
library(airway)
library(DESeq2)
library(biomaRt)

# Get `airway` data
data(airway)
cts <- assay(airway, "counts")
anno <- colData(airway) %>%
  as.data.frame()
anno$dex <- relevel(anno$dex, "untrt")
# DESeq
dds <- DESeqDataSetFromMatrix(
  countData = cts[,rownames(anno)],
  colData = anno,
  design = ~ dex)
dds <- DESeq(dds)
res <- results(dds, contrast = c("dex", "trt", "untrt"))
deseq_res_df <- Rubrary::output_DESeq(res)
# BioMart
mart <- biomaRt::useDataset(
  dataset = "hsapiens_gene_ensembl",
  mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL"))
conv <- Rubrary::convert_genes(
  genes = deseq_res_df$gene,
  from_to = c("ensembl_gene_id", "hgnc_symbol"),
  mart = mart, table = T) %>%
  filter(hgnc_symbol != "") %>% # Exclude unmatched HGNC symbols
  distinct(hgnc_symbol, .keep_all = T)
deseq_res_hgnc <- dplyr::left_join(
  conv, deseq_res_df, by = join_by(ensembl_gene_id == gene)) %>%
  arrange(desc(sign_log_p))
# PC genes
airway_deseq_res <- Rubrary::filter_genes(
  df = deseq_res_hgnc,
  gene_col = "hgnc_symbol",
  genes_filt = Rubrary::get_PC_genes(mart)
)

usethis::use_data(airway_deseq_res, overwrite = TRUE)
