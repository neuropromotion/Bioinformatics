library(biomaRt)

ensg_to_hgnc <- function(ensg_list) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  res <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
        filters = "ensembl_gene_id", 
        values = unique(ensg_list[grepl("^ENSG", ensg_list)]), 
        mart = mart)
  res$hgnc_symbol
