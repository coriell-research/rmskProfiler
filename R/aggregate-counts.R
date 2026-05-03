#' Aggregate counts to gene and TE level
#'
#' This function will sum counts to the gene level for all transcripts, and to the
#' specified level for TE-loci. A new SummarizedExperiment object is returned with
#' aggregated assays. rowData is added indicating the number of features summed ('Freq')
#' and a 'feature_length' column. For genes, 'feature_length' is the sum of reduced
#' exon widths. For TEs, it is the sum of the length of all TE-loci that are members
#' of the aggregated group, or the locus length if te_level is "locus".
#'
#' @param x SummarizedExperiment object produced by \code{importQuants()}
#' @param resource_dir Path to the rmskProfiler resources directory
#' @param te_level One of "subfamily" (default), "locus", "family", or "class".
#' "locus" keeps TE-loci independent, while others aggregate based on classification.
#'
#' @return SummarizedExperiment
#' @export
#'
#' @examples
#' \dontrun{
#' se <- importQuants("quants", resource_dir = "hg38-resources")
#'
#' # Aggregate transcripts to genes, keep TE loci intact
#' aggregated_locus <- aggregateCounts(se, resource_dir = "hg38-resources", te_level = "locus")
#'
#' # Aggregate transcripts to genes, TEs to subfamily
#' aggregated_sub <- aggregateCounts(se, resource_dir = "hg38-resources", te_level = "subfamily")
#' }
aggregateCounts <- function(
  x,
  resource_dir,
  te_level = c("subfamily", "family", "class", "locus")
) {
  te_level <- match.arg(te_level)
  rd_orig <- SummarizedExperiment::rowData(x)

  te_grouping <- switch(
    te_level,
    locus = rownames(x),
    subfamily = stringi::stri_c(
      rd_orig$Subfamily,
      rd_orig$Family,
      rd_orig$Class,
      sep = ":"
    ),
    family = stringi::stri_c(rd_orig$Family, rd_orig$Class, sep = ":"),
    class = rd_orig$Class
  )

  feature_id <- data.table::fcoalesce(
    rd_orig$gene_id,
    te_grouping
  )

  # Compute feature lengths
  resources <- list.files(resource_dir, full.names = TRUE)
  dbfile <- grep("annotation.txdb", resources, value = TRUE)
  txdb <- AnnotationDbi::loadDb(dbfile)

  # Gene lengths are reduced exon widths
  exons_by_gene <- suppressWarnings(GenomicFeatures::exonsBy(txdb, by = "gene"))
  reduced_exon_lengths <- sum(BiocGenerics::width(BiocGenerics::reduce(
    exons_by_gene
  )))
  names(reduced_exon_lengths) <- names(exons_by_gene)

  # TE lengths are sum of widths of ranges
  te_widths <- BiocGenerics::width(rd_orig$Ranges)
  sum_te_widths <- tapply(te_widths, feature_id, sum, na.rm = TRUE)
  sum_te_widths <- sum_te_widths[
    !names(sum_te_widths) %in% names(reduced_exon_lengths)
  ]

  final_widths <- c(reduced_exon_lengths, sum_te_widths)

  gene_meta <- rd_orig[
    startsWith(rownames(x), "ENS"),
    c("gene_id", "gene_name", "gene_type")
  ]
  gene_meta <- unique(data.table::as.data.table(data.frame(gene_meta)))

  rd_new <- data.table::as.data.table(data.frame(table(
    feature_id = feature_id
  )))
  rd_new <- data.table::merge.data.table(
    rd_new,
    gene_meta,
    by.x = "feature_id",
    by.y = "gene_id",
    all.x = TRUE
  )

  data.table::setDF(rd_new)
  rownames(rd_new) <- rd_new$feature_id
  rd_new$feature_length <- final_widths[rownames(rd_new)]

  assay_names <- names(SummarizedExperiment::assays(x))
  new_assays <- lapply(assay_names, function(an) {
    rowsum(SummarizedExperiment::assay(x, an), group = feature_id)
  })
  names(new_assays) <- assay_names

  result <- SummarizedExperiment::SummarizedExperiment(
    assays = new_assays,
    rowData = rd_new[rownames(new_assays[[1]]), ],
    colData = SummarizedExperiment::colData(x),
    metadata = S4Vectors::metadata(x)
  )

  return(result)
}
