#' Aggregate counts to gene/subfamily level
#'
#' This function will sum counts to the gene and the subfamily (default)
#' level for all transcripts and TE-loci. A new SummarizedExperiment object is
#' returned with three assays, 'counts', 'orig', and 'tpms' representing the aggregated
#' counts from the \code{catchSalmon()} down-scaled ('counts'), the original Salmon
#' counts ('orig'), and the aggregated TPM values ('tpms'), respectively. rowData is added to
#' the object indicating the number of TE-loci or transcripts that were summed for each resulting
#' feature ('Freq'). rowData of the aggregated object also contains a column called 'feature_length'.
#' For genes, 'feature_length' is the sum of the reduced exon widths. For TEs, 'feature_length'
#' is the sum of the length of all TE-loci that are members of the subfamily, family, or class.
#'
#' @param x SummarizedExperiment object produced by \code{importQuants()}
#' @param resource_dir Path to the rmskProfiler resources directory
#' @param level One of "subfamily" (default), "family", or "class" indicating
#' the level of classification to sum TE-loci to.
#'
#' @return SummarizedExperiment
#' @export
#'
#' @examples
#' \dontrun{
#'
#' se <- importQuants("quants", resource_dir = "hg38-resources")
#' aggregated <- aggregateCounts(se, resource_dir = "hg38-resources", level = "subfamily")
#'
#' }
aggregateCounts <- function(x, resource_dir, level = "subfamily") {
  # Create grouping variable for summarizing
  agg_level <- match.arg(level, choices = c("subfamily", "family", "class"))
  SummarizedExperiment::rowData(x)$repElem <-
    switch(
      agg_level,
      subfamily = stringr::str_c(
        SummarizedExperiment::rowData(x)$Subfamily,
        SummarizedExperiment::rowData(x)$Family,
        SummarizedExperiment::rowData(x)$Class,
        sep = ":"
      ),
      family = stringr::str_c(
        SummarizedExperiment::rowData(x)$Family,
        SummarizedExperiment::rowData(x)$Class,
        sep = ":"
      ),
      class = SummarizedExperiment::rowData(x)$Class
    )

  feature_id <- data.table::fcoalesce(
    SummarizedExperiment::rowData(x)$gene_id,
    SummarizedExperiment::rowData(x)$repElem
  )

  # Compute feature lengths
  # For genes, get the reduced exon lengths
  resources <- list.files(resource_dir, full.names = TRUE)
  dbfile <- grep("annotation.txdb", resources, value = TRUE)
  txdb <- AnnotationDbi::loadDb(dbfile)
  exons_by_gene <- suppressWarnings(GenomicFeatures::exonsBy(txdb, by = "gene"))
  reduced_exon_lengths <- sum(width(reduce(exons_by_gene)))
  names(reduced_exon_lengths) <- names(exons_by_gene)

  # For TEs, compute the sum of the widths of all loci
  feature_widths <- sum(width(SummarizedExperiment::rowData(x)$Ranges))
  sum_feature_widths <- tapply(feature_widths, feature_id, sum, na.rm = TRUE)
  sum_feature_widths <- sum_feature_widths[
    !startsWith(names(sum_feature_widths), "ENS")
  ]
  final_widths <- c(reduced_exon_lengths, sum_feature_widths)

  # Extract useful gene information
  df <- SummarizedExperiment::rowData(x)[
    startsWith(rownames(x), "ENS"),
    c("gene_id", "gene_name", "gene_type")
  ]
  df <- data.table::as.data.table(data.frame(df))
  df <- unique(df)

  # Combine the gene information with the length data
  rd <- data.table::as.data.table(data.frame(table(feature_id)))
  rd <- data.table::merge.data.table(
    rd,
    df,
    by.x = "feature_id",
    by.y = "gene_id",
    all.x = TRUE
  )
  data.table::setDF(rd, rownames = rd$feature_id)
  rd$feature_length <- final_widths[rownames(rd)]

  # Sum the assay data
  counts <- rowsum(SummarizedExperiment::assay(x, "counts"), group = feature_id)
  orig <- rowsum(SummarizedExperiment::assay(x, "orig"), group = feature_id)
  tpms <- rowsum(SummarizedExperiment::assay(x, "tpms"), group = feature_id)

  result <- SummarizedExperiment::SummarizedExperiment(
    assays = list("counts" = counts, "orig" = orig, "tpms" = tpms),
    rowData = rd[rownames(counts), ],
    colData = SummarizedExperiment::colData(x),
    metadata = S4Vectors::metadata(x)
  )

  return(result)
}
