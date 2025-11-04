#' Aggregate counts to gene/subfamily level
#'
#' This function will sum counts to the gene and the subfamily (default)
#' level for all transcripts and TE-loci. A new SummarizedExperiment object is
#' returned with two assays, 'counts' and 'orig' representing the aggregated
#' counts from the catchSalmon down-scaled and original Salmon counts matrices,
#' respectively. rowData is added to the object indicating the number of TE-loci
#' or transcripts that were summed for each resulting feature.
#'
#' @param x SummarizedExperiment object produced by importQuants()
#' @param level One of "subfamily" (default), "family", or "class" indicating
#' the level of classification to sum TE-loci to.
#'
#' @return SummarizedExperiment
#' @export
#'
#' @examples
#' \dontrun{
#'
#' se <- importQuants("quants", resources_dir = "hg38-resources")
#' aggregated <- aggregateCounts(se, level = "subfamily")
#'
#' }
aggregateCounts <- function(x, level = "subfamily") {
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

  rd <- data.table::as.data.table(data.frame(table(feature_id)))

  # Add on gene information from original se for easier downstream analysis
  df <- SummarizedExperiment::rowData(x)[
    startsWith(rownames(x), "ENS"),
    c("gene_id", "gene_name", "gene_type")
  ]
  df <- data.table::as.data.table(data.frame(df))
  df <- unique(df)
  rd <- data.table::merge.data.table(
    rd,
    df,
    by.x = "feature_id",
    by.y = "gene_id",
    all.x = TRUE
  )
  data.table::setDF(rd, rownames = rd$feature_id)

  counts <- rowsum(SummarizedExperiment::assay(x, "counts"), group = feature_id)
  orig <- rowsum(SummarizedExperiment::assay(x, "orig"), group = feature_id)

  result <- SummarizedExperiment::SummarizedExperiment(
    assays = list("counts" = counts, "orig" = orig),
    rowData = rd[rownames(counts), ],
    colData = SummarizedExperiment::colData(x),
    metadata = S4Vectors::metadata(x)
  )

  return(result)
}
