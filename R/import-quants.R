#' Import rmsk quants as a SummarizedExperiment object
#'
#' This function imports counts from the Salmon quant files adjusting for
#' overdispersion with \code{edgeR::catchSalmon()}. The resulting
#' \code{SummarizedExperiment} object contains rowData with transcript and TE
#' loci annotations and an additional column containing a \code{GRangesList} of
#' the locations for each of the TE loci and transcipts.
#'
#' @details
#' The SummarizedEperiment object contains two assays, 'counts' and 'orig.' The
#' 'counts' assay contains the counts \emph{after} adjusting for over dispersion.
#' This is the assay to be used for downstream differential expression analyses.
#' The orig' assay contains the counts before scaling. The metadata slot contains
#' the annotation data.frame returned from catchSalmon. The rowData slot of the
#' SummarizedExperiment contains columns with boolean values for each TE loci
#' indicating whether or not that loci has an overlap with a given feature. For
#' example, hasExonic==TRUE would indicate that loci has an overlap with some
#' exon on the same strand. hasUnstrandedExonic==TRUE would indicate that the
#' loci has an overlap with an exon on either strand. The rowData also contains
#' an additional column 'Ranges' that contains a GRangesList for every
#' transcript and TE hash location.
#'
#' @param quant_dir Path to the directories created by salmonQuant
#' @param resource_dir Path to the directory containing index generation resources
#' @param remove_zeros Should all zero rows be removed before returning object. Default TRUE
#'
#' @return SummarizedExperiment
#' @export
#'
importQuants <- function(quant_dir, resource_dir, remove_zeros = TRUE) {

  message("Importing quants with edgeR::catchSalmon...")
  paths <- list.dirs(path = quant_dir, full.names = TRUE, recursive = FALSE)
  catch <- edgeR::catchSalmon(paths, verbose = FALSE)
  counts <- catch$counts / catch$annotation$Overdispersion

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      counts = counts,
      orig = catch$counts
      ),
    metadata = list(annotation = catch$annotation)
  )
  colnames(se) <- basename(colnames(se))

  message("Reading in range information for transcripts and TE loci...")
  resources <- list.files(resource_dir, full.names = TRUE)
  gtf_file <- grep("annotation.gtf.gz", resources, value = TRUE)
  gtf <- rtracklayer::import(gtf_file)
  tx <- gtf[gtf$type == "transcript", ]
  tx_dt <- data.table::as.data.table(data.frame(tx))[, .(transcript_id, gene_id, gene_name, gene_type)]
  names(tx) <- tx$transcript_id
  tx <- as(tx, "GRangesList")

  # Get GRangesList for all loci
  rmsk_grlfile <- grep("rmsk-grl.rds", resources, value = TRUE)
  rmsk <- readRDS(rmsk_grlfile)

  # Align with features in SummarizedExperiment object and combine into one list
  grl <- c(tx, rmsk)
  grl <- grl[SummarizedExperiment::rownames(se)]

  message("Reading in TE loci and transcript annotations...")
  rdfile <- grep("rmsk-rowData.tsv.gz", resources, value = TRUE)
  rmsk_dt <- data.table::fread(rdfile)
  rd <- data.table::rbindlist(list(tx_dt, rmsk_dt), fill = TRUE)
  rd <- S4Vectors::DataFrame(rd)
  rownames(rd) <- c(tx_dt$transcript_id, rmsk_dt$Hash)
  rd <- rd[SummarizedExperiment::rownames(se), ]

  # Add Ranges as a new column in the rowData
  rd$Ranges <- grl
  SummarizedExperiment::rowData(se) <- rd

  if (isTRUE(remove_zeros)) {
    se <- se[rowSums(SummarizedExperiment::assay(se, "counts")) > 0, ]
  }

  message("Done.")
  return(se)
}
