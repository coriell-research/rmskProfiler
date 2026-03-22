#' Import rmsk quants as a SummarizedExperiment object
#'
#' This function imports counts from the Salmon quant files adjusting for
#' overdispersion with \code{edgeR::catchSalmon()}. The resulting
#' \code{SummarizedExperiment} object contains rowData with transcript and TE
#' loci annotations and an additional column containing a \code{GRangesList} of
#' the locations for each of the TE loci and transcipts.
#'
#' @details
#' The SummarizedExperiment object contains three assays, 'counts', 'orig.' and 'tpms'.
#' The 'counts' assay contains the counts \emph{after} adjusting for over dispersion.
#' This is the assay to be used for downstream differential expression analyses.
#' The 'orig' assay contains the counts before down-scaling. The 'tpms' assay
#' contains TPM values from the Salmon quant files. The metadata slot
#' contains the annotation data.frame returned from \code{catchSalmon()}. The rowData
#' slot of the SummarizedExperiment contains columns with boolean values for
#' each TE loci indicating whether or not that loci has an overlap with a given
#' feature. For example, hasExonic==TRUE would indicate that loci has an overlap
#' with some exon on the same strand. hasUnstrandedExonic==TRUE would indicate
#' that the loci has an overlap with an exon on either strand. The rowData also
#' contains an additional column 'Ranges' that contains a GRangesList for every
#' transcript and TE hash location.
#'
#' The rmskProfiler index contains transcripts from the "patch_hapl_scaff"
#' annotations. By default, \code{importQuants()} will exclude these transcripts
#' from the resulting object. The exclusion criteria a transcript where any of
#' its ranges does not lie within a standard chromosome will be excluded.
#'
#' @param quant_dir Path to the directories created by \code{salmonQuant()} or Salmon
#' @param resource_dir Path to the directory containing index generation resources
#' @param remove_zeros Should all zero rows be removed before returning object. Default TRUE
#' @param std_chromosomes Should only transcripts from standard chromosomes be
#' kept in the resulting object. Default TRUE.
#'
#' @return SummarizedExperiment
#' @export
#'
importQuants <- function(
  quant_dir,
  resource_dir,
  remove_zeros = TRUE,
  std_chromosomes = TRUE
) {
  message("Importing quants with edgeR::catchSalmon...")
  paths <- list.dirs(path = quant_dir, full.names = TRUE, recursive = FALSE)
  catch <- .catchSalmon2(paths, verbose = FALSE)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      counts = catch$counts / catch$annotation$Overdispersion,
      orig = catch$counts,
      tpms = catch$tpm
    ),
    metadata = list(annotation = catch$annotation)
  )
  colnames(se) <- basename(colnames(se))

  message("Reading in annotation information for transcripts and TE loci...")
  resources <- list.files(resource_dir, full.names = TRUE)
  rdfile <- grep("rmsk-rowData.rds", resources, value = TRUE)
  rd <- readRDS(rdfile)
  SummarizedExperiment::rowData(se) <- rd[SummarizedExperiment::rownames(se), ]

  if (isTRUE(remove_zeros)) {
    message("Removing any all zero rows...")
    se <- se[rowSums(SummarizedExperiment::assay(se, "counts")) > 0, ]
  }

  if (isTRUE(std_chromosomes)) {
    message("Dropping transcripts from non-standard chromosomes...")
    std_chrom <- paste0("chr", c(1:22, "X", "Y"))

    # Checking ranges only has to be done for transcript features
    tx_ranges <- SummarizedExperiment::rowData(se)[
      !is.na(SummarizedExperiment::rowData(se)$transcript_id),
    ]$Ranges

    # unlist works for tx but NOT all features since hashes have multiple chroms
    tx_chroms <- unlist(S4Vectors::runValue(GenomicRanges::seqnames(tx_ranges)))
    drop_tx <- names(tx_chroms)[which(!tx_chroms %in% std_chrom)]
    se <- se[setdiff(rownames(se), drop_tx), ]
  }

  message("Done.")
  return(se)
}
