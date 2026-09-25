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
#' Quants must be generated with Gibbs sampling (at least 2 Gibbs samples per
#' sample) since the Gibbs samples are used to estimate overdispersion.
#' \code{importQuants()} will throw an error listing any samples quantified
#' without them.
#'
#' The rmskProfiler index contains transcripts from the "patch_hapl_scaff"
#' annotations. By default, \code{importQuants()} will exclude these transcripts
#' from the resulting object. The exclusion criteria a transcript where any of
#' its ranges does not lie within a standard chromosome will be excluded.
#'
#' @param quant_dir Path to the directories created by \code{salmonQuant()} or Salmon
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus). Must match
#' the species used to generate the index.
#' @param exclude A character vector specifying which elements (repClass) were
#' excluded when generating the index. Must match the value passed to
#' \code{generateIndex()}. Default "Simple_repeat", "Low_complexity",
#' "Satellite", "RNA", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", and "Unknown".
#' @param min_len Minimum sequence length of a record used when generating the
#' index. Must match the value passed to \code{generateIndex()}. Default 32.
#' @param remove_zeros Should all zero rows be removed before returning object. Default TRUE
#' @param std_chromosomes Should only transcripts from standard chromosomes be
#' kept in the resulting object. Default TRUE.
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#'
#' @return SummarizedExperiment
#' @export
#'
importQuants <- function(
  quant_dir,
  species = c("Hs", "Mm"),
  exclude = .DEFAULT_EXCLUDE,
  min_len = 32,
  remove_zeros = TRUE,
  std_chromosomes = TRUE,
  cache = NULL
) {
  species <- match.arg(species)
  bfc <- .getCache(cache)
  rdfile <- .getResource(
    bfc,
    species,
    "rmsk-rowData.rds",
    exclude,
    min_len,
    hint = paste(
      "Check that species, exclude, and min_len match the values used to",
      "generate the index, or run createAnnotation() first."
    )
  )

  message("Importing quants...")
  paths <- list.dirs(path = quant_dir, full.names = TRUE, recursive = FALSE)

  # Overdispersion can only be estimated from Gibbs samples, so check for them
  # before reading in any quants
  meta_files <- file.path(paths, "aux_info", "meta_info.json")
  n_gibbs <- vapply(
    meta_files[file.exists(meta_files)],
    function(f) {
      n <- jsonlite::fromJSON(f)$num_bootstraps
      if (is.null(n)) 0L else as.integer(n)
    },
    integer(1)
  )
  no_gibbs <- n_gibbs < 2L
  if (any(no_gibbs)) {
    stop(
      "The following samples were quantified with fewer than 2 Gibbs samples: ",
      paste(basename(dirname(dirname(names(n_gibbs)[no_gibbs]))), collapse = ", "),
      ". Gibbs samples are required to estimate overdispersion. Re-run Salmon ",
      "with --numGibbsSamples (e.g. salmonQuant(n_gibbs = 30)).",
      call. = FALSE
    )
  }

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
