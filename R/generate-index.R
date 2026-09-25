#' Generate a Salmon index and annotations for TE quantification
#'
#' This function is a wrapper around several other functions for downloading and
#' generating all of the files needed to create a Salmon index in order to
#' quantify transcripts and TEs at the loci level.
#'
#' @details
#' This function creates a Salmon index from all GENCODE transcripts + selected
#' RepeatMasker elements + a genomic decoy. It does so by calling the following
#' functions:
#'
#' \itemize{
#'  \item{"downloadResources()"}{Downloads all neccessary files for downstream processing}
#'  \item{"rmskToBed()"}{Extracts and filters RepeatMasker to BED file}
#'  \item{"extractUniqueSeqs()"}{Extracts all unique RepeatMasker sequences from the genome}
#'  \item{"createAnnotation()"}{Creates annotations of all TE loci with respect to annotations in GTF}
#'  \item{"createGentrome()"}{Creates the files needed for decoy-aware Salmon index and calls 'Salmon index'}
#' }
#'
#' @param index_dir Directory in which to create the Salmon index. All other
#' resources are saved to the rmskProfiler cache.
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus) designating which
#' species to download resources for
#' @param check_integrity TRUE/FALSE, if TRUE check the md5sums of the GENCODE files after downloading
#' @param exclude A character vector specifying which elements (repClass) to exclude from the
#' resulting BED file. Default "Simple_repeat", "Low_complexity", "Satellite",
#' "RNA", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", and "Unknown".
#' @param min_len Minimum sequence length of a record. Default 32. Records
#' shorter than this length are excluded from the resulting BED file and index.
#' @param create_index TRUE/FALSE Create salmon index after generating resources? Default TRUE.
#' This assumes that "salmon" is available on your PATH
#' @param threads Number of threads to use for salmon index generation. Default 1
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#'
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' # Run pipeline for downloading and creating Salmon index and annotations
#' generateIndex(index_dir = "hg38-salmon_index", species = "Hs")
#' }
generateIndex <- function(
  index_dir,
  species = c("Hs", "Mm"),
  check_integrity = TRUE,
  exclude = .DEFAULT_EXCLUDE,
  min_len = 32,
  create_index = TRUE,
  threads = 1,
  cache = NULL
) {
  species <- match.arg(species)
  bfc <- .getCache(cache)

  message("Downloading resources --------------------------------------")
  downloadResources(
    species = species,
    check_integrity = check_integrity,
    cache = bfc
  )
  message("Converting rmsk ranges to BED ------------------------------")
  rmskToBed(
    species = species,
    exclude = exclude,
    min_len = min_len,
    cache = bfc
  )
  message("Extracting unique rmsk sequences from genome ---------------")
  extractUniqueSeqs(
    species = species,
    exclude = exclude,
    min_len = min_len,
    cache = bfc
  )
  message("Annotating unique sequences with genomic features ----------")
  createAnnotation(
    species = species,
    exclude = exclude,
    min_len = min_len,
    cache = bfc
  )
  message("Creating gentrome for Salmon index generation --------------")
  createGentrome(
    species = species,
    index_dir = index_dir,
    exclude = exclude,
    min_len = min_len,
    create_index = create_index,
    threads = threads,
    cache = bfc
  )

  return(invisible(NULL))
}
