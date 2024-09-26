#' Generate a Salmon index and annotations for TE quantification
#'
#' This function is a wrapper around several other functions for downloading and
#' generating all of the files needed to create a Salmon index in order to
#' quantify transcripts and TEs at the loci level. This function requires
#' calls to Python programs so before running ensure that you have created and
#' activated the neccessary "r-rmskProfiler" environment with
#' \code{install_rmskProfiler()} and \code{reticulate::use_condaenv("r-rmskProfiler", ...)}
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
#'  \item{"createGentrome()"}{Creates the files needed for decoy-aware Salmon index ans calls 'Salmon index'}
#' }
#'
#' @param out_dir Directory in which to save all resources including final generated Salmon index
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus) designating which
#' species to download resources for
#' @param check_integrity TRUE/FALSE, if TRUE check the md5sums of the GENCODE files
#' @param exclude A character vector specifying which elements to exclude from the
#' resulting BEd file. Default "Simple_repeat", "Low_complexity", "Satellite",
#' "RNA", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", and "Unknown".
#' @param min_len Minimum sequence length of a record. Default 31. Records
#' shorter than this length are excluded from the resulting BED file and index.
#' @param create_index TRUE/FALSE Create salmon index after generating resources? Default TRUE.
#' This assumes that "salmon" is available on your PATH
#' @param threads Number of threads to use for salmon index generation. Default 1
#'
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Ensure Python dependencies are installed
#' install_rmskProfiler()
#'
#' # Activate the Python environment
#' reticulate::use_condaenv("r-rmskProfiler")
#'
#' # Run pipeline for downloading and creating Salmon index and annotations
#' generateIndex(out_dir = "rmsk-resources")
#' }
generateIndex <- function(out_dir, species = c("Hs", "Mm"), check_integrity = TRUE,
                          exclude = c("Simple_repeat", "Low_complexity",
                                      "Satellite", "RNA", "rRNA", "snRNA",
                                      "scRNA", "srpRNA", "tRNA", "Unknown"),
                          min_len = 31, create_index = TRUE, threads = 1) {

  message("Downloading resources ----------")
  downloadResources(out_dir = out_dir, species = species, check_integrity = check_integrity)
  message("Converting rmsk out to BED ----------")
  rmskToBed(resource_dir = out_dir, exclude = exclude, min_len = min_len)
  message("Extracting unique rmsk sequences from genome ----------")
  extractUniqueSeqs(resource_dir = out_dir)
  message("Annotating unique sequences with genomic features ----------")
  createAnnotation(resource_dir = out_dir)
  message("Creating gentrome for Salmon index generation ----------")
  createGentrome(resource_dir = out_dir, create_index = create_index, threads = threads)

  return(invisible(NULL))
}
