#' Extract the contents of the RepeatMasker out file into BED format
#'
#' This function queries AnnotationHub for UCSC RepeatMasker annotations and removes all records
#' shorter than 31 bp, records derived from the following elements
#' "Simple_repeat", "Low_complexity", "Satellite", "RNA", "rRNA", "snRNA", "scRNA", "srpRNA",
#' "tRNA", "Unknown" and records not located in standard chromosomes. The extracted ranges are
#' then saved to a BED file ("rmsk.bed") in the rmskProfiler cache for downstream processing.
#'
#' @details
#' The RepeatMasker annotations used are:
#'
#' AH111333 : UCSC RepeatMasker annotations (Oct2022) for Human (hg38)
#' AH99012 : UCSC RepeatMasker annotations (Apr2021) for Mouse (mm10)
#'
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus)
#' @param exclude A character vector specifying which elements to exclude from the
#' resulting BEd file. Default "Simple_repeat", "Low_complexity", "Satellite",
#' "RNA", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", and "Unknown".
#' @param min_len Minimum sequence length of a record. Default 32. Records
#' shorter than this length are excluded from the resulting BED file.
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' rmskToBed(species = "Hs")
#' }
#'
rmskToBed <- function(
  species = c("Hs", "Mm"),
  exclude = .DEFAULT_EXCLUDE,
  min_len = 32,
  cache = NULL
) {
  species <- match.arg(species)

  ah <- AnnotationHub::AnnotationHub()
  if (species == "Hs") {
    # UCSC RepeatMasker annotations (Oct2022) for Human (hg38)
    gr <- ah[["AH111333"]]
  } else {
    # UCSC RepeatMasker annotations (Apr2021) for Mouse (mm10)
    gr <- ah[["AH99012"]]
  }

  # Standard chromosomes and any features not excluded
  gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- gr[!(gr$repClass %in% exclude) & GenomicRanges::width(gr) >= min_len]
  gr$name <- paste(gr$repName, gr$repFamily, gr$repClass, sep = ":")

  bfc <- .getCache(cache)
  outfile <- .newResource(bfc, species, "rmsk.bed", ".bed", exclude, min_len)
  rtracklayer::export.bed(gr, outfile)

  return(invisible(NULL))
}
