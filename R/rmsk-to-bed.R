#' Extract the contents of the RepeatMasker out file into BED format
#'
#' This function queries AnnotationHub for UCSC RepeatMasker annotations and removes all records
#' shorter than 31 bp, records derived from the following elements
#' "Simple_repeat", "Low_complexity", "Satellite", "RNA", "rRNA", "snRNA", "scRNA", "srpRNA",
#' "tRNA", "Unknown" and records not located in standard chromosomes. The extracted ranges are
#' then saved to a BED file ("rmsk.bed") for downstream processing.
#'
#' @details
#' The RepeatMasker annotations used are:
#'
#' AH111333 : UCSC RepeatMasker annotations (Oct2022) for Human (hg38)
#' AH99012 : UCSC RepeatMasker annotations (Apr2021) for Mouse (mm10)
#'
#' @param resource_dir Path to the directory where gentrome resources were
#' downloaded. this should be the same path specified by downloadResources().
#' @param exclude A character vector specifying which elements to exclude from the
#' resulting BEd file. Default "Simple_repeat", "Low_complexity", "Satellite",
#' "RNA", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", and "Unknown".
#' @param min_len Minimum sequence length of a record. Default 32. Records
#' shorter than this length are excluded from the resulting BED file.
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' rmskToBed(resource_dir = "/path/to/rmsk-resources")
#' }
#'
rmskToBed <- function(
  resource_dir,
  species = c("Hs", "Mm"),
  exclude = c(
    "Simple_repeat",
    "Low_complexity",
    "Satellite",
    "RNA",
    "rRNA",
    "snRNA",
    "scRNA",
    "srpRNA",
    "tRNA",
    "Unknown"
  ),
  min_len = 32
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
  gr <- gr[!(gr$repClass %in% exclude) & width(gr) >= min_len]
  gr$name <- paste(gr$repName, gr$repFamily, gr$repClass, sep = ":")

  outfile <- file.path(resource_dir, "rmsk.bed")
  rtracklayer::export.bed(gr, outfile)
}
