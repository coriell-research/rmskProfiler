#' Extract the contents of the RepeatMasker out file into BED format
#'
#' This function extracts the contents of the RepeatMasker out file into a BED
#' file. Existing BED files for these tracks exist however this function allows
#' for arbitrary filtering out of certain records by feature name or feature
#' length. The default behavior removes all records shorter than 31 bp and
#' records derived from the following elements "Simple_repeat", "Low_complexity",
#' "Satellite", "RNA", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", "Unknown".
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
