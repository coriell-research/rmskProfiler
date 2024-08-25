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
#' @param min_len Minimum sequence length of a record. Default 31. Records
#' shorter than this length are excluded from the resulting BED file.
#' @param overwrite Overwrite existing BED file if found in resource directory.
#' Default FALSE
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' rmskToBed(resource_dir = "/path/to/rmsk-resources")
#' }
rmskToBed <- function(resource_dir, exclude = c(
                       "Simple_repeat", "Low_complexity",
                       "Satellite", "RNA", "rRNA", "snRNA", "scRNA", "srpRNA",
                       "tRNA", "Unknown"
                     ),
                     min_len = 31,
                     overwrite = FALSE
                     ) {

  rmsk_file <- list.files(resource_dir, pattern = "*.out.gz", full.names = TRUE)
  if (length(rmsk_file) != 1) {
    stop("rmsk.out.fa.gz file not found in given directory. Check that the file exists")
  }

  out_file <- gsub(".gz", ".bed", rmsk_file)
  if (file.exists(out_file)) {
    if (isFALSE(overwrite)) {
      message(out_file, " exists! Skipping BED file generation.")
      return(invisible(NULL))
    }
  }

  message("Extracting contents of ", rmsk_file, " to a BED file...")
  tryCatch(
    rmsk_profiler$rmsk2bed(rmsk_file, exclude, min_len),
    warning = function(w) print(w),
    error = function(e) print(e)
    )
  message("BED file generation complete.")

  return(invisible(NULL))
}
