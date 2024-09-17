#' Quantify TEs and transcripts using Salmon with Gibbs sampling
#'
#' This function provides a simple wrapper around a system call to
#' \code{salmon quant ...}. The default arguments assume that you have paired
#' end reads (which you should if quantifying TEs is a priority). The defaults
#' DO NOT set any other flags for Salmon other than automatically detecting the
#' library type and setting the number of Gibbs samples. You can (and should)
#' pass additional arguments as character strings with the desired Salmon flags.
#' See the example below for how to do this.
#'
#' @param fq1 vector of file paths to fastq read 1 files
#' @param fq2 vector of file paths to fastq read 2 files
#' @param sample_names character vector of samples names matching each pair of fastq files
#' @param resource_dir path to the rmskProfiler resource directory containing an rmsk.salmon_index directory
#' @param out_dir path to save the quant directories for each sample. This will
#' be the parent directory, samples are saved in subdirectories like out_dir/<sample_name>_quants
#' @param nGibbs integer number of Gibbs samples to perform. Default 30. Published
#' work suggests this number can be set to around 180 / length(samples) but 30
#' ensures a good coverage with minimal cost in terms of speed.
#' @param ... Additional arguments passed to Salmon as character strings, e.g.
#' "--gcBias", "--seqBias" "--threads 12", see examples below
#'
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Specify paths to fastq files and names of the samples
#' fq1 <- c("sample1.R1.fq.gz", "sample2.R1.fq.gz", "sample3.R1.fq.gz")
#' fq2 <- c("sample1.R2.fq.gz", "sample2.R2.fq.gz", "sample3.R2.fq.gz")
#' sample_names <- c("sample1", "sample2", "sample3")
#'
#' # Perform with all default settings -- probably not recommended
#' salmonQuant(fq1, fq2, sample_names, resource_dir = "hg38-resources", out_dir = "quants")
#'
#' # More often though we will want to pass additional arguments to Salmon
#' # We can do so by providing additoinal flags as a character strings
#' salmonQuant(
#'   fq1 = fq1,
#'   fq2 = fq2,
#'   sample_names = sample_names,
#'   resource_dir = "hg38-resources",
#'   out_dir = "quants",
#'   "--gcBias",
#'   "--seqBias",
#'   "--posBias",
#'   "--threads 12"
#'   )
#' }
salmonQuant <- function(fq1, fq2, sample_names, resource_dir, out_dir,
                        nGibbs = 30, ...) {

  stopifnot("fq1, fq2, and sample_names differ in length" = (length(fq1) == length(fq2)) == (length(fq1) == length(sample_names)))

  idx <- file.path(resource_dir, "rmsk.salmon_index")
  stopifnot("rmsk.salmon_index does not exist in resource directory!" = dir.exists(idx))

  out_dirs <- file.path(out_dir, paste0(sample_names, "_quants"))

  dots <- list(...)
  more_args <- unlist(dots)

  for (i in seq_along(fq1)) {
    system2("salmon",
            args = c(
              "quant",
              "--libType", "A",
              "--mates1", fq1[i],
              "--mates2", fq2[i],
              "--output", out_dirs[i],
              "--index", idx,
              "--numGibbsSamples", nGibbs,
              more_args
            ))
  }

  return(invisible(NULL))
}

