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
#' @param sample_names character vector of samples names matching each pair of
#' fastq files
#' @param index_dir path to the Salmon index created by \code{generateIndex()}
#' or \code{createGentrome()}
#' @param out_dir path to save the quant directories for each sample. This will
#' be the parent directory, samples are saved in subdirectories like
#' out_dir/<sample_name>_quants
#' @param n_gibbs integer number of Gibbs samples to perform. Default 30.
#' Published
#' work suggests this number can be set to around 180 / length(samples) but 30
#' ensures a good coverage with minimal cost in terms of speed.
#' @param ... Additional arguments passed to Salmon as character strings, e.g.
#' "--gcBias", "--seqBias" "--threads 12", see examples below
#'
#' @return Invisibly \code{NULL}. Called for its side effect of writing the
#' Salmon quant directories for each sample to \code{out_dir}.
#' @export
#'
#' @examples
#' # Example reads simulated from the example region
#' fq_dir <- system.file("extdata", "fastq", package = "rmskProfiler")
#' list.files(fq_dir)
#'
#' \dontrun{
#' # Quantify the example reads against the index built from the example data.
#' # Salmon must be available on your PATH. The decoys are all of the primary
#' # chromosomes of the genome, so building the index takes ~1 h
#' bfc <- exampleCache(build = TRUE)
#' index_dir <- tempfile("index")
#' createGentrome(species = "Hs", index_dir = index_dir, cache = bfc)
#'
#' # Additional arguments are passed to Salmon as character strings
#' out_dir <- tempfile("quants")
#' salmonQuant(
#'     fq1 = file.path(fq_dir, "ctrl_1_R1.fq.gz"),
#'     fq2 = file.path(fq_dir, "ctrl_1_R2.fq.gz"),
#'     sample_names = "ctrl_1",
#'     index_dir = index_dir,
#'     out_dir = out_dir,
#'     n_gibbs = 2,
#'     "--threads 2"
#' )
#' list.files(out_dir)
#' }
salmonQuant <- function(
    fq1,
    fq2,
    sample_names,
    index_dir,
    out_dir,
    n_gibbs = 30,
    ...
) {
    .checkSalmon()
    stopifnot(
        "fq1, fq2, and sample_names differ in length" = length(fq1) ==
            length(fq2) &&
            length(fq1) == length(sample_names)
    )

    stopifnot("index_dir does not exist!" = dir.exists(index_dir))

    out_dirs <- file.path(out_dir, paste0(sample_names, "_quants"))

    dots <- list(...)
    more_args <- unlist(dots)

    for (i in seq_along(fq1)) {
        system2(
            "salmon",
            args = c(
                "quant",
                "--libType",
                "A",
                "--mates1",
                fq1[i],
                "--mates2",
                fq2[i],
                "--output",
                out_dirs[i],
                "--index",
                index_dir,
                "--numGibbsSamples",
                n_gibbs,
                more_args
            )
        )
    }

    return(invisible(NULL))
}


#' Stop if Salmon is not available on the PATH
#'
#' @param hint Message appended to the error
#'
#' @return Invisibly \code{NULL}. Called for its side effect of stopping with
#' an error if Salmon cannot be found.
#' @keywords internal
.checkSalmon <- function(hint = "") {
    if (!nzchar(Sys.which("salmon"))) {
        stop(
            "Salmon was not found on the PATH. Install Salmon ",
            "(https://salmon.readthedocs.io) and add it to the PATH. ",
            hint,
            call. = FALSE
        )
    }

    return(invisible(NULL))
}
