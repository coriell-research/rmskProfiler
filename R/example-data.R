#' Create a cache seeded with the example data
#'
#' Creates a new BiocFileCache in a temporary directory containing the example
#' data shipped with rmskProfiler in place of the GENCODE downloads. The example
#' data are a ~1 Mb region of hg38 chr22 (the rest of chr22 is masked with Ns)
#' with the GENCODE v48 annotations and transcripts of the genes in that region.
#' The cache can be passed to the \code{cache} argument of any rmskProfiler
#' function with \code{species = "Hs"}.
#'
#' @details
#' The example files are found in \code{system.file("extdata", package =
#' "rmskProfiler")}. The UCSC RepeatMasker records for the region are in
#' "rmsk-chr22.rds" and can be passed to \code{rmskToBed()}. Salmon quants for
#' 6 simulated samples quantified against the index built from the example data
#' are returned by \code{exampleQuants()}. See "inst/scripts/" for how these
#' files were created.
#'
#' @param build TRUE/FALSE. If TRUE, also run \code{rmskToBed()},
#' \code{extractUniqueSeqs()}, and \code{createAnnotation()} with the default
#' settings so that the cache contains everything needed by
#' \code{importQuants()}. Default FALSE.
#'
#' @return BiocFileCache
#' @export
#'
#' @examples
#' bfc <- exampleCache()
#' listResources(cache = bfc)
exampleCache <- function(build = FALSE) {
    bfc <- BiocFileCache::BiocFileCache(
        tempfile("rmskProfiler-example"),
        ask = FALSE
    )
    files <- c(
        genome = "genome-chr22.fa.gz",
        gtf = "gencode-chr22.gtf.gz",
        transcripts = "transcripts-chr22.fa.gz"
    )
    for (type in names(files)) {
        path <- system.file(
            "extdata",
            files[[type]],
            package = "rmskProfiler",
            mustWork = TRUE
        )
        added <- BiocFileCache::bfcadd(
            bfc,
            .resourceRname("Hs", type),
            fpath = path,
            action = "copy"
        )
        .recordResource(bfc, names(added), "Hs")
    }

    if (isTRUE(build)) {
        rmsk <- readRDS(system.file(
            "extdata",
            "rmsk-chr22.rds",
            package = "rmskProfiler",
            mustWork = TRUE
        ))
        rmskToBed(species = "Hs", cache = bfc, rmsk = rmsk)
        extractUniqueSeqs(species = "Hs", cache = bfc)
        createAnnotation(species = "Hs", cache = bfc)
    }

    return(bfc)
}


#' Path to the example Salmon quants
#'
#' Returns the path to a directory containing Salmon quants for 6 simulated
#' samples (3 per group) quantified with 10 Gibbs samples against the index
#' built from the example data. See \code{exampleCache()} for details on the
#' example data.
#'
#' @return Path to the directory containing one Salmon quant directory per
#' sample. It can be passed to the \code{quant_dir} argument of
#' \code{importQuants()}.
#' @export
#'
#' @examples
#' list.files(exampleQuants())
exampleQuants <- function() {
    system.file("extdata", "quants", package = "rmskProfiler", mustWork = TRUE)
}
