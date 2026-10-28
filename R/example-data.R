#' Create a cache seeded with the example data
#'
#' Creates a new BiocFileCache in a temporary directory containing the example
#' GTF shipped with rmskProfiler in place of the GENCODE download. The example
#' data are a ~1 Mb region of hg38 chr22 with the GENCODE v48 annotations of the
#' genes in that region. The cache can be passed to the \code{cache} argument of
#' any rmskProfiler function with \code{species = "Hs"}. Pass
#' \code{exampleGenome()} to the \code{genome} argument of functions that
#' extract sequences.
#'
#' @details
#' The example files are found in \code{system.file("extdata", package =
#' "rmskProfiler")}. The UCSC RepeatMasker records for the region are in
#' "rmsk-chr22.rds" and can be passed to the \code{rmsk} argument of
#' \code{rmskToBed()}, \code{extractUniqueSeqs()}, or \code{generateIndex()}.
#' Salmon quants for 6 simulated samples quantified against the index built
#' from the example data are returned by \code{exampleQuants()}. See
#' "inst/scripts/" for how these files were created.
#'
#' @param build TRUE/FALSE. If TRUE, also run \code{extractUniqueSeqs()} and
#' \code{createAnnotation()} on the example RepeatMasker records with the
#' default settings so that the cache contains everything needed by
#' \code{importQuants()}. The sequences are extracted from
#' \code{exampleGenome()}. Default FALSE.
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
    gtf <- system.file(
        "extdata",
        "gencode-chr22.gtf.gz",
        package = "rmskProfiler",
        mustWork = TRUE
    )
    added <- BiocFileCache::bfcadd(
        bfc,
        .resourceRname("Hs", "gtf"),
        fpath = gtf,
        action = "copy"
    )
    .recordResource(bfc, names(added), "Hs")

    if (isTRUE(build)) {
        rmsk <- readRDS(system.file(
            "extdata",
            "rmsk-chr22.rds",
            package = "rmskProfiler",
            mustWork = TRUE
        ))
        extractUniqueSeqs(
            species = "Hs",
            cache = bfc,
            genome = exampleGenome(),
            rmsk = rmsk
        )
        createAnnotation(species = "Hs", cache = bfc)
    }

    return(bfc)
}


#' Path to the example genome
#'
#' Returns the path to the example genome: hg38 chr22 from the start of the
#' chromosome to the end of the example region (chr22:37,000,000-37,999,999),
#' with every base before the region replaced by N. It can be passed to the
#' \code{genome} argument of \code{generateIndex()},
#' \code{extractUniqueSeqs()}, and \code{createGentrome()} together with the
#' cache returned by \code{exampleCache()}.
#'
#' @return Path to a gzipped FASTA file
#' @export
#'
#' @examples
#' Biostrings::fasta.seqlengths(exampleGenome())
exampleGenome <- function() {
    system.file(
        "extdata",
        "genome-chr22.fa.gz",
        package = "rmskProfiler",
        mustWork = TRUE
    )
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
