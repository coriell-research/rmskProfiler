#' Extract the contents of the RepeatMasker out file into BED format
#'
#' This function queries AnnotationHub for UCSC RepeatMasker annotations and
#' removes all records shorter than \code{min_len}, records whose repClass is
#' in \code{exclude}, and records not located in standard chromosomes. The
#' remaining ranges are then saved to a BED file ("rmsk.bed") in the
#' rmskProfiler cache for downstream processing.
#'
#' @details
#' The RepeatMasker annotations used are:
#'
#' AH111333 : UCSC RepeatMasker annotations (Oct2022) for Human (hg38)
#' AH99012 : UCSC RepeatMasker annotations (Apr2021) for Mouse (mm10)
#'
#' Other RepeatMasker annotations, or a subset of them, can be used instead by
#' passing them as a GRanges to \code{rmsk}.
#'
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus)
#' @param exclude A character vector specifying which elements (repClass) to
#' exclude from the resulting BED file. Default "Simple_repeat",
#' "Low_complexity", "Satellite", "RNA", "rRNA", "snRNA", "scRNA", "srpRNA",
#' "tRNA", and "Unknown".
#' @param min_len Minimum sequence length of a record. Default 32. Records
#' shorter than this length are excluded from the resulting BED file.
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#' @param rmsk NULL (default) to use the AnnotationHub RepeatMasker annotations
#' for \code{species}, or a GRanges of RepeatMasker records with metadata
#' columns repName, repFamily, and repClass.
#' @return Invisibly \code{NULL}. Called for its side effect of writing the
#' filtered RepeatMasker records (rmsk.bed) to the cache.
#' @export
#'
#' @examples
#' # Use the RepeatMasker records shipped with the example data. With
#' # rmsk = NULL, the full UCSC RepeatMasker annotations are retrieved from
#' # AnnotationHub.
#' bfc <- exampleCache()
#' rmsk <- readRDS(system.file("extdata", "rmsk-chr22.rds",
#'     package = "rmskProfiler"
#' ))
#' rmskToBed(species = "Hs", cache = bfc, rmsk = rmsk)
#' listResources(cache = bfc)
rmskToBed <- function(
    species = c("Hs", "Mm"),
    exclude = .DEFAULT_EXCLUDE,
    min_len = 32,
    cache = NULL,
    rmsk = NULL
) {
    species <- match.arg(species)

    if (is.null(rmsk)) {
        ah <- AnnotationHub::AnnotationHub()
        if (species == "Hs") {
            # UCSC RepeatMasker annotations (Oct2022) for Human (hg38)
            gr <- ah[["AH111333"]]
        } else {
            # UCSC RepeatMasker annotations (Apr2021) for Mouse (mm10)
            gr <- ah[["AH99012"]]
        }
    } else {
        stopifnot(
            "rmsk must be a GRanges" = methods::is(rmsk, "GRanges"),
            "rmsk must have repName, repFamily, and repClass columns" = all(
                c("repName", "repFamily", "repClass") %in%
                    names(GenomicRanges::mcols(rmsk))
            )
        )
        gr <- rmsk
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
