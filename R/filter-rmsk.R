#' Filter the UCSC RepeatMasker annotations
#'
#' This function queries AnnotationHub for UCSC RepeatMasker annotations and
#' removes all records shorter than \code{min_len}, records whose repClass is
#' in \code{exclude}, and records not located in standard chromosomes. The
#' remaining records are returned as a GRanges. \code{extractUniqueSeqs()}
#' calls this function with its own \code{exclude} and \code{min_len}, so use
#' it on its own to inspect the records that a build with these settings will
#' include.
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
#' The records can be written to a BED file with
#' \code{rtracklayer::export.bed()}.
#'
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus)
#' @param exclude A character vector specifying which elements (repClass) to
#' exclude. Default "Simple_repeat", "Low_complexity", "Satellite", "RNA",
#' "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", and "Unknown".
#' @param min_len Minimum sequence length of a record. Default 32. Records
#' shorter than this length are excluded.
#' @param rmsk NULL (default) to use the AnnotationHub RepeatMasker annotations
#' for \code{species}, or a GRanges of RepeatMasker records with metadata
#' columns repName, repFamily, and repClass.
#' @return GRanges of the filtered RepeatMasker records with an additional
#' metadata column, name, of the form "repName:repFamily:repClass".
#' @export
#'
#' @examples
#' # Use the RepeatMasker records shipped with the example data. With
#' # rmsk = NULL, the full UCSC RepeatMasker annotations are retrieved from
#' # AnnotationHub.
#' rmsk <- readRDS(system.file("extdata", "rmsk-chr22.rds",
#'     package = "rmskProfiler"
#' ))
#' filtered <- filterRmsk(species = "Hs", rmsk = rmsk)
#' filtered
#'
#' # Keep only LINEs and SINEs of at least 100 bp
#' classes <- unique(rmsk$repClass)
#' filterRmsk(
#'     species = "Hs",
#'     exclude = setdiff(classes, c("LINE", "SINE")),
#'     min_len = 100,
#'     rmsk = rmsk
#' )
filterRmsk <- function(
    species = c("Hs", "Mm"),
    exclude = .DEFAULT_EXCLUDE,
    min_len = 32,
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

    return(gr)
}
