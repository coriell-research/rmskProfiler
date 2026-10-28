#' Generate a Salmon index and annotations for TE quantification
#'
#' This function is a wrapper around several other functions for downloading and
#' generating all of the files needed to create a Salmon index in order to
#' quantify transcripts and TEs at the loci level.
#'
#' @details
#' This function creates a Salmon index from all GENCODE transcripts + selected
#' RepeatMasker elements + a genomic decoy. It does so by calling the following
#' functions:
#'
#' \describe{
#'   \item{\code{downloadResources()}}{Downloads the GENCODE GTF}
#'   \item{\code{extractUniqueSeqs()}}{Filters the RepeatMasker records with
#'   \code{filterRmsk()}, extracts their sequences from the genome, and finds
#'   the loci with identical sequences}
#'   \item{\code{createAnnotation()}}{Creates annotations of all TE loci with
#'   respect to annotations in GTF}
#'   \item{\code{createGentrome()}}{Creates the files needed for decoy-aware
#'   Salmon index and calls 'Salmon index'}
#' }
#'
#' The genome sequence is taken from the BSgenome package for the species,
#' BSgenome.Hsapiens.UCSC.hg38 (GRCh38) or BSgenome.Mmusculus.UCSC.mm10
#' (GRCm38). The transcript sequences are extracted from the genome using the
#' GENCODE GTF.
#'
#' @param index_dir Directory in which to create the Salmon index. The build
#' resources are saved to the rmskProfiler cache, see \code{listResources()}.
#' While the index is built, the gentrome is written to a temporary file next
#' to \code{index_dir}, which needs ~5.5 GB of free space for human.
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus) designating
#' which
#' species to download resources for
#' @param check_integrity TRUE/FALSE, if TRUE check the md5sum of the GENCODE
#' GTF after downloading
#' @param exclude A character vector specifying which elements (repClass) to
#' exclude from the index. Default "Simple_repeat", "Low_complexity",
#' "Satellite", "RNA", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", and
#' "Unknown".
#' @param min_len Minimum sequence length of a record. Default 32. Records
#' shorter than this length are excluded from the index.
#' @param create_index TRUE/FALSE Create salmon index after generating
#' resources? Default TRUE.
#' This assumes that "salmon" is available on your PATH
#' @param threads Number of threads to use for salmon index generation. Default
#' 1
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#' @param rmsk NULL (default) or a GRanges of RepeatMasker records passed to
#' \code{extractUniqueSeqs()}. Filter the records with \code{exclude} and
#' \code{min_len} rather than beforehand, see \code{extractUniqueSeqs()}.
#'
#' @return Invisibly \code{NULL}. Called for its side effect of creating the
#' build resources in the cache and either the Salmon index in
#' \code{index_dir} (\code{create_index = TRUE}) or the gentrome in the cache
#' (\code{create_index = FALSE}).
#' @export
#'
#' @examples
#' # The example cache contains the example GTF and the RepeatMasker records
#' # of the example region
#' bfc <- exampleCache()
#' rmsk <- readRDS(system.file("extdata", "rmsk-chr22.rds",
#'     package = "rmskProfiler"
#' ))
#'
#' \dontrun{
#' # Run the whole pipeline on the example data. The decoys are all of the
#' # primary chromosomes of the genome, so this requires Salmon and takes ~1 h
#' generateIndex(
#'     index_dir = tempfile("index"),
#'     species = "Hs",
#'     check_integrity = FALSE,
#'     cache = bfc,
#'     rmsk = rmsk
#' )
#' listResources(cache = bfc)
#'
#' # Build the full human index with the default cache
#' generateIndex("hg38-salmon_index", species = "Hs")
#' }
generateIndex <- function(
    index_dir,
    species = c("Hs", "Mm"),
    check_integrity = TRUE,
    exclude = .DEFAULT_EXCLUDE,
    min_len = 32,
    create_index = TRUE,
    threads = 1,
    cache = NULL,
    rmsk = NULL
) {
    species <- match.arg(species)
    if (isTRUE(create_index)) {
        if (missing(index_dir)) {
            stop("index_dir must be provided when create_index = TRUE")
        }
        .checkSalmon("Or set create_index = FALSE to skip building the index.")
    }
    bfc <- .getCache(cache)

    message("Downloading the GTF ----------------------------------------")
    downloadResources(
        species = species,
        check_integrity = check_integrity,
        cache = bfc
    )
    message("Extracting unique rmsk sequences from genome ---------------")
    extractUniqueSeqs(
        species = species,
        exclude = exclude,
        min_len = min_len,
        cache = bfc,
        rmsk = rmsk
    )
    message("Annotating unique sequences with genomic features ----------")
    createAnnotation(
        species = species,
        exclude = exclude,
        min_len = min_len,
        cache = bfc
    )
    message("Creating gentrome for Salmon index generation --------------")
    createGentrome(
        species = species,
        index_dir = index_dir,
        exclude = exclude,
        min_len = min_len,
        create_index = create_index,
        threads = threads,
        cache = bfc
    )

    return(invisible(NULL))
}
