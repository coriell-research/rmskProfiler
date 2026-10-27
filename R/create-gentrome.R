#' Create the gentrome fasta and Salmon index
#'
#' This function creates the gentrome (transcripts + unique rmsk sequences +
#' decoy) fasta and (optionally) the Salmon index. If index generation is
#' desired then the function assumes that Salmon is on your PATH. This function
#' will save a gentrome fasta file and decoys.txt file to the rmskProfiler
#' cache and (optionally) create a Salmon index from these files in
#' \code{index_dir}.
#'
#' The transcript sequences are extracted from the genome using the exons of
#' each transcript in the GENCODE GTF. These are identical to the GENCODE
#' transcript sequences. The primary chromosomes of the genome (chr1-22, X,
#' and Y) are used as the decoys.
#'
#' @inheritParams extractUniqueSeqs
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus)
#' @param index_dir Path to the directory where the Salmon index will be
#' created.
#' Only used if \code{create_index = TRUE}.
#' @param exclude A character vector specifying which elements (repClass) were
#' excluded by \code{rmskToBed()}. Default "Simple_repeat", "Low_complexity",
#' "Satellite", "RNA", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", and
#' "Unknown".
#' @param min_len Minimum sequence length of a record used by
#' \code{rmskToBed()}.
#' Default 32.
#' @param create_index Create salmon index after generating resources? Default
#' TRUE.
#' This assumes that "salmon" is available on your PATH
#' @param threads Number of threads to use for salmon index generation. Default
#' 1
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#'
#' @return Invisibly \code{NULL}. Called for its side effect of writing the
#' gentrome fasta and decoys file to the cache and, if
#' \code{create_index = TRUE}, creating the Salmon index in \code{index_dir}.
#' @export
#'
#' @examples
#' # Create the gentrome and decoys for the example data
#' bfc <- exampleCache(build = TRUE)
#' createGentrome(
#'     species = "Hs",
#'     create_index = FALSE,
#'     cache = bfc,
#'     genome = exampleGenome()
#' )
#' getResource("decoys.txt", species = "Hs", cache = bfc)
#'
#' # Also create the Salmon index if Salmon is available
#' if (nzchar(Sys.which("salmon"))) {
#'     index_dir <- tempfile("index")
#'     createGentrome(
#'         species = "Hs",
#'         index_dir = index_dir,
#'         cache = bfc,
#'         genome = exampleGenome()
#'     )
#' }
createGentrome <- function(
    species = c("Hs", "Mm"),
    index_dir,
    exclude = .DEFAULT_EXCLUDE,
    min_len = 32,
    create_index = TRUE,
    threads = 1,
    cache = NULL,
    genome = NULL
) {
    species <- match.arg(species)
    if (isTRUE(create_index)) {
        if (missing(index_dir)) {
            stop("index_dir must be provided when create_index = TRUE")
        }
        .checkSalmon("Or set create_index = FALSE to skip building the index.")
    }

    bfc <- .getCache(cache)
    gtf_file <- .getResource(
        bfc,
        species,
        "gtf",
        hint = "Run downloadResources() first."
    )
    rmsk_fa <- .getResource(
        bfc,
        species,
        "rmsk-unique.fa.gz",
        exclude,
        min_len,
        hint = paste(
            "Run extractUniqueSeqs() with the same species, exclude, and",
            "min_len first."
        )
    )
    genome <- .getGenome(species, genome)
    message("Using genome: ", genome$label)
    rmsk_genome <- .resourceGenome(
        bfc,
        species,
        "rmsk-unique.fa.gz",
        exclude,
        min_len
    )
    if (!is.na(rmsk_genome) && rmsk_genome != genome$label) {
        warning(
            "The unique RepeatMasker sequences were extracted from a ",
            "different genome (",
            rmsk_genome,
            "). Run extractUniqueSeqs() with the same genome first.",
            call. = FALSE
        )
    }

    # Gentrome generation ------------------------------------------------------
    # Transcripts are written in the order of the GTF, the same order as the
    # GENCODE transcript sequences. Salmon keeps the first of any duplicate
    # sequences so the order determines which duplicate is kept
    message("Extracting transcript sequences from the genome...")
    tx_order <- rtracklayer::import(
        gtf_file,
        format = "gtf",
        feature.type = "transcript",
        colnames = "transcript_id"
    )$transcript_id
    tx_seqs <- .transcriptSeqs(.getTxDb(species, bfc), genome, tx_order)

    message("Reading in unique RepeatMasker instances fasta...")
    rmsk_seqs <- Biostrings::readDNAStringSet(rmsk_fa, format = "fasta")

    gentrome_fa <- .newResource(
        bfc,
        species,
        "rmsk-gentrome.fa.gz",
        ".fa.gz",
        exclude,
        min_len,
        genome = genome$label
    )
    message(
        "Writing out gentrome to ",
        gentrome_fa,
        "... (this may take some time)"
    )
    Biostrings::writeXStringSet(
        c(tx_seqs, rmsk_seqs),
        filepath = gentrome_fa,
        compress = TRUE
    )
    rm(tx_seqs, rmsk_seqs)

    # The primary chromosomes (without chrM) are appended one at a time as the
    # decoys
    chrom_len <- .genomeLengths(genome)
    decoys <- grep(.PRIMARY_CHROMS, names(chrom_len), value = TRUE)
    decoys <- setdiff(decoys, "chrM")
    for (chrom in decoys) {
        chrom_seq <- .genomeSeqs(
            genome,
            GenomicRanges::GRanges(
                chrom,
                IRanges::IRanges(1L, chrom_len[[chrom]])
            )
        )
        names(chrom_seq) <- chrom
        Biostrings::writeXStringSet(
            chrom_seq,
            filepath = gentrome_fa,
            append = TRUE,
            compress = TRUE
        )
    }

    # Decoy generation ---------------------------------------------------------
    decoy_file <- .newResource(
        bfc,
        species,
        "decoys.txt",
        ".txt",
        exclude,
        min_len,
        genome = genome$label
    )
    message("Writing out decoys to ", decoy_file)
    utils::write.table(
        decoys,
        file = decoy_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )

    if (isTRUE(create_index)) {
        message("Creating salmon index...")
        status <- system2(
            "salmon",
            args = c(
                "index",
                "-t",
                shQuote(gentrome_fa),
                "-d",
                shQuote(decoy_file),
                "-p",
                threads,
                "-k",
                31,
                "-i",
                shQuote(index_dir),
                "--gencode",
                "--no-clip"
            )
        )
        if (status != 0L) {
            stop(
                "Salmon index generation failed with exit status ",
                status,
                ". Check the Salmon logs.",
                call. = FALSE
            )
        }
    }

    return(invisible(NULL))
}
