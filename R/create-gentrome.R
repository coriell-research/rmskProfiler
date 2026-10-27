#' Create the gentrome fasta and Salmon index
#'
#' This function creates the gentrome (transcripts + unique rmsk sequences +
#' decoy) fasta and (optionally) the Salmon index. If index generation is
#' desired then the function assumes that Salmon is on your PATH.
#'
#' The transcript sequences are extracted from the genome using the exons of
#' each transcript in the GENCODE GTF. These are identical to the GENCODE
#' transcript sequences. The unique RepeatMasker sequences are extracted from
#' the genome using the duplicate information created by
#' \code{extractUniqueSeqs()}. The primary chromosomes of the genome (chr1-22,
#' X, and Y) are used as the decoys and their names are saved to decoys.txt in
#' the rmskProfiler cache.
#'
#' With \code{create_index = TRUE}, the gentrome is written uncompressed to a
#' temporary file next to \code{index_dir}, which needs ~5.5 GB of free space
#' for human, and removed once the index is built. With
#' \code{create_index = FALSE}, the gentrome is saved to the rmskProfiler cache
#' as rmsk-gentrome.fa.gz so that the index can be built outside of R.
#'
#' @inheritParams extractUniqueSeqs
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus)
#' @param index_dir Path to the directory where the Salmon index will be
#' created.
#' Only used if \code{create_index = TRUE}.
#' @param exclude A character vector specifying which elements (repClass) were
#' excluded by \code{extractUniqueSeqs()}. Default "Simple_repeat",
#' "Low_complexity", "Satellite", "RNA", "rRNA", "snRNA", "scRNA", "srpRNA",
#' "tRNA", and "Unknown".
#' @param min_len Minimum sequence length of a record used by
#' \code{extractUniqueSeqs()}. Default 32.
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
#' decoys file to the cache and either creating the Salmon index in
#' \code{index_dir} (\code{create_index = TRUE}) or writing the gentrome fasta
#' to the cache (\code{create_index = FALSE}).
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
#' getResource("rmsk-gentrome.fa.gz", species = "Hs", cache = bfc)
#' getResource("decoys.txt", species = "Hs", cache = bfc)
#'
#' # Create the Salmon index instead if Salmon is available
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
    info_rds <- .getResource(
        bfc,
        species,
        "rmsk-duplicateInfo.rds",
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
        "rmsk-duplicateInfo.rds",
        exclude,
        min_len
    )
    if (!is.na(rmsk_genome) && rmsk_genome != genome$label) {
        warning(
            "The unique RepeatMasker sequences were found in a different ",
            "genome (",
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

    message("Extracting unique RepeatMasker sequences from the genome...")
    rmsk_seqs <- .uniqueSeqs(data.table::setDT(readRDS(info_rds)), genome)

    chrom_len <- .genomeLengths(genome)
    decoys <- grep(.PRIMARY_CHROMS, names(chrom_len), value = TRUE)
    decoys <- setdiff(decoys, "chrM")
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
    writeLines(decoys, decoy_file)

    if (isTRUE(create_index)) {
        # Salmon only reads the gentrome once, so write it uncompressed (much
        # faster than gzip) next to the index and remove it afterwards
        index_parent <- dirname(normalizePath(index_dir, mustWork = FALSE))
        dir.create(index_parent, recursive = TRUE, showWarnings = FALSE)
        gentrome_fa <- tempfile(
            "rmsk-gentrome",
            tmpdir = index_parent,
            fileext = ".fa"
        )
        on.exit(unlink(gentrome_fa), add = TRUE)
    } else {
        gentrome_fa <- .newResource(
            bfc,
            species,
            "rmsk-gentrome.fa.gz",
            ".fa.gz",
            exclude,
            min_len,
            genome = genome$label
        )
    }
    message(
        "Writing out gentrome to ",
        gentrome_fa,
        "... (this may take some time)"
    )
    .writeGentrome(
        gentrome_fa,
        c(tx_seqs, rmsk_seqs),
        genome,
        decoys,
        compress = !isTRUE(create_index)
    )
    rm(tx_seqs, rmsk_seqs)

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


#' Write the gentrome fasta
#'
#' The decoy chromosomes are extracted and appended one at a time to limit
#' memory use.
#'
#' @param path Path of the fasta file
#' @param seqs DNAStringSet of transcript and unique RepeatMasker sequences
#' @param genome Genome returned by \code{.getGenome()}
#' @param decoys Names of the decoy chromosomes
#' @param compress TRUE/FALSE, gzip the fasta file?
#'
#' @return NULL
#' @keywords internal
.writeGentrome <- function(path, seqs, genome, decoys, compress) {
    Biostrings::writeXStringSet(seqs, filepath = path, compress = compress)
    chrom_len <- .genomeLengths(genome)
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
            filepath = path,
            append = TRUE,
            compress = compress
        )
    }

    return(invisible(NULL))
}
