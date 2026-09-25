#' Create the gentrome fasta and Salmon index
#'
#' This function creates the gentrome (transcripts + unique rmsk sequences +
#' decoy) fasta and (optionally) the Salmon index. If index generation is
#' desired then the function assumes that Salmon is on your PATH. This function
#' will save a gentrome fasta file and decoys.txt file to the rmskProfiler
#' cache and (optionally) create a Salmon index from these files in
#' \code{index_dir}.
#'
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus)
#' @param index_dir Path to the directory where the Salmon index will be created.
#' Only used if \code{create_index = TRUE}.
#' @param exclude A character vector specifying which elements (repClass) were
#' excluded by \code{rmskToBed()}. Default "Simple_repeat", "Low_complexity",
#' "Satellite", "RNA", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", and "Unknown".
#' @param min_len Minimum sequence length of a record used by \code{rmskToBed()}.
#' Default 32.
#' @param create_index Create salmon index after generating resources? Default TRUE.
#' This assumes that "salmon" is available on your PATH
#' @param threads Number of threads to use for salmon index generation. Default 1
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#'
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' createGentrome(species = "Hs", index_dir = "/path/to/rmsk-salmon_index")
#' }
createGentrome <- function(
  species = c("Hs", "Mm"),
  index_dir,
  exclude = .DEFAULT_EXCLUDE,
  min_len = 32,
  create_index = TRUE,
  threads = 1,
  cache = NULL
) {
  species <- match.arg(species)
  if (isTRUE(create_index) && missing(index_dir)) {
    stop("index_dir must be provided when create_index = TRUE")
  }

  bfc <- .getCache(cache)
  key <- .settingsKey(exclude, min_len)
  genome_fa <- .getResource(
    bfc,
    .gencodeRname(species, "genome"),
    hint = "Run downloadResources() first."
  )
  tx_fa <- .getResource(
    bfc,
    .gencodeRname(species, "transcripts"),
    hint = "Run downloadResources() first."
  )
  rmsk_fa <- .getResource(
    bfc,
    .rname(species, "rmsk-unique.fa.gz", key),
    hint = "Run extractUniqueSeqs() with the same species, exclude, and min_len first."
  )

  # Gentrome generation ---------------------------------------------------------
  message("Reading in genome fasta...")
  genome_seqs <- Biostrings::readDNAStringSet(genome_fa, format = "fasta")

  # Subset for only the primary chromosomes
  genome_seqs <- genome_seqs[grepl("chr[0-9]+|chr[XY]", names(genome_seqs))]

  # Fix the names of the DNAStringSet (they import as "chr1 1", "chr2 2", etc.)
  names(genome_seqs) <- regmatches(
    names(genome_seqs),
    regexpr("chr[0-9]+|chr[XY]", names(genome_seqs))
  )

  message("Reading in transcripts fasta...")
  tx_seqs <- Biostrings::readDNAStringSet(tx_fa, format = "fasta")

  message("Reading in unique RepeatMasker instances fasta...")
  rmsk_seqs <- Biostrings::readDNAStringSet(rmsk_fa, format = "fasta")

  # Create the gentrome from combined seqs and write out
  gentrome <- c(tx_seqs, rmsk_seqs, genome_seqs)
  gentrome_fa <- .newResource(
    bfc,
    .rname(species, "rmsk-gentrome.fa.gz", key),
    ext = ".fa.gz"
  )

  message(
    "Writing out gentrome to ",
    gentrome_fa,
    "... (this may take some time)"
  )
  Biostrings::writeXStringSet(gentrome, filepath = gentrome_fa, compress = TRUE)

  # Decoy generation -------------------------------------------------------------
  # Get the names of the genome fasta headers for the decoys file
  decoys <- names(genome_seqs)
  decoy_file <- .newResource(
    bfc,
    .rname(species, "decoys.txt", key),
    ext = ".txt"
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
    cmd <- paste(
      "salmon index -t",
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

    tryCatch(
      system(cmd),
      warning = function(w) print(w),
      error = function(e) {
        stop("An error occurred during index generation! Check Salmon logs")
      }
    )
  }

  return(invisible(NULL))
}
