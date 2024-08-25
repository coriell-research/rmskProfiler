#' Create the gentrome fasta and Salmon index
#'
#' This function creates the gentrome (transcripts + unique rmsk sequences +
#' decoy) fasta and (optionally) the Salmon index. If index generation is
#' desired then the function assumes that salmon is on your PATH. This function
#' will create a gentrome.fa file, decoys.txt file, and (optionally) a salmon
#' index from these files in the resource directory.
#'
#' @param resource_dir Path to the directory containing index generation resources.
#' Output is saved to this location.
#' @param create_index Create salmon index after generating resources? Default TRUE.
#' This assumes that "salmon" is available on your PATH
#' @param threads Number of threads to use for salmon index generation. Default 1
#'
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' createGentrome(resource_dir = "/path/to/rmsk-resources")
#' }
createGentrome <- function(resource_dir, create_index = TRUE, threads = 1) {
  resources <- list.files(resource_dir, full.names = TRUE)
  genome_fa <- grep("primary_assembly.genome.fa.gz", resources, value = TRUE)
  tx_fa <- grep("transcripts.fa.gz", resources, value = TRUE)
  rmsk_fa <- grep("rmsk-unique.fa", resources, value = TRUE)

  if (length(genome_fa) != 1) {
    stop("<>.primary_assembly.genome.fa.gz file not found in given directory. Check that the file exists")
  }
  if (length(tx_fa) != 1) {
    stop("<>.transcripts.fa.gz file not found in given directory. Check that the file exists")
  }
  if (length(rmsk_fa) != 1) {
    stop("rmsk-unique.fa file not found in given directory. Check that the file exists")
  }

  # Gentrome generation ---------------------------------------------------------
  message("Reading in genome fasta...")
  genome_seqs <- Biostrings::readDNAStringSet(genome_fa, format = "fasta")

  # Subset for only the primary chromosomes
  genome_seqs <- genome_seqs[grepl("chr[0-9]+|chr[XY]", names(genome_seqs))]

  # Fix the names of the DNAStringSet (they import as "chr1 1", "chr2 2", etc.)
  names(genome_seqs) <- regmatches(names(genome_seqs), regexpr("chr[0-9]+|chr[XY]", names(genome_seqs)))

  message("Reading in transcripts fasta...")
  tx_seqs <- Biostrings::readDNAStringSet(tx_fa, format = "fasta")

  message("Reading in unique RepeatMasker instances fasta...")
  rmsk_seqs <- Biostrings::readDNAStringSet(rmsk_fa, format = "fasta")

  # Create the gentrome from combined seqs and write out
  gentrome <- c(tx_seqs, rmsk_seqs, genome_seqs)
  gentrome_fa <- file.path(resource_dir, "rmsk-gentrome.fa.gz")

  message("Writing out gentrome to ", gentrome_fa, "... (this may take some time)")
  Biostrings::writeXStringSet(gentrome, filepath = gentrome_fa, compress = TRUE)

  # Decoy generation -------------------------------------------------------------
  # Get the names of the genome fasta headers for the decoys file
  decoys <- names(genome_seqs)
  decoy_file <- file.path(resource_dir, "decoys.txt")
  message("Writing out decoys to ", decoy_file)
  write.table(
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
      "salmon index -t", gentrome_fa,
      "-d", decoy_file,
      "-p", threads,
      "-i", file.path(resource_dir, "rmsk.salmon_index"),
      "--gencode"
    )

    tryCatch(
      system(cmd),
      warning = function(w) print(w),
      error = function(e) stop("An error occurred during index generation! Check Salmon logs")
    )
  }

  return(invisible(NULL))
}
