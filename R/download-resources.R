#' Attempt to download a file from the internet
#'
#' @param url URL
#' @param dest filepath of the destination file
#' @param ... additional args not used
#'
#' @return NULL
.tryDownload <- function(url, dest, ...) {
  tryCatch(
    {
      curl::curl_download(url, dest, mode = "wb", ...)
    },
    warning = function(w) {
      print(w)
      message("A warning occurred downloading ", dest)
      message("Removing ", dest, " just to be safe. Please try again.")

      if (file.exists(dest)) {
        file.remove(dest)
      }
    },
    error = function(e) {
      print(e)
      message("An error occurred downloading ", dest)
      message("Removing ", dest, ". Please try again.")

      if (file.exists(dest)) {
        file.remove(dest)
      }
    }
  )
}


#' Download files needed for index generation
#'
#' This function will attempt to download all of the necessary resources for
#' generating the rmsk-gentrome index. For humans, it will download the RepeatMasker
#' hg38.fa.out file from RepeatMasker.org, GENCODE v36 transcript sequences,
#' primary assembly, and annotation GTF. For mouse, it will download the
#' RepeatMasker mm10.fa.out file, GENCODE M25 transcript sequences, primary
#' assembly, and annotation GTF. If any of these file names already exist in the
#' out_dir they will be skipped.
#'
#' @param out_dir Directory to save files to. If it does not exist it will be created.
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus) designating which
#' species to download resources for
#' @param check_integrity TRUE/FALSE, if TRUE check the md5sums of the GENCODE files
#'
#' @return NULL
#' @export
#' @examples
#' \dontrun{
#' downloadResources(out_dir = "/path/to/rmsk-resources")
#' }
downloadResources <- function(out_dir, species = c("Hs", "Mm"), check_integrity = TRUE) {
  species <- match.arg(species)

  urls <- c(
    "http://repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz",
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz",
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.transcripts.fa.gz",
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz"
  )
  fnames <- basename(urls)

  # Hashsums only for GENCODE files - annotation, transcripts, assembly
  md5sums <- c("c03931958d4572148650d62eb6dec41a", "d9046028f532a5f42e6af438a7330c34", "e7d5fc50346e2d6dfd2861db31871dfa")

  if (species == "Mm") {
    urls <- c(
      "http://repeatmasker.org/genomes/mm10/RepeatMasker-rm405-db20140131/mm10.fa.out.gz",
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz",
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz",
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz"
    )
    fnames <- basename(urls)
    md5sums <- c("0c38fc4ccbc731a2708fc91e7f1c2efd", "a821c0dde39c48b9d2c4b48d36b0180c", "3bc591be24b77f710b6ba5d41022fc5a")
  }

  if (!dir.exists(out_dir)) {
    message(out_dir, " does not exist. Creating.")
    dir.create(out_dir, recursive = TRUE)
  }
  outfiles <- file.path(out_dir, fnames)

  for (i in seq_along(urls)) {
    if (file.exists(outfiles[i])) {
      message(outfiles[i], " already exists in ", out_dir, ". Skipping.")
      next
    }
    message("Attempting to download ", fnames[i], "...")
    .tryDownload(urls[i], outfiles[i])
  }

  if (isTRUE(check_integrity)) {
    message("Checking file integrity of downloaded files...")
    badfile <- md5sums != as.vector(tools::md5sum(outfiles[2:4]))
    if (any(badfile)) {
      msg <- paste(outfiles[2:4][which(badfile)], "Did not download properly. Remove this file and retry.")
      stop(msg)
    }
    message("Success!")
  }

  return(invisible(NULL))
}
