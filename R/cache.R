#' Default RepeatMasker classes excluded from the index
#' @keywords internal
.DEFAULT_EXCLUDE <- c(
  "Simple_repeat",
  "Low_complexity",
  "Satellite",
  "RNA",
  "rRNA",
  "snRNA",
  "scRNA",
  "srpRNA",
  "tRNA",
  "Unknown"
)


#' Get the BiocFileCache used to store rmskProfiler resources
#'
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' If NULL, the default rmskProfiler cache located at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")} is used.
#'
#' @return BiocFileCache
#' @keywords internal
.getCache <- function(cache = NULL) {
  if (methods::is(cache, "BiocFileCache")) {
    return(cache)
  }
  if (is.null(cache)) {
    cache <- tools::R_user_dir("rmskProfiler", which = "cache")
  }
  BiocFileCache::BiocFileCache(cache, ask = FALSE)
}


#' GENCODE resources for a given species
#'
#' @param species Either "Hs" or "Mm"
#'
#' @return data.frame with columns type, url, file, and md5
#' @keywords internal
.gencodeResources <- function(species) {
  if (species == "Hs") {
    # https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/MD5SUMS
    base_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/"
    files <- c(
      "gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz",
      "gencode.v48.transcripts.fa.gz",
      "GRCh38.primary_assembly.genome.fa.gz"
    )
    md5sums <- c(
      "f7ffc813464f52e428c116bc3b83dce1",
      "e4a4d396cca5dd6d0889248b9e93b42a",
      "42e38e8dd5027dd2ae8aeb8f3a990d07"
    )
  } else {
    base_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/"
    files <- c(
      "gencode.vM25.annotation.gtf.gz",
      "gencode.vM25.transcripts.fa.gz",
      "GRCm38.primary_assembly.genome.fa.gz"
    )
    md5sums <- c(
      "0c38fc4ccbc731a2708fc91e7f1c2efd",
      "a821c0dde39c48b9d2c4b48d36b0180c",
      "3bc591be24b77f710b6ba5d41022fc5a"
    )
  }

  data.frame(
    type = c("gtf", "transcripts", "genome"),
    url = paste0(base_url, files),
    file = files,
    md5 = md5sums,
    row.names = c("gtf", "transcripts", "genome")
  )
}


#' Create a short key identifying the settings used to build TE resources
#'
#' @param exclude Character vector of excluded RepeatMasker classes
#' @param min_len Minimum TE record length
#'
#' @return character(1)
#' @keywords internal
.settingsKey <- function(exclude, min_len) {
  settings <- paste(c(sort(unique(exclude)), min_len), collapse = ",")
  substr(unname(tools::md5sum(bytes = charToRaw(settings))), 1, 10)
}


#' Create the name of a cached resource
#'
#' Resources that depend only on the species are named like
#' "rmskProfiler:Hs:<file>". Resources that also depend on the exclude and
#' min_len settings include the settings key, "rmskProfiler:Hs:<key>:<file>".
#'
#' @param species Either "Hs" or "Mm"
#' @param file Name(s) of the resource file(s)
#' @param key Optional settings key created by \code{.settingsKey()}
#'
#' @return character vector the same length as file
#' @keywords internal
.rname <- function(species, file, key = NULL) {
  prefix <- paste(c("rmskProfiler", species, key), collapse = ":")
  paste(prefix, file, sep = ":")
}


#' Name of a cached GENCODE resource
#'
#' @param species Either "Hs" or "Mm"
#' @param type One of "gtf", "transcripts", or "genome"
#'
#' @return character(1)
#' @keywords internal
.gencodeRname <- function(species, type) {
  .rname(species, .gencodeResources(species)[type, "file"])
}


#' Get the path of an existing resource in the cache
#'
#' @param bfc BiocFileCache
#' @param rname Name of the resource
#' @param hint Message appended to the error if the resource is missing
#'
#' @return Path to the cached file
#' @keywords internal
.getResource <- function(bfc, rname, hint = "") {
  hit <- BiocFileCache::bfcquery(bfc, rname, field = "rname", exact = TRUE)
  if (nrow(hit) > 0L) {
    path <- unname(BiocFileCache::bfcrpath(bfc, rids = hit$rid[1L]))
    if (file.exists(path)) {
      return(path)
    }
  }
  stop(
    "Resource '",
    rname,
    "' not found in cache ",
    BiocFileCache::bfccache(bfc),
    ". ",
    hint,
    call. = FALSE
  )
}


#' Get a path in the cache to write a generated resource to
#'
#' If the resource already exists in the cache its path is returned so that it
#' is overwritten. Otherwise a new cache entry is created.
#'
#' @param bfc BiocFileCache
#' @param rname Name of the resource
#' @param ext File extension of the resource, e.g. ".bed"
#'
#' @return Path to the cached file
#' @keywords internal
.newResource <- function(bfc, rname, ext) {
  hit <- BiocFileCache::bfcquery(bfc, rname, field = "rname", exact = TRUE)
  if (nrow(hit) > 0L) {
    return(unname(BiocFileCache::bfcrpath(bfc, rids = hit$rid[1L])))
  }
  unname(BiocFileCache::bfcnew(bfc, rname, ext = ext))
}
