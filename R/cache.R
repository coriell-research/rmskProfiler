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
#' The genome and transcript sequences are not downloaded. They come from the
#' BSgenome package for the species, see \code{.getGenome()}.
#'
#' @param species Either "Hs" or "Mm"
#'
#' @return data.frame with columns type, url, file, and md5
#' @keywords internal
.gencodeResources <- function(species) {
    if (species == "Hs") {
        # MD5SUMS file in base_url
        base_url <- paste0(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/",
            "Gencode_human/release_48/"
        )
        file <- "gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz"
        md5sum <- "f7ffc813464f52e428c116bc3b83dce1"
    } else {
        base_url <- paste0(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/",
            "Gencode_mouse/release_M25/"
        )
        file <- "gencode.vM25.annotation.gtf.gz"
        md5sum <- "0c38fc4ccbc731a2708fc91e7f1c2efd"
    }

    data.frame(
        type = "gtf",
        url = paste0(base_url, file),
        file = file,
        md5 = md5sum,
        row.names = "gtf"
    )
}


#' GENCODE files downloaded by earlier versions of rmskProfiler
#'
#' These are no longer used but are listed by \code{listResources()} and
#' removed by \code{removeResources(downloads = TRUE)}.
#'
#' @param species Either "Hs" or "Mm"
#'
#' @return Named character vector of file names
#' @keywords internal
.legacyDownloads <- function(species) {
    if (species == "Hs") {
        c(
            transcripts = "gencode.v48.transcripts.fa.gz",
            genome = "GRCh38.primary_assembly.genome.fa.gz"
        )
    } else {
        c(
            transcripts = "gencode.vM25.transcripts.fa.gz",
            genome = "GRCm38.primary_assembly.genome.fa.gz"
        )
    }
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


#' Files generated for each combination of species and exclude/min_len settings
#' @keywords internal
.BUILD_FILES <- c(
    "rmsk.bed",
    "rmsk-duplicateInfo.rds",
    "rmsk-unique.fa.gz",
    "rmsk-rowData.rds",
    "rmsk-gentrome.fa.gz",
    "decoys.txt"
)


#' Name of the BiocFileCache metadata table holding build information
#' @keywords internal
.META_TABLE <- "rmskProfiler"


#' Name of a cached resource
#'
#' @param species Either "Hs" or "Mm"
#' @param resource "gtf" for the GENCODE GTF, "annotation.txdb", or one of the
#' files in \code{.BUILD_FILES}
#' @param exclude,min_len Settings used to build the resource. Required for
#' files in \code{.BUILD_FILES}, ignored otherwise.
#'
#' @return character(1)
#' @keywords internal
.resourceRname <- function(species, resource, exclude = NULL, min_len = NULL) {
    gencode <- .gencodeResources(species)
    if (resource %in% gencode$type) {
        return(.rname(species, gencode[resource, "file"]))
    }
    if (resource %in% .BUILD_FILES) {
        stopifnot(!is.null(exclude), !is.null(min_len))
        return(.rname(species, resource, .settingsKey(exclude, min_len)))
    }
    .rname(species, resource)
}


#' Look up the path of a resource in the cache
#'
#' @param bfc BiocFileCache
#' @param rname Name of the resource
#'
#' @return Named path to the cached file (named by rid), or NULL if the
#' resource is not in the cache or its file does not exist
#' @keywords internal
.findResource <- function(bfc, rname) {
    hit <- BiocFileCache::bfcquery(bfc, rname, field = "rname", exact = TRUE)
    if (nrow(hit) == 0L) {
        return(NULL)
    }
    path <- BiocFileCache::bfcrpath(bfc, rids = hit$rid[1L])
    if (!file.exists(path)) {
        return(NULL)
    }
    path
}


#' Get the path of an existing resource in the cache
#'
#' @inheritParams .resourceRname
#' @param bfc BiocFileCache
#' @param hint Message appended to the error if the resource is missing
#'
#' @return Path to the cached file
#' @keywords internal
.getResource <- function(
    bfc,
    species,
    resource,
    exclude = NULL,
    min_len = NULL,
    hint = ""
) {
    rname <- .resourceRname(species, resource, exclude, min_len)
    path <- .findResource(bfc, rname)
    if (is.null(path)) {
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
    unname(path)
}


#' Get a path in the cache to write a generated resource to
#'
#' If the resource already exists in the cache its path is returned so that it
#' is overwritten. Otherwise a new cache entry is created. In both cases the
#' build information for the resource is recorded in the cache metadata.
#'
#' @inheritParams .resourceRname
#' @param bfc BiocFileCache
#' @param ext File extension of the resource, e.g. ".bed"
#' @param genome Label of the genome the resource was built from, or NULL
#'
#' @return Path to the cached file
#' @keywords internal
.newResource <- function(
    bfc,
    species,
    resource,
    ext,
    exclude = NULL,
    min_len = NULL,
    genome = NULL
) {
    rname <- .resourceRname(species, resource, exclude, min_len)
    hit <- BiocFileCache::bfcquery(bfc, rname, field = "rname", exact = TRUE)
    if (nrow(hit) > 0L) {
        path <- BiocFileCache::bfcrpath(bfc, rids = hit$rid[1L])
    } else {
        path <- BiocFileCache::bfcnew(bfc, rname, ext = ext)
    }
    .recordResource(bfc, names(path), species, exclude, min_len, genome)

    unname(path)
}


#' Get the label of the genome a cached resource was built from
#'
#' @inheritParams .resourceRname
#' @param bfc BiocFileCache
#'
#' @return character(1), NA if no genome was recorded for the resource
#' @keywords internal
.resourceGenome <- function(
    bfc,
    species,
    resource,
    exclude = NULL,
    min_len = NULL
) {
    rname <- .resourceRname(species, resource, exclude, min_len)
    hit <- BiocFileCache::bfcquery(bfc, rname, field = "rname", exact = TRUE)
    meta <- .readMeta(bfc)
    meta$genome[match(hit$rid[1L], meta$rid)]
}


#' Read the rmskProfiler metadata table from the cache
#'
#' Rows for resources that are no longer in the cache are dropped.
#'
#' @param bfc BiocFileCache
#'
#' @return data.frame with columns rid, species, build, exclude, min_len, and
#' genome
#' @keywords internal
.readMeta <- function(bfc) {
    empty <- data.frame(
        rid = character(),
        species = character(),
        build = character(),
        exclude = character(),
        min_len = integer(),
        genome = character()
    )
    if (!.META_TABLE %in% BiocFileCache::bfcmetalist(bfc)) {
        return(empty)
    }
    meta <- BiocFileCache::bfcmeta(bfc, .META_TABLE)
    # Caches created by earlier versions have no genome column
    if (!"genome" %in% names(meta)) {
        meta$genome <- rep(NA_character_, nrow(meta))
    }
    meta[meta$rid %in% BiocFileCache::bfcrid(bfc), names(empty)]
}


#' Record build information for a resource in the cache metadata
#'
#' Any existing record for the resource is replaced.
#'
#' @param bfc BiocFileCache
#' @param rid Resource id
#' @param species Either "Hs" or "Mm"
#' @param exclude,min_len Settings used to build the resource, or NULL for
#' resources that do not depend on them
#' @param genome Label of the genome the resource was built from, or NULL for
#' resources that do not depend on it
#'
#' @return NULL
#' @keywords internal
.recordResource <- function(
    bfc,
    rid,
    species,
    exclude = NULL,
    min_len = NULL,
    genome = NULL
) {
    has_settings <- !is.null(exclude) && !is.null(min_len)
    record <- data.frame(
        rid = rid,
        species = species,
        build = if (has_settings) {
            .settingsKey(exclude, min_len)
        } else {
            NA_character_
        },
        exclude = if (has_settings) {
            paste(sort(unique(exclude)), collapse = ",")
        } else {
            NA_character_
        },
        min_len = if (has_settings) as.integer(min_len) else NA_integer_,
        genome = if (is.null(genome)) NA_character_ else genome
    )
    meta <- .readMeta(bfc)
    .writeMeta(bfc, rbind(meta[meta$rid != rid, ], record))
}


#' Write the rmskProfiler metadata table to the cache
#'
#' @param bfc BiocFileCache
#' @param meta data.frame returned by \code{.readMeta()}
#'
#' @return NULL
#' @keywords internal
.writeMeta <- function(bfc, meta) {
    if (nrow(meta) > 0L) {
        BiocFileCache::bfcmeta(bfc, .META_TABLE, overwrite = TRUE) <- meta
    } else if (.META_TABLE %in% BiocFileCache::bfcmetalist(bfc)) {
        BiocFileCache::bfcmetaremove(bfc, .META_TABLE)
    }

    return(invisible(NULL))
}
