#' List resources stored in the rmskProfiler cache
#'
#' Lists all of the downloaded and generated resources that rmskProfiler has
#' stored in the cache along with the settings used to build them.
#'
#' @details
#' Resources are one of three types:
#'
#' \describe{
#'  \item{download}{The GENCODE GTF downloaded by \code{downloadResources()}.
#'  GENCODE genome and transcript sequences downloaded by earlier versions of
#'  rmskProfiler are also listed as downloads. They are no longer used and can
#'  be removed with \code{removeResources(downloads = TRUE)}.}
#'  \item{annotation}{The TxDb created from the GENCODE GTF}
#'  \item{build}{Files generated for a specific combination of \code{exclude}
#'  and \code{min_len} settings. The \code{build} column contains a short key
#'  identifying these settings. The \code{genome} column contains the genome
#'  the TE and gentrome sequences were extracted from.}
#' }
#'
#' The Salmon index is not stored in the cache and is not listed.
#'
#' @param species NULL (default) to list resources for all species, or either
#' "Hs" (Homo sapiens) or "Mm" (Mus musculus)
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#'
#' @return data.frame with columns species, type, build, exclude, min_len,
#' genome, resource, size_mb, created, and path. Values in the resource column
#' can be passed to \code{getResource()}, except for the "genome" and
#' "transcripts" downloads of earlier versions.
#' @export
#'
#' @examples
#' bfc <- exampleCache(build = TRUE)
#' listResources(cache = bfc)
#' listResources(species = "Hs", cache = bfc)
listResources <- function(species = NULL, cache = NULL) {
    if (!is.null(species)) {
        species <- match.arg(species, c("Hs", "Mm"))
    }
    bfc <- .getCache(cache)

    info <- as.data.frame(BiocFileCache::bfcinfo(bfc))
    info <- info[
        startsWith(info$rname, "rmskProfiler:"),
        c("rid", "rname", "rpath", "create_time")
    ]

    # Resource names are "rmskProfiler:<species>[:<build>]:<file>"
    parts <- strsplit(info$rname, ":", fixed = TRUE)
    info$species <- vapply(parts, `[`, character(1), 2L)
    info$build <- vapply(
        parts,
        function(x) if (length(x) == 4L) x[3L] else NA_character_,
        character(1)
    )
    file <- vapply(parts, function(x) x[length(x)], character(1))

    downloads <- unlist(lapply(c("Hs", "Mm"), function(sp) {
        gtf <- .gencodeResources(sp)
        c(stats::setNames(gtf$file, gtf$type), .legacyDownloads(sp))
    }))
    is_download <- file %in% downloads
    info$resource <- file
    info$resource[is_download] <- names(downloads)[match(
        file[is_download],
        downloads
    )]
    info$type <- ifelse(
        is_download,
        "download",
        ifelse(is.na(info$build), "annotation", "build")
    )

    meta <- .readMeta(bfc)
    info <- merge(
        info,
        meta[, c("rid", "exclude", "min_len", "genome")],
        by = "rid",
        all.x = TRUE
    )
    info$size_mb <- round(file.size(info$rpath) / 1e6, 1)

    if (!is.null(species)) {
        info <- info[info$species == species, ]
    }
    info <- info[
        order(
            info$species,
            match(info$type, c("download", "annotation", "build")),
            info$build,
            info$resource
        ),
    ]

    result <- data.frame(
        species = info$species,
        type = info$type,
        build = info$build,
        exclude = info$exclude,
        min_len = info$min_len,
        genome = info$genome,
        resource = info$resource,
        size_mb = info$size_mb,
        created = info$create_time,
        path = info$rpath
    )

    return(result)
}


#' Get the path to a resource in the rmskProfiler cache
#'
#' Returns the local path to a downloaded or generated resource so that it can
#' be used outside of rmskProfiler.
#'
#' @param resource Name of the resource. One of "gtf" for the GENCODE GTF,
#' "annotation.txdb" for the TxDb, or one of "rmsk.bed",
#' "rmsk-duplicateInfo.rds", "rmsk-unique.fa.gz", "rmsk-rowData.rds",
#' "rmsk-gentrome.fa.gz", or "decoys.txt" for files generated for a specific
#' build.
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus)
#' @param exclude A character vector specifying which elements (repClass) were
#' excluded when building the resource. Only used for build resources. Default
#' "Simple_repeat", "Low_complexity", "Satellite", "RNA", "rRNA", "snRNA",
#' "scRNA", "srpRNA", "tRNA", and "Unknown".
#' @param min_len Minimum sequence length of a record used when building the
#' resource. Only used for build resources. Default 32.
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#'
#' @return Path to the cached file
#' @export
#'
#' @examples
#' bfc <- exampleCache(build = TRUE)
#'
#' # Path to the GTF used to build the index
#' getResource("gtf", species = "Hs", cache = bfc)
#'
#' # Read in the TE and transcript annotations
#' rd <- readRDS(getResource("rmsk-rowData.rds", species = "Hs", cache = bfc))
#' head(rd)
getResource <- function(
    resource,
    species = c("Hs", "Mm"),
    exclude = .DEFAULT_EXCLUDE,
    min_len = 32,
    cache = NULL
) {
    if (resource %in% c("genome", "transcripts")) {
        stop(
            "The GENCODE ",
            resource,
            " sequences are no longer downloaded or used. The genome is ",
            "provided by the BSgenome package for the species and the ",
            "transcript sequences are extracted from it. Use ",
            "removeResources(downloads = TRUE) to remove files downloaded by ",
            "earlier versions.",
            call. = FALSE
        )
    }
    resource <- match.arg(
        resource,
        c("gtf", "annotation.txdb", .BUILD_FILES)
    )
    species <- match.arg(species)

    hint <- if (resource == "gtf") {
        "Run downloadResources() first."
    } else if (resource == "annotation.txdb") {
        "Run createAnnotation() first."
    } else {
        paste(
            "Check that exclude and min_len match the values used to build",
            "the index, or run generateIndex() first. Use listResources() to",
            "see the available builds."
        )
    }

    .getResource(.getCache(cache), species, resource, exclude, min_len, hint)
}


#' Remove resources from the rmskProfiler cache
#'
#' Removes all files generated for a build (a combination of species,
#' \code{exclude}, and \code{min_len} settings) from the cache. Optionally, the
#' GENCODE downloads and the TxDb created from them can also be removed,
#' including the GENCODE genome and transcript sequences downloaded by earlier
#' versions of rmskProfiler. The
#' Salmon index is not stored in the cache and is not removed.
#'
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus)
#' @param exclude A character vector specifying which elements (repClass) were
#' excluded when building the resources to remove. Default "Simple_repeat",
#' "Low_complexity", "Satellite", "RNA", "rRNA", "snRNA", "scRNA", "srpRNA",
#' "tRNA", and "Unknown".
#' @param min_len Minimum sequence length of a record used when building the
#' resources to remove. Default 32.
#' @param downloads TRUE/FALSE, also remove the GENCODE downloads and the TxDb
#' for the species? Default FALSE.
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#'
#' @return Invisibly, the names of the removed resources
#' @export
#'
#' @examples
#' bfc <- exampleCache(build = TRUE)
#'
#' # Remove the files generated for the default build
#' removeResources(species = "Hs", cache = bfc)
#'
#' # Remove everything, including downloads
#' removeResources(species = "Hs", downloads = TRUE, cache = bfc)
#' listResources(cache = bfc)
removeResources <- function(
    species = c("Hs", "Mm"),
    exclude = .DEFAULT_EXCLUDE,
    min_len = 32,
    downloads = FALSE,
    cache = NULL
) {
    species <- match.arg(species)
    bfc <- .getCache(cache)

    resources <- .BUILD_FILES
    if (isTRUE(downloads)) {
        resources <- c(
            resources,
            .gencodeResources(species)$type,
            "annotation.txdb"
        )
    }
    rnames <- vapply(
        resources,
        .resourceRname,
        character(1),
        species = species,
        exclude = exclude,
        min_len = min_len,
        USE.NAMES = FALSE
    )
    if (isTRUE(downloads)) {
        rnames <- c(rnames, .rname(species, .legacyDownloads(species)))
    }

    info <- BiocFileCache::bfcinfo(bfc)
    hits <- info$rname %in% rnames
    if (!any(hits)) {
        message("No matching resources found in the cache.")
        return(invisible(character()))
    }

    removed <- info$rname[hits]
    BiocFileCache::bfcremove(bfc, info$rid[hits])
    .writeMeta(bfc, .readMeta(bfc))
    message("Removed ", length(removed), " resources from the cache.")

    return(invisible(removed))
}
