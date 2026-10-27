#' Download files needed for index generation
#'
#' This function will download the GENCODE annotation GTF needed for generating
#' the rmsk-gentrome index into the rmskProfiler cache: GENCODE v48
#' (chr_patch_hapl_scaff) for human and GENCODE M25 for mouse. Files that
#' already exist in the cache are not downloaded again. The genome and
#' transcript sequences are not downloaded. They are taken from the BSgenome
#' package for the species, see \code{generateIndex()}.
#'
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus) designating
#' which
#' species to download resources for
#' @param check_integrity TRUE/FALSE, if TRUE check the md5sum of the GENCODE
#' GTF
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#'
#' @return Invisibly \code{NULL}. Called for its side effect of downloading
#' the GENCODE GTF into the cache.
#' @export
#' @examples
#' # The example cache already contains an (example) GENCODE GTF, so nothing
#' # is downloaded. With the default cache, downloadResources(species = "Hs")
#' # downloads the GENCODE v48 GTF.
#' bfc <- exampleCache()
#' downloadResources(species = "Hs", check_integrity = FALSE, cache = bfc)
downloadResources <- function(
    species = c("Hs", "Mm"),
    check_integrity = TRUE,
    cache = NULL
) {
    species <- match.arg(species)
    bfc <- .getCache(cache)
    res <- .gencodeResources(species)
    rnames <- vapply(
        res$type,
        .resourceRname,
        character(1),
        species = species,
        USE.NAMES = FALSE
    )

    paths <- character(nrow(res))
    for (i in seq_len(nrow(res))) {
        path <- .findResource(bfc, rnames[i])
        if (!is.null(path)) {
            message(res$file[i], " already exists in the cache. Skipping.")
        } else {
            # Remove any entry whose file has gone missing before downloading
            # again
            hit <- BiocFileCache::bfcquery(
                bfc,
                rnames[i],
                field = "rname",
                exact = TRUE
            )
            if (nrow(hit) > 0L) {
                BiocFileCache::bfcremove(bfc, hit$rid)
            }
            message("Attempting to download ", res$file[i], "...")
            path <- BiocFileCache::bfcadd(
                bfc,
                rnames[i],
                fpath = res$url[i],
                rtype = "web"
            )
        }
        .recordResource(bfc, names(path), species)
        paths[i] <- unname(path)
    }

    if (isTRUE(check_integrity)) {
        message("Checking file integrity of downloaded files...")
        badfile <- res$md5 != unname(tools::md5sum(paths))
        if (any(badfile)) {
            for (rname in rnames[badfile]) {
                hit <- BiocFileCache::bfcquery(
                    bfc,
                    rname,
                    field = "rname",
                    exact = TRUE
                )
                BiocFileCache::bfcremove(bfc, hit$rid)
            }
            stop(
                paste(res$file[badfile], collapse = ", "),
                " did not download properly and were removed from the cache. ",
                "Please retry."
            )
        }
        message("Success!")
    }

    return(invisible(NULL))
}
