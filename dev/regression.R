# Regression test against validated rmskProfiler resources
#
# Rebuilds every build found in a reference rmskProfiler cache with the current
# code and checks that the results are identical to the reference:
#
# - the GTF is the GENCODE release the package downloads (md5)
# - rmsk-duplicateInfo.rds and rmsk-rowData.rds are identical()
# - the TxDb has the same transcripts, exons, CDS, and genes
# - decoys.txt has the same lines
# - the gentrome has the md5 recorded in GENTROME_MD5 below. The reference cache
#   does not keep the gentrome (it is removed when the index is built), so
#   createGentrome() is run without Salmon and the md5 of the uncompressed
#   gentrome is compared instead.
#
# This is not part of the package or its tests. It takes ~45 minutes per human
# build and needs:
# - the BSgenome version recorded in the reference (see listResources())
# - AnnotationHub access to AH111333 (hg38) / AH99012 (mm10). Set
#   ANNOTATION_HUB_LOCAL=TRUE to use the local AnnotationHub cache offline
# - testthat, to run createGentrome() without Salmon
#
# Run from the package root with:
#   Rscript dev/regression.R <reference_cache_dir> [<work_dir>]
# e.g.
#   Rscript dev/regression.R /mnt/data/references/GRCh38/rmskProfiler/resources
#
# The reference cache is only read. The builds are written to work_dir
# (default: a temporary directory). The script exits with status 1 if any check
# fails.
#
# To record the gentrome of a newly validated build, run the script against its
# cache and add the md5 it reports as "not recorded" to GENTROME_MD5.

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

# md5 of the uncompressed gentrome of validated builds
GENTROME_MD5 <- data.frame(
    species = "Hs",
    build = "e2d0ff3d9d",
    genome = "BSgenome.Hsapiens.UCSC.hg38 1.4.5",
    gtf = "gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz",
    md5 = "e1e39e90f8b783328d188e918961f1ea"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || !dir.exists(args[1])) {
    stop("Usage: Rscript dev/regression.R <reference_cache_dir> [<work_dir>]")
}
reference_dir <- normalizePath(args[1])
work_dir <- if (length(args) > 1L) {
    args[2]
} else {
    tempfile("rmskProfiler-regression")
}
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

# Results -----------------------------------------------------------------------
results <- data.frame(
    check = character(),
    status = character(),
    detail = character()
)
record <- function(check, status, detail = "") {
    message(sprintf("[%s] %s %s", status, check, detail))
    results[nrow(results) + 1L, ] <<- list(check, status, detail)
}
compare <- function(check, new, reference) {
    if (identical(new, reference)) {
        record(check, "PASS")
    } else {
        detail <- paste(
            utils::head(all.equal(new, reference), 5),
            collapse = "; "
        )
        record(check, "FAIL", detail)
    }
}

# Reference resources -----------------------------------------------------------
# BiocFileCache needs to write a lock file, so open a copy of the reference
# database. Paths are rebased onto the reference directory
db_dir <- tempfile("reference-db")
dir.create(db_dir)
file.copy(file.path(reference_dir, "BiocFileCache.sqlite"), db_dir)
reference <- listResources(cache = .getCache(db_dir))
reference$path <- file.path(reference_dir, basename(reference$path))
message("Reference resources in ", reference_dir, ":")
print(reference[, c("species", "type", "build", "genome", "resource")])

refPath <- function(species, resource, build = NA) {
    hit <- reference$species == species & reference$resource == resource
    if (!is.na(build)) {
        hit <- hit & reference$build %in% build
    }
    reference$path[hit][1]
}

# Rebuild and compare -----------------------------------------------------------
for (species in unique(reference$species[reference$resource == "gtf"])) {
    gencode <- .gencodeResources(species)
    gtf <- refPath(species, "gtf")
    if (tools::md5sum(gtf) == gencode$md5) {
        record(paste(species, "gtf"), "PASS")
    } else {
        record(
            paste(species, "gtf"),
            "FAIL",
            paste("not", gencode$file, "which the package downloads")
        )
    }

    bfc <- .getCache(file.path(work_dir, species))
    if (is.null(.findResource(bfc, .resourceRname(species, "gtf")))) {
        added <- BiocFileCache::bfcadd(
            bfc,
            .resourceRname(species, "gtf"),
            fpath = gtf,
            action = "copy"
        )
        .recordResource(bfc, names(added), species)
    }
    genome_label <- .genomeLabel(.getGenome(species))

    builds <- unique(reference[
        reference$species == species & reference$type == "build",
        c("build", "exclude", "min_len", "genome")
    ])
    for (i in seq_len(nrow(builds))) {
        build <- builds$build[i]
        exclude <- strsplit(builds$exclude[i], ",", fixed = TRUE)[[1]]
        min_len <- builds$min_len[i]
        label <- paste(species, build)
        message(
            "\nBuild ",
            label,
            " (min_len ",
            min_len,
            ", exclude ",
            builds$exclude[i],
            ")"
        )

        if (.settingsKey(exclude, min_len) != build) {
            record(
                label,
                "FAIL",
                "exclude and min_len do not give the build key"
            )
            next
        }
        if (!identical(builds$genome[i], genome_label)) {
            record(
                label,
                "FAIL",
                paste0(
                    "built from ",
                    builds$genome[i],
                    " but ",
                    genome_label,
                    " is installed"
                )
            )
            next
        }

        tryCatch(
            {
                extractUniqueSeqs(
                    species = species,
                    exclude = exclude,
                    min_len = min_len,
                    cache = bfc
                )
                compare(
                    paste(label, "rmsk-duplicateInfo.rds"),
                    readRDS(getResource(
                        "rmsk-duplicateInfo.rds",
                        species,
                        exclude,
                        min_len,
                        cache = bfc
                    )),
                    readRDS(refPath(species, "rmsk-duplicateInfo.rds", build))
                )

                createAnnotation(
                    species = species,
                    exclude = exclude,
                    min_len = min_len,
                    cache = bfc
                )
                compare(
                    paste(label, "rmsk-rowData.rds"),
                    readRDS(getResource(
                        "rmsk-rowData.rds",
                        species,
                        exclude,
                        min_len,
                        cache = bfc
                    )),
                    readRDS(refPath(species, "rmsk-rowData.rds", build))
                )

                # Run createGentrome() without Salmon and keep the md5 of the
                # gentrome it would index
                gentrome_md5 <- NA_character_
                testthat::with_mocked_bindings(
                    createGentrome(
                        species = species,
                        index_dir = file.path(work_dir, "index"),
                        exclude = exclude,
                        min_len = min_len,
                        cache = bfc
                    ),
                    .checkSalmon = function(...) invisible(NULL),
                    .salmonIndex = function(gentrome_fa, ...) {
                        gentrome_md5 <<- unname(tools::md5sum(gentrome_fa))
                        0L
                    },
                    .package = "rmskProfiler"
                )
                compare(
                    paste(label, "decoys.txt"),
                    readLines(getResource(
                        "decoys.txt",
                        species,
                        exclude,
                        min_len,
                        cache = bfc
                    )),
                    readLines(refPath(species, "decoys.txt", build))
                )
                expected <- GENTROME_MD5$md5[
                    GENTROME_MD5$species == species &
                        GENTROME_MD5$build == build &
                        GENTROME_MD5$genome == genome_label &
                        GENTROME_MD5$gtf == gencode$file
                ]
                if (length(expected) == 0L) {
                    record(
                        paste(label, "gentrome"),
                        "UNCHECKED",
                        paste(
                            "md5",
                            gentrome_md5,
                            "not recorded, add it to GENTROME_MD5 once validated"
                        )
                    )
                } else if (identical(gentrome_md5, expected)) {
                    record(paste(label, "gentrome"), "PASS")
                } else {
                    record(
                        paste(label, "gentrome"),
                        "FAIL",
                        paste("md5", gentrome_md5, "expected", expected)
                    )
                }
            },
            error = function(e) {
                record(label, "FAIL", conditionMessage(e))
            }
        )
    }

    # The TxDb is created by the first createAnnotation() of the species
    txdb_file <- .findResource(bfc, .resourceRname(species, "annotation.txdb"))
    if (!is.null(txdb_file)) {
        txdb <- suppressPackageStartupMessages(AnnotationDbi::loadDb(txdb_file))
        txdb_ref <- AnnotationDbi::loadDb(refPath(species, "annotation.txdb"))
        for (f in c("transcripts", "exons", "cds", "genes")) {
            extract <- getExportedValue("GenomicFeatures", f)
            compare(
                paste(species, "annotation.txdb", f),
                as.data.frame(suppressMessages(extract(txdb))),
                as.data.frame(suppressMessages(extract(txdb_ref)))
            )
        }
    }
}

# Summary -----------------------------------------------------------------------
message("\nSummary:")
print(results, right = FALSE)
if (any(results$status == "FAIL")) {
    quit(status = 1L)
}
