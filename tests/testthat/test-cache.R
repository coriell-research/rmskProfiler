test_that("settings keys do not depend on the order of exclude", {
    expect_identical(
        .settingsKey(c("a", "b"), 32),
        .settingsKey(c("b", "a", "a"), 32)
    )
    expect_false(.settingsKey("a", 32) == .settingsKey("a", 50))
})

test_that("resource names identify species and build settings", {
    expect_identical(
        .resourceRname("Hs", "gtf"),
        "rmskProfiler:Hs:gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz"
    )
    expect_identical(
        .resourceRname("Mm", "decoys.txt", "Satellite", 32),
        paste0("rmskProfiler:Mm:", .settingsKey("Satellite", 32), ":decoys.txt")
    )
    expect_error(.resourceRname("Hs", "decoys.txt"))
})

test_that("exampleCache() contains the example GTF", {
    bfc <- exampleCache()
    res <- listResources(cache = bfc)
    expect_identical(res$resource, "gtf")
    expect_true(all(res$type == "download"))
    expect_true(file.exists(getResource("gtf", species = "Hs", cache = bfc)))
    expect_error(
        getResource("decoys.txt", species = "Hs", cache = bfc),
        "not found in cache"
    )
})

test_that("rmskToBed() validates rmsk", {
    expect_error(
        rmskToBed(species = "Hs", rmsk = data.frame()),
        "rmsk must be a GRanges"
    )
    expect_error(
        rmskToBed(
            species = "Hs",
            rmsk = GenomicRanges::GRanges("chr22:1-100")
        ),
        "repName, repFamily, and repClass"
    )
})

test_that("removeResources() removes build files and optionally downloads", {
    bfc <- suppressMessages(exampleCache(build = TRUE))
    expect_true(
        "rmsk-rowData.rds" %in% listResources(cache = bfc)$resource
    )

    suppressMessages(removeResources(species = "Hs", cache = bfc))
    res <- listResources(cache = bfc)
    expect_setequal(res$resource, c("gtf", "annotation.txdb"))

    suppressMessages(
        removeResources(species = "Hs", downloads = TRUE, cache = bfc)
    )
    expect_identical(nrow(listResources(cache = bfc)), 0L)
})

test_that("genome and transcript downloads of earlier versions are removed", {
    bfc <- exampleCache()
    for (file in .legacyDownloads("Hs")) {
        BiocFileCache::bfcadd(
            bfc,
            .rname("Hs", file),
            fpath = exampleGenome(),
            action = "copy"
        )
    }
    res <- listResources(cache = bfc)
    expect_setequal(res$resource, c("gtf", "transcripts", "genome"))
    expect_true(all(res$type == "download"))

    expect_error(
        getResource("genome", species = "Hs", cache = bfc),
        "no longer downloaded.*removeResources\\(downloads = TRUE\\)"
    )

    suppressMessages(
        removeResources(species = "Hs", downloads = TRUE, cache = bfc)
    )
    expect_identical(nrow(listResources(cache = bfc)), 0L)
})

test_that("caches without a genome column in the metadata can be read", {
    bfc <- exampleCache()
    meta <- .readMeta(bfc)
    BiocFileCache::bfcmeta(bfc, .META_TABLE, overwrite = TRUE) <-
        meta[, setdiff(names(meta), "genome")]
    expect_identical(.readMeta(bfc)$genome, NA_character_)
    expect_true(is.na(listResources(cache = bfc)$genome))
})

test_that("build files of earlier versions are listed and removed", {
    bfc <- suppressMessages(exampleCache(build = TRUE))
    for (file in .LEGACY_BUILD_FILES) {
        BiocFileCache::bfcadd(
            bfc,
            .resourceRname("Hs", file, .DEFAULT_EXCLUDE, 32),
            fpath = exampleGenome(),
            action = "copy"
        )
    }
    res <- listResources(cache = bfc)
    expect_true(all(.LEGACY_BUILD_FILES %in% res$resource))
    expect_error(
        getResource("rmsk.bed", species = "Hs", cache = bfc),
        "no longer created or used"
    )

    suppressMessages(removeResources(species = "Hs", cache = bfc))
    expect_setequal(
        listResources(cache = bfc)$resource,
        c("gtf", "annotation.txdb")
    )
})
