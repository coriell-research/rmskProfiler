test_that("settings keys do not depend on the order of exclude", {
    expect_identical(
        .settingsKey(c("a", "b"), 32),
        .settingsKey(c("b", "a", "a"), 32)
    )
    expect_false(.settingsKey("a", 32) == .settingsKey("a", 50))
})

test_that("resource names identify species and build settings", {
    expect_identical(
        .resourceRname("Hs", "genome"),
        "rmskProfiler:Hs:GRCh38.primary_assembly.genome.fa.gz"
    )
    expect_identical(
        .resourceRname("Mm", "rmsk.bed", "Satellite", 32),
        paste0("rmskProfiler:Mm:", .settingsKey("Satellite", 32), ":rmsk.bed")
    )
    expect_error(.resourceRname("Hs", "rmsk.bed"))
})

test_that("exampleCache() contains the example downloads", {
    bfc <- exampleCache()
    res <- listResources(cache = bfc)
    expect_setequal(res$resource, c("genome", "gtf", "transcripts"))
    expect_true(all(res$type == "download"))
    expect_true(file.exists(getResource("gtf", species = "Hs", cache = bfc)))
    expect_error(
        getResource("rmsk.bed", species = "Hs", cache = bfc),
        "not found in cache"
    )
})

test_that("rmskToBed() validates rmsk", {
    bfc <- exampleCache()
    expect_error(
        rmskToBed(species = "Hs", cache = bfc, rmsk = data.frame()),
        "rmsk must be a GRanges"
    )
    expect_error(
        rmskToBed(
            species = "Hs",
            cache = bfc,
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
    expect_setequal(
        res$resource,
        c("genome", "gtf", "transcripts", "annotation.txdb")
    )

    suppressMessages(
        removeResources(species = "Hs", downloads = TRUE, cache = bfc)
    )
    expect_identical(nrow(listResources(cache = bfc)), 0L)
})
