# Run code with an empty PATH so that Salmon cannot be found
withoutSalmon <- function(code) {
    old <- Sys.getenv("PATH")
    on.exit(Sys.setenv(PATH = old))
    Sys.setenv(PATH = "")
    code
}

test_that("functions that run Salmon fail before doing any work without it", {
    bfc <- exampleCache()
    before <- listResources(cache = bfc)

    expect_error(
        withoutSalmon(generateIndex(tempfile(), species = "Hs", cache = bfc)),
        "Salmon was not found on the PATH.*create_index = FALSE"
    )
    expect_error(
        withoutSalmon(createGentrome("Hs", tempfile(), cache = bfc)),
        "Salmon was not found on the PATH.*create_index = FALSE"
    )
    expect_error(
        withoutSalmon(salmonQuant(
            "r1.fq",
            "r2.fq",
            "s1",
            tempdir(),
            tempdir()
        )),
        "Salmon was not found on the PATH"
    )
    expect_identical(listResources(cache = bfc)$resource, before$resource)
})

test_that("Salmon is not required when the index is not created", {
    local_small_decoy()
    bfc <- exampleCache()
    rmsk <- readRDS(system.file(
        "extdata",
        "rmsk-chr22.rds",
        package = "rmskProfiler"
    ))
    expect_no_error(withoutSalmon(suppressMessages(generateIndex(
        species = "Hs",
        check_integrity = FALSE,
        create_index = FALSE,
        cache = bfc,
        rmsk = rmsk
    ))))
    expect_true("decoys.txt" %in% listResources(cache = bfc)$resource)
})

test_that("index_dir is required to create the index", {
    expect_error(
        generateIndex(species = "Hs", cache = exampleCache()),
        "index_dir must be provided"
    )
    expect_error(
        createGentrome(species = "Hs", cache = exampleCache()),
        "index_dir must be provided"
    )
})
