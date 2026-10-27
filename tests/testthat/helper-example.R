# Build the example cache once and share it between tests
.example <- new.env()

builtCache <- function() {
    if (is.null(.example$bfc)) {
        .example$bfc <- suppressMessages(exampleCache(build = TRUE))
    }
    .example$bfc
}

exampleSE <- function() {
    if (is.null(.example$se)) {
        .example$se <- suppressMessages(
            importQuants(exampleQuants(), species = "Hs", cache = builtCache())
        )
    }
    .example$se
}
