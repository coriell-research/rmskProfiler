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

# Use chrM (16.6 kb) as the only decoy so that createGentrome() writes and
# indexes a small gentrome instead of every primary chromosome of the genome
local_small_decoy <- function(env = parent.frame()) {
    local_mocked_bindings(.decoys = function(genome) "chrM", .env = env)
}
