writeToyInputs <- function(bed_lines) {
    genome_fa <- tempfile(fileext = ".fa")
    # 1-4 AAAA, 5-8 CCCC, 9-12 GGGG, 13-16 TTTT, 17-24 ACGTACGT, 25-32 AAAACCCC
    Biostrings::writeXStringSet(
        Biostrings::DNAStringSet(c(
            "chrA 1" = "AAAACCCCGGGGTTTTACGTACGTAAAACCCC"
        )),
        genome_fa
    )
    bed <- tempfile(fileext = ".bed")
    writeLines(bed_lines, bed)
    list(bed = bed, genome = genome_fa)
}

test_that("sequences are hashed, reverse complemented, and deduplicated", {
    files <- writeToyInputs(c(
        "chrA\t0\t8\tL1a:L1:LINE\t0\t+",
        "chrA\t8\t16\tAluY:Alu:SINE\t0\t-",
        "chrA\t16\t24\tMER1:hAT:DNA\t0\t+",
        "chrA\t24\t32\tL1a:L1:LINE\t0\t+"
    ))
    res <- suppressMessages(.extractUniqueRecords(files$bed, files$genome))
    sha1 <- function(x) digest::digest(x, algo = "sha1", serialize = FALSE)

    # GGGGTTTT on the minus strand is AAAACCCC, a duplicate of the first locus
    expect_identical(
        as.character(res$seqs),
        c("AAAACCCC", "ACGTACGT") |>
            stats::setNames(c(sha1("AAAACCCC"), sha1("ACGTACGT")))
    )
    expect_identical(
        res$info$Hash,
        c(
            sha1("AAAACCCC"),
            sha1("AAAACCCC"),
            sha1("ACGTACGT"),
            sha1("AAAACCCC")
        )
    )
    expect_named(
        res$info,
        c("Hash", "RepName", "seqnames", "strand", "start", "end")
    )
    expect_identical(res$info$start, c(0L, 8L, 16L, 24L))
    expect_identical(res$info$end, c(8L, 16L, 24L, 32L))
    expect_identical(res$info$strand, c("+", "-", "+", "+"))
})

test_that("records outside of the genome are skipped with a warning", {
    files <- writeToyInputs(c(
        "chrA\t0\t8\tL1a:L1:LINE\t0\t+",
        "chrA\t28\t40\tL1a:L1:LINE\t0\t+",
        "chrB\t0\t8\tL1a:L1:LINE\t0\t+"
    ))
    expect_warning(
        res <- suppressMessages(.extractUniqueRecords(files$bed, files$genome)),
        "2 records not contained in the genome were skipped"
    )
    expect_identical(nrow(res$info), 1L)
    expect_length(res$seqs, 1L)
})
