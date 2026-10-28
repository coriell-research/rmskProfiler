# RepeatMasker records as returned by filterRmsk() from BED-style lines
bedToRecords <- function(bed_lines) {
    bed <- data.table::fread(
        text = bed_lines,
        header = FALSE,
        col.names = c("chrom", "start", "end", "name", "score", "strand")
    )
    gr <- GenomicRanges::GRanges(
        bed$chrom,
        IRanges::IRanges(bed$start + 1L, bed$end),
        strand = ifelse(bed$strand == ".", "*", bed$strand)
    )
    gr$name <- bed$name
    gr
}

toyInputs <- function(bed_lines) {
    # 1-4 AAAA, 5-8 CCCC, 9-12 GGGG, 13-16 TTTT, 17-24 ACGTACGT, 25-32 AAAACCCC
    genome <- Biostrings::DNAStringSet(c(
        chrA = "AAAACCCCGGGGTTTTACGTACGTAAAACCCC"
    ))
    list(bed = bedToRecords(bed_lines), genome = genome)
}

test_that("sequences are hashed, reverse complemented, and deduplicated", {
    files <- toyInputs(c(
        "chrA\t0\t8\tL1a:L1:LINE\t0\t+",
        "chrA\t8\t16\tAluY:Alu:SINE\t0\t-",
        "chrA\t16\t24\tMER1:hAT:DNA\t0\t+",
        "chrA\t24\t32\tL1a:L1:LINE\t0\t+"
    ))
    res <- suppressMessages(.extractUniqueRecords(
        files$bed,
        files$genome
    ))
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
    files <- toyInputs(c(
        "chrA\t0\t8\tL1a:L1:LINE\t0\t+",
        "chrA\t28\t40\tL1a:L1:LINE\t0\t+",
        "chrB\t0\t8\tL1a:L1:LINE\t0\t+"
    ))
    expect_warning(
        res <- suppressMessages(.extractUniqueRecords(
            files$bed,
            files$genome
        )),
        "2 records not contained in the genome were skipped"
    )
    expect_identical(nrow(res$info), 1L)
    expect_length(res$seqs, 1L)
})

test_that("unstranded records are extracted from the plus strand", {
    files <- toyInputs(c(
        "chrA\t0\t8\tL1a:L1:LINE\t0\t.",
        "chrA\t8\t16\tAluY:Alu:SINE\t0\t-"
    ))
    res <- suppressMessages(
        .extractUniqueRecords(files$bed, files$genome)
    )
    expect_identical(unname(as.character(res$seqs)), "AAAACCCC")
    expect_identical(res$info$strand, c(".", "-"))
})

test_that("unstranded records work with the BSgenome", {
    bed <- bedToRecords(c(
        "chr22\t37000000\t37000010\tL1a:L1:LINE\t0\t.",
        "chr22\t37000000\t37000010\tL1a:L1:LINE\t0\t+",
        "chr22\t37000020\t37000030\tL1a:L1:LINE\t0\t-"
    ))
    res <- suppressMessages(.extractUniqueRecords(bed, .getGenome("Hs")))
    expect_identical(res$info$Hash[1], res$info$Hash[2])
    expect_length(res$seqs, 2L)
})

test_that(".uniqueSeqs() rebuilds the unique sequences from the loci", {
    files <- toyInputs(c(
        "chrA\t0\t8\tL1a:L1:LINE\t0\t+",
        "chrA\t8\t16\tAluY:Alu:SINE\t0\t-",
        "chrA\t16\t24\tMER1:hAT:DNA\t0\t.",
        "chrA\t24\t32\tL1a:L1:LINE\t0\t+"
    ))
    genome <- files$genome
    res <- suppressMessages(.extractUniqueRecords(files$bed, genome))
    expect_identical(
        as.character(.uniqueSeqs(res$info, genome)),
        as.character(res$seqs)
    )
})

test_that("filterRmsk() filters and names the RepeatMasker records", {
    rmsk <- readRDS(system.file(
        "extdata",
        "rmsk-chr22.rds",
        package = "rmskProfiler"
    ))
    gr <- filterRmsk(species = "Hs", rmsk = rmsk)
    expect_s4_class(gr, "GRanges")
    expect_false(any(gr$repClass %in% .DEFAULT_EXCLUDE))
    expect_true(all(GenomicRanges::width(gr) >= 32))
    expect_identical(
        gr$name,
        paste(gr$repName, gr$repFamily, gr$repClass, sep = ":")
    )
    expect_identical(
        filterRmsk(species = "Hs", min_len = 100, rmsk = gr),
        gr[GenomicRanges::width(gr) >= 100]
    )
})
