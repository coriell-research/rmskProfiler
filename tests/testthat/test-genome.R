toyGenome <- function() {
    Biostrings::DNAStringSet(c(
        "chrA 1" = "AAAACCCCGGGGTTTTACGTACGTAAAACCCC",
        "chrB 2" = "ACGTACGTAC"
    ))
}

test_that(".getGenome() accepts a DNAStringSet, FASTA paths, and FaFile", {
    gr <- GenomicRanges::GRanges(c("chrA:9-16:-", "chrA:1-4:+", "chrB:3-6:+"))
    expected <- c("AAAACCCC", "AAAA", "GTAC")

    from_set <- .getGenome("Hs", toyGenome())
    expect_s3_class(from_set, "rmskGenome")
    expect_match(from_set$label, "^custom:[0-9a-f]{10}$")
    expect_identical(as.character(.genomeSeqs(from_set, gr)), expected)
    expect_identical(
        .genomeLengths(from_set),
        c(chrA = 32L, chrB = 10L)
    )

    gz <- tempfile(fileext = ".fa.gz")
    Biostrings::writeXStringSet(toyGenome(), gz, compress = TRUE)
    from_gz <- .getGenome("Hs", gz)
    # The same genome gets the same label whatever the input type
    expect_identical(from_gz$label, from_set$label)
    expect_identical(as.character(.genomeSeqs(from_gz, gr)), expected)

    # An indexed FASTA is used in place
    fa <- tempfile(fileext = ".fa")
    Biostrings::writeXStringSet(toyGenome(), fa)
    Rsamtools::indexFa(fa)
    from_fa <- .getGenome("Hs", fa)
    expect_identical(BiocGenerics::path(from_fa$seqs), fa)
    expect_identical(as.character(.genomeSeqs(from_fa, gr)), expected)

    # An unindexed FaFile is indexed
    fa2 <- tempfile(fileext = ".fa")
    Biostrings::writeXStringSet(
        stats::setNames(toyGenome(), c("chrA", "chrB")),
        fa2
    )
    from_fafile <- .getGenome("Hs", Rsamtools::FaFile(fa2))
    expect_true(file.exists(paste0(fa2, ".fai")))
    expect_identical(as.character(.genomeSeqs(from_fafile, gr)), expected)

    expect_identical(from_fa$label, from_set$label)
    expect_identical(from_fafile$label, from_set$label)
    expect_false(
        .getGenome("Hs", toyGenome()[1])$label == from_set$label
    )

    # A FASTA path is only copied and indexed once per session
    expect_identical(
        BiocGenerics::path(.getGenome("Hs", gz)$seqs),
        BiocGenerics::path(from_gz$seqs)
    )

    # A resolved genome is returned as is
    expect_identical(.getGenome("Hs", from_set), from_set)
})

test_that(".getGenome() errors on invalid input", {
    expect_error(.getGenome("Hs", tempfile()), "Genome FASTA file not found")
    expect_error(.getGenome("Hs", 1L), "genome must be NULL, a BSgenome")
})

test_that(".getGenome() asks for the BSgenome package if it is missing", {
    local_mocked_bindings(.isInstalled = function(pkg) FALSE)
    expect_error(
        .getGenome("Hs"),
        "BSgenome.Hsapiens.UCSC.hg38.*BiocManager::install"
    )
    expect_error(.getGenome("Mm"), "BSgenome.Mmusculus.UCSC.mm10")
})

test_that("GENCODE sequence names are matched to UCSC names", {
    target <- c(
        "chr1",
        "chr14_GL000009v2_random",
        "chr6_GL000250v2_alt",
        "chr1_KN196472v1_fix",
        "chrUn_KI270302v1",
        "chr1_GL455991_alt"
    )
    expect_identical(
        .matchSeqnames(
            c(
                "chr1",
                "GL000009.2",
                "GL000250.2",
                "KN196472.1",
                "KI270302.1",
                "GL455991.1",
                "GL000009.1",
                "chr2"
            ),
            target
        ),
        c(target, NA, NA)
    )
})

test_that("transcripts extracted from the genome match GENCODE", {
    txdb <- suppressWarnings(suppressMessages(txdbmaker::makeTxDbFromGFF(
        system.file(
            "extdata",
            "gencode-chr22.gtf.gz",
            package = "rmskProfiler"
        ),
        format = "gtf"
    )))
    gencode <- Biostrings::readDNAStringSet(test_path(
        "fixtures",
        "transcripts-chr22.fa.gz"
    ))
    names(gencode) <- sub("\\|.*", "", names(gencode))

    tx <- .transcriptSeqs(
        txdb,
        .getGenome("Hs", exampleGenome()),
        tx_order = names(gencode)
    )
    expect_identical(as.character(tx), as.character(gencode))
})

test_that("transcripts on sequences missing from the genome are skipped", {
    txdb <- suppressWarnings(suppressMessages(txdbmaker::makeTxDbFromGFF(
        system.file(
            "extdata",
            "gencode-chr22.gtf.gz",
            package = "rmskProfiler"
        ),
        format = "gtf"
    )))
    renamed <- stats::setNames(toyGenome(), c("chr1", "chr2"))
    expect_error(
        .transcriptSeqs(txdb, .getGenome("Hs", renamed)),
        "No transcripts found on the sequences of the genome"
    )
})

test_that("the example genome matches BSgenome.Hsapiens.UCSC.hg38", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
    region <- GenomicRanges::GRanges("chr22:37000000-37999999")
    bsgenome <- .getGenome("Hs")
    expect_match(bsgenome$label, "^BSgenome.Hsapiens.UCSC.hg38 ")
    expect_identical(
        .genomeSeqs(bsgenome, region),
        .genomeSeqs(.getGenome("Hs", exampleGenome()), region)
    )
})
