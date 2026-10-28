test_that(".getGenome() returns the BSgenome of the species", {
    hs <- .getGenome("Hs")
    expect_s4_class(hs, "BSgenome")
    expect_match(
        .genomeLabel(hs),
        "^BSgenome\\.Hsapiens\\.UCSC\\.hg38 [0-9.]+$"
    )
    expect_match(
        .genomeLabel(.getGenome("Mm")),
        "^BSgenome\\.Mmusculus\\.UCSC\\.mm10 [0-9.]+$"
    )
})

test_that("the primary chromosomes are the decoys", {
    expect_identical(
        .decoys(.getGenome("Hs")),
        paste0("chr", c(1:22, "X", "Y"))
    )
    expect_identical(
        .decoys(.getGenome("Mm")),
        paste0("chr", c(1:19, "X", "Y"))
    )
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
        .getGenome("Hs"),
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
    toy <- Biostrings::DNAStringSet(c(chr1 = "AAAACCCC", chr2 = "ACGTACGT"))
    expect_error(
        .transcriptSeqs(txdb, toy),
        "No transcripts found on the sequences of the genome"
    )
})
