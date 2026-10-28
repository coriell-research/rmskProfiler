# One gene on the plus strand of chr1 (1-based coordinates):
#   promoter 8001-10200 (2000 bp upstream, 200 bp downstream of the TSS)
#   exon 1   10001-11000, 5' UTR 10001-10500, CDS 10501-11000
#   intron   11001-12000
#   exon 2   12001-13000, CDS 12001-12500, 3' UTR 12501-13000
toyGTF <- function() {
    attrs <- paste0(
        'gene_id "G1"; transcript_id "T1"; gene_name "GENE1"; ',
        'gene_type "protein_coding";'
    )
    records <- data.frame(
        feature = c("gene", "transcript", "exon", "exon", "CDS", "CDS"),
        start = c(10001, 10001, 10001, 12001, 10501, 12001),
        end = c(13000, 13000, 11000, 13000, 11000, 12500),
        frame = c(".", ".", ".", ".", "0", "1")
    )
    gtf <- tempfile(fileext = ".gtf")
    writeLines(
        paste(
            "chr1",
            "TEST",
            records$feature,
            records$start,
            records$end,
            ".",
            "+",
            records$frame,
            attrs,
            sep = "\t"
        ),
        gtf
    )
    gtf
}

# TE loci as saved by extractUniqueSeqs(): 0-based start, "." if unstranded
toyDuplicateInfo <- function() {
    data.table::data.table(
        Hash = c(
            "exon",
            "intron",
            "antisense",
            "utr3",
            "promoter",
            "intergenic",
            "multi",
            "multi",
            "unstranded"
        ),
        RepName = c(
            "AluY:Alu:SINE",
            "L1a:L1:LINE",
            "MER1:hAT:DNA",
            "AluY:Alu:SINE",
            "MIR:MIR:SINE",
            "L2:L2:LINE",
            "L1a:L1:LINE",
            "L1b:L1:LINE",
            "L2:L2:LINE"
        ),
        seqnames = "chr1",
        strand = c("+", "+", "-", "+", "+", "+", "+", "+", "."),
        start = c(10600, 11400, 11400, 12600, 9000, 50000, 11600, 60000, 11800),
        end = c(10700, 11500, 11500, 12700, 9100, 50100, 11700, 60100, 11900)
    )
}

# Build the toy annotation once and share it between tests
.toy <- new.env()

toyAnnotation <- function() {
    if (!is.null(.toy$rd)) {
        return(.toy$rd)
    }
    bfc <- BiocFileCache::BiocFileCache(tempfile("toy-cache"), ask = FALSE)
    added <- BiocFileCache::bfcadd(
        bfc,
        .resourceRname("Hs", "gtf"),
        fpath = toyGTF(),
        action = "copy"
    )
    .recordResource(bfc, names(added), "Hs")
    info <- .newResource(
        bfc,
        "Hs",
        "rmsk-duplicateInfo.rds",
        ".rds",
        .DEFAULT_EXCLUDE,
        32
    )
    saveRDS(toyDuplicateInfo(), info)
    suppressWarnings(suppressMessages(
        createAnnotation(species = "Hs", cache = bfc)
    ))
    .toy$rd <- readRDS(
        getResource("rmsk-rowData.rds", species = "Hs", cache = bfc)
    )
    .toy$rd
}

test_that("TE loci are annotated with the features they overlap", {
    rd <- toyAnnotation()
    te <- rd[!is.na(rd$Hash), ]
    flag <- function(column) stats::setNames(te[[column]], te$Hash)
    expected <- function(...) {
        hashes <- c(...)
        stats::setNames(te$Hash %in% hashes, te$Hash)
    }

    # Stranded flags require the TE and the feature to be on the same strand
    expect_identical(flag("hasExonic"), expected("exon", "utr3"))
    expect_identical(flag("hasIntronic"), expected("intron", "multi"))
    expect_identical(flag("has3UTR"), expected("utr3"))
    expect_identical(flag("has5UTR"), expected())
    expect_identical(flag("hasPromoter"), expected("promoter"))
    # Promoters do not count as genic
    expect_identical(
        flag("hasIntergenic"),
        expected("antisense", "promoter", "intergenic", "unstranded")
    )

    # Unstranded flags ignore the strand of the TE
    expect_identical(
        flag("hasUnstrandedIntronic"),
        expected("intron", "antisense", "multi", "unstranded")
    )
    expect_identical(flag("hasUnstrandedExonic"), expected("exon", "utr3"))
    expect_identical(
        flag("hasUnstrandedIntergenic"),
        expected("promoter", "intergenic")
    )

    features <- stats::setNames(
        strsplit(te$overlapping_features, ";", fixed = TRUE),
        te$Hash
    )
    expect_setequal(features$utr3, c("exon:G1", "3utr:G1"))
    expect_setequal(features$promoter, "promoter:G1")
    expect_length(features$intergenic, 0L)
})

test_that("loci with the same sequence are collapsed into one row", {
    rd <- toyAnnotation()
    te <- rd[!is.na(rd$Hash), ]
    expect_setequal(rownames(te), unique(toyDuplicateInfo()$Hash))

    multi <- te["multi", ]
    expect_identical(multi$N_Loci, 2L)
    expect_identical(multi$Subfamily, "L1a,L1b")
    expect_identical(multi$Family, "L1")
    expect_identical(multi$Class, "LINE")
    expect_identical(multi$N_Subfamily, 2L)
    # Ranges are 1-based and keep the strand of each locus
    expect_identical(
        as.character(multi$Ranges[[1]]),
        c("chr1:11601-11700:+", "chr1:60001-60100:+")
    )
    expect_identical(
        as.character(te["unstranded", ]$Ranges[[1]]),
        "chr1:11801-11900"
    )
    expect_identical(
        as.character(GenomicRanges::strand(te["unstranded", ]$Ranges[[1]])),
        "*"
    )

    expect_identical(te["exon", ]$N_Loci, 1L)
    expect_identical(te["exon", ]$Subfamily, "AluY")
    expect_identical(te["exon", ]$Family, "Alu")
    expect_identical(te["exon", ]$Class, "SINE")
})

test_that("transcripts are annotated with their gene", {
    rd <- toyAnnotation()
    tx <- rd["T1", ]
    expect_identical(tx$transcript_id, "T1")
    expect_identical(tx$gene_id, "G1")
    expect_identical(tx$gene_name, "GENE1")
    expect_identical(tx$gene_type, "protein_coding")
    expect_true(is.na(tx$Hash))
    expect_identical(as.character(tx$Ranges[[1]]), "chr1:10001-13000:+")
})
