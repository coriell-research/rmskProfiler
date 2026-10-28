test_that("the example build contains annotations for transcripts and TEs", {
    bfc <- builtCache()
    rd <- readRDS(getResource("rmsk-rowData.rds", species = "Hs", cache = bfc))
    expect_s4_class(rd, "DataFrame")
    expect_true(all(
        c(
            "transcript_id",
            "gene_id",
            "Hash",
            "Subfamily",
            "N_Loci",
            "hasExonic",
            "hasIntronic",
            "hasIntergenic",
            "Ranges"
        ) %in%
            colnames(rd)
    ))
    expect_true(any(!is.na(rd$transcript_id)))
    expect_true(any(!is.na(rd$Hash)))
    expect_identical(rownames(rd)[!is.na(rd$Hash)], rd$Hash[!is.na(rd$Hash)])

    # Every unique sequence is annotated
    info <- readRDS(
        getResource("rmsk-duplicateInfo.rds", species = "Hs", cache = bfc)
    )
    expect_setequal(unique(info$Hash), rd$Hash[!is.na(rd$Hash)])
    expect_false(any(
        c("rmsk.bed", "rmsk-unique.fa.gz") %in%
            listResources(cache = bfc)$resource
    ))
})

test_that("createGentrome() writes the gentrome and decoys", {
    local_chr22_decoy()
    bfc <- builtCache()
    suppressMessages(createGentrome(
        species = "Hs",
        create_index = FALSE,
        cache = bfc
    ))
    decoys <- readLines(getResource("decoys.txt", species = "Hs", cache = bfc))
    expect_identical(decoys, "chr22")

    # Transcripts in GTF order, then the unique TE sequences, then the decoys
    gentrome <- Biostrings::readDNAStringSet(
        getResource("rmsk-gentrome.fa.gz", species = "Hs", cache = bfc)
    )
    tx <- Biostrings::readDNAStringSet(test_path(
        "fixtures",
        "transcripts-chr22.fa.gz"
    ))
    info <- readRDS(
        getResource("rmsk-duplicateInfo.rds", species = "Hs", cache = bfc)
    )
    expect_identical(
        names(gentrome),
        c(sub("\\|.*", "", names(tx)), unique(info$Hash), "chr22")
    )
    expect_identical(
        as.character(gentrome[seq_along(tx)]),
        stats::setNames(as.character(tx), sub("\\|.*", "", names(tx)))
    )
    # Every build file records the genome it was created from
    res <- listResources(cache = bfc)
    expect_true(all(
        res$genome[res$type == "build"] == .genomeLabel(.getGenome("Hs"))
    ))
})

test_that("createGentrome() gives Salmon the same gentrome uncompressed", {
    local_chr22_decoy()
    bfc <- suppressMessages(exampleCache(build = TRUE))
    suppressMessages(createGentrome(
        species = "Hs",
        create_index = FALSE,
        cache = bfc
    ))
    cached <- Biostrings::readDNAStringSet(
        getResource("rmsk-gentrome.fa.gz", species = "Hs", cache = bfc)
    )

    # Record what Salmon is given instead of running it
    salmon <- new.env()
    local_mocked_bindings(
        .checkSalmon = function(...) invisible(NULL),
        .salmonIndex = function(gentrome_fa, decoy_file, index_dir, threads) {
            salmon$gentrome_fa <- gentrome_fa
            salmon$gzipped <- identical(
                readBin(gentrome_fa, "raw", 2L),
                as.raw(c(0x1f, 0x8b))
            )
            salmon$gentrome <- Biostrings::readDNAStringSet(gentrome_fa)
            salmon$decoys <- readLines(decoy_file)
            salmon$index_dir <- index_dir
            salmon$threads <- threads
            0L
        }
    )
    parent <- tempfile("index-parent")
    index_dir <- file.path(parent, "index")
    suppressMessages(createGentrome(
        species = "Hs",
        index_dir = index_dir,
        threads = 3,
        cache = bfc
    ))

    expect_false(salmon$gzipped)
    expect_identical(salmon$gentrome, cached)
    expect_identical(salmon$decoys, "chr22")
    expect_identical(salmon$index_dir, index_dir)
    expect_identical(salmon$threads, 3)
    # The gentrome is written next to the index and removed afterwards
    expect_identical(
        dirname(salmon$gentrome_fa),
        normalizePath(parent)
    )
    expect_false(file.exists(salmon$gentrome_fa))
    # The cached gentrome no longer matches the index and is removed
    res <- listResources(cache = bfc)$resource
    expect_false("rmsk-gentrome.fa.gz" %in% res)
    expect_true("decoys.txt" %in% res)

    local_mocked_bindings(.salmonIndex = function(...) 1L)
    expect_error(
        suppressMessages(createGentrome(
            species = "Hs",
            index_dir = index_dir,
            cache = bfc
        )),
        "Salmon index generation failed with exit status 1"
    )
})

test_that("the temporary gentrome is written next to the index if possible", {
    parent <- tempfile("index-parent")
    expect_identical(
        .gentromeDir(file.path(parent, "index")),
        normalizePath(parent)
    )

    skip_on_os("windows")
    skip_if(identical(Sys.info()[["user"]], "root"), "root can write anywhere")
    read_only <- tempfile("read-only")
    dir.create(read_only)
    Sys.chmod(read_only, "555")
    on.exit(Sys.chmod(read_only, "755"), add = TRUE)
    expect_identical(.gentromeDir(file.path(read_only, "index")), tempdir())
})

test_that("createGentrome() warns if the TEs come from a different genome", {
    local_chr22_decoy()
    bfc <- suppressMessages(exampleCache(build = TRUE))

    # The genome the TEs were extracted from does not warn
    expect_no_warning(suppressMessages(createGentrome(
        species = "Hs",
        create_index = FALSE,
        cache = bfc
    )))

    # A different version of the BSgenome warns
    rid <- names(.findResource(
        bfc,
        .resourceRname("Hs", "rmsk-duplicateInfo.rds", .DEFAULT_EXCLUDE, 32)
    ))
    .recordResource(
        bfc,
        rid,
        "Hs",
        .DEFAULT_EXCLUDE,
        32,
        genome = "BSgenome.Hsapiens.UCSC.hg38 0.0.0"
    )
    expect_warning(
        suppressMessages(createGentrome(
            species = "Hs",
            create_index = FALSE,
            cache = bfc
        )),
        "different genome"
    )
})

test_that("importQuants() returns a SummarizedExperiment", {
    se <- exampleSE()
    expect_s4_class(se, "SummarizedExperiment")
    expect_identical(
        SummarizedExperiment::assayNames(se),
        c("counts", "orig", "tpms")
    )
    expect_identical(ncol(se), 6L)
    expect_setequal(
        colnames(se),
        basename(list.dirs(
            exampleQuants(),
            recursive = FALSE
        ))
    )
    expect_true(all(rowSums(SummarizedExperiment::assay(se, "counts")) > 0))
})

test_that("importQuants() requires Gibbs samples", {
    quant_dir <- tempfile()
    dir.create(quant_dir)
    file.copy(
        list.dirs(exampleQuants(), recursive = FALSE)[1],
        quant_dir,
        recursive = TRUE
    )
    meta_file <- list.files(
        quant_dir,
        "meta_info.json",
        recursive = TRUE,
        full.names = TRUE
    )
    meta <- jsonlite::fromJSON(meta_file)
    meta$num_bootstraps <- 0L
    jsonlite::write_json(meta, meta_file, auto_unbox = TRUE)

    expect_error(
        suppressMessages(
            importQuants(quant_dir, species = "Hs", cache = builtCache())
        ),
        "fewer than 2 Gibbs samples"
    )
})

test_that("aggregateCounts() aggregates to genes and TE levels", {
    se <- exampleSE()
    rd <- SummarizedExperiment::rowData(se)
    is_te <- !is.na(rd$Hash)
    for (level in c("locus", "subfamily", "family", "class")) {
        agg <- suppressMessages(
            aggregateCounts(se, "Hs", te_level = level, cache = builtCache())
        )
        expect_s4_class(agg, "SummarizedExperiment")
        expect_identical(ncol(agg), ncol(se))
        expect_equal(
            sum(SummarizedExperiment::assay(agg, "counts")),
            sum(SummarizedExperiment::assay(se, "counts"))
        )
        expect_true(all(
            c("Freq", "feature_length") %in%
                colnames(SummarizedExperiment::rowData(agg))
        ))
    }
    n_sub <- length(unique(paste(rd$Subfamily, rd$Family, rd$Class)[is_te]))
    agg <- suppressMessages(
        aggregateCounts(se, "Hs", te_level = "subfamily", cache = builtCache())
    )
    n_genes <- length(unique(rd$gene_id[!is_te]))
    expect_identical(nrow(agg), n_sub + n_genes)
})

test_that("estimateReadThrough() fits the example data", {
    se <- exampleSE()
    group <- factor(sub("_.*", "", colnames(se)))
    rt <- suppressMessages(estimateReadThrough(se, group = group))
    expect_s3_class(rt, "data.frame")
    expect_gt(nrow(rt), 0L)
    expect_true(all(
        c(
            "Hash",
            "host_id",
            "efficiency",
            "rt_fraction",
            "padj_readthrough",
            "padj_autonomous"
        ) %in%
            colnames(rt)
    ))
    expect_true(all(rt$Hash %in% rownames(se)))
    expect_identical(ncol(attr(rt, "readthrough")), ncol(se))
    expect_identical(ncol(attr(rt, "design")), 2L)
})

test_that("importQuants() reads the Salmon counts and scales them", {
    se <- exampleSE()
    quant <- utils::read.delim(
        file.path(exampleQuants(), "ctrl_1", "quant.sf"),
        row.names = 1
    )
    expect_equal(
        SummarizedExperiment::assay(se, "orig")[, "ctrl_1"],
        stats::setNames(quant[rownames(se), "NumReads"], rownames(se))
    )
    expect_equal(
        SummarizedExperiment::assay(se, "tpms")[, "ctrl_1"],
        stats::setNames(quant[rownames(se), "TPM"], rownames(se))
    )
    # Counts are divided by the overdispersion of each target
    overdispersion <- S4Vectors::metadata(se)$annotation[
        rownames(se),
        "Overdispersion"
    ]
    expect_equal(
        SummarizedExperiment::assay(se, "counts"),
        SummarizedExperiment::assay(se, "orig") / overdispersion
    )

    rd <- readRDS(getResource(
        "rmsk-rowData.rds",
        species = "Hs",
        cache = builtCache()
    ))
    expect_identical(SummarizedExperiment::rowData(se), rd[rownames(se), ])

    # Without removing all zero rows, every target is kept
    all_rows <- suppressMessages(importQuants(
        exampleQuants(),
        species = "Hs",
        remove_zeros = FALSE,
        cache = builtCache()
    ))
    expect_setequal(rownames(all_rows), rownames(quant))
})

test_that("aggregateCounts() sums counts and lengths of each feature", {
    se <- exampleSE()
    rd <- SummarizedExperiment::rowData(se)
    counts <- SummarizedExperiment::assay(se, "counts")
    agg <- suppressMessages(
        aggregateCounts(se, "Hs", te_level = "subfamily", cache = builtCache())
    )
    agg_counts <- SummarizedExperiment::assay(agg, "counts")
    agg_rd <- SummarizedExperiment::rowData(agg)

    # A gene is the sum of its transcripts, and its length the width of its
    # reduced exons
    tx_per_gene <- table(rd$gene_id)
    gene <- names(tx_per_gene)[tx_per_gene > 1L][1]
    in_gene <- which(rd$gene_id == gene)
    expect_equal(agg_counts[gene, ], colSums(counts[in_gene, ]))
    expect_identical(agg_rd[gene, "Freq"], length(in_gene))
    expect_identical(agg_rd[gene, "gene_name"], unique(rd$gene_name[in_gene]))
    gtf <- rtracklayer::import(system.file(
        "extdata",
        "gencode-chr22.gtf.gz",
        package = "rmskProfiler"
    ))
    exons <- gtf[gtf$type == "exon" & gtf$gene_id == gene]
    expect_equal(
        agg_rd[gene, "feature_length"],
        sum(GenomicRanges::width(GenomicRanges::reduce(exons)))
    )

    # A subfamily is the sum of its unique sequences, and its length the total
    # width of their loci
    is_te <- !is.na(rd$Hash)
    key <- paste(rd$Subfamily, rd$Family, rd$Class, sep = ":")
    hashes_per_key <- table(key[is_te])
    subfamily <- names(hashes_per_key)[hashes_per_key > 1L][1]
    in_subfamily <- which(is_te & key == subfamily)
    expect_equal(agg_counts[subfamily, ], colSums(counts[in_subfamily, ]))
    expect_identical(agg_rd[subfamily, "Freq"], length(in_subfamily))
    expect_equal(
        agg_rd[subfamily, "feature_length"],
        sum(GenomicRanges::width(unlist(rd$Ranges[in_subfamily])))
    )

    # At the locus level the TEs are not aggregated
    by_locus <- suppressMessages(
        aggregateCounts(se, "Hs", te_level = "locus", cache = builtCache())
    )
    te_hashes <- rownames(se)[is_te]
    expect_equal(
        SummarizedExperiment::assay(by_locus, "counts")[te_hashes, ],
        counts[te_hashes, ]
    )
})
