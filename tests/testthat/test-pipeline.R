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
    bfc <- builtCache()
    suppressMessages(createGentrome(
        species = "Hs",
        create_index = FALSE,
        cache = bfc,
        genome = exampleGenome()
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
        res$genome[res$type == "build"] ==
            .getGenome("Hs", exampleGenome())$label
    ))
})

test_that("createGentrome() builds the index from a temporary gentrome", {
    skip_if(!nzchar(Sys.which("salmon")), "Salmon is not available")
    bfc <- suppressMessages(exampleCache(build = TRUE))
    # A gentrome cached by an earlier build is removed
    suppressMessages(createGentrome(
        species = "Hs",
        create_index = FALSE,
        cache = bfc,
        genome = exampleGenome()
    ))
    parent <- tempfile("index-parent")
    index_dir <- file.path(parent, "index")
    suppressMessages(createGentrome(
        species = "Hs",
        index_dir = index_dir,
        cache = bfc,
        genome = exampleGenome()
    ))
    expect_true(file.exists(file.path(index_dir, "info.json")))
    # The gentrome is removed and not saved to the cache
    expect_identical(list.files(parent), "index")
    res <- listResources(cache = bfc)$resource
    expect_true("decoys.txt" %in% res)
    expect_false("rmsk-gentrome.fa.gz" %in% res)
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
    bfc <- builtCache()
    genome <- Biostrings::readDNAStringSet(exampleGenome())

    # The same genome given as a DNAStringSet does not warn
    expect_no_warning(suppressMessages(createGentrome(
        species = "Hs",
        create_index = FALSE,
        cache = bfc,
        genome = genome
    )))

    other <- c(genome, Biostrings::DNAStringSet(c(chrZ = "ACGT")))
    expect_warning(
        suppressMessages(createGentrome(
            species = "Hs",
            create_index = FALSE,
            cache = bfc,
            genome = other
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
