# Read the Gibbs samples of one sample as a targets x samples matrix
readGibbs <- function(path) {
    meta <- jsonlite::fromJSON(file.path(path, "aux_info", "meta_info.json"))
    n_targets <- if (is.null(meta$num_targets)) {
        meta$num_valid_targets
    } else {
        meta$num_targets
    }
    con <- gzcon(file(
        file.path(path, "aux_info", "bootstrap", "bootstraps.gz"),
        open = "rb"
    ))
    on.exit(close(con))
    boot <- readBin(
        con,
        what = "double",
        n = n_targets * meta$num_bootstraps
    )
    dim(boot) <- c(n_targets, meta$num_bootstraps)
    rownames(boot) <- utils::read.delim(file.path(path, "quant.sf"))$Name
    boot
}

quantPaths <- function() {
    list.dirs(exampleQuants(), full.names = TRUE, recursive = FALSE)
}

test_that("importQuants() imports effective lengths", {
    se <- exampleSE()
    quant <- utils::read.delim(
        file.path(exampleQuants(), "ctrl_1", "quant.sf"),
        row.names = 1
    )
    expect_equal(
        SummarizedExperiment::assay(se, "length")[, "ctrl_1"],
        stats::setNames(quant[rownames(se), "EffectiveLength"], rownames(se))
    )
})

test_that("importQuants() estimates overdispersion of aggregated features", {
    se <- exampleSE()
    od <- S4Vectors::metadata(se)$overdispersion
    expect_named(od, c("locus", "subfamily", "family", "class"))
    rd <- SummarizedExperiment::rowData(se)
    for (level in names(od)) {
        expect_true(all(.featureId(rd, level) %in% names(od[[level]])))
        expect_true(all(od[[level]] >= 1))
        expect_true(is.finite(attr(od[[level]], "prior")))
    }
})

test_that(".catchSalmon2() group overdispersion matches target overdispersion", {
    paths <- quantPaths()
    ids <- utils::read.delim(file.path(paths[1], "quant.sf"))$Name
    catch <- .catchSalmon2(
        paths,
        verbose = FALSE,
        groups = list(identity = ids)
    )
    expect_equal(
        unname(catch$group.overdispersion$identity[ids]),
        catch$annotation[ids, "Overdispersion"]
    )
    expect_equal(
        attr(catch$group.overdispersion$identity, "prior"),
        catch$overdispersion.prior
    )
})

test_that("gene overdispersion is estimated from summed Gibbs samples", {
    se <- exampleSE()
    rd <- SummarizedExperiment::rowData(se)
    od <- S4Vectors::metadata(se)$overdispersion$subfamily

    # A gene with several transcripts, all with counts
    tx_per_gene <- table(rd$gene_id)
    gene <- names(tx_per_gene)[tx_per_gene > 2L][1]
    tx <- rownames(rd)[which(rd$gene_id == gene)]

    # Every transcript of the gene in the index (zero-count ones included)
    rd_all <- readRDS(getResource(
        "rmsk-rowData.rds",
        species = "Hs",
        cache = builtCache()
    ))
    tx_all <- rownames(rd_all)[which(rd_all$gene_id == gene)]

    od_sum <- 0
    df <- 0L
    for (path in quantPaths()) {
        boot <- readGibbs(path)
        summed <- colSums(boot[
            intersect(tx_all, rownames(boot)),
            ,
            drop = FALSE
        ])
        m <- mean(summed)
        if (m > 0) {
            od_sum <- od_sum + sum((summed - m)^2) / m
            df <- df + length(summed) - 1L
        }
    }
    prior <- attr(od, "prior")
    expected <- max((3 * prior + df * (od_sum / df)) / (3 + df), 1)
    expect_equal(unname(od[gene]), expected)
})

test_that("aggregateCounts() returns tximport-style gene lengths", {
    se <- exampleSE()
    rd <- SummarizedExperiment::rowData(se)
    tpms <- SummarizedExperiment::assay(se, "tpms")
    eff_length <- SummarizedExperiment::assay(se, "length")
    agg <- suppressMessages(
        aggregateCounts(se, "Hs", te_level = "subfamily", cache = builtCache())
    )
    agg_length <- SummarizedExperiment::assay(agg, "length")

    tx_per_gene <- table(rd$gene_id)
    gene <- names(tx_per_gene)[tx_per_gene > 1L][1]
    in_gene <- which(rd$gene_id == gene)
    expect_equal(
        agg_length[gene, ],
        colSums(tpms[in_gene, ] * eff_length[in_gene, ]) /
            colSums(tpms[in_gene, ])
    )
    expect_equal(
        SummarizedExperiment::assay(agg, "abundance")[gene, ],
        colSums(tpms[in_gene, ])
    )

    # TE lengths are the same in every sample: the summed effective lengths
    # averaged over samples, so that they are not corrected for length
    is_te <- !is.na(rd$Hash)
    key <- paste(rd$Subfamily, rd$Family, rd$Class, sep = ":")
    subfamily <- names(which(table(key[is_te]) > 1L))[1]
    in_subfamily <- which(is_te & key == subfamily)
    expect_equal(
        unname(agg_length[subfamily, ]),
        rep(mean(colSums(eff_length[in_subfamily, ])), ncol(se))
    )
    expect_identical(
        SummarizedExperiment::rowData(agg)[subfamily, "n_loci"],
        sum(rd$N_Loci[in_subfamily])
    )
    expect_equal(
        SummarizedExperiment::rowData(agg)[subfamily, "Overdispersion"],
        unname(S4Vectors::metadata(se)$overdispersion$subfamily[subfamily])
    )
    expect_true(all(is.na(SummarizedExperiment::rowData(agg)[gene, "n_loci"])))
})

test_that("aggregateCounts() fills lengths of genes without abundance", {
    se <- exampleSE()
    rd <- SummarizedExperiment::rowData(se)
    tx_per_gene <- table(rd$gene_id)
    genes <- names(tx_per_gene)[tx_per_gene > 1L][1:2]
    in_one <- which(rd$gene_id == genes[1])
    in_all <- which(rd$gene_id == genes[2])

    tpms <- SummarizedExperiment::assay(se, "tpms")
    tpms[in_one, 1] <- 0
    tpms[in_all, ] <- 0
    SummarizedExperiment::assay(se, "tpms") <- tpms
    eff_length <- SummarizedExperiment::assay(se, "length")

    agg <- suppressMessages(
        aggregateCounts(se, "Hs", te_level = "class", cache = builtCache())
    )
    agg_length <- SummarizedExperiment::assay(agg, "length")

    # Geometric mean of the other samples
    expect_equal(
        unname(agg_length[genes[1], 1]),
        exp(mean(log(agg_length[genes[1], -1])))
    )
    # Mean length of the transcripts in all samples
    expect_equal(
        unname(agg_length[genes[2], ]),
        rep(mean(rowMeans(eff_length[in_all, ])), ncol(se))
    )
})

test_that("aggregateCounts() generates gene counts from abundance", {
    se <- exampleSE()
    no <- suppressMessages(
        aggregateCounts(se, "Hs", te_level = "family", cache = builtCache())
    )
    is_gene <- rownames(no) %in% SummarizedExperiment::rowData(se)$gene_id
    raw <- SummarizedExperiment::assay(no, "counts")
    abundance <- SummarizedExperiment::assay(no, "abundance")
    gene_length <- SummarizedExperiment::assay(no, "length")

    for (method in c("scaledTPM", "lengthScaledTPM")) {
        agg <- suppressMessages(aggregateCounts(
            se,
            "Hs",
            te_level = "family",
            counts_from_abundance = method,
            cache = builtCache()
        ))
        counts <- SummarizedExperiment::assay(agg, "counts")
        expect_identical(S4Vectors::metadata(agg)$countsFromAbundance, method)

        # Library size of the genes is kept and TEs are unchanged
        expect_equal(colSums(counts[is_gene, ]), colSums(raw[is_gene, ]))
        expect_identical(counts[!is_gene, ], raw[!is_gene, ])

        expected <- abundance[is_gene, ]
        if (method == "lengthScaledTPM") {
            expected <- expected * rowMeans(gene_length[is_gene, ])
        }
        expected <- t(t(expected) * colSums(raw[is_gene, ]) / colSums(expected))
        expect_equal(counts[is_gene, ], expected)
    }
})

test_that("aggregateCounts() checks for lengths and overdispersions", {
    se <- exampleSE()
    no_length <- se
    SummarizedExperiment::assay(no_length, "length") <- NULL
    expect_error(
        aggregateCounts(no_length, "Hs", cache = builtCache()),
        "Re-run importQuants"
    )

    no_od <- se
    S4Vectors::metadata(no_od)$overdispersion <- NULL
    expect_message(
        agg <- aggregateCounts(no_od, "Hs", cache = builtCache()),
        "Overdispersions of the aggregated features are not available"
    )
    expect_false(
        "Overdispersion" %in% colnames(SummarizedExperiment::rowData(agg))
    )
})

test_that("gene overdispersion is smaller than that of its transcripts", {
    # Ambiguity between isoforms cancels when they are summed
    se <- exampleSE()
    rd <- SummarizedExperiment::rowData(se)
    od <- S4Vectors::metadata(se)$overdispersion$subfamily
    tx_od <- S4Vectors::metadata(se)$annotation[rownames(se), "Overdispersion"]
    w <- rowSums(SummarizedExperiment::assay(se, "orig"))
    is_tx <- !is.na(rd$gene_id)
    tx_per_gene <- table(rd$gene_id)
    genes <- names(tx_per_gene)[tx_per_gene > 1L]
    member_od <- vapply(
        genes,
        function(g) {
            i <- which(is_tx & rd$gene_id == g)
            sum(tx_od[i] * w[i]) / sum(w[i])
        },
        numeric(1)
    )
    expect_lt(stats::median(od[genes]), stats::median(member_od))
})

test_that("aggregateCounts() output works with edgeR::DGEListFromTximeta()", {
    skip_if_not_installed("edgeR", "4.10.0")
    se <- exampleSE()
    agg <- suppressMessages(
        aggregateCounts(se, "Hs", te_level = "subfamily", cache = builtCache())
    )
    y <- edgeR::DGEListFromTximeta(agg, divide = TRUE)
    expect_identical(dim(y$offset.prior), dim(agg))
    expect_equal(
        unname(y$counts),
        unname(
            SummarizedExperiment::assay(agg, "counts") /
                SummarizedExperiment::rowData(agg)$Overdispersion
        )
    )

    # Genes get the tximport length offset, TEs are simple sums with no offset
    is_te <- !rownames(agg) %in% SummarizedExperiment::rowData(se)$gene_id
    expect_true(all(y$offset.prior[is_te, ] == 0))
    expect_true(any(y$offset.prior[!is_te, ] != 0))
    expect_equal(
        unname(y$counts[is_te, ]),
        unname(
            rowsum(
                SummarizedExperiment::assay(se, "orig"),
                .featureId(SummarizedExperiment::rowData(se), "subfamily")
            )[rownames(agg)[is_te], ] /
                SummarizedExperiment::rowData(agg)$Overdispersion[is_te]
        )
    )

    # The list interface gives the same DGEList
    txi <- list(
        abundance = SummarizedExperiment::assay(agg, "abundance"),
        counts = SummarizedExperiment::assay(agg, "counts"),
        length = SummarizedExperiment::assay(agg, "length"),
        countsFromAbundance = S4Vectors::metadata(agg)$countsFromAbundance
    )
    y2 <- edgeR::DGEListFromTximport(
        txi,
        genes = as.data.frame(SummarizedExperiment::rowData(agg)),
        divide = TRUE
    )
    expect_identical(y2$counts, y$counts)
    expect_identical(y2$offset.prior, y$offset.prior)
})

test_that("TE counts are not scaled with counts from abundance", {
    skip_if_not_installed("edgeR", "4.10.0")
    se <- exampleSE()
    no <- suppressMessages(
        aggregateCounts(se, "Hs", te_level = "subfamily", cache = builtCache())
    )
    agg <- suppressMessages(aggregateCounts(
        se,
        "Hs",
        te_level = "subfamily",
        counts_from_abundance = "lengthScaledTPM",
        cache = builtCache()
    ))
    is_te <- !rownames(agg) %in% SummarizedExperiment::rowData(se)$gene_id
    y <- edgeR::DGEListFromTximeta(agg, divide = TRUE)
    expect_null(y$offset.prior)
    expect_identical(
        SummarizedExperiment::assay(agg, "counts")[is_te, ],
        SummarizedExperiment::assay(no, "counts")[is_te, ]
    )
    expect_identical(
        SummarizedExperiment::assay(agg, "length")[is_te, ],
        SummarizedExperiment::assay(no, "length")[is_te, ]
    )
})

# Reference results from tximport and edgeR ------------------------------------
# tximport and edgeR are treated as the truth. They are run on the Salmon output
# of the example quants, with a tx2gene table that maps every target that
# importQuants() keeps to its aggregated feature.

# quant.sf files of the example quants, named by sample
quantFiles <- function() {
    paths <- quantPaths()
    stats::setNames(file.path(paths, "quant.sf"), basename(paths))
}

# tx2gene table mapping targets to the features of te_level. With genes_only,
# TE-loci are left out so that tximport returns genes only
referenceTx2Gene <- function(te_level, genes_only = FALSE) {
    ids <- utils::read.delim(quantFiles()[1])$Name
    rd <- readRDS(getResource(
        "rmsk-rowData.rds",
        species = "Hs",
        cache = builtCache()
    ))[ids, ]
    feature <- .featureId(rd, te_level)
    keep <- !ids %in% .nonStandardTranscripts(rd)
    if (genes_only) {
        keep <- keep & !is.na(rd$gene_id)
    }
    data.frame(tx = ids[keep], gene = feature[keep])
}

test_that("aggregateCounts() gene assays match tximport", {
    skip_if_not_installed("tximport")
    se <- exampleSE()
    tx2gene <- referenceTx2Gene("subfamily", genes_only = TRUE)

    for (method in c("no", "scaledTPM", "lengthScaledTPM")) {
        txi <- suppressMessages(tximport::tximport(
            quantFiles(),
            type = "salmon",
            tx2gene = tx2gene,
            countsFromAbundance = method
        ))
        agg <- suppressMessages(aggregateCounts(
            se,
            "Hs",
            te_level = "subfamily",
            counts_from_abundance = method,
            cache = builtCache()
        ))
        is_gene <- rownames(agg) %in% tx2gene$gene
        genes <- rownames(agg)[is_gene]

        # Genes left out by importQuants() have no counts
        expect_true(all(genes %in% rownames(txi$counts)))
        dropped <- setdiff(rownames(txi$counts), genes)
        expect_true(all(txi$counts[dropped, ] == 0))

        for (a in c("counts", "abundance", "length")) {
            expect_equal(
                SummarizedExperiment::assay(agg, a)[genes, ],
                txi[[a]][genes, ],
                label = paste0(a, " (", method, ")")
            )
        }
        expect_identical(
            S4Vectors::metadata(agg)$countsFromAbundance,
            txi$countsFromAbundance
        )
    }
})

test_that("divided counts and offsets match edgeR::DGEListFromTximport()", {
    skip_if_not_installed("tximport")
    skip_if_not_installed("edgeR", "4.10.0")
    se <- exampleSE()

    for (level in .TE_LEVELS) {
        # tximport sums the Gibbs samples over each feature and edgeR estimates
        # the overdispersion of the feature from them
        txi <- suppressMessages(tximport::tximport(
            quantFiles(),
            type = "salmon",
            tx2gene = referenceTx2Gene(level),
            dropInfReps = FALSE
        ))
        expected <- edgeR::DGEListFromTximport(txi, divide = TRUE)

        agg <- suppressMessages(
            aggregateCounts(se, "Hs", te_level = level, cache = builtCache())
        )
        y <- edgeR::DGEListFromTximeta(agg, divide = TRUE)
        features <- rownames(agg)
        is_te <- !features %in% SummarizedExperiment::rowData(se)$gene_id
        expect_true(any(is_te) && any(!is_te))
        expect_true(all(features %in% rownames(expected)))

        expect_equal(
            y$genes$Overdispersion,
            expected$genes[features, "Overdispersion"],
            label = paste("Overdispersion of", level)
        )
        expect_equal(
            y$counts[is_te, ],
            expected$counts[features[is_te], ],
            label = paste("divided TE counts of", level)
        )
        expect_equal(
            y$counts[!is_te, ],
            expected$counts[features[!is_te], ],
            label = paste("divided gene counts of", level)
        )
        expect_equal(
            y$offset.prior[!is_te, ],
            expected$offset.prior[features[!is_te], ],
            label = paste("gene offsets of", level)
        )
    }
})

test_that("edgeR::DGEListFromTximeta() gives a DGEList ready for DE testing", {
    skip_if_not_installed("edgeR", "4.10.0")
    se <- exampleSE()
    agg <- suppressMessages(
        aggregateCounts(se, "Hs", te_level = "subfamily", cache = builtCache())
    )
    rd <- SummarizedExperiment::rowData(agg)
    is_te <- !rownames(agg) %in% SummarizedExperiment::rowData(se)$gene_id
    group <- factor(sub("_[0-9]+$", "", colnames(agg)))

    # edgeR only re-estimates overdispersions from inferential replicates, so
    # the 'Overdispersion' column must be used as is
    expect_false(any(grepl("infRep", SummarizedExperiment::assayNames(agg))))

    y <- edgeR::DGEListFromTximeta(agg, group = group, divide = TRUE)
    expect_s4_class(y, "DGEList")
    expect_no_error(edgeR:::validDGEList(y))

    # Counts are the summed counts divided by the overdispersions
    expected <- SummarizedExperiment::assay(agg, "counts") / rd$Overdispersion
    expect_identical(dimnames(y$counts), dimnames(agg))
    expect_equal(y$counts, expected)
    expect_true(all(is.finite(y$counts) & y$counts >= 0))
    expect_true(y$divided.counts)
    expect_identical(y$tximport.counts, "raw")

    # Samples: library sizes of the divided counts and the given groups
    expect_identical(rownames(y$samples), colnames(agg))
    expect_equal(y$samples$lib.size, unname(colSums(expected)))
    expect_true(all(y$samples$norm.factors == 1))
    expect_identical(y$samples$group, group)

    # Gene annotation: the rowData of agg plus edgeR's length summaries
    expect_identical(rownames(y$genes), rownames(agg))
    expect_true(all(colnames(rd) %in% colnames(y$genes)))
    expect_true(all(c("AveLength", "Max2MinLength") %in% colnames(y$genes)))
    expect_equal(y$genes$Overdispersion, rd$Overdispersion)
    expect_true(all(is.finite(y$genes$Overdispersion)))
    expect_true(all(y$genes$Overdispersion >= 1))
    expect_true(all(y$genes$Max2MinLength[is_te] == 1))

    # The offset: centred log gene lengths, zero for TEs
    expect_identical(dimnames(y$offset.prior), dimnames(agg))
    expect_true(all(is.finite(y$offset.prior)))
    expect_equal(unname(rowMeans(y$offset.prior)), rep(0, nrow(agg)))
    expect_true(all(y$offset.prior[is_te, ] == 0))
    expect_true(any(y$offset.prior[!is_te, ] != 0))

    # Filtering keeps the offset aligned with the counts
    keep <- edgeR::filterByExpr(y)
    expect_true(any(keep))
    y <- y[keep, , keep.lib.sizes = FALSE]
    expect_identical(dimnames(y$offset.prior), dimnames(y$counts))
    expect_equal(y$samples$lib.size, unname(colSums(y$counts)))

    # The offset is combined with the library sizes and normalization factors
    y <- edgeR::normLibSizes(y)
    expect_false(all(y$samples$norm.factors == 1))
    offset <- edgeR::getOffset(y)
    expect_equal(
        offset,
        t(t(y$offset.prior) + log(y$samples$lib.size * y$samples$norm.factors))
    )

    # Downstream functions run and use the offset
    design <- stats::model.matrix(~group)
    y <- edgeR::estimateDisp(y, design)
    expect_true(is.finite(y$common.dispersion))
    fit <- edgeR::glmQLFit(y, design)
    expect_equal(unname(as.matrix(fit$offset)), unname(offset))
    res <- edgeR::glmQLFTest(fit, coef = 2)
    expect_identical(rownames(res$table), rownames(y))
    expect_true(all(res$table$PValue >= 0 & res$table$PValue <= 1))
    expect_s4_class(edgeR::topTags(res), "TopTags")
    # CPMs use the offset as the effective library size of each count
    expect_equal(edgeR::cpm(y), y$counts / exp(offset) * 1e6)
    expect_s4_class(edgeR::voomLmFit(y, design), "MArrayLM")
})
