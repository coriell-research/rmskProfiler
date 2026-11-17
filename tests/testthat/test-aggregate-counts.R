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
