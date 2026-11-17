#' Aggregate counts to gene and TE level
#'
#' This function sums transcripts to genes, and TE-loci (unique sequences) to
#' the specified level. A new SummarizedExperiment object is returned with
#' gene-level and TE-level assays laid out like the output of
#' \code{tximeta::summarizeToGene()}, so that it can be passed directly to
#' \code{edgeR::DGEListFromTximeta()}.
#'
#' @details
#' The returned object has three assays:
#' \describe{
#'   \item{counts}{The estimated counts (the 'orig' assay of \code{x}) summed
#'   for each feature. These are \emph{not} divided by the overdispersion. For
#'   genes they can instead be generated from abundances, see
#'   \code{counts_from_abundance}.}
#'   \item{abundance}{The TPMs summed for each feature.}
#'   \item{length}{For genes, the average effective length of their
#'   transcripts weighted by their abundance in each sample, as in tximport.
#'   For genes with no abundance in a sample, the geometric mean of their
#'   lengths in the other samples is used, or, when they have no abundance in
#'   any sample, the mean length of their transcripts. For TEs, the sum of the
#'   effective lengths of the TE-loci in each feature, averaged over samples
#'   so that it is the same in every sample. TEs are not corrected for
#'   length, and the offset that edgeR forms from this assay is zero for
#'   them.}
#' }
#'
#' The rowData contains the number of features summed ('Freq'), the gene name
#' and type, the 'Overdispersion' of each feature, 'feature_length', and
#' 'n_loci'. 'feature_length' is, for genes, the sum of reduced exon widths
#' and, for TEs, the sum of the widths of all genomic loci of the TE-loci in
#' the feature. 'n_loci' is the number of genomic loci of the TEs in each
#' feature (NA for genes). Both describe only the TE-loci present in \code{x}.
#'
#' \strong{Overdispersion.} The overdispersion of each aggregated feature is
#' estimated by \code{importQuants()} from the Gibbs samples summed over the
#' transcripts or TE-loci in the feature. This is the overdispersion of the
#' summed counts themselves. Read-to-target ambiguity among the members of a
#' feature (e.g. isoforms of a gene, or near-identical copies of a TE
#' subfamily) cancels when they are summed, while ambiguity between features
#' (e.g. between related TE subfamilies, or between TE-loci and the
#' transcripts of genes) does not. The overdispersion of a gene is therefore
#' usually much smaller than that of its transcripts, while the overdispersion
#' of a TE feature can be smaller or larger than that of its TE-loci. For this
#' reason the counts are summed from the raw estimated counts, rather than
#' from the 'counts' assay of \code{x} (which is divided by the
#' overdispersion of each transcript or TE-locus). See the sections below for
#' using them in differential expression analyses. The overdispersions are
#' estimated for the features of the object returned by \code{importQuants()}
#' and are not updated if \code{x} is subset afterwards.
#'
#' \strong{Counts from abundance.} Summed counts are biased when differential
#' transcript usage changes the average length of a gene's transcripts
#' (Soneson et al. 2015). With \code{counts_from_abundance = "no"} (default)
#' the summed counts are returned and the 'length' assay can be used as an
#' offset (\code{edgeR::DGEListFromTximeta()} does this automatically). With
#' \code{"scaledTPM"}, gene counts are the gene TPMs scaled up to the library
#' size of the gene counts, and with \code{"lengthScaledTPM"}, they are the
#' gene TPMs multiplied by the average length of the gene over samples, then
#' scaled up to the library size of the gene counts. These methods are applied
#' to genes only. TEs are always the sum of their estimated counts: TE signal
#' is mostly fragments of longer transcripts, so counts (which measure
#' transcribed bases) are more meaningful than abundances (which count
#' fragments equally regardless of their length). TE counts are always
#' simple sums of the estimated counts, with no correction for length,
#' whatever \code{counts_from_abundance} is.
#'
#' Note that TE-loci overlapping exons, or in transcribed introns, may count
#' reads from their host genes. See \code{estimateReadThrough()}.
#'
#' @section Recommended workflow for differential expression:
#' The recommended workflow for genes and TE subfamilies follows the methods
#' of the tximport and edgeR developers:
#' \enumerate{
#'   \item Quantify with Gibbs sampling (\code{salmonQuant()}) and import with
#'   \code{importQuants()}.
#'   \item Aggregate with \code{aggregateCounts(te_level = "subfamily")},
#'   keeping the default \code{counts_from_abundance = "no"}.
#'   \item Create a \code{DGEList} with
#'   \code{edgeR::DGEListFromTximeta(se, divide = TRUE)} (edgeR >= 4.10). This
#'   divides the counts by the 'Overdispersion' column of the rowData
#'   (Baldoni et al. 2024) and adds the centred log of the 'length' assay as
#'   \code{offset.prior}. Together these are tximport's "original counts and
#'   offset" method for genes. The offset is zero for TEs, which stay simple
#'   sums. \code{edgeR::DGEListFromTximport()} gives the same result from a
#'   list of the assays (see the example below).
#'   \item Analyse as usual. The offset is combined with the library sizes and
#'   normalization factors by \code{glmQLFit()}, \code{estimateDisp()},
#'   \code{voomLmFit()}, and \code{cpm()}.
#' }
#' \preformatted{
#' y <- edgeR::DGEListFromTximeta(se, group = group, divide = TRUE)
#' keep <- edgeR::filterByExpr(y)
#' y <- y[keep, , keep.lib.sizes = FALSE]
#' y <- edgeR::normLibSizes(y)
#' fit <- edgeR::glmQLFit(y, design)
#' res <- edgeR::glmQLFTest(fit, coef = 2)
#'
#' # Equivalent, with the list interface of edgeR::DGEListFromTximport()
#' txi <- list(
#'     abundance = assay(se, "abundance"),
#'     counts = assay(se, "counts"),
#'     length = assay(se, "length"),
#'     countsFromAbundance = metadata(se)$countsFromAbundance
#' )
#' y <- edgeR::DGEListFromTximport(
#'     txi,
#'     genes = as.data.frame(rowData(se)),
#'     group = group,
#'     divide = TRUE
#' )
#' }
#' The aggregated counts should not be analysed without the offset (e.g. with
#' \code{edgeR::SE2DGEList()}). The tximport authors advise against summed
#' gene counts without an offset, because they do not correct for
#' differential transcript usage.
#'
#' The alternative tximport method, "bias corrected counts without an offset",
#' is \code{counts_from_abundance = "lengthScaledTPM"} (or
#' \code{"scaledTPM"}). \code{edgeR::DGEListFromTximeta()} then adds no offset
#' for any row. TE counts and lengths are the same as with \code{"no"}. For
#' 3' tagged RNA-seq, which has no length bias,
#' tximport recommends summed counts without an offset: remove it with
#' \code{y$offset.prior <- NULL}.
#'
#' @section Relation to tximport and edgeR:
#' For genes, the 'counts', 'abundance', and 'length' assays (including the
#' lengths of genes without abundance) and the \code{"scaledTPM"} and
#' \code{"lengthScaledTPM"} counts are computed as by
#' \code{tximport::summarizeToGene()}. The overdispersions are computed as by
#' \code{edgeR::DGEListFromTximport()} from inferential replicates summed to
#' the gene level by tximport. Three choices go beyond the tximport and edgeR
#' documentation, which does not deal with TEs:
#' \itemize{
#'   \item \emph{Dividing gene counts.} tximport and edgeR show divided counts
#'   for transcript-level analyses, and their gene-level examples use the
#'   default \code{divide = FALSE}. The rmskProfiler index also contains TE
#'   sequences, so reads can be ambiguous between a gene and TEs, and that
#'   ambiguity remains after summing. Dividing has no effect where the
#'   overdispersion is 1, so \code{divide = TRUE} is recommended for genes and
#'   TEs alike.
#'   \item \emph{TE lengths.} TE groups are not given tximport's abundance
#'   weighted length, which estimates a number of transcript molecules. TE
#'   reads mostly come from fragments of longer transcripts, so the summed
#'   counts of a TE group measure its transcribed bases. TE groups are
#'   therefore simple sums with no length correction, as transcripts and
#'   TE-loci are in edgeR's transcript-level workflow. A change in which loci
#'   of a group are expressed counts as a change in expression of the group.
#'   \item \emph{Moderation prior.} The overdispersions of genes and TE groups
#'   are moderated towards a prior estimated from both together, as edgeR does
#'   for all rows of an object. On data where the prior is above 1, gene
#'   overdispersions therefore differ slightly from those of edgeR run on
#'   genes alone.
#' }
#'
#' @section Earlier versions of edgeR:
#' Divide the counts and set the offset with \code{edgeR::scaleOffset()}.
#' This gives the same results as \code{edgeR::DGEListFromTximeta()}. The
#' offset is only needed with \code{counts_from_abundance = "no"}:
#' \preformatted{
#' counts <- assay(se, "counts") / rowData(se)$Overdispersion
#' y <- edgeR::DGEList(counts, group = group)
#' keep <- edgeR::filterByExpr(y)
#' y <- edgeR::normLibSizes(y[keep, , keep.lib.sizes = FALSE])
#' len <- assay(se, "length")[keep, ]
#' offset <- log(len / exp(rowMeans(log(len)))) +
#'     rep(log(y$samples$lib.size * y$samples$norm.factors), each = sum(keep))
#' y <- edgeR::scaleOffset(y, offset)
#' }
#'
#' @references
#' Soneson C, Love MI, Robinson MD (2015). Differential analyses for RNA-seq:
#' transcript-level estimates improve gene-level inferences. F1000Research,
#' 4:1521.
#'
#' Baldoni PL, Chen Y, Hediyeh-zadeh S, Liao Y, Dong X, Ritchie ME, Shi W,
#' Smyth GK (2024). Dividing out quantification uncertainty allows efficient
#' assessment of differential transcript expression with edgeR. Nucleic Acids
#' Research, 52(3):e13.
#'
#' @param x SummarizedExperiment object produced by \code{importQuants()}
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus). Must match
#' the species used to generate the index.
#' @param te_level One of "subfamily" (default), "locus", "family", or "class".
#' "locus" keeps TE-loci independent, while others aggregate based on
#' classification.
#' @param counts_from_abundance One of "no" (default), "scaledTPM", or
#' "lengthScaledTPM". How gene counts are computed. See Details.
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#'
#' @return SummarizedExperiment
#' @export
#'
#' @examples
#' bfc <- exampleCache(build = TRUE)
#' se <- importQuants(exampleQuants(), species = "Hs", cache = bfc)
#'
#' # Aggregate transcripts to genes, keep TE loci intact
#' locus <- aggregateCounts(se, "Hs", te_level = "locus", cache = bfc)
#'
#' # Aggregate transcripts to genes, TEs to subfamily
#' subfamily <- aggregateCounts(se, "Hs", te_level = "subfamily", cache = bfc)
#' subfamily
#'
#' # Generate gene counts from abundances
#' scaled <- aggregateCounts(
#'     se,
#'     "Hs",
#'     counts_from_abundance = "lengthScaledTPM",
#'     cache = bfc
#' )
aggregateCounts <- function(
    x,
    species = c("Hs", "Mm"),
    te_level = c("subfamily", "family", "class", "locus"),
    counts_from_abundance = c("no", "scaledTPM", "lengthScaledTPM"),
    cache = NULL
) {
    species <- match.arg(species)
    te_level <- match.arg(te_level)
    counts_from_abundance <- match.arg(counts_from_abundance)

    required <- c("orig", "tpms", "length")
    missing_assays <- setdiff(
        required,
        SummarizedExperiment::assayNames(x)
    )
    if (length(missing_assays)) {
        stop(
            "x is missing the assay(s): ",
            paste(missing_assays, collapse = ", "),
            ". Re-run importQuants() to create them.",
            call. = FALSE
        )
    }

    rd_orig <- SummarizedExperiment::rowData(x)
    feature_id <- .featureId(rd_orig, te_level)
    gene_ids <- unique(feature_id[!is.na(rd_orig$gene_id)])

    # Summed counts and abundances ---------------------------------------------
    orig <- SummarizedExperiment::assay(x, "orig")
    tpms <- SummarizedExperiment::assay(x, "tpms")
    eff_length <- SummarizedExperiment::assay(x, "length")

    counts <- rowsum(orig, group = feature_id)
    abundance <- rowsum(tpms, group = feature_id)
    is_gene <- rownames(counts) %in% gene_ids

    # Lengths: genes are the abundance weighted average of their transcripts in
    # each sample, as in tximport. TEs are simple sums and are not corrected
    # for length, so their length is the same in every sample (the summed
    # effective lengths averaged over samples). The offset that edgeR forms
    # from the centred log lengths is then zero for TEs.
    te_length <- rowMeans(rowsum(eff_length, group = feature_id))
    lengths_mat <- matrix(
        te_length,
        nrow = length(te_length),
        ncol = ncol(eff_length),
        dimnames = list(names(te_length), colnames(eff_length))
    )
    if (any(is_gene)) {
        lengths_mat[is_gene, ] <- .geneLength(
            tpms,
            eff_length,
            feature_id,
            rownames(counts)[is_gene]
        )
    }

    if (counts_from_abundance != "no" && any(is_gene)) {
        counts[is_gene, ] <- .countsFromAbundance(
            counts[is_gene, , drop = FALSE],
            abundance[is_gene, , drop = FALSE],
            lengths_mat[is_gene, , drop = FALSE],
            counts_from_abundance
        )
    }

    # Compute feature lengths
    txdb <- .getTxDb(species, .getCache(cache))

    # Gene lengths are reduced exon widths
    exons_by_gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
    reduced_exon_lengths <- sum(GenomicRanges::width(GenomicRanges::reduce(
        exons_by_gene
    )))
    names(reduced_exon_lengths) <- names(exons_by_gene)

    # TE lengths are sum of widths of ranges
    te_widths <- sum(GenomicRanges::width(rd_orig$Ranges))
    sum_te_widths <- tapply(te_widths, feature_id, sum, na.rm = TRUE)
    sum_te_widths <- sum_te_widths[
        !names(sum_te_widths) %in% names(reduced_exon_lengths)
    ]

    final_widths <- c(reduced_exon_lengths, sum_te_widths)

    # Number of genomic loci of the TEs in each feature
    is_te <- is.na(rd_orig$gene_id)
    n_loci <- tapply(rd_orig$N_Loci[is_te], feature_id[is_te], sum)

    gene_meta <- rd_orig[
        startsWith(rownames(x), "ENS"),
        c("gene_id", "gene_name", "gene_type")
    ]
    gene_meta <- unique(data.table::as.data.table(data.frame(gene_meta)))

    rd_new <- data.table::as.data.table(data.frame(table(
        feature_id = feature_id
    )))
    rd_new <- data.table::merge.data.table(
        rd_new,
        gene_meta,
        by.x = "feature_id",
        by.y = "gene_id",
        all.x = TRUE
    )

    data.table::setDF(rd_new, rownames = rd_new$feature_id)
    overdispersion <- S4Vectors::metadata(x)$overdispersion[[te_level]]
    # Without overdispersions the column is left out, so that
    # edgeR::DGEListFromTximeta(divide = TRUE) does not divide by NA
    if (is.null(overdispersion)) {
        message(
            "Overdispersions of the aggregated features are not available, ",
            "so the 'Overdispersion' column is not added. Re-run ",
            "importQuants() to estimate them."
        )
    } else {
        rd_new$Overdispersion <- unname(overdispersion[rownames(rd_new)])
    }
    rd_new$feature_length <- final_widths[rownames(rd_new)]
    rd_new$n_loci <- as.integer(n_loci[rownames(rd_new)])

    metadata <- S4Vectors::metadata(x)
    metadata$overdispersion <- NULL
    metadata$overdispersion.prior <- if (is.null(overdispersion)) {
        NA_real_
    } else {
        attr(overdispersion, "prior")
    }
    metadata$countsFromAbundance <- counts_from_abundance

    result <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            counts = counts,
            abundance = abundance,
            length = lengths_mat
        ),
        rowData = rd_new[rownames(counts), ],
        colData = SummarizedExperiment::colData(x),
        metadata = metadata
    )

    return(result)
}


# Levels that TE-loci can be aggregated to
.TE_LEVELS <- c("locus", "subfamily", "family", "class")


#' Feature IDs used to aggregate transcripts and TE-loci
#'
#' Transcripts are aggregated to their gene_id. TE-loci are kept as is
#' ("locus") or aggregated to their RepeatMasker classification.
#'
#' @param rd rowData of transcripts and TE-loci
#' @param te_level One of "locus", "subfamily", "family", or "class"
#'
#' @return Character vector of the feature ID of each row of \code{rd}
#' @keywords internal
.featureId <- function(rd, te_level) {
    te_grouping <- switch(
        te_level,
        locus = rownames(rd),
        subfamily = stringi::stri_c(
            rd$Subfamily,
            rd$Family,
            rd$Class,
            sep = ":"
        ),
        family = stringi::stri_c(rd$Family, rd$Class, sep = ":"),
        class = rd$Class
    )

    data.table::fcoalesce(rd$gene_id, te_grouping)
}


#' Abundance weighted gene lengths
#'
#' The average effective length of the transcripts of each gene, weighted by
#' their abundance in each sample. Genes without abundance in a sample are
#' given the geometric mean of their lengths in the other samples, and genes
#' without abundance in any sample the mean length of their transcripts.
#'
#' @param tpms,eff_length Matrices of transcript (and TE) TPMs and effective
#' lengths
#' @param feature_id Feature ID of each row
#' @param genes IDs of the genes to compute lengths for
#'
#' @return Matrix of gene lengths with rows in the order of \code{genes}
#' @keywords internal
.geneLength <- function(tpms, eff_length, feature_id, genes) {
    is_tx <- feature_id %in% genes
    tx_gene <- feature_id[is_tx]
    tx_tpms <- tpms[is_tx, , drop = FALSE]
    tx_length <- eff_length[is_tx, , drop = FALSE]

    weighted <- rowsum(tx_tpms * tx_length, tx_gene)[genes, , drop = FALSE]
    total <- rowsum(tx_tpms, tx_gene)[genes, , drop = FALSE]
    gene_length <- weighted / total

    no_tpm <- !is.finite(gene_length)
    if (any(no_tpm)) {
        # Mean length of the transcripts of each gene
        mean_length <- tapply(rowMeans(tx_length), tx_gene, mean)[genes]

        # Geometric mean of the other samples
        log_length <- log(gene_length)
        log_length[no_tpm] <- NA_real_
        geo_mean <- exp(rowMeans(log_length, na.rm = TRUE))
        fill <- ifelse(is.nan(geo_mean), mean_length, geo_mean)

        fill_mat <- matrix(fill, nrow(gene_length), ncol(gene_length))
        gene_length[no_tpm] <- fill_mat[no_tpm]
    }

    gene_length
}


#' Generate counts from abundance
#'
#' Gene counts proportional to their TPMs ("scaledTPM") or to their TPMs
#' multiplied by their average length over samples ("lengthScaledTPM"),
#' scaled so that the total counts of each sample are unchanged.
#'
#' @param counts,abundance,gene_length Matrices of gene counts, abundances, and
#' lengths
#' @param method One of "scaledTPM" or "lengthScaledTPM"
#'
#' @return Matrix of counts
#' @keywords internal
.countsFromAbundance <- function(counts, abundance, gene_length, method) {
    new_counts <- abundance
    if (method == "lengthScaledTPM") {
        new_counts <- new_counts * rowMeans(gene_length)
    }
    scale <- colSums(counts) / colSums(new_counts)
    sweep(new_counts, 2L, scale, "*")
}
