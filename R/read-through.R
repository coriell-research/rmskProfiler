#' Estimate read-through transcription of TE loci (experimental)
#'
#' TE loci located within the intron of a transcript, or downstream of its
#' 3' end, on the same strand may be quantified from reads derived from the
#' unspliced (pre-mRNA) or read-through transcript of the host rather than from
#' autonomous TE expression. This function estimates, for each TE, the number
#' of reads expected from host transcription alone and compares it to the
#' observed count.
#'
#' @details
#' \strong{Pairing TE loci with hosts.} Hosts are the transcripts in \code{x}
#' or genes (the span of their transcripts in \code{x}), with locations taken
#' from the \code{Ranges} column of the rowData. TE sequences that overlap an
#' exon on the same strand (\code{hasExonic == TRUE}) are excluded, since their
#' counts are confounded with mature mRNA. A remaining TE locus is
#' \emph{intronic} to a host if it lies entirely within the host's span on the
#' same strand, and \emph{downstream} of a host if it lies entirely within
#' \code{downstream} bp of the host's 3' end on the same strand. Transcripts
#' removed by \code{importQuants(remove_zeros = TRUE)} had no counts and so
#' cannot contribute read-through.
#'
#' \strong{Read-through rate.} Pre-mRNA reads outside of TE loci are not
#' quantified by the rmskProfiler index, so the read-through rate is estimated
#' from the TE loci themselves. The rate is the number of reads per bp of TE
#' per unit of host TPM, estimated separately for each sample so that
#' differences in sequencing depth cancel. It is pooled over \emph{reference}
#' TE loci: TE sequences with a single locus that is paired with a single host.
#' Only hosts with a TPM of at least \code{min_host_tpm} in a sample contribute
#' to rate estimates in that sample. For a group of reference loci the pooled
#' rate is
#'
#' \deqn{rate = \sum count / \sum (width \times TPM_{host})}
#'
#' which uses zero counts as evidence of a low rate rather than discarding them.
#'
#' For intronic TE loci, each host's rate is estimated from the other reference
#' loci in the same host (the TE locus being scored is left out) and shrunk
#' towards the genome-wide intronic rate using a gamma-Poisson model. The
#' amount of shrinkage is estimated from the variability of host rates across
#' the data, so hosts with many reference reads get rates close to their own,
#' while hosts with few get rates close to the genome-wide rate. The
#' \code{host_weight} column gives the share of the rate that comes from the
#' host itself. For downstream TE loci, read-through decays with distance from
#' the 3' end, so the rate is pooled over downstream reference loci within the
#' same distance bin. Bins are the \code{downstream_bins} quantiles of the
#' reference distances.
#'
#' \strong{Expected counts.} The expected read-through count for a TE locus in
#' sample j is \eqn{rate_j \times width \times TPM_{host,j}}. Expected counts are
#' summed over all hosts of a TE locus (e.g. overlapping genes) and over all
#' loci of a TE sequence.
#'
#' \strong{Per-TE model.} The expected counts assume every TE captures
#' read-through at the pooled rate, but TE loci differ in how efficiently they
#' do so, so the level of a TE's counts cannot separate read-through from
#' autonomous expression. The pattern across samples can: read-through counts
#' rise and fall with the host, autonomous counts do not. Each TE is fitted
#' with
#'
#' \deqn{count_j \sim NB(a \times expected_j + depth_j \times exp(X_j \beta))}
#'
#' with \eqn{a \ge 0}, where \eqn{depth_j} is the total count of sample j
#' relative to the mean over samples and the negative binomial dispersion is
#' shared by all TEs. The autonomous term is a log-linear model on the design
#' matrix \eqn{X}, interpreted as in edgeR or limma. By default \eqn{X} is an
#' intercept, so autonomous expression is constant across samples apart from
#' depth. Supply \code{group} or \code{design} to let it differ, for example
#' between treatments.
#'
#' The model can only separate the two components if host expression varies
#' across samples in a way that the design does not explain, since the design
#' can otherwise absorb the host's variation. TEs are only fitted if the
#' coefficient of variation of \eqn{expected_j / depth_j} after regressing on
#' the design is at least \code{min_host_cv}, and they have a mean count of at
#' least \code{min_mean_count}. The columns are
#' \itemize{
#'  \item{\code{host_cv}: coefficient of variation of expected read-through
#'  relative to depth across samples that is not explained by the design. Low
#'  values mean read-through and autonomous expression are hard to tell
#'  apart.}
#'  \item{\code{mean_observed}, \code{mean_expected}: the observed counts and
#'  the counts expected at the pooled read-through rate, averaged over
#'  samples.}
#'  \item{\code{efficiency}: \eqn{a}, the read-through efficiency of the TE
#'  relative to the pooled rate.}
#'  \item{\code{mean_readthrough}, \code{mean_autonomous}: the fitted
#'  read-through (\eqn{a \times expected}) and autonomous counts, averaged over
#'  samples.}
#'  \item{\code{rt_fraction}: the fraction of fitted counts from
#'  read-through.}
#'  \item{\code{cor}: Spearman correlation of observed and expected counts
#'  across samples.}
#'  \item{\code{pvalue_readthrough}, \code{padj_readthrough}: likelihood ratio
#'  test of whether the TE has a read-through component (\eqn{a > 0}), adjusted
#'  using the Benjamini-Hochberg method.}
#'  \item{\code{pvalue_autonomous}, \code{padj_autonomous}: likelihood ratio
#'  test of whether the TE has an autonomous component, comparing the full
#'  model to the model with read-through only.}
#'  \item{\code{zeros_observed}, \code{zeros_expected}: the number of samples
#'  with a count of zero, observed and predicted by the fitted model.}
#'  \item{\code{mapping_overdispersion}: the overdispersion of the TE's counts
#'  from mapping ambiguity, estimated from Salmon's Gibbs samples by
#'  \code{importQuants()}. The \code{counts} assay is already divided by it.}
#' }
#'
#' This function is experimental and its interface may change.
#'
#' @param x SummarizedExperiment object produced by \code{importQuants()}
#' @param level One of "gene" (default) or "transcript". The host features to
#' pair TE loci with. Gene-level TPMs are the sum of transcript TPMs.
#' @param downstream Maximum distance in bp downstream of the host 3' end for a
#' TE locus to be considered downstream. Default 10000. Set to 0 to only
#' consider intronic TE loci.
#' @param downstream_bins Number of distance bins used to estimate downstream
#' read-through rates. Default 5.
#' @param min_host_tpm Minimum host TPM in a sample for its TE loci to be used
#' to estimate read-through rates in that sample. Default 1.
#' @param min_mean_count Minimum mean observed count across samples for a TE
#' to be fitted. Default 5.
#' @param min_host_cv Minimum coefficient of variation of expected read-through
#' relative to depth across samples for a TE to be fitted. Default 0.1.
#' @param unique_only TRUE/FALSE, only consider TE sequences with a single locus
#' in the genome? Default TRUE. TE sequences with multiple loci are never used
#' as references.
#' @param design Optional design matrix for the autonomous component, with one
#' row per sample, e.g. from \code{model.matrix()}. Default NULL, an intercept
#' only.
#' @param group Optional factor with one value per sample giving the sample
#' groups, used to create the design \code{model.matrix(~group)}. Only one of
#' \code{design} and \code{group} can be given.
#'
#' @return data.frame with one row per TE with at least one host, containing
#' the TE hash, classification, and location, the primary host (the host
#' contributing the most expected counts) with its position ("intron" or
#' "downstream"), distance from its 3' end (0 for intronic TEs), and
#' host_weight (NA for downstream TEs), the number of hosts and all host IDs,
#' and the columns described in Details. Matrices of the expected counts at the
#' pooled rate and of the fitted read-through counts for each TE and sample are
#' attached as the "expected" and "readthrough" attributes, along with the
#' sample depth factors ("depth") and the dispersion ("dispersion"). The fitted
#' autonomous counts for each TE and sample, the autonomous coefficients for
#' each TE, and the design are attached as the "autonomous", "coefficients",
#' and "design" attributes.
#' @export
#'
#' @examples
#' \dontrun{
#' se <- importQuants("quants", species = "Hs")
#' rt <- estimateReadThrough(se)
#'
#' # Allow autonomous expression to differ between treatments
#' rt <- estimateReadThrough(se, group = se$treatment)
#'
#' # Or use any design matrix, as for edgeR or limma
#' design <- model.matrix(~ treatment + batch, data = colData(se))
#' rt <- estimateReadThrough(se, design = design)
#'
#' # TEs whose counts are mostly explained by read-through
#' subset(rt, padj_readthrough < 0.05 & rt_fraction > 0.8)
#'
#' # TEs with autonomous expression in addition to any read-through
#' subset(rt, padj_autonomous < 0.05)
#'
#' # Fitted read-through counts per sample
#' readthrough <- attr(rt, "readthrough")
#'
#' # The Hash column matches rownames(se) so it can be used to remove TEs whose
#' # counts are mostly explained by read-through before downstream analysis
#' drop <- rt$Hash[which(rt$padj_readthrough < 0.05 & rt$rt_fraction > 0.8)]
#' filtered <- se[!rownames(se) %in% drop, ]
#' }
estimateReadThrough <- function(
  x,
  level = c("gene", "transcript"),
  downstream = 10000,
  downstream_bins = 5,
  min_host_tpm = 1,
  min_mean_count = 5,
  min_host_cv = 0.1,
  unique_only = TRUE,
  design = NULL,
  group = NULL
) {
  level <- match.arg(level)
  design <- .autonomousDesign(design, group, ncol(x))

  rd <- SummarizedExperiment::rowData(x)
  counts <- SummarizedExperiment::assay(x, "counts")
  orig <- SummarizedExperiment::assay(x, "orig")
  tpms <- SummarizedExperiment::assay(x, "tpms")
  is_tx <- !is.na(rd$transcript_id)
  is_te <- !is.na(rd$Hash) & !rd$hasExonic
  if (isTRUE(unique_only)) {
    is_te <- is_te & rd$N_Loci == 1L
  }

  # TE loci ----------------------------------------------------------------------
  te_loci <- unlist(rd$Ranges[is_te])
  te_hash <- names(te_loci)
  names(te_loci) <- NULL
  te_width <- GenomicRanges::width(te_loci)

  # Salmon quantifies a TE from fragments that fit within it, so expected
  # counts scale with its effective length rather than its width
  te_ids <- unique(te_hash)
  eff_length <- .effectiveLengths(
    orig[te_ids, , drop = FALSE],
    tpms[te_ids, , drop = FALSE],
    te_width[match(te_ids, te_hash)]
  )

  # Host spans and expression ------------------------------------------------------
  tx_ids <- rownames(x)[is_tx]
  tx_span <- unlist(rd$Ranges[is_tx])

  if (level == "transcript") {
    host_span <- tx_span
    host_tpm <- tpms[tx_ids, , drop = FALSE]
    host_name <- stats::setNames(rd$gene_name[is_tx], tx_ids)
  } else {
    gene_ids <- rd$gene_id[is_tx]
    gene_span <- range(S4Vectors::splitAsList(tx_span, gene_ids))

    # Genes whose transcripts span more than one chromosome or strand are dropped
    host_span <- unlist(gene_span[lengths(gene_span) == 1L])
    host_tpm <- rowsum(tpms[tx_ids, , drop = FALSE], gene_ids)
    host_name <- stats::setNames(
      rd$gene_name[is_tx],
      gene_ids
    )[!duplicated(gene_ids)]
  }

  # Pair TE loci with hosts ----------------------------------------------------------
  # TE loci overlapping a same-strand exon were removed above, so any TE locus
  # within a same-strand host span lies in one of its introns
  message("Finding intronic and downstream TE loci...")
  hits <- GenomicRanges::findOverlaps(te_loci, host_span, type = "within")
  pairs <- data.frame(
    locus = S4Vectors::queryHits(hits),
    host = S4Vectors::subjectHits(hits),
    position = rep("intron", length(hits)),
    distance = rep(0L, length(hits))
  )
  if (downstream > 0) {
    pairs <- rbind(pairs, .downstreamPairs(te_loci, host_span, downstream))
  }
  if (nrow(pairs) == 0L) {
    message("No intronic or downstream TE loci found.")
    return(data.frame())
  }

  pair_hash <- te_hash[pairs$locus]
  pair_host <- names(host_span)[pairs$host]
  host_expr <- host_tpm[pair_host, , drop = FALSE]
  pair_counts <- counts[pair_hash, , drop = FALSE]
  exposure <- eff_length[pair_hash, , drop = FALSE] * host_expr

  # Reference pairs: single-locus TE sequences paired with a single host, used
  # in samples where the host is expressed
  pairs_per_locus <- tabulate(pairs$locus, nbins = length(te_loci))
  is_ref <- rd[pair_hash, "N_Loci"] == 1L & pairs_per_locus[pairs$locus] == 1L
  use <- is_ref & host_expr >= min_host_tpm
  ref_counts <- pair_counts * use
  ref_exposure <- exposure * use

  # Read-through rates for each pair and sample --------------------------------------
  # Rates are a level for each host (intronic) or distance bin (downstream),
  # constant across samples, multiplied by robust sample factors
  message("Estimating read-through rates for ", nrow(pairs), " TE-host pairs...")
  sample_rate <- .sampleRates(ref_counts, ref_exposure, pair_host)
  rate <- matrix(NA_real_, nrow(pairs), ncol(x))
  host_weight <- rep(NA_real_, nrow(pairs))
  host_dispersion <- NA_real_

  is_intron <- pairs$position == "intron"
  if (any(is_intron)) {
    intron <- .intronRates(
      ref_counts[is_intron, , drop = FALSE],
      ref_exposure[is_intron, , drop = FALSE],
      pair_host[is_intron],
      sample_rate
    )
    rate[is_intron, ] <- intron$rate
    host_weight[is_intron] <- intron$weight
    host_dispersion <- intron$dispersion
  }

  is_down <- !is_intron
  if (any(is_down)) {
    rate[is_down, ] <- .downstreamRates(
      ref_counts[is_down, , drop = FALSE],
      ref_exposure[is_down, , drop = FALSE],
      pairs$distance[is_down],
      is_ref[is_down],
      downstream_bins,
      sample_rate
    )
  }

  # Expected counts ----------------------------------------------------------------------
  expected_pairs <- rate * exposure
  expected <- rowsum(expected_pairs, pair_hash)
  colnames(expected) <- colnames(x)

  hash_ids <- rownames(expected)
  observed <- counts[hash_ids, , drop = FALSE]
  obs_total <- rowSums(observed)
  exp_total <- rowSums(expected)

  # Primary host: the pair contributing the most expected counts to each TE
  pair_total <- rowSums(expected_pairs)
  ord <- order(pair_hash, -pair_total, na.last = TRUE)
  primary <- ord[!duplicated(pair_hash[ord])]
  primary <- primary[match(hash_ids, pair_hash[primary])]
  host_list <- tapply(pair_host, pair_hash, function(h) {
    stringi::stri_flatten(unique(h), collapse = ";")
  })

  # Per-TE read-through model ------------------------------------------------------------
  # count_ij ~ NB(a_i * expected_ij + depth_j * exp(X_j beta_i)), a_i >= 0. a_i is
  # the read-through efficiency of the TE relative to the pooled rate and the
  # second term its host-independent (autonomous) expression, modelled on the
  # design as in edgeR
  message("Fitting read-through and autonomous components...")
  depth <- colSums(counts) / mean(colSums(counts))
  depth_mat <- matrix(depth, nrow(observed), ncol(observed), byrow = TRUE)

  # Host variation relative to depth that is not explained by the design, which
  # is needed to separate the two components
  shape <- expected / depth_mat
  hat <- design %*% solve(crossprod(design), t(design))
  resid <- shape - shape %*% hat
  host_cv <- sqrt(rowSums(resid^2) / (ncol(x) - ncol(design))) / rowMeans(shape)

  fitted <- rowMeans(observed) >= min_mean_count &
    rowSums(expected) > 0 &
    is.finite(host_cv) &
    host_cv >= min_host_cv
  fit <- .fitReadThrough(
    observed[fitted, , drop = FALSE],
    expected[fitted, , drop = FALSE],
    depth,
    design
  )

  n_te <- length(hash_ids)
  a <- rt_fraction <- p_rt <- p_auto <- zeros_expected <- rep(NA_real_, n_te)
  a[fitted] <- fit$a
  readthrough <- expected * a
  autonomous <- matrix(
    NA_real_,
    n_te,
    ncol(x),
    dimnames = list(hash_ids, colnames(x))
  )
  autonomous[fitted, ] <- fit$autonomous
  coefficients <- matrix(
    NA_real_,
    n_te,
    ncol(design),
    dimnames = list(hash_ids, colnames(design))
  )
  coefficients[fitted, ] <- fit$beta
  rt_fraction[fitted] <- rowSums(readthrough[fitted, , drop = FALSE]) /
    rowSums(fit$mu)
  p_rt[fitted] <- fit$pvalue_readthrough
  p_auto[fitted] <- fit$pvalue_autonomous
  zeros_expected[fitted] <- fit$zeros_expected
  zeros_observed <- rowSums(round(observed) == 0)

  r <- rep(NA_real_, n_te)
  if (ncol(x) >= 3L) {
    r[fitted] <- .rowCor(
      observed[fitted, , drop = FALSE],
      expected[fitted, , drop = FALSE],
      method = "spearman"
    )
  }

  # Mapping ambiguity from Salmon's Gibbs samples (see importQuants())
  ann <- S4Vectors::metadata(x)$annotation
  mapping_od <- if (!is.null(ann$Overdispersion)) {
    ann[hash_ids, "Overdispersion"]
  } else {
    rep(NA_real_, n_te)
  }

  primary_locus <- te_loci[pairs$locus[primary]]
  primary_host <- pair_host[primary]
  result <- data.frame(
    Hash = hash_ids,
    Subfamily = rd[hash_ids, "Subfamily"],
    Family = rd[hash_ids, "Family"],
    Class = rd[hash_ids, "Class"],
    N_Loci = rd[hash_ids, "N_Loci"],
    seqnames = as.character(GenomicRanges::seqnames(primary_locus)),
    start = GenomicRanges::start(primary_locus),
    end = GenomicRanges::end(primary_locus),
    strand = as.character(GenomicRanges::strand(primary_locus)),
    host_id = primary_host,
    host_name = unname(host_name[primary_host]),
    position = pairs$position[primary],
    distance = pairs$distance[primary],
    host_weight = host_weight[primary],
    n_hosts = as.integer(lengths(strsplit(host_list[hash_ids], ";"))),
    hosts = unname(host_list[hash_ids]),
    host_cv = host_cv,
    mean_observed = rowMeans(observed),
    mean_expected = rowMeans(expected),
    efficiency = a,
    mean_readthrough = rowMeans(readthrough),
    mean_autonomous = rowMeans(autonomous),
    rt_fraction = rt_fraction,
    cor = r,
    pvalue_readthrough = p_rt,
    padj_readthrough = stats::p.adjust(p_rt, method = "BH"),
    pvalue_autonomous = p_auto,
    padj_autonomous = stats::p.adjust(p_auto, method = "BH"),
    zeros_observed = zeros_observed,
    zeros_expected = zeros_expected,
    mapping_overdispersion = mapping_od
  )

  ord <- order(result$padj_readthrough, -result$rt_fraction, na.last = TRUE)
  result <- result[ord, ]
  rownames(result) <- NULL
  attr(result, "expected") <- expected[result$Hash, , drop = FALSE]
  attr(result, "readthrough") <- readthrough[result$Hash, , drop = FALSE]
  attr(result, "autonomous") <- autonomous[result$Hash, , drop = FALSE]
  attr(result, "coefficients") <- coefficients[result$Hash, , drop = FALSE]
  attr(result, "design") <- design
  attr(result, "depth") <- depth
  attr(result, "dispersion") <- fit$dispersion
  attr(result, "host_dispersion") <- host_dispersion

  return(result)
}


#' Effective lengths of TEs from counts and TPMs
#'
#' Salmon TPMs are counts divided by effective length and rescaled in each
#' sample, so the ratio of counts to TPM is the effective length multiplied by
#' a constant for each sample. The constant cancels when read-through rates are
#' estimated separately for each sample. In samples where a TE has no counts,
#' its effective length is taken from its other samples. TEs with no counts in
#' any sample are given the median ratio of effective length to width of TEs
#' with a similar width.
#'
#' @param orig,tpms Matrices of Salmon counts and TPMs for TEs
#' @param width Width of each TE
#'
#' @return Matrix of effective lengths on the scale of each sample
#' @keywords internal
.effectiveLengths <- function(orig, tpms, width) {
  log_eff <- log(orig / tpms)
  log_eff[!(orig > 0 & tpms > 0)] <- NA_real_

  te_mean <- rowMeans(log_eff, na.rm = TRUE)
  sample_scale <- apply(log_eff - te_mean, 2, stats::median, na.rm = TRUE)

  # TEs without counts in any sample: use TEs of similar width
  no_counts <- !is.finite(te_mean)
  if (any(no_counts)) {
    breaks <- unique(stats::quantile(
      width[!no_counts],
      probs = seq(0, 1, length.out = 21)
    ))
    bin <- findInterval(width, breaks, all.inside = TRUE)
    log_ratio <- tapply(
      te_mean[!no_counts] - log(width[!no_counts]),
      bin[!no_counts],
      stats::median
    )
    te_mean[no_counts] <- log(width[no_counts]) +
      log_ratio[as.character(bin[no_counts])]
  }

  fill <- outer(te_mean, sample_scale, "+")
  missing <- is.na(log_eff)
  log_eff[missing] <- fill[missing]

  exp(log_eff)
}


#' Robust sample factors for read-through rates
#'
#' Estimates the read-through rate in each sample up to a constant, in the
#' same way as DESeq2 size factors: each host's rate in each sample is divided
#' by its geometric mean over samples, and the factor for a sample is the
#' median of these ratios over hosts. Using the median means a minority of
#' hosts containing TEs with autonomous expression cannot bias the factors.
#' Only hosts with reference counts in every sample and at least 50 reference
#' counts in total are used. If there are fewer than 10 such hosts, the rate
#' pooled over all reference loci in each sample is used instead.
#'
#' @param counts,exposure Matrices of reference counts and exposures
#' (effective length x host TPM), zero for pairs or samples not used as
#' references
#' @param host Host of each pair
#'
#' @return Numeric vector of sample factors
#' @keywords internal
.sampleRates <- function(counts, exposure, host) {
  pooled <- colSums(counts) / colSums(exposure)
  if (any(!is.finite(pooled) | pooled <= 0)) {
    stop(
      "Not enough reference TE loci in expressed hosts to estimate ",
      "read-through rates in every sample."
    )
  }

  host_counts <- rowsum(counts, host)
  host_exposure <- rowsum(exposure, host)
  keep <- rowSums(host_counts) >= 50 &
    apply(host_counts > 0 & host_exposure > 0, 1, all)
  if (sum(keep) < 10L) {
    return(pooled)
  }

  log_rate <- log(host_counts[keep, , drop = FALSE] /
    host_exposure[keep, , drop = FALSE])
  factors <- exp(apply(log_rate - rowMeans(log_rate), 2, stats::median))

  factors
}


#' Read-through rates for intronic TE loci
#'
#' The rate for a pair in sample j is the sample factor for sample j
#' multiplied by an intronic level for its host that is constant across
#' samples. The host level is the genome-wide intronic level multiplied by the
#' host's relative rate, which is pooled over its other reference loci (the TE
#' locus being scored is left out) and all samples, and shrunk towards 1 using
#' a gamma-Poisson model. The prior is a gamma distribution with mean 1 and
#' squared coefficient of variation equal to the variance of host relative
#' rates beyond Poisson noise, estimated by the method of moments. Keeping
#' host rates constant across samples avoids adding sampling noise from a few
#' reference loci to every sample.
#'
#' @param counts,exposure Matrices of reference counts and exposures
#' (effective length x host TPM) for intronic pairs, zero for pairs or samples
#' not used as references
#' @param host Host of each pair
#' @param sample_rate Sample factors from \code{.sampleRates()}
#'
#' @return List with elements rate (matrix), weight (the share of the host
#' relative rate from the host itself for each pair), and dispersion (the
#' squared coefficient of variation of host relative rates)
#' @keywords internal
.intronRates <- function(counts, exposure, host, sample_rate) {
  sample_mat <- matrix(sample_rate, nrow(counts), ncol(counts), byrow = TRUE)
  level <- sum(counts) / sum(exposure * sample_mat)
  base_rate <- sample_mat * level

  # Counts and counts expected at the genome-wide rate, summed over samples
  pair_counts <- rowSums(counts)
  pair_expected <- rowSums(exposure * base_rate)
  host_counts <- rowsum(pair_counts, host)[, 1]
  host_expected <- rowsum(pair_expected, host)[, 1]

  # Variance of host relative rates beyond Poisson noise, from hosts with at
  # least 5 expected reads
  keep <- host_expected >= 5
  dispersion <- if (sum(keep) >= 2L) {
    r <- host_counts[keep] / host_expected[keep]
    max(0, mean((r - 1)^2 - 1 / host_expected[keep]))
  } else {
    0
  }

  if (dispersion <= 0) {
    # No detectable variation between hosts, use the genome-wide rate
    return(list(
      rate = base_rate,
      weight = rep(0, nrow(counts)),
      dispersion = 0
    ))
  }

  # Leave out each pair's own counts, then shrink towards 1 with a gamma prior
  # with shape and rate alpha
  loo_counts <- pmax(host_counts[host] - pair_counts, 0)
  loo_expected <- pmax(host_expected[host] - pair_expected, 0)
  alpha <- 1 / dispersion
  relative <- (loo_counts + alpha) / (loo_expected + alpha)

  list(
    rate = base_rate * relative,
    weight = loo_expected / (loo_expected + alpha),
    dispersion = dispersion
  )
}


#' Read-through rates for downstream TE loci
#'
#' The rate for a pair in sample j is the sample factor for sample j
#' multiplied by a level for its distance bin that is constant across samples,
#' pooled over the downstream reference loci in the bin. Bins without
#' reference exposure use the level over all downstream reference loci.
#'
#' @param counts,exposure Matrices of reference counts and exposures
#' (effective length x host TPM) for downstream pairs, zero for pairs or
#' samples not used as references
#' @param distance Distance of each pair from the host 3' end
#' @param is_ref Logical, is each pair a reference pair
#' @param n_bins Number of distance bins
#' @param sample_rate Sample factors from \code{.sampleRates()}
#'
#' @return Matrix of rates
#' @keywords internal
.downstreamRates <- function(
  counts,
  exposure,
  distance,
  is_ref,
  n_bins,
  sample_rate
) {
  sample_mat <- matrix(sample_rate, nrow(counts), ncol(counts), byrow = TRUE)
  pair_counts <- rowSums(counts)
  pair_expected <- rowSums(exposure * sample_mat)
  overall <- sum(pair_counts) / sum(pair_expected)
  if (!is.finite(overall)) {
    stop(
      "Not enough downstream reference TE loci in expressed hosts to estimate ",
      "read-through rates."
    )
  }

  breaks <- unique(stats::quantile(
    distance[is_ref],
    probs = seq(0, 1, length.out = n_bins + 1)
  ))
  bin <- findInterval(distance, breaks, all.inside = TRUE)

  bin_level <- rowsum(pair_counts, bin)[, 1] / rowsum(pair_expected, bin)[, 1]
  level <- unname(bin_level[as.character(bin)])
  level[!is.finite(level)] <- overall

  sample_mat * level
}


#' Design matrix for the autonomous component
#'
#' @param design NULL or a design matrix with one row per sample
#' @param group NULL or a factor with one value per sample
#' @param n_samples Number of samples
#'
#' @return Design matrix. An intercept-only design if both design and group
#' are NULL.
#' @keywords internal
.autonomousDesign <- function(design, group, n_samples) {
  if (!is.null(design) && !is.null(group)) {
    stop("Only one of design and group can be given.")
  }
  if (!is.null(group)) {
    group <- as.factor(group)
    if (length(group) != n_samples) {
      stop("group must have one value for each sample.")
    }
    design <- stats::model.matrix(~group)
  }
  if (is.null(design)) {
    design <- matrix(1, n_samples, 1, dimnames = list(NULL, "(Intercept)"))
  }
  design <- as.matrix(design)
  if (nrow(design) != n_samples) {
    stop("design must have one row for each sample.")
  }
  if (qr(design)$rank < ncol(design)) {
    stop("design is not of full rank.")
  }
  if (n_samples - ncol(design) < 2L) {
    stop(
      "At least 2 more samples than design coefficients are needed to ",
      "separate read-through from autonomous expression."
    )
  }

  design
}


#' Fit the per-TE read-through model
#'
#' Fits count_ij ~ NB(a_i * expected_ij + depth_j * exp(X_j beta_i)) with
#' a_i >= 0 for every TE, along with the models with only read-through (no
#' autonomous term) and only autonomous expression (a_i = 0), using a common
#' negative binomial dispersion. Each component is tested with a likelihood
#' ratio test against the model without it. Because the null models are on the
#' boundary of the parameter space, the read-through test uses a 50:50 mixture
#' of chi-squared distributions with 0 and 1 degrees of freedom and the
#' autonomous test a binomial mixture of chi-squared distributions with 0 to p
#' degrees of freedom, where p is the number of design coefficients.
#'
#' @param y Matrix of observed counts
#' @param e Matrix of expected read-through counts
#' @param depth Sample depth factors
#' @param design Design matrix for the autonomous component
#' @param n_rounds Number of rounds alternating between fitting the models and
#' estimating the dispersion
#'
#' @return List with a, beta, autonomous (fitted autonomous counts), mu
#' (fitted means), pvalue_readthrough, pvalue_autonomous, zeros_expected, and
#' dispersion
#' @keywords internal
.fitReadThrough <- function(y, e, depth, design, n_rounds = 3L) {
  n_coef <- ncol(design)
  dispersion <- 0.1
  for (round in seq_len(n_rounds + 1L)) {
    rt_only <- .fitReadThroughOnly(y, e, dispersion)
    auto_only <- .fitAutonomousOnly(y, depth, design, dispersion)
    full <- .fitFull(y, e, depth, design, dispersion, rt_only, auto_only)
    if (round <= n_rounds) {
      dispersion <- .momentDispersion(y, full$mu, n_params = n_coef + 1L)
    }
  }

  stat_rt <- 2 * (full$loglik - auto_only$loglik)
  stat_auto <- 2 * (full$loglik - rt_only$loglik)

  list(
    a = full$a,
    beta = full$beta,
    autonomous = full$mu - full$a * e,
    mu = full$mu,
    pvalue_readthrough = .chiBarSquare(stat_rt, 1L),
    pvalue_autonomous = .chiBarSquare(stat_auto, n_coef),
    zeros_expected = rowSums((1 + dispersion * full$mu)^(-1 / dispersion)),
    dispersion = dispersion
  )
}


#' Fit the read-through only model
#'
#' Fits count_ij ~ NB(a_i * expected_ij) with a_i >= 0 by iteratively
#' reweighted least squares.
#'
#' @param y,e Matrices of observed and expected counts
#' @param dispersion Negative binomial dispersion
#' @param n_iter Number of iterations
#'
#' @return List with a, mu, and loglik
#' @keywords internal
.fitReadThroughOnly <- function(y, e, dispersion, n_iter = 25L) {
  a <- rowSums(y) / rowSums(e)
  for (iter in seq_len(n_iter)) {
    mu <- pmax(a * e, 1e-8)
    w <- 1 / (mu + dispersion * mu^2)
    a <- rowSums(w * e * y) / rowSums(w * e^2)
  }
  mu <- pmax(a * e, 1e-8)

  list(a = a, mu = mu, loglik = .nbLogLik(y, mu, dispersion))
}


#' Fit the autonomous only model
#'
#' Fits count_ij ~ NB(depth_j * exp(X_j beta_i)), a negative binomial GLM with
#' a log link and log(depth) offset, by iteratively reweighted least squares.
#'
#' @param y Matrix of observed counts
#' @param depth Sample depth factors
#' @param design Design matrix
#' @param dispersion Negative binomial dispersion
#' @param n_iter Number of iterations
#'
#' @return List with beta, mu, and loglik
#' @keywords internal
.fitAutonomousOnly <- function(y, depth, design, dispersion, n_iter = 25L) {
  offset <- matrix(log(depth), nrow(y), ncol(y), byrow = TRUE)
  n_coef <- ncol(design)

  # Start from the least squares fit of log counts
  z <- log(y + 0.5) - offset
  beta <- t(solve(crossprod(design), crossprod(design, t(z))))
  for (iter in seq_len(n_iter)) {
    eta <- pmin(beta %*% t(design) + offset, 50)
    mu <- exp(eta)
    w <- mu / (1 + dispersion * mu)
    z <- eta - offset + (y - mu) / mu
    info <- .crossWeighted(lapply(seq_len(n_coef), function(k) {
      matrix(design[, k], nrow(y), ncol(y), byrow = TRUE)
    }), w)
    score <- vapply(
      seq_len(n_coef),
      function(k) rowSums(w * z * rep(design[, k], each = nrow(y))),
      numeric(nrow(y))
    )
    beta <- pmax(.batchSolve(info, matrix(score, nrow(y))), -30)
  }
  mu <- pmax(exp(pmin(beta %*% t(design) + offset, 50)), 1e-8)

  list(beta = beta, mu = mu, loglik = .nbLogLik(y, mu, dispersion))
}


#' Fit the full read-through and autonomous model
#'
#' Fits count_ij ~ NB(a_i * expected_ij + depth_j * exp(X_j beta_i)) with
#' a_i >= 0 by Fisher scoring with step halving, starting from the read-through
#' only and autonomous only fits. The fit for each TE is never worse than
#' either of these, which are limits of the full model.
#'
#' @param y,e Matrices of observed and expected counts
#' @param depth Sample depth factors
#' @param design Design matrix
#' @param dispersion Negative binomial dispersion
#' @param rt_only,auto_only Fits from \code{.fitReadThroughOnly()} and
#' \code{.fitAutonomousOnly()}
#' @param n_iter Number of iterations
#'
#' @return List with a, beta, mu, and loglik
#' @keywords internal
.fitFull <- function(
  y,
  e,
  depth,
  design,
  dispersion,
  rt_only,
  auto_only,
  n_iter = 50L
) {
  n <- nrow(y)
  n_coef <- ncol(design)
  depth_mat <- matrix(depth, n, ncol(y), byrow = TRUE)
  mean_fn <- function(a, beta, rows) {
    auto <- depth_mat[rows, , drop = FALSE] *
      exp(pmin(beta %*% t(design), 50))
    list(mu = pmax(a * e[rows, , drop = FALSE] + auto, 1e-8), auto = auto)
  }

  # Start halfway between the two single-component fits. Halving the
  # autonomous term shifts the linear predictor by log(0.5), which is possible
  # when the design contains a constant
  a <- rt_only$a / 2
  beta <- auto_only$beta
  const <- qr.coef(qr(design), rep(1, nrow(design)))
  if (max(abs(design %*% const - 1)) < 1e-8) {
    beta <- beta + log(0.5) * matrix(const, n, n_coef, byrow = TRUE)
  }
  m <- mean_fn(a, beta, seq_len(n))
  ll <- .nbLogLik(y, m$mu, dispersion)

  # Only TEs that have not converged are updated in each iteration
  active <- seq_len(n)
  for (iter in seq_len(n_iter)) {
    n_active <- length(active)
    y_a <- y[active, , drop = FALSE]
    mu_a <- m$mu[active, , drop = FALSE]
    auto_a <- m$auto[active, , drop = FALSE]
    a_a <- a[active]
    beta_a <- beta[active, , drop = FALSE]
    ll_a <- ll[active]

    w <- 1 / (mu_a + dispersion * mu_a^2)
    jac <- c(
      list(e[active, , drop = FALSE]),
      lapply(seq_len(n_coef), function(k) {
        auto_a * matrix(design[, k], n_active, ncol(y), byrow = TRUE)
      })
    )
    info <- .crossWeighted(jac, w)
    resid <- w * (y_a - mu_a)
    score <- vapply(jac, function(j) rowSums(j * resid), numeric(n_active))
    step <- .batchSolve(info, matrix(score, n_active))

    # Step halving to ensure the log-likelihood does not decrease
    scale <- rep(1, n_active)
    todo <- rep(TRUE, n_active)
    new_a <- a_a
    new_beta <- beta_a
    new_ll <- ll_a
    new_mu <- mu_a
    new_auto <- auto_a
    for (half in seq_len(10L)) {
      cand_a <- pmax(a_a + scale * step[, 1], 0)
      cand_beta <- pmax(beta_a + scale * step[, -1, drop = FALSE], -30)
      cand_m <- mean_fn(cand_a, cand_beta, active)
      cand_ll <- .nbLogLik(y_a, cand_m$mu, dispersion)
      accept <- todo & cand_ll >= ll_a - 1e-10
      new_a[accept] <- cand_a[accept]
      new_beta[accept, ] <- cand_beta[accept, ]
      new_ll[accept] <- cand_ll[accept]
      new_mu[accept, ] <- cand_m$mu[accept, ]
      new_auto[accept, ] <- cand_m$auto[accept, ]
      todo <- todo & !accept
      if (!any(todo)) {
        break
      }
      scale[todo] <- scale[todo] / 2
    }

    a[active] <- new_a
    beta[active, ] <- new_beta
    ll[active] <- new_ll
    m$mu[active, ] <- new_mu
    m$auto[active, ] <- new_auto
    active <- active[abs(new_ll - ll_a) >= 1e-8]
    if (length(active) == 0L) {
      break
    }
  }

  # The single-component fits are limits of the full model
  use_rt <- rt_only$loglik > ll
  a[use_rt] <- rt_only$a[use_rt]
  beta[use_rt, ] <- -30
  m$mu[use_rt, ] <- rt_only$mu[use_rt, ]
  ll[use_rt] <- rt_only$loglik[use_rt]

  use_auto <- auto_only$loglik > ll
  a[use_auto] <- 0
  beta[use_auto, ] <- auto_only$beta[use_auto, ]
  m$mu[use_auto, ] <- auto_only$mu[use_auto, ]
  ll[use_auto] <- auto_only$loglik[use_auto]

  list(a = a, beta = beta, mu = m$mu, loglik = ll)
}


#' Weighted cross products of a list of matrices for every row
#'
#' @param mats List of q matrices with the same dimensions
#' @param w Matrix of weights with the same dimensions
#'
#' @return Array with dimensions nrow x q x q where element \code{[i, u, v]} is
#' \code{sum_j w[i, j] * mats[[u]][i, j] * mats[[v]][i, j]}
#' @keywords internal
.crossWeighted <- function(mats, w) {
  q <- length(mats)
  out <- array(0, c(nrow(w), q, q))
  for (u in seq_len(q)) {
    for (v in seq_len(u)) {
      out[, u, v] <- out[, v, u] <- rowSums(w * mats[[u]] * mats[[v]])
    }
  }

  out
}


#' Solve many small symmetric positive definite systems
#'
#' Solves A_i x_i = b_i for every row i using a Cholesky decomposition,
#' vectorized over rows. A small ridge is added to the diagonal for stability.
#'
#' @param A Array with dimensions n x q x q
#' @param b Matrix with dimensions n x q
#'
#' @return Matrix with dimensions n x q
#' @keywords internal
.batchSolve <- function(A, b) {
  n <- dim(A)[1]
  q <- dim(A)[2]
  ridge <- 1e-10 * pmax(apply(A, 1, function(m) max(diag(m))), 1)
  L <- array(0, c(n, q, q))
  for (j in seq_len(q)) {
    prev <- seq_len(j - 1L)
    d <- A[, j, j] + ridge
    if (j > 1L) {
      d <- d - rowSums(matrix(L[, j, prev], n)^2)
    }
    L[, j, j] <- sqrt(pmax(d, 1e-12))
    for (i in seq_len(q)[-seq_len(j)]) {
      s <- A[, i, j]
      if (j > 1L) {
        s <- s - rowSums(matrix(L[, i, prev], n) * matrix(L[, j, prev], n))
      }
      L[, i, j] <- s / L[, j, j]
    }
  }

  z <- matrix(0, n, q)
  for (i in seq_len(q)) {
    s <- b[, i]
    if (i > 1L) {
      prev <- seq_len(i - 1L)
      s <- s - rowSums(matrix(L[, i, prev], n) * matrix(z[, prev], n))
    }
    z[, i] <- s / L[, i, i]
  }
  x <- matrix(0, n, q)
  for (i in rev(seq_len(q))) {
    s <- z[, i]
    if (i < q) {
      nxt <- seq(i + 1L, q)
      s <- s - rowSums(matrix(L[, nxt, i], n) * matrix(x[, nxt], n))
    }
    x[, i] <- s / L[, i, i]
  }

  x
}


#' P-values for a likelihood ratio test with parameters on the boundary
#'
#' @param stat Likelihood ratio statistics
#' @param df Number of parameters tested
#'
#' @return p-values from a binomial mixture of chi-squared distributions with 0
#' to df degrees of freedom
#' @keywords internal
.chiBarSquare <- function(stat, df) {
  stat <- pmax(stat, 0)
  weights <- stats::dbinom(0:df, df, 0.5)
  p <- vapply(
    seq_len(df),
    function(k) weights[k + 1L] * stats::pchisq(stat, k, lower.tail = FALSE),
    numeric(length(stat))
  )
  p <- matrix(p, length(stat))
  ifelse(stat > 1e-8, rowSums(p), 1)
}


#' Moment estimate of a common negative binomial dispersion
#'
#' @param y Matrix of observed counts
#' @param mu Matrix of fitted means
#' @param n_params Number of fitted parameters per row
#'
#' @return Median over rows with a mean fitted count of at least 1 of the
#' per-row moment estimates, where Var(y) = mu + dispersion * mu^2
#' @keywords internal
.momentDispersion <- function(y, mu, n_params) {
  keep <- rowMeans(mu) >= 1
  if (sum(keep) < 10L) {
    warning(
      "Fewer than 10 TEs with enough counts to estimate the dispersion. ",
      "Using a dispersion of 0.1."
    )
    return(0.1)
  }
  n <- ncol(y)
  y <- y[keep, , drop = FALSE]
  mu <- mu[keep, , drop = FALSE]
  per_te <- (rowSums((y - mu)^2) * n / (n - n_params) - rowSums(mu)) /
    rowSums(mu^2)

  max(stats::median(per_te), 1e-4)
}


#' Negative binomial log-likelihood summed over the columns of each row
#'
#' Uses the gamma function form so that non-integer counts are allowed.
#'
#' @param y Matrix of observed counts
#' @param mu Matrix of means
#' @param dispersion Negative binomial dispersion
#'
#' @return Numeric vector of log-likelihoods
#' @keywords internal
.nbLogLik <- function(y, mu, dispersion) {
  size <- 1 / dispersion
  ll <- lgamma(y + size) -
    lgamma(size) -
    lgamma(y + 1) +
    size * log(size / (size + mu)) +
    y * log(mu / (size + mu))

  rowSums(ll)
}


#' Find TE loci downstream of the 3' end of same-strand hosts
#'
#' @param te GRanges of TE loci
#' @param span GRanges of host spans
#' @param downstream Maximum distance in bp downstream of the host 3' end
#'
#' @return data.frame with columns locus, host, position, and distance
#' @keywords internal
.downstreamPairs <- function(te, span, downstream) {
  region <- suppressWarnings(GenomicRanges::trim(
    GenomicRanges::flank(span, width = downstream, start = FALSE)
  ))
  hits <- GenomicRanges::findOverlaps(te, region, type = "within")
  locus <- S4Vectors::queryHits(hits)
  host <- S4Vectors::subjectHits(hits)

  data.frame(
    locus = locus,
    host = host,
    position = rep("downstream", length(hits)),
    distance = GenomicRanges::distance(te[locus], span[host])
  )
}


#' Row-wise correlation between two matrices of the same dimensions
#'
#' @param x,y Numeric matrices
#' @param method One of "spearman" or "pearson"
#'
#' @return Numeric vector of correlations, NA where a row has zero variance
#' @keywords internal
.rowCor <- function(x, y, method = c("spearman", "pearson")) {
  method <- match.arg(method)
  if (nrow(x) == 0L) {
    return(numeric())
  }
  if (method == "spearman") {
    x <- matrix(t(apply(x, 1, rank)), nrow = nrow(x))
    y <- matrix(t(apply(y, 1, rank)), nrow = nrow(y))
  }
  x <- x - rowMeans(x)
  y <- y - rowMeans(y)
  r <- rowSums(x * y) / sqrt(rowSums(x^2) * rowSums(y^2))
  r[!is.finite(r)] <- NA_real_

  return(r)
}
