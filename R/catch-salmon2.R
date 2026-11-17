#' Import transcript-level counts with offset
#'
#' @details
#' This function is mostly lifted from:
#'
#' https://github.com/plbaldoni/TranscriptDE-code/blob/main/code/pkg/R/utils.R#L52
#'
#' and combined with existing catchSalmon code. I only used data.table::fread in
#' place of readr so
#' as to not add more dependencies.
#'
#' When \code{groups} is given, the Gibbs samples of the targets in each group
#' are also summed and the overdispersion of each group is estimated from the
#' summed samples, in the same way as for single targets. Read-to-target
#' ambiguity between the members of a group cancels when they are summed,
#' while ambiguity between groups does not.
#'
#' @param paths character vector giving paths to the sample-specific directories
#' created by a kallisto or
#' Salmon. Each entry corresponds to one RNA-seq sample.
#' @param verbose logical If TRUE, progress information will be sent to standard
#' output as each sample is processed.
#' @param groups NULL or a named list of character vectors, each giving the
#' group of every target in the order of the quant.sf files. Targets with an
#' NA group are left out of that grouping.
#'
#' @return A list with the counts, TPMs, and effective lengths of each target,
#' an annotation data.frame with the overdispersion of each target, the
#' overdispersion prior, the resample type of each sample, and (when
#' \code{groups} is given) a list with the overdispersion of each group for
#' each grouping.
#' @keywords internal
.catchSalmon2 <- function(paths, verbose = TRUE, groups = NULL) {
    NSamples <- length(paths)

    OK <- requireNamespace("jsonlite", quietly = TRUE)
    if (!OK) {
        stop(
            "jsonlite package required but is not installed (or can't be ",
            "loaded)"
        )
    }

    # Accumulators for the overdispersion of each group, filled from the
    # first sample with Gibbs samples
    GroupOverDisp <- GroupDF <- NULL

    ResampleType <- rep_len("", NSamples)
    for (j in 1L:NSamples) {
        if (verbose) {
            message("Reading ", paths[j], ", ", appendLF = FALSE)
        }
        MetaFile <- file.path(paths[j], "aux_info", "meta_info.json")
        QuantFile <- file.path(paths[j], "quant.sf")
        BootFile <- file.path(
            paths[j],
            "aux_info",
            "bootstrap",
            "bootstraps.gz"
        )

        if (!file.exists(QuantFile)) {
            stop("quant.sf file not found at specified path")
        }

        Meta <- jsonlite::fromJSON(MetaFile)
        NTx <- Meta$num_targets

        if (is.null(NTx)) {
            NTx <- Meta$num_valid_targets
        }

        if (is.null(NTx)) {
            stop("Can't find number of targets")
        }

        NBoot <- Meta$num_bootstraps
        if (is.null(NBoot)) {
            stop("Can't find number of bootstraps")
        }

        Type <- Meta$samp_type
        if (is.null(Type)) {
            ResampleType[j] <- "bootstrap"
        } else {
            ResampleType[j] <- Type
        }

        if (verbose) {
            message(NTx, " transcripts, ", NBoot, " ", Type, " samples")
        }

        if (j == 1L) {
            Counts <- matrix(0, NTx, NSamples)
            TPM <- matrix(0, NTx, NSamples)
            Length <- matrix(0, NTx, NSamples)
            DF <- rep_len(0L, NTx)
            OverDisp <- rep_len(0, NTx)
            Quant1 <- data.table::fread(
                QuantFile,
                sep = "\t",
                select = c(
                    "Name" = "character",
                    "EffectiveLength" = "numeric",
                    "TPM" = "numeric",
                    "NumReads" = "numeric"
                ),
                showProgress = FALSE
            )
            Counts[, 1L] <- Quant1$NumReads
            TPM[, 1L] <- Quant1$TPM
            Length[, 1L] <- Quant1$EffectiveLength
            if (!is.null(groups)) {
                stopifnot(all(lengths(groups) == NTx))
            }
        } else {
            Quant <- data.table::fread(
                QuantFile,
                sep = "\t",
                select = c(
                    "EffectiveLength" = "numeric",
                    "TPM" = "numeric",
                    "NumReads" = "numeric"
                ),
                showProgress = FALSE
            )
            Counts[, j] <- Quant$NumReads
            TPM[, j] <- Quant$TPM
            Length[, j] <- Quant$EffectiveLength
        }

        if (NBoot > 0L) {
            BootFileCon <- gzcon(file(BootFile, open = "rb"))
            Boot <- readBin(
                BootFileCon,
                what = "double",
                n = NTx * NBoot
            )
            close(BootFileCon)
            dim(Boot) <- c(NTx, NBoot)
            M <- rowMeans(Boot)
            i <- (M > 0)
            OverDisp[i] <- OverDisp[i] +
                rowSums(
                    (Boot[i, ] -
                        M[i])^2
                ) /
                    M[i]
            DF[i] <- DF[i] + NBoot - 1L

            # Sum the Gibbs samples of the targets in each group. rowsum()
            # sorts the groups, so they are in the same order in every sample
            if (!is.null(groups)) {
                if (is.null(GroupOverDisp)) {
                    GroupOverDisp <- GroupDF <- vector("list", length(groups))
                    names(GroupOverDisp) <- names(GroupDF) <- names(groups)
                }
                for (g in names(groups)) {
                    keep <- !is.na(groups[[g]])
                    B <- rowsum(Boot[keep, , drop = FALSE], groups[[g]][keep])
                    if (is.null(GroupOverDisp[[g]])) {
                        GroupOverDisp[[g]] <- stats::setNames(
                            rep_len(0, nrow(B)),
                            rownames(B)
                        )
                        GroupDF[[g]] <- stats::setNames(
                            rep_len(0L, nrow(B)),
                            rownames(B)
                        )
                    }
                    MB <- rowMeans(B)
                    ib <- (MB > 0)
                    GroupOverDisp[[g]][ib] <- GroupOverDisp[[g]][ib] +
                        rowSums((B[ib, , drop = FALSE] - MB[ib])^2) / MB[ib]
                    GroupDF[[g]][ib] <- GroupDF[[g]][ib] + NBoot - 1L
                }
            }
        }
    }

    OverDisp <- .moderateOverdispersion(OverDisp, DF)
    OverDispPrior <- attr(OverDisp, "prior")
    attr(OverDisp, "prior") <- NULL

    GroupOverDispersion <- NULL
    if (!is.null(groups)) {
        GroupOverDispersion <- lapply(names(groups), function(g) {
            if (is.null(GroupOverDisp[[g]])) {
                ids <- sort(unique(groups[[g]][!is.na(groups[[g]])]))
                od <- stats::setNames(rep_len(0, length(ids)), ids)
                df <- stats::setNames(rep_len(0L, length(ids)), ids)
            } else {
                od <- GroupOverDisp[[g]]
                df <- GroupDF[[g]]
            }
            .moderateOverdispersion(od, df)
        })
        names(GroupOverDispersion) <- names(groups)
    }

    Quant1 <- as.data.frame(Quant1, stringsAsFactors = FALSE)
    dimnames(Counts) <- list(Quant1$Name, paths)
    dimnames(TPM) <- list(Quant1$Name, paths)
    dimnames(Length) <- list(Quant1$Name, paths)
    row.names(Quant1) <- Quant1$Name
    Quant1$Name <- NULL
    Quant1$TPM <- Quant1$NumReads <- Quant1$EffectiveLength <- NULL
    Quant1$Overdispersion <- OverDisp

    list(
        counts = Counts,
        tpm = TPM,
        length = Length,
        annotation = Quant1,
        overdispersion.prior = OverDispPrior,
        resample.type = ResampleType,
        group.overdispersion = GroupOverDispersion
    )
}


#' Moderate overdispersion estimates
#'
#' Takes the sums of the scaled squared deviations of the inferential
#' replicates (Gibbs samples) about their means, and their degrees of freedom,
#' accumulated over samples, and returns overdispersion estimates moderated
#' towards a common prior, as in \code{edgeR::catchSalmon()}. Features without
#' any degrees of freedom are given the prior.
#'
#' @param od_sum Numeric vector of accumulated sums of squares divided by means
#' @param df Integer vector of accumulated degrees of freedom
#'
#' @return Numeric vector of overdispersions (keeping the names of
#' \code{od_sum}) with the prior as attribute "prior".
#' @keywords internal
.moderateOverdispersion <- function(od_sum, df) {
    od <- od_sum
    i <- (df > 0L)
    if (sum(i) > 0L) {
        od[i] <- od[i] / df[i]
        df_median <- stats::median(df[i])
        df_prior <- 3
        prior <- stats::median(od[i]) /
            stats::qf(0.5, df1 = df_median, df2 = df_prior)
        if (prior < 1) {
            prior <- 1
        }
        od[i] <- (df_prior * prior + df[i] * od[i]) / (df_prior + df[i])
        od <- pmax(od, 1)
        od[!i] <- prior
    } else {
        od[] <- NA_real_
        prior <- NA_real_
    }
    attr(od, "prior") <- prior
    od
}
