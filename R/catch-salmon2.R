#' Import transcript-level counts with offset
#'
#' @details
#' This function is lifted from:
#'
#' https://github.com/plbaldoni/TranscriptDE-code/blob/main/code/pkg/R/utils.R#L52
#'
#' It is analogous to catchSalmon but imports TPM values as well
#' @param character vector giving paths to the sample-specific directories created by a kallisto or
#' Salmon. Each entry corresponds to one RNA-seq sample.
#' @param logical. If TRUE, progress information will be sent to standard output as each sample is processed.
#' @keywords internal
.catchSalmon2 <- function(paths, verbose = TRUE) {
  NSamples <- length(paths)

  OK <- requireNamespace("jsonlite", quietly = TRUE)
  if (!OK) {
    stop("jsonlite package required but is not installed (or can't be loaded)")
  }

  OK <- requireNamespace("readr", quietly = TRUE)
  if (!OK) {
    stop("readr package required but is not installed (or can't be loaded)")
  }

  ResampleType <- rep_len("", NSamples)
  for (j in 1L:NSamples) {
    if (verbose) {
      cat("Reading ", paths[j], ", ", sep = "")
    }
    MetaFile <- file.path(paths[j], "aux_info", "meta_info.json")
    QuantFile <- file.path(paths[j], "quant.sf")
    BootFile <- file.path(paths[j], "aux_info", "bootstrap", "bootstraps.gz")

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
    if (is.null(ResampleType)) {
      Type <- "bootstrap"
    } else {
      ResampleType[j] <- Type
    }

    if (verbose) {
      cat(NTx, "transcripts,", NBoot, Type, "samples\n")
    }

    if (j == 1L) {
      Counts <- matrix(0, NTx, NSamples)
      TPM <- matrix(0, NTx, NSamples)
      DF <- rep_len(0L, NTx)
      OverDisp <- rep_len(0, NTx)
      Quant1 <- suppressWarnings(readr::read_tsv(
        QuantFile,
        col_types = "cdd_d",
        progress = FALSE
      ))
      Counts[, 1L] <- Quant1$NumReads
      TPM[, 1L] <- Quant1$TPM
    } else {
      Quant <- suppressWarnings(readr::read_tsv(
        QuantFile,
        col_types = "____d",
        progress = FALSE
      ))
      Counts[, j] <- Quant$NumReads
      TPM[, j] <- Quant$TPM
    }

    if (NBoot > 0L) {
      BootFileCon <- gzcon(file(BootFile, open = "rb"))
      Boot <- readBin(
        BootFileCon,
        what = "double",
        n = NTx *
          NBoot
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
    }
  }

  i <- (DF > 0L)
  if (sum(i) > 0L) {
    OverDisp[i] <- OverDisp[i] / DF[i]
    DFMedian <- median(DF[i])
    DFPrior <- 3
    OverDispPrior <- median(OverDisp[i]) /
      qf(0.5, df1 = DFMedian, df2 = DFPrior)
    if (OverDispPrior < 1) {
      OverDispPrior <- 1
    }
    OverDisp[i] <- (DFPrior * OverDispPrior + DF[i] * OverDisp[i]) /
      (DFPrior +
        DF[i])
    OverDisp <- pmax(OverDisp, 1)
    OverDisp[!i] <- OverDispPrior
  } else {
    OverDisp[] <- NA_real_
    OverDispPrior <- NA_real_
  }

  Quant1 <- as.data.frame(Quant1, stringsAsFactors = FALSE)
  dimnames(Counts) <- list(Quant1$Name, paths)
  dimnames(TPM) <- list(Quant1$Name, paths)
  row.names(Quant1) <- Quant1$Name
  Quant1$Name <- NULL
  Quant1$TPM <- Quant1$NumReads <- NULL
  Quant1$Overdispersion <- OverDisp

  list(
    counts = Counts,
    tpm = TPM,
    annotation = Quant1,
    overdispersion.prior = OverDispPrior,
    resample.type = ResampleType
  )
}
