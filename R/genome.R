#' Get the BSgenome of a species
#'
#' @param species Either "Hs" or "Mm"
#'
#' @return BSgenome.Hsapiens.UCSC.hg38 or BSgenome.Mmusculus.UCSC.mm10
#' @keywords internal
.getGenome <- function(species) {
    if (species == "Hs") {
        BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    } else {
        BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    }
}


#' Label of a genome recorded in the cache metadata
#'
#' @param genome BSgenome returned by \code{.getGenome()}
#'
#' @return The package name and version of the BSgenome, e.g.
#' "BSgenome.Hsapiens.UCSC.hg38 1.4.5"
#' @keywords internal
.genomeLabel <- function(genome) {
    pkg <- methods::slot(genome, "pkgname")
    paste(pkg, utils::packageVersion(pkg))
}


#' Names of the decoy chromosomes
#'
#' The primary chromosomes of the genome (e.g. chr1-22, X, and Y for human)
#' are used as the decoys.
#'
#' @param genome BSgenome returned by \code{.getGenome()}
#'
#' @return Character vector of chromosome names in the order of the genome
#' @keywords internal
.decoys <- function(genome) {
    grep(
        "^chr([0-9]+|[XY])$",
        names(GenomeInfoDb::seqlengths(genome)),
        value = TRUE
    )
}


#' Match GENCODE sequence names to the names of a UCSC genome
#'
#' GENCODE names scaffolds, patches, and haplotypes by their GenBank accession
#' (e.g. "GL000009.2") while UCSC genomes encode the accession in the name
#' (e.g. "chr14_GL000009v2_random", or "chr1_GL455991_alt" without a version
#' for mm10).
#'
#' @param x Sequence names to match, e.g. the seqlevels of a TxDb
#' @param target Sequence names of the genome
#'
#' @return Character vector the same length as \code{x} with the matching name
#' in \code{target}, or NA if there is none
#' @keywords internal
.matchSeqnames <- function(x, target) {
    out <- target[match(x, target)]

    # Accessions encoded in UCSC names, "GL000009v2" -> "GL000009.2"
    scaffold <- grepl("_", target, fixed = TRUE)
    acc <- sub("^chr[^_]*_([^_]+).*$", "\\1", target[scaffold])
    acc <- sub("v([0-9]+)$", ".\\1", acc)
    todo <- is.na(out)
    out[todo] <- target[scaffold][match(x[todo], acc)]

    # UCSC names without a version match the accession without a version
    todo <- is.na(out)
    no_version <- !grepl(".", acc, fixed = TRUE)
    out[todo] <- target[scaffold][no_version][match(
        sub("\\.[0-9]+$", "", x[todo]),
        acc[no_version]
    )]

    out
}


#' Extract the sequences of the transcripts in a TxDb from a genome
#'
#' Sequence names of the TxDb are matched to the genome with
#' \code{.matchSeqnames()}. Transcripts with exons on sequences that are not in
#' the genome are dropped with a message.
#'
#' @param txdb TxDb
#' @param genome BSgenome returned by \code{.getGenome()}
#' @param tx_order Optional character vector of transcript ids giving the
#' order of the returned sequences. Transcripts not in \code{tx_order} are
#' placed after those that are.
#'
#' @return DNAStringSet named by transcript id
#' @keywords internal
.transcriptSeqs <- function(txdb, genome, tx_order = NULL) {
    exons <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
    seqmap <- .matchSeqnames(
        GenomeInfoDb::seqlevels(exons),
        names(GenomeInfoDb::seqlengths(genome))
    )
    names(seqmap) <- GenomeInfoDb::seqlevels(exons)
    exon_seq <- as.integer(GenomicRanges::seqnames(unlist(
        exons,
        use.names = FALSE
    )))
    missing <- unique(rep(names(exons), lengths(exons))[
        is.na(seqmap[exon_seq])
    ])
    if (length(missing) == length(exons)) {
        stop(
            "No transcripts found on the sequences of the genome. Check that ",
            "the GTF uses UCSC chromosome names (e.g. 'chr1').",
            call. = FALSE
        )
    }
    if (length(missing) > 0L) {
        message(
            length(missing),
            " transcripts on sequences not in the genome were skipped"
        )
        exons <- exons[!names(exons) %in% missing]
    }
    exons <- GenomeInfoDb::keepSeqlevels(
        exons,
        GenomeInfoDb::seqlevelsInUse(exons),
        pruning.mode = "coarse"
    )
    GenomeInfoDb::seqlevels(exons) <- unname(seqmap[
        GenomeInfoDb::seqlevels(exons)
    ])
    if (!is.null(tx_order)) {
        exons <- exons[order(match(names(exons), tx_order))]
    }

    GenomicFeatures::extractTranscriptSeqs(genome, exons)
}
