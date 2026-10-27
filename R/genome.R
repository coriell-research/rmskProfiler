#' BSgenome packages providing the default genome for each species
#' @keywords internal
.BSGENOME <- c(
    Hs = "BSgenome.Hsapiens.UCSC.hg38",
    Mm = "BSgenome.Mmusculus.UCSC.mm10"
)


#' Indexed copies of genome FASTA files made in this session, by file
#' @keywords internal
.GENOME_COPIES <- new.env(parent = emptyenv())


#' Chromosomes used for transcripts (with chrM) and decoys (without chrM)
#' @keywords internal
.PRIMARY_CHROMS <- "^chr([0-9]+|[XYM])$"


#' Get the genome sequence used to build the index
#'
#' @param species Either "Hs" or "Mm"
#' @param genome NULL to use the BSgenome package for the species, a BSgenome,
#' a DNAStringSet, an Rsamtools FaFile, or a path to a (gzipped) FASTA file. A
#' genome returned by this function is returned as is.
#'
#' @return An object of class "rmskGenome", a list with elements \code{seqs},
#' a BSgenome or FaFile supporting \code{Biostrings::getSeq()}, and
#' \code{label}, a description of the genome recorded in the cache metadata.
#' The label of a BSgenome is its package name and version. Other genomes are
#' labelled "custom" followed by a key computed from their sequence names and
#' lengths, so the same genome gets the same label whatever the input type.
#' @keywords internal
.getGenome <- function(species, genome = NULL) {
    if (inherits(genome, "rmskGenome")) {
        return(genome)
    }

    label <- NULL
    if (is.null(genome)) {
        pkg <- .BSGENOME[[species]]
        if (!.isInstalled(pkg)) {
            stop(
                pkg,
                " is required to use the default genome for species '",
                species,
                "'. Install it with BiocManager::install(\"",
                pkg,
                "\") or provide a genome with the genome argument.",
                call. = FALSE
            )
        }
        genome <- getExportedValue(pkg, pkg)
    }

    if (methods::is(genome, "BSgenome")) {
        pkg <- methods::slot(genome, "pkgname")
        label <- paste(pkg, utils::packageVersion(pkg))
    } else if (is.character(genome) && length(genome) == 1L) {
        if (!file.exists(genome)) {
            stop("Genome FASTA file not found: ", genome, call. = FALSE)
        }
        # Use an indexed FASTA directly, otherwise read it in and index a copy
        # once per session
        if (file.exists(paste0(genome, ".fai"))) {
            genome <- Rsamtools::FaFile(genome)
        } else {
            info <- file.info(genome)
            key <- paste(normalizePath(genome), info$size, info$mtime)
            copy <- .GENOME_COPIES[[key]]
            if (is.null(copy) || !file.exists(paste0(copy, ".fai"))) {
                copy <- .indexedCopy(Biostrings::readDNAStringSet(genome))
                assign(key, copy, envir = .GENOME_COPIES)
            }
            genome <- Rsamtools::FaFile(copy)
        }
    }

    if (methods::is(genome, "DNAStringSet")) {
        genome <- Rsamtools::FaFile(.indexedCopy(genome))
    }

    if (methods::is(genome, "FaFile")) {
        if (!file.exists(paste0(BiocGenerics::path(genome), ".fai"))) {
            Rsamtools::indexFa(genome)
        }
    } else if (!methods::is(genome, "BSgenome")) {
        stop(
            "genome must be NULL, a BSgenome, a DNAStringSet, an FaFile, or ",
            "the path to a FASTA file",
            call. = FALSE
        )
    }

    genome <- structure(
        list(seqs = genome, label = label),
        class = "rmskGenome"
    )
    if (is.null(genome$label)) {
        lengths <- .genomeLengths(genome)
        key <- paste(names(lengths), lengths, sep = ":", collapse = ",")
        genome$label <- paste0(
            "custom:",
            substr(unname(tools::md5sum(bytes = charToRaw(key))), 1, 10)
        )
    }
    genome
}


#' Write an indexed FASTA copy of a DNAStringSet to a temporary file
#'
#' @param seqs DNAStringSet. Names are trimmed at the first space.
#'
#' @return Path to the FASTA file
#' @keywords internal
.indexedCopy <- function(seqs) {
    names(seqs) <- sub(" .*", "", names(seqs))
    fa <- tempfile("genome", fileext = ".fa")
    Biostrings::writeXStringSet(seqs, fa)
    Rsamtools::indexFa(fa)
    fa
}


#' Check that a package is installed
#'
#' @param pkg Name of the package
#'
#' @return TRUE/FALSE
#' @keywords internal
.isInstalled <- function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
}


#' Lengths of the sequences in a genome
#'
#' @param genome Genome returned by \code{.getGenome()}
#'
#' @return Named integer vector
#' @keywords internal
.genomeLengths <- function(genome) {
    GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(genome$seqs))
}


#' Extract the sequences of ranges from a genome
#'
#' Sequences of ranges on the minus strand are reverse complemented.
#'
#' @param genome Genome returned by \code{.getGenome()}
#' @param gr GRanges contained in the genome
#'
#' @return Unnamed DNAStringSet
#' @keywords internal
.genomeSeqs <- function(genome, gr) {
    seqs <- Biostrings::getSeq(genome$seqs, gr)
    names(seqs) <- NULL
    seqs
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
#' @param genome Genome returned by \code{.getGenome()}
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
        names(.genomeLengths(genome))
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
            "the genome uses the same chromosome names as the GTF ",
            "(e.g. 'chr1').",
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

    GenomicFeatures::extractTranscriptSeqs(genome$seqs, exons)
}
