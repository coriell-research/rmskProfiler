#' Extract unique RepeatMasker sequences and their loci from the genome
#'
#' @param rmsk GRanges of filtered RepeatMasker records returned by
#' \code{rmskToBed()}
#' @param genome Genome returned by \code{.getGenome()}
#' @param chunk_size Number of sequences hashed at a time
#'
#' @return list with elements \code{seqs}, a DNAStringSet of unique sequences
#' named by their SHA1 hash, and \code{info}, a data.table mapping each hash to
#' the RepeatMasker loci with that sequence (0-based start)
#' @import data.table
#' @keywords internal
.extractUniqueRecords <- function(rmsk, genome, chunk_size = 1e5) {
    # One row per record as in a BED file: 0-based start and "." for records
    # without a strand
    strand <- as.character(GenomicRanges::strand(rmsk))
    bed <- data.table::data.table(
        seqnames = as.character(GenomicRanges::seqnames(rmsk)),
        start = GenomicRanges::start(rmsk) - 1L,
        end = GenomicRanges::end(rmsk),
        RepName = rmsk$name,
        strand = ifelse(strand == "*", ".", strand)
    )

    # Like bedtools getfasta, skip records not contained in the genome
    chrom_len <- .genomeLengths(genome)[bed$seqnames]
    in_bounds <- !is.na(chrom_len) & bed$end <= chrom_len
    if (!all(in_bounds)) {
        warning(
            sum(!in_bounds),
            " records not contained in the genome were skipped",
            call. = FALSE
        )
        bed <- bed[in_bounds]
    }

    # Unstranded records are extracted from the plus strand
    message("Extracting RepeatMasker sequences from the genome...")
    seqs <- .genomeSeqs(
        genome,
        GenomicRanges::GRanges(
            bed$seqnames,
            IRanges::IRanges(bed$start + 1L, bed$end),
            strand = ifelse(bed$strand == "-", "-", "+")
        )
    )

    message("Hashing RepeatMasker sequences and checking for duplicates...")
    sha1 <- digest::getVDigest("sha1")
    chunks <- split(seq_along(seqs), ceiling(seq_along(seqs) / chunk_size))
    bed[,
        Hash := unlist(
            lapply(chunks, function(i) {
                sha1(as.character(seqs[i]), serialize = FALSE)
            }),
            use.names = FALSE
        )
    ]

    first <- !duplicated(bed$Hash)
    uniq_seqs <- seqs[first]
    names(uniq_seqs) <- bed$Hash[first]

    info <- bed[, .(Hash, RepName, seqnames, strand, start, end)]

    return(list(seqs = uniq_seqs, info = info))
}


#' Extract the unique RepeatMasker sequences from the genome
#'
#' The sequence of each hash is extracted from its first locus, which gives the
#' same sequences in the same order as \code{.extractUniqueRecords()}.
#'
#' @param info data.table of duplicate information created by
#' \code{.extractUniqueRecords()}
#' @param genome Genome returned by \code{.getGenome()}
#'
#' @return DNAStringSet of unique sequences named by their SHA1 hash
#' @keywords internal
.uniqueSeqs <- function(info, genome) {
    first <- info[!duplicated(info$Hash)]
    seqs <- .genomeSeqs(
        genome,
        GenomicRanges::GRanges(
            first$seqnames,
            IRanges::IRanges(first$start + 1L, first$end),
            strand = ifelse(first$strand == "-", "-", "+")
        )
    )
    names(seqs) <- first$Hash
    seqs
}


#' Find the unique RepeatMasker sequences and their loci
#'
#' This function will extract the sequences of the RepeatMasker records
#' returned by \code{rmskToBed()} from the genome, hash them, and save a
#' data.table, rmsk-duplicateInfo.rds, to the rmskProfiler cache which maps the
#' SHA1 hash of each sequence to the RepeatMasker loci with that sequence. Loci
#' with identical sequences share a hash and are quantified as a single target.
#' The unique sequences themselves are not saved: \code{createGentrome()}
#' extracts them from the genome using the first locus of each hash.
#'
#' The sequences are extracted from the BSgenome package for the species,
#' BSgenome.Hsapiens.UCSC.hg38 (GRCh38) or BSgenome.Mmusculus.UCSC.mm10
#' (GRCm38), unless another genome is given with \code{genome}.
#'
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus)
#' @param exclude A character vector specifying which elements (repClass) to
#' exclude. Default "Simple_repeat", "Low_complexity", "Satellite", "RNA",
#' "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", and "Unknown".
#' @param min_len Minimum sequence length of a record. Default 32.
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#' @param genome NULL (default) to use the BSgenome package for the species, or
#' a BSgenome, DNAStringSet, Rsamtools FaFile, or path to a (gzipped) FASTA
#' file. The sequence names must match the RepeatMasker records and GTF
#' (e.g. "chr1").
#' @param rmsk NULL (default) to use the AnnotationHub RepeatMasker annotations
#' for \code{species}, or a GRanges of RepeatMasker records with metadata
#' columns repName, repFamily, and repClass, such as the GRanges returned by
#' \code{rmskToBed()}. The records are filtered with \code{rmskToBed()} using
#' \code{exclude} and \code{min_len}.
#'
#' @return Invisibly \code{NULL}. Called for its side effect of writing the
#' duplicate information (rmsk-duplicateInfo.rds) to the cache.
#' @export
#'
#' @examples
#' bfc <- exampleCache()
#' rmsk <- readRDS(system.file("extdata", "rmsk-chr22.rds",
#'     package = "rmskProfiler"
#' ))
#' extractUniqueSeqs(
#'     species = "Hs",
#'     cache = bfc,
#'     genome = exampleGenome(),
#'     rmsk = rmsk
#' )
#'
#' # Loci with identical sequences share a hash
#' info <- readRDS(
#'     getResource("rmsk-duplicateInfo.rds", species = "Hs", cache = bfc)
#' )
#' head(info)
extractUniqueSeqs <- function(
    species = c("Hs", "Mm"),
    exclude = .DEFAULT_EXCLUDE,
    min_len = 32,
    cache = NULL,
    genome = NULL,
    rmsk = NULL
) {
    species <- match.arg(species)

    bfc <- .getCache(cache)
    genome <- .getGenome(species, genome)
    message("Filtering RepeatMasker records...")
    rmsk <- rmskToBed(
        species = species,
        exclude = exclude,
        min_len = min_len,
        rmsk = rmsk
    )

    message("Using genome: ", genome$label)
    res <- .extractUniqueRecords(rmsk, genome)

    out_info <- .newResource(
        bfc,
        species,
        "rmsk-duplicateInfo.rds",
        ".rds",
        exclude,
        min_len,
        genome = genome$label
    )
    message("Writing duplicate information to: ", out_info)
    saveRDS(res$info, out_info)

    return(invisible(NULL))
}
