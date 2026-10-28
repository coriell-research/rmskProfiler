# Create the example data in inst/extdata
#
# The example data are a 1 Mb region of hg38 chr22 (chr22:37,000,000-37,999,999)
# with the GENCODE v48 annotations of the genes contained in it, the UCSC
# RepeatMasker records in it, and Salmon quants of simulated reads against the
# rmskProfiler index built from these files. The GENCODE transcript sequences
# of these genes are written to tests/testthat/fixtures.
#
# Requirements:
# - the GENCODE v48 GTF in the rmskProfiler cache, i.e.
#   downloadResources(species = "Hs")
# - BSgenome.Hsapiens.UCSC.hg38
# - AnnotationHub access to AH111333 (UCSC RepeatMasker, hg38)
# - Salmon (>= 1.10.0) on PATH
#
# Run from the package root with:
#   Rscript inst/scripts/make-example-data.R

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(Biostrings)
})
devtools::load_all(".", quiet = TRUE)
set.seed(20261027)

out_dir <- file.path("inst", "extdata")
region <- GRanges("chr22", IRanges(37e6, 38e6 - 1))

# Genome ----------------------------------------------------------------------
# chr22 up to the end of the region with every base before the region masked,
# used in place of the BSgenome to build the example index with chr22 as the
# only decoy. Coordinates are unchanged so the GTF and RepeatMasker records need
# no shifting
message("Creating genome...")
genome <- .getGenome("Hs")
genome_label <- .genomeLabel(genome)
chr22 <- genome[["chr22"]]
masked <- DNAStringSet(paste0(
    strrep("N", start(region) - 1L),
    as.character(subseq(chr22, start(region), end(region)))
))
names(masked) <- "chr22"
chr22 <- masked[[1]]

# GTF -------------------------------------------------------------------------
message("Creating GTF...")
gtf <- rtracklayer::import(getResource("gtf", species = "Hs"), which = region)
genes <- subsetByOverlaps(gtf[gtf$type == "gene"], region, type = "within")
gtf <- gtf[gtf$gene_id %in% genes$gene_id]
gtf_file <- tempfile(fileext = ".gtf")
rtracklayer::export(gtf, gtf_file, format = "gtf")
con <- gzfile(file.path(out_dir, "gencode-chr22.gtf.gz"), "w")
writeLines(readLines(gtf_file), con)
close(con)

# Transcripts -----------------------------------------------------------------
# Extracted from the genome as in createGentrome(), in the order of the GTF
message("Creating transcripts...")
txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(
    file.path(out_dir, "gencode-chr22.gtf.gz"),
    format = "gtf"
))
tx_ids <- unique(gtf$transcript_id[gtf$type == "transcript"])
tx_seqs <- .transcriptSeqs(txdb, genome, tx_order = tx_ids)
writeXStringSet(
    tx_seqs,
    file.path("tests", "testthat", "fixtures", "transcripts-chr22.fa.gz"),
    compress = TRUE,
    width = 60
)

# RepeatMasker ----------------------------------------------------------------
message("Creating RepeatMasker records...")
ah <- AnnotationHub::AnnotationHub()
rmsk <- subsetByOverlaps(ah[["AH111333"]], region, type = "within")
GenomeInfoDb::seqlevels(rmsk) <- "chr22"
saveRDS(rmsk, file.path(out_dir, "rmsk-chr22.rds"), compress = "xz")

# Index -----------------------------------------------------------------------
message("Building the rmskProfiler index...")
work_dir <- tempfile("rmskProfiler-example-")
dir.create(work_dir)
index_dir <- file.path(work_dir, "index")
bfc <- exampleCache()
testthat::with_mocked_bindings(
    generateIndex(
        index_dir = index_dir,
        species = "Hs",
        check_integrity = FALSE,
        threads = 4,
        cache = bfc,
        rmsk = rmsk
    ),
    .getGenome = function(species) masked,
    .genomeLabel = function(genome) genome_label,
    .package = "rmskProfiler"
)

# Simulate reads --------------------------------------------------------------
# Each sample has reads from mature transcripts, from pre-mRNA (gene bodies plus
# a run-on downstream of the 3' end, at a rate that varies between samples) so
# that intronic and downstream TEs receive read-through, and from a set of TE
# loci that are autonomously expressed in the "trt" group only
message("Simulating reads...")
rd <- readRDS(getResource("rmsk-rowData.rds", species = "Hs", cache = bfc))
te_loci <- unlist(rd$Ranges[!is.na(rd$Hash)])

gene_span <- genes
run_on <- pmin(round(stats::rexp(length(gene_span), 1 / 5000)), 15000L)
gene_span <- resize(
    gene_span,
    width(gene_span) + run_on,
    fix = "start"
)
gene_span <- trim(gene_span)
gene_span <- gene_span[end(gene_span) <= end(region)]
pre_mrna <- DNAStringSet(Views(chr22, ranges(gene_span)))
minus <- as.logical(strand(gene_span) == "-")
pre_mrna[minus] <- reverseComplement(pre_mrna[minus])
names(pre_mrna) <- gene_span$gene_id

auto_te <- te_loci[sample(
    which(
        !overlapsAny(te_loci, genes, ignore.strand = TRUE) &
            width(te_loci) >= 300
    ),
    20
)]
auto_seqs <- DNAStringSet(Views(chr22, ranges(auto_te)))

tx_gene <- gtf$gene_id[match(
    sub("\\|.*", "", names(tx_seqs)),
    gtf$transcript_id
)]
tx_level <- stats::rlnorm(length(tx_seqs), meanlog = 2, sdlog = 1.5)
de <- sample(length(tx_seqs), round(0.1 * length(tx_seqs)))
gene_level <- tapply(tx_level, tx_gene, sum)[names(pre_mrna)]
gene_level[is.na(gene_level)] <- 0

simulate_fragments <- function(seqs, n) {
    keep <- width(seqs) >= 150L & n > 0
    seqs <- seqs[keep]
    n <- n[keep]
    idx <- rep(seq_along(seqs), n)
    flen <- pmin(
        pmax(round(stats::rnorm(length(idx), 250, 25)), 150L),
        width(seqs)[idx]
    )
    start <- floor(stats::runif(length(idx)) * (width(seqs)[idx] - flen + 1)) +
        1
    subseq(seqs[idx], start = start, width = flen)
}

write_reads <- function(frags, prefix) {
    r1 <- subseq(frags, start = 1L, width = 100L)
    r2 <- reverseComplement(subseq(frags, end = width(frags), width = 100L))
    flip <- stats::runif(length(frags)) < 0.5
    tmp <- r1[flip]
    r1[flip] <- r2[flip]
    r2[flip] <- tmp
    ids <- BStringSet(sprintf("read%d", seq_along(frags)))
    qual <- BStringSet(rep(strrep("I", 100L), length(frags)))
    for (mate in c("1", "2")) {
        reads <- if (mate == "1") r1 else r2
        ShortRead::writeFastq(
            ShortRead::ShortReadQ(reads, qual, ids),
            paste0(prefix, "_R", mate, ".fq.gz"),
            compress = TRUE
        )
    }
}

sample_names <- c(paste0("ctrl_", 1:3), paste0("trt_", 1:3))
n_pairs <- 5e5
for (s in sample_names) {
    is_trt <- startsWith(s, "trt")
    depth <- stats::runif(1, 0.8, 1.2)
    rt_rate <- stats::runif(1, 0.05, 0.3)

    level <- tx_level * stats::rlnorm(length(tx_level), 0, 0.2)
    if (is_trt) {
        level[de] <- level[de] * 3
    }
    tx_w <- level * width(tx_seqs)
    pre_w <- rt_rate * gene_level * width(pre_mrna) / 10
    auto_w <- if (is_trt) {
        rep(200, length(auto_seqs))
    } else {
        rep(0, length(auto_seqs))
    }
    auto_w <- auto_w * width(auto_seqs)

    total <- round(n_pairs * depth)
    w <- c(tx_w, pre_w, auto_w)
    n <- as.vector(stats::rmultinom(1, total, w / sum(w)))
    n_tx <- n[seq_along(tx_seqs)]
    n_pre <- n[length(tx_seqs) + seq_along(pre_mrna)]
    n_auto <- n[length(tx_seqs) + length(pre_mrna) + seq_along(auto_seqs)]

    frags <- c(
        simulate_fragments(tx_seqs, n_tx),
        simulate_fragments(pre_mrna, n_pre),
        simulate_fragments(auto_seqs, n_auto)
    )
    frags <- frags[sample(length(frags))]
    write_reads(frags, file.path(work_dir, s))
}

# Keep the first 2000 read pairs of one sample for the salmonQuant() example
dir.create(file.path(out_dir, "fastq"), showWarnings = FALSE)
for (mate in c("1", "2")) {
    reads <- ShortRead::readFastq(
        file.path(work_dir, paste0("ctrl_1_R", mate, ".fq.gz"))
    )
    out_fq <- file.path(out_dir, "fastq", paste0("ctrl_1_R", mate, ".fq.gz"))
    unlink(out_fq)
    ShortRead::writeFastq(reads[seq_len(2000L)], out_fq, compress = TRUE)
}

# Quantify --------------------------------------------------------------------
message("Quantifying with Salmon...")
quant_dir <- file.path(work_dir, "quants")
salmonQuant(
    fq1 = file.path(work_dir, paste0(sample_names, "_R1.fq.gz")),
    fq2 = file.path(work_dir, paste0(sample_names, "_R2.fq.gz")),
    sample_names = sample_names,
    index_dir = index_dir,
    out_dir = quant_dir,
    n_gibbs = 10,
    "--threads 4"
)

# Keep only the files needed by importQuants()
keep <- c(
    "quant.sf",
    "cmd_info.json",
    "lib_format_counts.json",
    file.path("aux_info", "meta_info.json"),
    file.path("aux_info", "bootstrap", "bootstraps.gz")
)
unlink(file.path(out_dir, "quants"), recursive = TRUE)
for (s in sample_names) {
    from <- file.path(quant_dir, paste0(s, "_quants"))
    to <- file.path(out_dir, "quants", s)
    for (f in keep) {
        dir.create(
            dirname(file.path(to, f)),
            recursive = TRUE,
            showWarnings = FALSE
        )
        file.copy(file.path(from, f), file.path(to, f))
    }
}

message("Done.")
