#' Load the TxDb for a species from the cache, creating it if needed
#'
#' @param species Either "Hs" or "Mm"
#' @param bfc BiocFileCache
#'
#' @return TxDb
#' @keywords internal
.getTxDb <- function(species, bfc) {
    dbfile <- .findResource(bfc, .resourceRname(species, "annotation.txdb"))
    if (!is.null(dbfile)) {
        message("Loading TxDb...")
        # loadDb() loads packages that print startup messages
        return(suppressPackageStartupMessages(AnnotationDbi::loadDb(dbfile)))
    }

    gtf_file <- .getResource(
        bfc,
        species,
        "gtf",
        hint = "Run downloadResources() first."
    )
    organism <- "Homo sapiens"
    taxid <- 9606
    if (species == "Mm") {
        organism <- "Mus musculus"
        taxid <- 10090
    }

    message("Creating TxDb from GTF...")
    # makeTxDbFromGFF() warns about GENCODE records it does not import
    txdb <- suppressWarnings(
        txdbmaker::makeTxDbFromGFF(
            file = gtf_file,
            format = "gtf",
            organism = organism,
            taxonomyId = taxid,
            dataSource = .gencodeResources(species)["gtf", "url"]
        )
    )
    AnnotationDbi::saveDb(
        txdb,
        .newResource(bfc, species, "annotation.txdb", ".txdb")
    )

    return(txdb)
}


#' Find overlaps of TE loci regions with genomic features
#'
#' Finds any TE that overlaps with an exon by transcript, and intron by
#' transcript, a 3 prime UTR by transcript, a 5 prime UTR by transcript, or
#' a promoter by gene.
#'
#' @param x A GRanges object of TE loci
#' @param txdb TxDb to create annotation from
#'
#' @return List of hash vectors overlapping genomic features
.getHashOverlaps <- function(x, txdb) {
    message("Preparing reduced genomic features...")

    gr_exons <- unlist(GenomicFeatures::exonsBy(txdb, by = "gene"))
    GenomicRanges::mcols(gr_exons)$feature_id <- names(gr_exons)
    GenomicRanges::mcols(gr_exons)$type <- "exon"

    gr_promoters <- GenomicRanges::promoters(GenomicFeatures::genes(txdb))
    GenomicRanges::mcols(gr_promoters)$feature_id <- names(gr_promoters)
    GenomicRanges::mcols(gr_promoters)$type <- "promoter"

    # Create a Transcript-to-Gene map to standardize Introns/UTRs
    message("Mapping transcripts to genes...")
    # select() reports how keys map to columns, e.g. "1:1 mapping"
    tx2gene <- suppressMessages(
        AnnotationDbi::select(
            txdb,
            keys = AnnotationDbi::keys(txdb, "TXNAME"),
            columns = "GENEID",
            keytype = "TXNAME"
        )
    )
    tx2gene_map <- stats::setNames(tx2gene$GENEID, tx2gene$TXNAME)

    gr_introns <- unlist(GenomicFeatures::intronsByTranscript(
        txdb,
        use.names = TRUE
    ))
    GenomicRanges::mcols(gr_introns)$feature_id <- tx2gene_map[names(
        gr_introns
    )]
    GenomicRanges::mcols(gr_introns)$type <- "intron"

    gr_3utr <- unlist(GenomicFeatures::threeUTRsByTranscript(
        txdb,
        use.names = TRUE
    ))
    GenomicRanges::mcols(gr_3utr)$feature_id <- tx2gene_map[names(gr_3utr)]
    GenomicRanges::mcols(gr_3utr)$type <- "3utr"

    gr_5utr <- unlist(GenomicFeatures::fiveUTRsByTranscript(
        txdb,
        use.names = TRUE
    ))
    GenomicRanges::mcols(gr_5utr)$feature_id <- tx2gene_map[names(gr_5utr)]
    GenomicRanges::mcols(gr_5utr)$type <- "5utr"

    cols_to_keep <- c("feature_id", "type")
    GenomicRanges::mcols(gr_exons) <- GenomicRanges::mcols(gr_exons)[,
        cols_to_keep
    ]
    GenomicRanges::mcols(gr_promoters) <- GenomicRanges::mcols(gr_promoters)[,
        cols_to_keep
    ]
    GenomicRanges::mcols(gr_introns) <- GenomicRanges::mcols(gr_introns)[,
        cols_to_keep
    ]
    GenomicRanges::mcols(gr_3utr) <- GenomicRanges::mcols(gr_3utr)[,
        cols_to_keep
    ]
    GenomicRanges::mcols(gr_5utr) <- GenomicRanges::mcols(gr_5utr)[,
        cols_to_keep
    ]

    gr_features <- c(gr_exons, gr_introns, gr_promoters, gr_3utr, gr_5utr)

    message("Finding overlaps...")
    hits <- GenomicRanges::findOverlaps(x, gr_features, ignore.strand = TRUE)

    dt_hits <- data.table(
        query_idx = S4Vectors::queryHits(hits),
        subject_idx = S4Vectors::subjectHits(hits)
    )

    dt_hits[, Hash := GenomicRanges::mcols(x)$Hash[query_idx]]
    dt_hits[, q_strand := as.character(GenomicRanges::strand(x)[query_idx])]
    dt_hits[, feature := GenomicRanges::mcols(gr_features)$type[subject_idx]]
    dt_hits[,
        feature_id := GenomicRanges::mcols(gr_features)$feature_id[subject_idx]
    ]
    dt_hits[,
        s_strand := as.character(GenomicRanges::strand(gr_features)[
            subject_idx
        ])
    ]

    dt_hits[, strand_match := (q_strand == s_strand)]

    dt_hits[,
        feature_label := data.table::fifelse(
            is.na(feature_id),
            NA_character_,
            paste0(feature, ":", feature_id)
        )
    ]

    # Aggregate features by Hash
    mapping_dt <- dt_hits[,
        .(
            overlapping_features = stringi::stri_flatten(
                unique(feature_label[!is.na(feature_label)]),
                collapse = ";"
            )
        ),
        by = Hash
    ]

    getHashes <- function(dt, feat, stranded = FALSE) {
        if (isTRUE(stranded)) {
            unique(dt[feature == feat & strand_match == TRUE, Hash])
        } else {
            unique(dt[feature == feat, Hash])
        }
    }

    result <- list(
        hash_in_exon = getHashes(dt_hits, "exon", TRUE),
        hash_in_intron = getHashes(dt_hits, "intron", TRUE),
        hash_in_promoter = getHashes(dt_hits, "promoter", TRUE),
        hash_in_3utr = getHashes(dt_hits, "3utr", TRUE),
        hash_in_5utr = getHashes(dt_hits, "5utr", TRUE),

        u_hash_in_exon = getHashes(dt_hits, "exon", FALSE),
        u_hash_in_intron = getHashes(dt_hits, "intron", FALSE),
        u_hash_in_promoter = getHashes(dt_hits, "promoter", FALSE),
        u_hash_in_3utr = getHashes(dt_hits, "3utr", FALSE),
        u_hash_in_5utr = getHashes(dt_hits, "5utr", FALSE),

        mapping_dt = mapping_dt
    )

    return(result)
}


#' Generate annotation resources for annotating TE loci and genes
#'
#' This function produces multiple annotation files used to annotate downstream
#' counts generated after quantification with Salmon. It creates a GRangesList
#' object for each hash and a tab-delimited file containing annotation
#' information for each hash with indicator values for the types of genic
#' features they overlap. Overlap annotations are generated with respect to the
#' downloaded gencode GTF file.
#'
#' @param species Either "Hs" (Homo sapiens) or "Mm" (Mus musculus)
#' @param exclude A character vector specifying which elements (repClass) were
#' excluded by \code{rmskToBed()}. Default "Simple_repeat", "Low_complexity",
#' "Satellite", "RNA", "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA", and
#' "Unknown".
#' @param min_len Minimum sequence length of a record used by
#' \code{rmskToBed()}.
#' Default 32.
#' @param cache NULL, a path to a cache directory, or a BiocFileCache object.
#' Default NULL uses the rmskProfiler cache at
#' \code{tools::R_user_dir("rmskProfiler", which = "cache")}.
#'
#' @return Invisibly \code{NULL}. Called for its side effect of writing the
#' annotation (rmsk-rowData.rds) to the cache.
#' @import data.table
#' @export
#'
#' @examples
#' # Build the annotation for the example data
#' bfc <- exampleCache()
#' rmsk <- readRDS(system.file("extdata", "rmsk-chr22.rds",
#'     package = "rmskProfiler"
#' ))
#' rmskToBed(species = "Hs", cache = bfc, rmsk = rmsk)
#' extractUniqueSeqs(species = "Hs", cache = bfc)
#' createAnnotation(species = "Hs", cache = bfc)
#'
#' rd <- readRDS(getResource("rmsk-rowData.rds", species = "Hs", cache = bfc))
#' rd[!is.na(rd$Hash), c("Subfamily", "N_Loci", "hasIntronic")]
createAnnotation <- function(
    species = c("Hs", "Mm"),
    exclude = .DEFAULT_EXCLUDE,
    min_len = 32,
    cache = NULL
) {
    species <- match.arg(species)
    bfc <- .getCache(cache)
    info_rds <- .getResource(
        bfc,
        species,
        "rmsk-duplicateInfo.rds",
        exclude,
        min_len,
        hint = paste(
            "Run extractUniqueSeqs() with the same species, exclude, and",
            "min_len first."
        )
    )
    gtf_file <- .getResource(
        bfc,
        species,
        "gtf",
        hint = "Run downloadResources() first."
    )

    message("Reading in the duplicate information...")
    dt <- data.table::setDT(readRDS(info_rds))

    message("Creating a GRangesList for all TE ranges...")
    gr <- GenomicRanges::makeGRangesFromDataFrame(
        dt,
        starts.in.df.are.0based = TRUE,
        keep.extra.columns = TRUE
    )
    grl <- S4Vectors::splitAsList(gr, gr$Hash)

    message("Computing overlaps of TE-loci with transcript annotations...")
    txdb <- .getTxDb(species, bfc)
    ov <- .getHashOverlaps(gr, txdb)

    message("Getting all unique hash-element pairs...")
    hash_dt <- dt[, .(N_Loci = .N), by = .(Hash, RepName)]
    hash_dt[,
        c("Subfamily", "Family", "Class") := data.table::tstrsplit(
            RepName,
            ":",
            fixed = TRUE
        )
    ]

    message("Collapsing hash-level information into rowData...")
    by_hash <- hash_dt[,
        .(
            Class = stringi::stri_flatten(unique(Class), collapse = ","),
            Family = stringi::stri_flatten(unique(Family), collapse = ","),
            Subfamily = stringi::stri_flatten(
                unique(Subfamily),
                collapse = ","
            ),
            N_Loci = sum(N_Loci),
            N_Class = length(unique(Class)),
            N_Family = length(unique(Family)),
            N_Subfamily = length(unique(Subfamily))
        ),
        by = Hash
    ]

    by_hash[, `:=`(
        hasPromoter = Hash %chin% ov$hash_in_promoter,
        hasExonic = Hash %chin% ov$hash_in_exon,
        hasIntronic = Hash %chin% ov$hash_in_intron,
        has3UTR = Hash %chin% ov$hash_in_3utr,
        has5UTR = Hash %chin% ov$hash_in_5utr,
        hasUnstrandedPromoter = Hash %chin% ov$u_hash_in_promoter,
        hasUnstrandedExonic = Hash %chin% ov$u_hash_in_exon,
        hasUnstrandedIntronic = Hash %chin% ov$u_hash_in_intron,
        hasUnstranded3UTR = Hash %chin% ov$u_hash_in_3utr,
        hasUnstranded5UTR = Hash %chin% ov$u_hash_in_5utr
    )][, `:=`(
        hasIntergenic = (!hasExonic & !hasIntronic & !has3UTR & !has5UTR),
        hasUnstrandedIntergenic = (!hasUnstrandedExonic &
            !hasUnstrandedIntronic &
            !hasUnstranded3UTR &
            !hasUnstranded5UTR)
    )]

    message("Merging exact feature overlaps...")
    by_hash <- merge(by_hash, ov$mapping_dt, by = "Hash", all.x = TRUE)
    by_hash[is.na(overlapping_features), overlapping_features := ""]

    message("Reading in transcript annotations...")
    gtf <- rtracklayer::import(gtf_file)
    tx <- gtf[gtf$type == "transcript", ]
    tx_dt <- data.table::as.data.table(data.frame(tx))[, .(
        transcript_id,
        gene_id,
        gene_name,
        gene_type
    )]
    names(tx) <- tx$transcript_id
    tx <- methods::as(tx, "GRangesList")
    rmsk_grl <- c(tx, grl)

    # Strip redundant metadata and names
    unlisted_rmsk <- unlist(rmsk_grl, use.names = FALSE)
    GenomicRanges::mcols(unlisted_rmsk) <- NULL
    rmsk_grl <- BiocGenerics::relist(unlisted_rmsk, rmsk_grl)

    # Combine annotation DataFrames
    rd <- data.table::rbindlist(list(tx_dt, by_hash), fill = TRUE)
    rd <- S4Vectors::DataFrame(rd)
    rownames(rd) <- c(tx_dt$transcript_id, by_hash$Hash)
    rd$Ranges <- rmsk_grl[rownames(rd)]

    rd_file <- .newResource(
        bfc,
        species,
        "rmsk-rowData.rds",
        ".rds",
        exclude,
        min_len
    )
    message("Writing out rowData to: ", rd_file)
    saveRDS(rd, rd_file)
    message("Done.")

    return(invisible(NULL))
}
