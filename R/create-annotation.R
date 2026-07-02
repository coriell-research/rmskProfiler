#' Extract duplicate info from json file as a data.table
#'
#' @param jsonfile JSON file containing information about rmsk sequences
#' @import data.table
#' @return data.table
#'
.dupInfoToDT <- function(jsonfile) {
  message(
    "Reading in the duplicate information from ",
    basename(jsonfile),
    "..."
  )
  info <- jsonlite::fromJSON(jsonfile)

  message("Creating data.table from json information...")
  dt <- data.table::data.table(Hash = names(info), Instance = info)

  message("Unnesting json information...")
  dt <- dt[, .(Instance = as.character(unlist(Instance))), by = Hash]

  message("Extracting repetitive element names from location strings...")
  dt[,
    c("RepName", "Location") := data.table::tstrsplit(
      Instance,
      "::",
      fixed = TRUE
    )
  ]

  message("Extracting position information...")
  dt[,
    c("seqnames", "position") := data.table::tstrsplit(
      Location,
      ":",
      fixed = TRUE
    )
  ]
  dt[,
    strand := data.table::fifelse(
      stringi::stri_detect(position, regex = "\\(\\+\\)$"),
      "+",
      "-"
    )
  ]
  dt[, position := stringi::stri_replace(position, "", regex = "\\(.\\)")][,
    c("start", "end") := data.table::tstrsplit(position, "-", fixed = TRUE)
  ][,
    `:=`(start = as.integer(start), end = as.integer(end))
  ]
  dt[, `:=`(Location = NULL, position = NULL, Instance = NULL)]

  return(dt)
}


#' Find overlaps of TE loci regions with genomic features
#'
#' Finds any TE that overlaps with an exon by transcript, and intron by
#' transcript, a 3 prime UTR by transcript, a 5 prime UTR by transcript, or
#' a promoter by gene.
#'
#' @param x A GRanges object of TE loci
#' @param gtffile Path to GTF file to create annotation from
#' @param resource_dir Path to the rmsk resource directory. TxDb will be saved here.
#'
#' @return List of hash vectors overlapping genomic features
.getHashOverlaps <- function(x, gtffile, resource_dir) {
  organism <- "Homo sapiens"
  taxid <- 9606
  data_source <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz"

  if (grepl("M25", gtffile)) {
    organism <- "Mus Musculus"
    taxid <- 10090
    data_source <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
  }

  dbfile <- gsub(".gtf.gz", ".txdb", gtffile, fixed = TRUE)

  if (!file.exists(dbfile)) {
    message("Creating TxDb from GTF...")
    txdb <- suppressWarnings(
      txdbmaker::makeTxDbFromGFF(
        file = gtffile,
        format = "gtf",
        organism = organism,
        taxonomyId = taxid,
        dataSource = data_source
      )
    )
    AnnotationDbi::saveDb(txdb, dbfile)
  } else {
    message("Loading TxDb...")
    txdb <- AnnotationDbi::loadDb(dbfile)
  }

  message("Preparing reduced genomic features...")

  gr_exons <- unlist(GenomicFeatures::exonsBy(txdb, by = "gene"))
  GenomicRanges::mcols(gr_exons)$feature_id <- names(gr_exons)
  GenomicRanges::mcols(gr_exons)$type <- "exon"

  gr_promoters <- GenomicRanges::promoters(GenomicFeatures::genes(txdb))
  GenomicRanges::mcols(gr_promoters)$feature_id <- names(gr_promoters)
  GenomicRanges::mcols(gr_promoters)$type <- "promoter"

  # Create a Transcript-to-Gene map to standardize Introns/UTRs
  message("Mapping transcripts to genes...")
  tx2gene <- suppressMessages(
    AnnotationDbi::select(
      txdb,
      keys = AnnotationDbi::keys(txdb, "TXNAME"),
      columns = "GENEID",
      keytype = "TXNAME"
    )
  )
  tx2gene_map <- setNames(tx2gene$GENEID, tx2gene$TXNAME)

  gr_introns <- unlist(GenomicFeatures::intronsByTranscript(
    txdb,
    use.names = TRUE
  ))
  GenomicRanges::mcols(gr_introns)$feature_id <- tx2gene_map[names(gr_introns)]
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
  GenomicRanges::mcols(gr_3utr) <- GenomicRanges::mcols(gr_3utr)[, cols_to_keep]
  GenomicRanges::mcols(gr_5utr) <- GenomicRanges::mcols(gr_5utr)[, cols_to_keep]

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
    s_strand := as.character(GenomicRanges::strand(gr_features)[subject_idx])
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
#' @param resource_dir Path to the directory containing index generation resources.
#' Output is saved to this location.
#' @param keep_ranges Should a GRangesList of each transcript/TE locus be saved in the annotation
#' object? default FALSE
#'
#' @return NULL
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' createAnnotation(resource_dir = "/path/to/rmsk-resources")
#' }
createAnnotation <- function(resource_dir, keep_ranges = FALSE) {
  resources <- list.files(resource_dir, full.names = TRUE)
  info_json <- grep("rmsk-duplicateInfo.json", resources, value = TRUE)
  gtf_file <- grep("annotation.gtf.gz", resources, value = TRUE)

  dt <- .dupInfoToDT(info_json)

  message("Creating a GRangesList for all TE ranges...")
  gr <- GenomicRanges::makeGRangesFromDataFrame(
    dt,
    starts.in.df.are.0based = TRUE,
    keep.extra.columns = TRUE
  )
  grl <- S4Vectors::splitAsList(gr, gr$Hash)

  message("Computing overlaps of TE-loci with transcript annotations...")
  ov <- .getHashOverlaps(gr, gtf_file, resource_dir)

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
      Subfamily = stringi::stri_flatten(unique(Subfamily), collapse = ","),
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
  tx <- as(tx, "GRangesList")
  rmsk_grl <- c(tx, grl)

  # Combine annotation DataFrames
  rd <- data.table::rbindlist(list(tx_dt, by_hash), fill = TRUE)
  rd <- S4Vectors::DataFrame(rd)
  rownames(rd) <- c(tx_dt$transcript_id, by_hash$Hash)

  if (isTRUE(keep_ranges)) {
    rd$Ranges <- rmsk_grl[rownames(rd)]
  }

  message(
    "Writing out rowData to: ",
    file.path(resource_dir, "rmsk-rowData.rds")
  )
  saveRDS(rd, file.path(resource_dir, "rmsk-rowData.rds"))
  message("Done.")

  return(invisible(NULL))
}
