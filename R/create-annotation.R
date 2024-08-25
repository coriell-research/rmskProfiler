#' Extract duplicate info from json file as a data.table
#'
#' @param jsonfile JSON file containing information about rmsk sequences
#' @import data.table
#' @return data.table
#'
.dupInfoToDT <- function(jsonfile) {
  message("Reading in the duplicate information from ", basename(jsonfile), "...")
  info <- jsonlite::fromJSON(jsonfile)

  message("Creating data.table from json information...")
  dt <- data.table::data.table(Hash = names(info), Instance = info)

  message("Unnesting json information...")
  dt <- dt[, .(Instance = as.character(unlist(Instance))), by = Hash]

  message("Extracting repetitive element names from location strings...")
  dt[, c("ID", "Location") := data.table::tstrsplit(Instance, "::", fixed = TRUE)]
  dt[, `:=`(
    RepID = stringi::stri_replace(ID, "", regex = "\\..*"),
    RepName = stringi::stri_replace(ID, "", regex = "^[0-9]+\\.")
  )]

  message("Extracting position information...")
  dt[, c("seqnames", "position") := data.table::tstrsplit(Location, ":", fixed = TRUE)]
  dt[, strand := data.table::fifelse(stringi::stri_detect(position, regex = "\\(+\\)"), "+", "-")]
  dt[, position := stringi::stri_replace(position, "", regex = "\\(.\\)")][,
       c("start", "end") := data.table::tstrsplit(position, "-", fixed = TRUE)][,
       `:=`(start = as.integer(start), end = as.integer(end))]
  dt[, `:=`(Location = NULL, position = NULL, Instance = NULL, ID = NULL)]

  return(dt)
}


#' Find overlaps of TE loci regions with genomic features
#'
#' Finds any TE that overlaps with an exon by transcript, and intron by
#' transcript, a 3 prime UTR by transcript, a 5 prime UTR by transcript, or
#' a promoter by gene.
#'
#' @param gtffile Path to GTF file to create annotation from
#' @param x A GRanges object of TE loci
#' @param stranded TRUE/FALSE should overlaps be computed with respect to strand
#'
#' @return List of hash vectors overlapping genomic features
#'
.getHashOverlaps <- function(x, ignore, gtffile) {
  message("Creating TxDb from GTF...")
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(gtffile))

  message("Extracting genomic regions from txdb...")
  exons_by_tx <- unlist(GenomicFeatures::exonsBy(txdb, by = "tx"))
  introns_by_tx <- unlist(GenomicFeatures::intronsByTranscript(txdb))
  promoters_by_gene <- GenomicRanges::promoters(GenomicFeatures::genes(txdb))
  threeUTR_by_tx <- unlist(GenomicFeatures::threeUTRsByTranscript(txdb))
  fiveUTR_by_tx <- unlist(GenomicFeatures::fiveUTRsByTranscript(txdb))

  message("Finding overlaps between TE loci and genomic features...")
  exon_hits <- GenomicRanges::findOverlaps(x, exons_by_tx, ignore.strand = !ignore)
  intron_hits <- GenomicRanges::findOverlaps(x, introns_by_tx, ignore.strand = !ignore)
  promoter_hits <- GenomicRanges::findOverlaps(x, promoters_by_gene, ignore.strand = !ignore)
  threeUTR_hits <- GenomicRanges::findOverlaps(x, threeUTR_by_tx, ignore.strand = !ignore)
  fiveUTR_hits <- GenomicRanges::findOverlaps(x, fiveUTR_by_tx, ignore.strand = !ignore)

  message("Collecting results...")
  hash_in_exon <- unique(x[S4Vectors::queryHits(exon_hits), ]$Hash)
  hash_in_intron <- unique(x[S4Vectors::queryHits(intron_hits), ]$Hash)
  hash_in_promoter <- unique(x[S4Vectors::queryHits(promoter_hits), ]$Hash)
  hash_in_3utr <- unique(x[S4Vectors::queryHits(threeUTR_hits), ]$Hash)
  hash_in_5utr <- unique(x[S4Vectors::queryHits(fiveUTR_hits), ]$Hash)

  result <- list(
     hash_in_exon = hash_in_exon,
     hash_in_intron = hash_in_intron,
     hash_in_promoter = hash_in_promoter,
     hash_in_3utr = hash_in_3utr,
     hash_in_5utr = hash_in_5utr
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
#' @param stranded Should the annotation be created with respect to strand? Default
#' TRUE. If stranded=FALSE then TE overlap annotation is performed without respect
#' to strand, i.e. findOverlaps(..., ignore.strand = TRUE).
#'
#' @return NULL
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' createAnnotation(resource_dir = "/path/to/rmsk-resources")
#' }
createAnnotation <- function(resource_dir, stranded = TRUE) {

  resources <- list.files(resource_dir, full.names = TRUE)
  info_json <- grep("rmsk-duplicateInfo.json", resources, value = TRUE)
  gtf_file <- grep("annotation.gtf.gz", resources, value = TRUE)

  dt <- .dupInfoToDT(info_json)

  message("Creating a GRangesList for all TE ranges...")
  gr <- GenomicRanges::makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE)
  grl <- S4Vectors::splitAsList(gr, gr$Hash)

  ov <- .getHashOverlaps(gr, stranded, gtf_file)

  message("Getting all unique hash-element pairs...")
  hash_dt <- dt[, .(N_Loci = .N), by = .(Hash, RepName)]
  hash_dt[, c("Class", "Family", "Subfamily") := data.table::tstrsplit(RepName, ".", fixed = TRUE)]

  message("Collapsing hash-level information into rowData...")
  by_hash <- hash_dt[, .(
    Class = stringi::stri_flatten(unique(Class), collapse = ","),
    Family = stringi::stri_flatten(unique(Family), collapse = ","),
    Subfamily = stringi::stri_flatten(unique(Subfamily), collapse = ","),
    N_Loci = sum(N_Loci),
    N_Class = length(unique(Class)),
    N_Family = length(unique(Family)),
    N_Subfamily = length(unique(Subfamily))
  ), by = Hash]

  by_hash[, `:=`(
    hasPromoter = Hash %chin% ov$hash_in_promoter,
    hasExonic = Hash %chin% ov$hash_in_exon,
    hasIntronic = Hash %chin% ov$hash_in_intron,
    has3UTR = Hash %chin% ov$hash_in_3utr,
    has5UTR = Hash %chin% ov$hash_in_5utr)][,
    hasIntergenic := (!hasExonic & !hasIntronic & !has3UTR & !has5UTR)]

  message("Writing out hash-level rowData to: ", file.path(resource_dir, "rmsk-rowData.tsv"))
  data.table::fwrite(by_hash, file.path(resource_dir, "rmsk-rowData.tsv"), sep = "\t")

  message("Writing TE GRangesList to: ", file.path(resource_dir, "rmsk-grl.rds"))
  saveRDS(grl, file.path(resource_dir, "rmsk-grl.rds"))
  message("Done.")

  return(invisible(NULL))
}
