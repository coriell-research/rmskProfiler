#' Generate annotation resources for annotating TE loci and genes
#'
#' This function produces multiple annotation files used to annotate downstream
#' counts generated after quantification with Salmon. It creates: GRangesList
#' object for each Hash value and all unique ranges for the hash and a
#' tab-delimited file of unique hash-by-overlap information. This function
#' generates overlap annotations with respect to the annotations provided by
#' annotatr::build_annotations() on either hg38_basicgenes or mm10_basicgenes,
#' depending on the resources detected.
#'
#' @param resource_dir Path to the directory containing index generation resources.
#' Output is saved to this location.
#' @param stranded Should the annotation be created with respect to strand? Default
#' TRUE. If stranded=FALSE then TE overlap annotation is performed on sequence
#' positions without respect to strand.
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
  genome_fa <- grep("primary_assembly.genome.fa.gz", resources, value = TRUE)
  build <- "hg38"
  if (grepl("GRCm38", genome_fa)) {
    build <- "mm10"
  }

  # Process duplicate info --------------------------------------------------

  message("Reading in the duplicate information from ", basename(info_json), "...")
  info <- jsonlite::fromJSON(info_json)

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

  if (isTRUE(stranded)) {
    data.table::setkey(dt, seqnames, strand, start, end)
  } else {
    data.table::setkey(dt, seqnames, start, end)
  }

  message("Creating a GRangesList for all TE ranges...")
  gr <- GenomicRanges::makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE)
  grl <- S4Vectors::splitAsList(gr, gr$Hash)

  # Create annotations ------------------------------------------------------

  message("Building genic annotation with `annotatr`...")
  if (build == "hg38") {
    anno <- suppressMessages(annotatr::build_annotations("hg38", "hg38_basicgenes"))
  } else if (build == "mm10") {
    anno <- suppressMessages(annotatr::build_annotations("mm10", "mm10_basicgenes"))
  } else {
    stop("build must be one of 'hg38' or 'mm10'")
  }

  # Coerce to data.table, key, and remove non-standard chromosomes
  anno <- data.table::as.data.table(data.frame(anno))

  if (isTRUE(stranded)) {
    data.table::setkey(anno, seqnames, strand, start, end)
  } else {
    data.table::setkey(anno, seqnames, start, end)
  }

  keep_chrom <- paste("chr", c(1:22, "X", "Y"), sep = "")
  anno <- anno[seqnames %in% keep_chrom]

  # Perform overlap ---------------------------------------------------------

  message("Overlapping TE locations with genic annotations...")
  ov <- data.table::foverlaps(dt, anno, type = "any")

  # If no overlaps with genic regions then type is intergenic
  ov[is.na(type), `:=`(
    type = stringi::stri_c(build, "_genes_intergenic"),
    start = 0L, end = 0L, width = 0L
  )]

  # Rowdata -----------------------------------------------------------------

  message("Getting all unique hash-element pairs...")
  hash_dt <- dt[, .(N_Loci = .N), by = .(Hash, RepName)]
  hash_dt[, c("Class", "Family", "Subfamily") := data.table::tstrsplit(RepName, ".", fixed = TRUE)]

  message("Collapsing hash-level information into rowData (this may take some time)...")
  by_hash <- hash_dt[, .(
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

  message("Annotating hashes with genomic location...")
  is_intergenic <- ov[stringi::stri_detect(type, fixed = "intergenic"), unique(Hash)]
  is_1to5kb <- ov[stringi::stri_detect(type, fixed = "1to5kb"), unique(Hash)]
  is_promoter <- ov[stringi::stri_detect(type, fixed = "promoters"), unique(Hash)]
  is_exonic <- ov[stringi::stri_detect(type, fixed = "exons"), unique(Hash)]
  is_intronic <- ov[stringi::stri_detect(type, fixed = "introns"), unique(Hash)]
  is_3utr <- ov[stringi::stri_detect(type, fixed = "3UTRs"), unique(Hash)]
  is_5utr <- ov[stringi::stri_detect(type, fixed = "5UTRs"), unique(Hash)]

  by_hash[, `:=`(
    hasIntergenic = Hash %chin% is_intergenic,
    has1to5kb = Hash %chin% is_1to5kb,
    hasPromoter = Hash %chin% is_promoter,
    hasExonic = Hash %chin% is_exonic,
    hasIntronic = Hash %chin% is_intronic,
    has3UTR = Hash %chin% is_3utr,
    has5UTR = Hash %chin% is_5utr
  )]

  message("Writing TE GRangesList to: ", file.path(resource_dir, "rmsk-grl.rds"))
  saveRDS(grl, file.path(resource_dir, "rmsk-grl.rds"))

  message("Writing out hash-level rowData to: ", file.path(resource_dir, "rmsk-rowData.tsv"))
  data.table::fwrite(by_hash, file.path(resource_dir, "rmsk-rowData.tsv"), sep = "\t")
  message("Done.")

  return(invisible(NULL))
}
