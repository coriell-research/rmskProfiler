# rmskProfiler 0.99.0

* The package is now licensed GPL (>= 2) instead of MIT, because it adapts code
  from edgeR (GPL >= 2) and tximport (LGPL >= 2). Copyright notices for the
  adapted code were added to the source files, and its authors are listed as
  contributors and copyright holders.
* `aggregateCounts()` gains a `counts_from_abundance` argument (`"no"`,
  `"scaledTPM"`, or `"lengthScaledTPM"`) to generate gene counts from
  abundances, as in tximport. The aggregated object now has the assays
  `counts` (summed original counts, not divided by overdispersion),
  `abundance`, and `length` (abundance weighted effective gene lengths; a
  length that is the same in every sample for TEs, which are simple sums
  with no length correction), so it can be passed to
  `edgeR::DGEListFromTximeta()` or `edgeR::DGEListFromTximport()`. Its `rowData` gains `Overdispersion`, the
  overdispersion of each aggregated feature estimated from summed Gibbs
  samples, and `n_loci`.
* `importQuants()` adds a `length` assay with Salmon's effective lengths and
  estimates the overdispersion of genes and TE subfamilies, families, and
  classes from the summed Gibbs samples, stored in
  `metadata(se)$overdispersion`.
* Initial commit with full R implementation and BiocFileCache resource management.
* Unique RepeatMasker sequences are extracted in R (Biostrings, digest)
  instead of Python. The duplicate information is stored as
  `rmsk-duplicateInfo.rds`.
* The genome sequence is taken from BSgenome.Hsapiens.UCSC.hg38 or
  BSgenome.Mmusculus.UCSC.mm10, which are now imported, instead of downloading
  the GENCODE genome, and the transcript sequences are extracted from it using
  the GENCODE GTF instead of downloading the GENCODE transcript sequences. Both
  are identical to the GENCODE files. `downloadResources()` now only downloads
  the GTF. `listResources()` reports the genome each build was extracted from,
  and `removeResources(downloads = TRUE)` removes genome and transcript files
  downloaded by earlier versions.
* Index generation writes fewer intermediate files, which saves ~8 minutes and
  ~2.3 GB of cache for human. `rmskToBed()` is renamed `filterRmsk()` and
  returns the filtered RepeatMasker records instead of writing rmsk.bed, and `extractUniqueSeqs()` takes
  RepeatMasker records with its new `rmsk` argument and filters them itself. rmsk-unique.fa.gz is no longer saved:
  `createGentrome()` extracts the unique TE sequences from the genome using
  rmsk-duplicateInfo.rds. When the index is built, the gentrome is written
  uncompressed to a temporary file and removed afterwards, along with any
  gentrome cached by an earlier build. It is only saved to the cache with
  `create_index = FALSE`. `createAnnotation()` imports only
  the transcript records of the GTF. The gentrome and Salmon index are
  unchanged.
* Downloaded and generated files are managed with BiocFileCache. Use
  `listResources()`, `getResource()`, and `removeResources()` to manage them.