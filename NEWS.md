# rmskProfiler 0.99.0

* Initial commit with full R implementation and BiocFileCache resource management.
* Unique RepeatMasker sequences are extracted in R (Biostrings, digest)
  instead of Python. The duplicate information is stored as
  `rmsk-duplicateInfo.rds`.
* The genome sequence is taken from BSgenome.Hsapiens.UCSC.hg38 or
  BSgenome.Mmusculus.UCSC.mm10 instead of downloading the GENCODE genome, and
  the transcript sequences are extracted from it using the GENCODE GTF instead
  of downloading the GENCODE transcript sequences. Both are identical to the
  GENCODE files. `downloadResources()` now only downloads the GTF.
  `generateIndex()`, `extractUniqueSeqs()`, and `createGentrome()` gain a
  `genome` argument to use a different genome, and `exampleGenome()` returns
  the example genome. `listResources()` reports the genome each build was
  extracted from, and `removeResources(downloads = TRUE)` removes genome and
  transcript files downloaded by earlier versions.
* Index generation writes fewer intermediate files, which saves ~8 minutes and
  ~2.3 GB of cache for human. `rmskToBed()` returns the filtered RepeatMasker
  records instead of writing rmsk.bed, and `extractUniqueSeqs()` takes them
  with its new `rmsk` argument. rmsk-unique.fa.gz is no longer saved:
  `createGentrome()` extracts the unique TE sequences from the genome using
  rmsk-duplicateInfo.rds. When the index is built, the gentrome is written
  uncompressed to a temporary file and removed afterwards, and it is only saved
  to the cache with `create_index = FALSE`. `createAnnotation()` imports only
  the transcript records of the GTF. The gentrome and Salmon index are
  unchanged.
* Downloaded and generated files are managed with BiocFileCache. Use
  `listResources()`, `getResource()`, and `removeResources()` to manage them.