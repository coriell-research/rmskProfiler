# rmskProfiler 0.99.0

* Submitted to Bioconductor.
* Unique RepeatMasker sequences are extracted in R (ShortRead, Biostrings,
  digest) instead of Python. The duplicate information is stored as
  `rmsk-duplicateInfo.rds`.
* Downloaded and generated files are managed with BiocFileCache. Use
  `listResources()`, `getResource()`, and `removeResources()` to manage them.
