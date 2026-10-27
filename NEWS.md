# rmskProfiler 0.99.0

* Initial commit with full R implementation and BiocFileCache resource management.
* Unique RepeatMasker sequences are extracted in R (ShortRead, Biostrings,
  digest) instead of Python. The duplicate information is stored as
  `rmsk-duplicateInfo.rds`.
* Downloaded and generated files are managed with BiocFileCache. Use
  `listResources()`, `getResource()`, and `removeResources()` to manage them.