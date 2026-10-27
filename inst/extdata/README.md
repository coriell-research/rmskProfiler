# rmskProfiler example data

All files were created by `inst/scripts/make-example-data.R` for the region
chr22:37,000,000-37,999,999 of the human reference genome (GRCh38/hg38).
`exampleCache()`, `exampleGenome()`, and `exampleQuants()` provide access to
them.

| File | Description | Source |
|---|---|---|
| `genome-chr22.fa.gz` | chr22 from the start of the chromosome to the end of the region, with every base before the region replaced by N | `BSgenome.Hsapiens.UCSC.hg38` (identical to GENCODE release 48, `GRCh38.primary_assembly.genome.fa.gz`) |
| `gencode-chr22.gtf.gz` | Annotations of the genes contained in the region | GENCODE release 48, `gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz` |
| `rmsk-chr22.rds` | GRanges of the RepeatMasker records contained in the region | AnnotationHub AH111333, UCSC RepeatMasker annotations (Oct 2022) for hg38 |
| `quants/` | Salmon quants with 10 Gibbs samples for 6 simulated samples (`ctrl_1`-`ctrl_3`, `trt_1`-`trt_3`) against the rmskProfiler index built from the files above | Simulated |
| `fastq/` | The first 2000 simulated read pairs of `ctrl_1`, used in the `salmonQuant()` example | Simulated |

The simulated reads come from mature transcripts, from pre-mRNA (gene bodies
plus a run-on downstream of the 3' end, at a rate that varies between samples),
and from 20 intergenic TE loci expressed in the `trt` samples only. 10% of the
transcripts are expressed 3-fold higher in the `trt` samples.

GENCODE data are freely available without restriction
(https://www.gencodegenes.org/pages/data_access.html). UCSC RepeatMasker
annotations are freely available for all uses
(https://genome.ucsc.edu/conditions.html).

The GENCODE transcript sequences of the genes in the region, used to test that
the transcript sequences extracted from the genome match GENCODE, are in
`tests/testthat/fixtures/transcripts-chr22.fa.gz`.
