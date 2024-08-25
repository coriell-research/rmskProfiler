
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rmskProfiler

<!-- badges: start -->
<!-- badges: end -->

This package provides an end-to-end solution for accurately quantifying
transposable elements from RNA-seq data at the loci-level and properly
importing into R for downstream analysis.

## Installation

You can install the development version of rmskProfiler from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("coriell-research/rmskProfiler")
```

### Python dependencies

This package also requires an installation of Python (either conda or
virtualenv). A helper function in this package exists to download the
necessary Python dependencies to a special environment called
“r-rmskProfiler” using `reticulate`. For example, on my machine which
uses mambaforge, I can install the env like so:

    install_rmskProfiler(
      method = "conda",
      conda = "/home/gennaro/mambaforge/condabin/conda",
      channel = "bioconda"
      )

Once installed, you should load subsequent runs of the package with:

``` r
library(rmskProfiler)
reticulate::use_condaenv("r-rmskProfiler")
```

### Salmon dependency

This package assumes that you have a recent version of
[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) installed
and available on your PATH. If not, please follow the latest
[installation](https://salmon.readthedocs.io/en/latest/building.html#binary-releases)
instructions before using this package.

## Usage

The package has three main components:

1.  Generate a `Salmon` index using unique RepeatMasker elements +
    transcripts + genomic decoy
2.  Quantify reads with Salmon using Gibbs sampling
3.  Import the Salmon quants as a `SummarizedExperiment` for downstream
    differential expression analysis

A full analysis pipeline for generation of the index, quantification,
and importing of quants looks like:

``` r
library(rmskProfiler)
reticulate::use_condaenv("r-rmskProfiler")


# Generate the Salmon Index for humans using 12 threads
generateIndex("rmsk-resources", species = "Hs", threads = 12)

# TODO: Perform quantification with Salmon on fastq files
fq1 <- c("sample1.R1.fq.gz", "sample2.R1.fq.gz", "sample3.R1.fq.gz")
fq2 <- c("sample1.R2.fq.gz", "sample2.R2.fq.gz", "sample3.R2.fq.gz")
sample_names <- c("sample1", "sample2", "sample3")
salmonQuant(fq1, fq2, sample_names, out_dir = "quants", gibbs = 30, threads = 12)

# TODO: Import quants as SummarizedExperiment
se <- importQuants("quants")

# TODO: Summarize counts to gene / subfamily level
summarized <- aggregateCounts(se)
```
