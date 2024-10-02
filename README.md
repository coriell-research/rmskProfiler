
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

**However**, if you do not wish to generate your own indexes and instead
use one of the [pre-built
versions](https://drive.google.com/drive/folders/1pvxQ9evNGOotktH6Kp2p44UAIU5mQw0U?usp=drive_link)
(hg38 and mm10 using function defaults), then you don’t have to worry
about Python dependencies since Python is only used during index
generation.

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

A full pipeline for generation of the index, quantification, and
importing of quants looks like:

``` r
library(rmskProfiler)
reticulate::use_condaenv("r-rmskProfiler")


# Generate the Salmon index for humans using 12 threads 
# and save to directory named 'hg38-resources'
generateIndex("hg38-resources", species = "Hs", threads = 12)

fq1 <- c("/path/to/sample1.R1.fq.gz", "/path/to/sample2.R1.fq.gz", "/path/to/sample3.R1.fq.gz")
fq2 <- c("/path/to/sample1.R2.fq.gz", "/path/to/sample2.R2.fq.gz", "/path/to/sample3.R2.fq.gz")
sample_names <- c("sample1", "sample2", "sample3")

# Perform quantification with Salmon on fastq files
salmonQuant(
  fq1 = fq1, 
  fq2 = fq2, 
  sample_names = sample_names, 
  resource_dir = "hg38-resources", 
  out_dir = "quants", 
  "--gcBias",                       # Additional arguments can be passed as character strings
  "--seqBias",
  "--posBias",
  "--threads 12"
  )

# Import the transcripts and TE loci counts as a SummarizedExperiment object
# rowData contains transcript and TE annotations and GRanges
se <- importQuants("quants", resources_dir = "hg38-resources")

# Proceed to downstream analysis
```
