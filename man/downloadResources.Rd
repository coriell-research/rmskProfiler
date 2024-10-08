% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download-resources.R
\name{downloadResources}
\alias{downloadResources}
\title{Download files needed for index generation}
\usage{
downloadResources(out_dir, species = c("Hs", "Mm"), check_integrity = TRUE)
}
\arguments{
\item{out_dir}{Directory to save files to. If it does not exist it will be created.}

\item{species}{Either "Hs" (Homo sapiens) or "Mm" (Mus musculus) designating which
species to download resources for}

\item{check_integrity}{TRUE/FALSE, if TRUE check the md5sums of the GENCODE files}
}
\description{
This function will attempt to download all of the necessary resources for
generating the rmsk-gentrome index. For humans, it will download the RepeatMasker
hg38.fa.out file from RepeatMasker.org, GENCODE v36 transcript sequences,
primary assembly, and annotation GTF. For mouse, it will download the
RepeatMasker mm10.fa.out file, GENCODE M25 transcript sequences, primary
assembly, and annotation GTF. If any of these file names already exist in the
out_dir they will be skipped.
}
\examples{
\dontrun{
downloadResources(out_dir = "/path/to/rmsk-resources")
}
}
