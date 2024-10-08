% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import-quants.R
\name{importQuants}
\alias{importQuants}
\title{Import rmsk quants as a SummarizedExperiment object}
\usage{
importQuants(quant_dir, resource_dir, remove_zeros = TRUE)
}
\arguments{
\item{quant_dir}{Path to the directories created by salmonQuant}

\item{resource_dir}{Path to the directory containing index generation resources}

\item{remove_zeros}{Should all zero rows be removed before returning object. Default TRUE}
}
\value{
SummarizedExperiment
}
\description{
This function imports counts from the Salmon quant files adjusting for
overdispersion with \code{edgeR::catchSalmon()}. The resulting
\code{SummarizedExperiment} object contains rowData with transcript and TE
loci annotations and an additional column containing a \code{GRangesList} of
the locations for each of the TE loci and transcipts.
}
\details{
The SummarizedEperiment object contains two assays, 'counts' and 'orig.' The
'counts' assay contains the counts \emph{after} adjusting for over dispersion.
This is the assay to be used for downstream differential expression analyses.
The orig' assay contains the counts before scaling. The metadata slot contains
the annotation data.frame returned from catchSalmon. The rowData slot of the
SummarizedExperiment contains columns with boolean values for each TE loci
indicating whether or not that loci has an overlap with a given feature. For
example, hasExonic==TRUE would indicate that loci has an overlap with some
exon on the same strand. hasUnstrandedExonic==TRUE would indicate that the
loci has an overlap with an exon on either strand. The rowData also contains
an additional column 'Ranges' that contains a GRangesList for every
transcript and TE hash location.
}
