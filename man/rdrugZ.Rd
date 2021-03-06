% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdrugZ.R
\name{rdrugZ}
\alias{rdrugZ}
\title{Main Function}
\usage{
rdrugZ(readfile, drugz_outfile, control_samples, drug_samples,
  fc_outfile = F, remove_genes = NULL, pseudocount = 5, minObs = 1,
  half_window_size = 500, index_column = 1, verbose = T)
}
\arguments{
\item{readfile}{an absolute path with filename to counts matrix data (string arrary)}

\item{drugz_outfile}{an absolute path with file name to save drugZ output (string arrary)}

\item{control_samples}{vector of column names used for controls (string vector, i.e. c('CTRL_repA', 'CTRL_repB'))}

\item{drug_samples}{vector of column names used for treatment conditions (string vector, i.e. c('Drug_repA', 'Drug_repB'))}

\item{fc_outfile}{logical T or F (default = F) to save foldchange results to disk}

\item{remove_genes}{vector of gene names (default = NULL) used to remove control genes (string vector, i.e. c('LacZ', 'luciferase', etc..))}

\item{pseudocount}{numeric value (default = 5) for the amount of pseudo counts to add to raw read counts}

\item{minObs}{numeric value (default = 1) for the minimum amount of allowable observations per gene}

\item{half_window_size}{numeric value (default = 500) for calculating empiracle bayes estimate of standard deviation iteratively scanning through a window for calculation}

\item{index_column}{numeric value (default = 1) to indicate which column contains sgRNA identifiers}

\item{verbose}{logical T or F (default = T) to print progress information to console}
}
\value{
Writes drugZ results to disk
}
\description{
rdrugZ
}
\details{
DRUGZ:  Identify drug-gene interactions in paired sample genomic perturbation screens
This script was adpated from the original drugZ python version developed in the Hart Lab.
}
\references{
Publication on drugZ can be found at: \url{https://www.biorxiv.org/content/biorxiv/early/2017/12/12/232736.full.pdf}

Original Python source code can be found at: \url{https://github.com/hart-lab/drugz}
}
\author{
drugZ algorithm was developed in the Traver Hart Lab, by Gang Wang \url{http://www.hart-lab.org}

Script conversion to R: Justin Sing \url{https://github.com/singjc}
}
