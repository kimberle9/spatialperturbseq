% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{dotplot_upregulated}
\alias{dotplot_upregulated}
\title{Dotplot of upregulated genes for each perturbation.}
\usage{
dotplot_upregulated(
  obj,
  barcodes_to_include,
  include_sgrna = FALSE,
  logfc.threshold = 0.25,
  min.pct = 0.2
)
}
\arguments{
\item{obj}{Seurat v5 object with annotate_gRNA ran}

\item{barcodes_to_include}{List of barcodes/gRNA from annotate_gRNA to be included}

\item{include_sgrna}{Boolean. If FALSE, do not show sgrna feature itself}

\item{logfc.threshold}{Minimum log fold change threshold between the 2 populations. Default 0.25.}

\item{min.pct}{Only consider genes with expressed in this minimum percentage of cells. Default 0.2.}
}
\value{
ggplot2 barplot
}
\description{
Dotplot of upregulated genes for each perturbation.
}
\author{
Kimberle Shen
}
