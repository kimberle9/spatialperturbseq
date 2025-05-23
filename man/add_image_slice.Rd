% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{add_image_slice}
\alias{add_image_slice}
\title{title Adds an image component to a Seurat object so that it is compatible with Seurat's spatial analysis functions. Adapted from: https://github.com/STOmics/STCellbin}
\usage{
add_image_slice(obj)
}
\arguments{
\item{obj}{Seurat v5 object generated by Stereopy's H5ad to RDS script}
}
\value{
Seurat object with image present
}
\description{
title Adds an image component to a Seurat object so that it is compatible with Seurat's spatial analysis functions. Adapted from: https://github.com/STOmics/STCellbin
}
\author{
Kimberle Shen
}
