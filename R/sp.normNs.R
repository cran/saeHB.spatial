#' Synthetic Data for Small Area Estimation under Spatial SAR Model and Normal Distribution with non-sampled area
#'
#' @description Synthetic data of 64 regions to simulate Small Area Estimation under Spatial SAR Model and Normal Distribution with non-sampled area using Hierarchical Bayesian Method
#'
#' This data contains \code{NA} values that indicates no sampled at one or more regions. It uses the \code{sp.norm} dataset with the direct estimators and the related variances of 5 regions are missing.
#'
#' @usage data(sp.normNs)
#'
#' @format A data frame with 64 observations on the following 4 variables:
#' \describe{
#'   \item{y}{Direct estimators for each region}
#'   \item{x1}{Auxiliary variable of x1}
#'   \item{x2}{Auxiliary variable of x2}
#'   \item{vardir}{Sampling variance of the direct estimators for each region}
#' }
#'
"sp.normNs"
