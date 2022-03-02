#' saeHB.spatial : Small Area Estimation under Spatial SAR Model using Hierarchical Bayesian Method
#'
#' Provides several functions and datasets for area level of Small Area Estimation under Spatial SAR Model using Hierarchical Bayesian (HB) Method. Model-based estimators include the HB estimators based on a Spatial Fay-Herriot model with univariate normal distribution for variable of interest.The 'rjags' package is employed to obtain parameter estimates. For the reference, see Rao and Molina (2015) <doi:10.1002/9781118735855>.
#'
#' @section Author(s):
#' Arina Mana Sikana, Azka Ubaidillah
#'
#' \strong{Maintaner}: Arina Mana Sikana \email{221810195@@stis.ac.id}
#'
#' @section Functions:
#' \describe{
#'   \item{\code{\link{spatial.normal}}}{This function gives small area estimator under Spatial SAR Model and is implemented to variable of interest (y) that assumed to be a Normal Distribution. The range of data is \eqn{(-\infty < y < \infty)}}.
#' }
#'
#' @section Reference:
#' \itemize{
#'    \item{Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New Jersey: John Wiley and Sons, Inc. <doi:10.1002/9781118735855>.}
#'    \item{J. Kubacki and A. Jedrzejczak. (2016). Small Area Estimation of Income Under Spatial SAR Model. Statistics in Transition New Series, Vol. 17, No. 3, pp. 365–390. <doi: 10.21307/stattrans-2016-028>.}
#'    \item{H. C. Chung and G. S. Datta. (2020). Bayesian Hierarchical Spatial Models for Small Area Estimation. Research Report Series. Washington, D.C.: U.S. Census Bureau.}
#' }
#'
#' @docType package
#' @name saeHB.spatial
#'
#' @import stringr
#' @import coda
#' @import rjags
#' @import stats
#' @import grDevices
#' @import graphics
#'
NULL
