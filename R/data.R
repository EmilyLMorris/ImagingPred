#' Simulated imaging data for example
#'
#' A dataset containing simulated imaging data, summarized
#' as one time series for each region of interest.
#'
#' @format A list of 11 data frames, each with 264 rows and 120 columns.
#' \describe{
#' Each row represents the time series for one region of interest in the brain,
#' a summary of multiple voxels.
#' }
"img_data"

#' Simulated covariate data for example
#'
#' X matrix of covariates used to generate the connectivity matrix for an example.
#'
#' @format A data frame with 11 rows (subjects) and 10 columns (covariates).
#' \describe{
#'   Each row contains the covariate values for one subject.
#' }
"X"

#' Beta used to simulate data for the example
#'
#' Matrix of beta values used to simulate the connectivity matrix
#' for multiple subjects for illustration of the package.
#'
#' @format A matrix of coefficient values used to simulate data, with 100 rows and
#' 10 columns.
#' \describe{
#'   Each row contains the coefficients for each covariate used to simulate connectivity
#'   data for a given edge in the network. Each row contains the value for the models for each
#'   edge.
#' }
"true_beta"
