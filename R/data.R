#' Simulated multi-state data from a 9-state model of the acute myeloid leukemia
#' (AML) disease pathway
#'
#'
#' @format ## `sim_data_aml`
#' A data frame with 500 rows and 12 columns:
#' \describe{
#'   \item{stop1, stop2, stop3, stop4}{Event times for nested series of competing risks}
#'   \item{CR1}{State of Complete Remission 1 (CR1)}
#'   \item{Death_noCR}{State of death without CR}
#'   \item{Relapse1}{State of first relapse}
#'   \item{NRM_CR}{State of non-relapse mortality (NRM) in CR}
#'   \item{CR2}{State of Complete Remission 2 (CR2)}
#'   \item{RM}{State of relapse mortality (RM)}
#'   \item{Relapse2}{State of second relapse}
#'   \item{Death_CR2}{State of death in CR2}
#'   \item{X1}{Binary covariate X1}
#'   \item{X2}{Binary covariate X2}
#'   ...
#' }
#' @source View data-raw/sim_data_aml.R
"who"
