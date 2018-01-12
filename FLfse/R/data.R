#' North Sea cod
#'
#' Assessment input data for North Sea cod as used by ICES WGNSSK 2017. The time
#' series covers the years 1963-2017, with catches up to and including 2016.
#'
#' @format An object of class \code{FLStock}
#'
#' @source \url{https://www.stockassessment.org}
"cod4_stk"

#' North Sea cod indices
#'
#' Survey indices for North Sea cod. The time
#' series covers the years 1963-2016. Two surveys are included: Q1 IBTS and
#' Q3 IBTS. The survey values are the values as used by WGNSSk 2017 and not the
#' raw values.
#'
#' @format An object of class \code{FLIndices} with elements:
#' \describe{
#'   \item{IBTS_Q1_gam}{An object of class \code{FLIndex}}
#'   \item{IBTS_Q3_gam}{An object of class \code{FLIndex}}
#' }
#'
#' @source \url{https://www.stockassessment.org}
"cod4_idx"

#' North Sea cod SAM configuration file
#'
#' Configuration file for the SAM stockassessment for North Sea cod as used by
#' WGNSSK 2017.
#'
#' @source \url{https://www.stockassessment.org}
"cod4_conf_sam"

#' Irish Sea plaice
#'
#' Assessment input data for Irish Sea plaice as used by ICES WGCSE 2017.
#' The time series covers the years 1963-2016.
#'
#' @format An object of class \code{FLStock}
"ple7a_stk"

#' Irish Sea plaice indices
#'
#' Survey indices for Irish Sea plaice as used by WGCSE 2017. Three survey
#' indices are included: "UK BT SURVEY (Sept) - all stations Rcode", "DARDS"
#' and "DARDA".
#'
#' @format An object of class \code{FLIndices} with elements:
#' \describe{
#'   \item{UK BT SURVEY (Sept) - all stations Rcode}{An object of class
#'     \code{FLIndex}}
#'   \item{DARDS}{An object of class \code{FLIndex}}
#'   \item{DARDA}{An object of class \code{FLIndex}}
#' }
"ple7a_idx"

#' Irish Sea plaice SAM configuration file
#'
#' Configuration file for the SAM stockassessment for Irish Sea plaice as used
#' by WGCSE 2017.
"ple7a_conf_sam"



