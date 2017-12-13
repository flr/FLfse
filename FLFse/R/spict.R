### ------------------------------------------------------------------------ ###
### run SPiCT with FLStock object ####
### ------------------------------------------------------------------------ ###

#' Run SPiCT with FLR objects.
#'
#' This function extracts the catch from an \code{FLStock} object and the survey
#' indices/index from a \code{FLIndices} or \code{FLIndex} object and runs a
#' SPiCT (Surplus Production in Continuous Time) stock assessment.
#'
#' The catch time series is obtained from the \code{catch} slot of \code{stk.}
#'
#' The survey index/indices are obtained from the \code{index} slot(s) of
#' \code{idx}. If the slot(s) contains an age age structure, the sum over all
#' ages is used.
#'
#' Additional configurations can be passed as a list to SPiCT with the
#' \code{conf} argument. They are passed directly to SPiCT
#' (\code{\link[spict]{fit.spict}}) without checking. Any configurations accepted by
#' (\code{\link[spict]{fit.spict}}) can be used.
#'
#' @section Warning:
#' This methods requires the \code{spict} package and all its dependencies to be
#'   installed. For details how to obtain \code{spict}, see
#'   \url{https://github.com/mawp/spict/}.
#'
#' @param stk Object of class \linkS4class{FLStock} with catch time series.
#' @param idx Object of class \linkS4class{FLIndices} or \linkS4class{FLIndex}
#'   object with survey index time series.
#' @param conf Optional configurations passed to SPiCT. Should be a list.
#'
#' @return An object of class \code{spictcls} with the model results.
#'
#' @examples
#' # fit SPiCT to Irish Sea plaice
#' fit <- FLR_SPiCT(stk = ple7a_stk, idx = ple7a_idx)
#' fit
#'
#' # pass additional configuration, set time step to 1 per year
#' conf <- list(dteuler = 1)
#' fit <- FLR_SPiCT(stk = ple7a_stk, idx = ple7a_idx, conf = conf)
#' fit
#'
#' @export

### create generic for FLR_SPiCT
setGeneric("FLR_SPiCT", function(stk, idx, conf = NULL) {
  standardGeneric("FLR_SPiCT")
})

### stk = FLStock, idx = FLIndices
#' @export
setMethod(f = "FLR_SPiCT",
          signature = signature(stk = "FLStock", idx = "FLIndices"),
          definition = function(stk, idx, conf = NULL) {

            FLR_SPiCT_run(stk = stk, idx = idx, conf = conf)

          })
### stk = FLStock, idx = FLIndex
#' @export
setMethod(f = "FLR_SPiCT",
          signature = signature(stk = "FLStock", idx = "FLIndex"),
          definition = function(stk, idx, conf = NULL) {

            ### coerce FLIndex into FLIndices
            idx <- FLIndices(idx)

            ### run SPiCT
            FLR_SPiCT_run(stk = stk, idx = idx, conf = conf)

          })


# fit <- FLR_SPiCT(stk = stk, idx = idx)
# fit2 <- FLR_SPiCT(stk = stk, idx = idx[[1]])

### load FLR objects and run SAM
FLR_SPiCT_run <- function(stk, idx, conf = NULL) {

  ### check if required package is available
  if (!requireNamespace("spict", quietly = TRUE)) {
    stop(paste("Package 'spict' needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }

  # ### check argument classes
  # if (!"FLStock" %in% is(stk)) {
  #   stop("FLR_SPICT(): argument stk has to be object of class FLStock!")
  # }
  # if (!"FLIndex" %in% is(idx) | !"FLIndices" %in% is(idx)) {
  #   stop("FLR_SPICT(): argument stk has to be object of class FLStock!")
  # }

  # ### if idx only one index, convert into FLIndices
  # if (is(idx, "FLIndex")) {
  #   idx <- FLIndices(idx)
  # }

  ### create SPiCT input object list
  inp <- list(obsC = NULL, timeC = NULL, obsI = NULL, timeI = NULL)

  ### extract indices
  ### first survey years
  inp$timeI <- lapply(idx, function(x) {

    ### extract years
    timeI <- dims(x)$minyear:dims(x)$maxyear
    ### set timing of survey as mean of start & end of survey
    idx_t <- mean(range(x)[c("startf", "endf")], na.rm = TRUE)
    if (!is.nan(idx_t)) timeI <- timeI + idx_t

    return(timeI)

  })

  ### then survey index values
  inp$obsI <- lapply(idx, function(x) {

    ### extract index
    obsI <- index(x)
    ### if > age, calc sum over all ages
    if (dim(obsI)[1] > 1) obsI <- quantSums(obsI)
    ### coerce into vector
    obsI <- c(obsI)

    return(obsI)

  })

  ### catch
  ### catch values
  inp$obsC <- c(catch(stk))
  ### catch timing
  inp$timeC <- as.numeric(dimnames(catch(stk))$year)

  ### additional arguments ...
  if (!is.null(conf)) {

    ### add to input object
    inp <- c(inp, conf[!names(conf) %in% names(inp)])

  }

  ### check if input complies with SPiCT's demands
  . <- capture.output(inp <- spict::check.inp(inp))

  ### run SPiCT
  fit <- spict::fit.spict(inp)

  return(fit)

}
