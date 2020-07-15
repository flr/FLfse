### suppress warning about variables used in foreach loop
### when checking package
utils::globalVariables(c("fit_i", "iter_i", "i"))
### ------------------------------------------------------------------------ ###
### run SAM with FLStock object ####
### ------------------------------------------------------------------------ ###

### convert matrix into FLQuant
matrix2FLQuant <- function(input) {
  ### get year & age range
  years <- dimnames(input)[[1]]
  ages <- dimnames(input)[[2]]
  ### create dimnames
  if (length(ages) == 1 & isTRUE(ages == -1)) {
    dims_qnt <- list(quant = "all", year = years)
  } else {
    dims_qnt <- list(age = ages, year = years)
  }
  ### create FLQuant
  FLQuant(t(input), dimnames = dims_qnt)
}

### function for running SAM
### runs inside FLR_SAM_run
FLR_SAM_iter <- function(stk, idx, conf, par_ini, i, NA_rm, conf_full, ...) {

  ### subset stock and index to current iter
  stk_i <- FLCore::iter(stk, i)
  idx_i <- lapply(idx, FLCore::iter, i)
  par_i <- par_ini[[i]]

  ### calculate landings fraction
  lf <- landings.n(stk_i) / catch.n(stk_i)
  ### assume only landings if NA
  lf[is.na(lf) | is.infinite(lf)] <- 1

  ### do some sanity checks ...
  ### fill empty landings/discards weights

  ### create SAM input objects
  dat_lst <- list(residual.fleet = catch.n(stk_i),
                  prop.mature = mat(stk_i),
                  stock.mean.weight = stock.wt(stk_i),
                  catch.mean.weight = catch.wt(stk_i),
                  dis.mean.weight = discards.wt(stk_i),
                  land.mean.weight = landings.wt(stk_i),
                  natural.mortality = m(stk_i),
                  prop.f = harvest.spwn(stk_i),
                  prop.m = m.spwn(stk_i),
                  land.frac = lf)

  ### reduce dimensions, keep only age & year
  dat_lst <- lapply(dat_lst, function(x) {
    x[, drop = TRUE]
  })

  ### transpose matrix (required for SAM)
  dat_lst <- lapply(dat_lst, t)

  ### remove trailing years which contain only NAs
  ### otherwise SAM will complain
  if (isTRUE(NA_rm)) {

    dat_lst <- lapply(dat_lst, function(x) {
      ### find rows/years with NAs
      years_NA <- apply(x, 1, function(y) {
        all(is.na(y))
      })
      ### check if last year contains only NAs, if yes, find length of sequence
      if (isTRUE(c(tail(years_NA, 1), use.names = FALSE))) {
        ### length of TRUE/FALSE sequences
        seq_lngth <- rle(years_NA)
        ### length of last TRUE sequence
        lngth_NA <- as.vector(tail(seq_lngth$lengths, 1))
        ### find positions for corresponding years
        pos_remove <- which(rownames(x) == tail(rownames(x), lngth_NA))
        x <- x[-pos_remove, ]
      }
      return(x)
    })
    ### trim years of land.frac to dimension of catch
    ### (after removing trailing NAs)
    ### otherwise, weird things happen when a forecast is performed...
    if (any(dim(dat_lst$residual.fleet) != dim(dat_lst$land.frac))) {
      dat_lst$land.frac <- dat_lst$land.frac[seq(nrow(dat_lst$residual.fleet)),]
    }

  }

  ### add recapture data, if they exists
  if (!is.null(attr(catch.n(stk_i), "recap"))) {
    recap <- attr(catch.n(stk_i), "recap")
    if (is.list(recap) & !is.data.frame(recap)) {
      if (isTRUE(length(recap) > 1)) {
        recap <- recap[[i]]
      }
    }
    dat_lst$recapture <- recap
  }

  ### catch number weight
  if (!is.null(attr(catch.n(stk_i), "weight"))) {
    cn_weight <- attr(catch.n(stk_i), "weight")
    if (isTRUE(dims(cn_weight)$iter > 1)) {
      cn_weight <- FLCore::iter(cn_weight, i)
    }
    cn_weight <- t(cn_weight[, drop = TRUE])
    attr(dat_lst$residual.fleet, "weight") <- cn_weight
  }

  ### extract indices
  dat_idx <- lapply(idx_i, function(x){
    ### extract index slot
    tmp <- index(x)
    ### get dims
    tmp_dim <- dim(tmp)
    tmp_dimnames <- dimnames(tmp)
    ### coerce into array, drop all dimensions apart from age/year
    tmp <- array(data = tmp, dim = tmp_dim[1:2], dimnames = tmp_dimnames[1:2])
    ### transpose
    tmp <- t(tmp)
    ### change age description for ssb/biomass surveys
    if (dim(tmp)[2] == 1) {
      if (!grepl(x = dimnames(tmp)[[2]],
                 pattern = "^[0-9]+$") | ### non numeric age
          is.na(dimnames(tmp)[[2]]) | ### NA as age
          dimnames(tmp)[[2]] == -1) {
        dimnames(tmp)[[2]] <- -1 ### recognized by SAM as SSB index
      }
    }

    ### add survey timing as attribute
    attr(tmp, "time") <- as.vector(range(x)[c("startf", "endf")])

    return(tmp)
  })
  ### add partial attribute "part" to surveys, e.g. used for herring
  part_n <- which(sapply(idx_i, type) == "partial")
  if (length(part_n) > 1) {
    for (idx_i_i in part_n) {
      attr(dat_idx[[idx_i_i]], "part") <- as.vector(which(part_n == idx_i_i))
    }
  }

  ### add indices to input list
  dat_lst$surveys <- dat_idx

  ### create SAM input object
  dat_sam <- do.call(stockassessment::setup.sam.data, dat_lst)

  ### create default configuration
  conf_sam <- stockassessment::defcon(dat_sam)

  ### fbar range
  if (all(!is.na(range(stk_i)[c("minfbar", "maxfbar")]))) {
    conf_sam$fbarRange <- range(stk_i)[c("minfbar", "maxfbar")]
  }

  ### insert configuration, if supplied to function
  if (!is.null(conf) & !isTRUE(conf_full)) {

    ### find slots that can be used
    conf_names <- intersect(names(conf), names(conf_sam))

    ### go through slot names
    if (length(conf_names) > 0) {

      for (y in conf_names) {

        ### workaround for keyParScaledYA
        ### default 0x0 matrix, needs to extended
        if (y == "keyParScaledYA") {

          conf_sam[[y]] <- conf[[y]]

          ### otherwise enter only values
        } else {

          ### position where values supplied (i.e. not NA)
          pos <- NULL ### reset from previous iteration
          pos <- which(!is.na(conf[[y]]))
          ### insert values
          conf_sam[[y]][pos] <- conf[[y]][pos]

        }

      }

    }

    ### otherwise use provided configuration without ANY checking
  } else if (isTRUE(conf_full)) {

    conf_sam <- conf

  }

  ### define parameters for SAM if none are passed through
  if (is.null(par_i)) {

    par_i <- stockassessment::defpar(dat_sam, conf_sam)

    ### check dimensions of parameters and adapt if neccessary
  } else {

    ### use stock weight to get year dimensions of data
    n_yrs <- dim(dat_sam$stockMeanWeight)[1]
    ### same for initial parameter values (data transposed ...)
    n_yrs_ini <- dim(par_i$logF)[2]

    ### if dimensions differ,
    ### remove redudant years at end of time series if more values provided
    ### or replicate values from last year if values missing

    ### adapt stock numbers and fishing mortality initial values
    if (n_yrs_ini < n_yrs) {

      ### find missing years
      n_missing <- n_yrs - n_yrs_ini
      ### extend
      par_i$logF <- par_i$logF[, c(seq(n_yrs_ini),
                                   rep(tail(n_yrs_ini, n_missing)))]
      par_i$logN <- par_i$logN[, c(seq(n_yrs_ini),
                                   rep(tail(n_yrs_ini, n_missing)))]

    } else if (n_yrs_ini > n_yrs) {

      par_i$logF <- par_i$logF[, seq(n_yrs)]
      par_i$logN <- par_i$logN[, seq(n_yrs)]

    }

  }

  ### arguments for SAM
  SAM_args <- list(data = dat_sam, conf = conf_sam,
                   parameters = par_i)
  ### additional configuration for recapture data
  if (!is.null(attr(catch.n(stk_i), "recap_conf"))) {
    SAM_args <- append(SAM_args, attr(catch.n(stk_i), "recap_conf"))
  }
  ### more arguments from "..."
  if (length(names(list(...))) > 0) {
    SAM_args <- append(SAM_args, list(...))
  }
  ### run SAM
  sam_msg <- capture.output(
    fit <- do.call(stockassessment::sam.fit, SAM_args)
  )

  ### save screen message(s) as attribute
  attr(x = fit, which = "messages") <- sam_msg

  ### return
  return(fit)

}

### load FLR objects and run SAM
### function not exported, use FLR_SAM instead
FLR_SAM_run <- function(stk, idx, conf = NULL,
                        conf_full = FALSE, ### use provided conf in full
                        force_list_output = FALSE,
                        DoParallel = FALSE, ### compute iterations in parallel
                        par_ini = NULL, ### initial parameter values
                        NA_rm = TRUE, ### remove trailing years with NAs
                        ... ### passed to sam.fit()
                        ) {

  ### check if required package is available
  if (!requireNamespace("stockassessment", quietly = TRUE)) {
    stop(paste("Package 'stockassessment' needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }

  ### check iter dimension and increase if necessary
  lengths <- c(dims(stk)$iter, sapply(lapply(idx, dim), "[", 6))
  if (length(unique(lengths)) > 1) {

    ### check if iterations can be propagated, i.e. if 1 or >1
    if (any(lengths[-which(lengths == max(lengths))] > 1)) {
      stop("incompatible number of iterations in stock and index!")
    }

    ### propagate stock
    if (lengths[1] < max(lengths)) {
      stk <- propagate(stk, max(lengths), fill.iter = TRUE)
    }

    ### propagate index/indices
    if (any(lengths[-1] < max(lengths))) {
      idx <- lapply(idx, propagate, iter = max(lengths), fill.iter = TRUE)
    }

  }

  ### number of iterations
  it <- dim(stk)[6]

  ### check dimensions of supplied parameters and propagate if neccessary
  if (length(par_ini) != it) {

    par_ini <- rep(list(par_ini), it)

  }

  ### set parallel or serial processing of iterations
  `%do_tmp%` <- ifelse(isTRUE(DoParallel), foreach::`%dopar%`, foreach::`%do%`)

  ### go through iterations
  res_iter <- foreach(i = seq(dims(stk)$iter), .errorhandling = "pass",
                      .packages = c("FLCore", "stockassessment")) %do_tmp% {

    FLR_SAM_iter(stk = stk, idx = idx, conf = conf, par_ini = par_ini, i = i,
                 NA_rm = NA_rm, conf_full = conf_full, ...)

  }

  ### set class to "sam_list"
  class(res_iter) <- "sam_list"

  ### if only 1 iteration, return remove list structure
  if (length(res_iter) == 1 & !isTRUE(force_list_output)) {

    return(res_iter[[1]])

  ### otherwise, return full list
  } else {

    return(res_iter)

  }

}


### ------------------------------------------------------------------------ ###
### Run SAM with FLR objects ####
### ------------------------------------------------------------------------ ###

#' Run SAM with FLR objects.
#'
#' This function runs a SAM assessment using FLR objects as input. Stock and
#' fishery data is extracted from an object of class \code{FLStock}, survey
#' indices from \code{FLIndex} or \code{FLIndices}. Additional model
#' configurations can be passed to SAM.
#'
#' Stock and fishery data is extracted from \code{stk}, survey indices from
#' \code{idx}.
#' Stock data used are:
#' \itemize{
#'   \item catch numbers at age (\code{catch.n} slot from \code{stk})
#'   \item maturity ogive (\code{mat} slot from \code{stk})
#'   \item stock weight at age (\code{stock.wt} slot from \code{stk})
#'   \item catch weight at age (\code{catch.wt} slot from \code{stk})
#'   \item discards weigth at age (\code{discards.wt} slot from \code{stk})
#'   \item landings weight at age (\code{landings.wt} slot from \code{stk})
#'   \item natural mortality at age (\code{m} slot from \code{stk})
#'   \item proportion of fishing mortality before spawning at age
#'     (\code{harvest.spwn} slot from \code{stk})
#'    \item proportion of natural mortality before spawning at age
#'      (\code{catch.n} slot from \code{stk})
#'    \item landing fraction at age (calculated as landings numbers at age
#'      divided by the catch numbers at age)
#' }
#'
#' Trailing years without data in \code{stk} are removed, unless turned
#' of by setting \code{NA_rm = FALSE}.
#'
#' Survey indices are extracted from \code{idx}, using the \code{index}
#' slot(s).
#' SSB indices can be used. To define an index as SSB index, the first (age)
#' dimension of the \code{index} slot of \code{idx} has to be of length 1 and
#' the name of this dimension can be either missing (\code{NA}), non-numeric
#' (e.g. "ssb") or \code{-1}.
#'
#' Additional configurations can be passed as a list to SAM with the
#' \code{conf} argument. If argument \code{conf_full} is set to \code{TRUE},
#' the configuration is passed straight on to SAM without any checking.
#' If argument \code{conf_full} is set to \code{FALSE} (default), then
#' \code{FLR_SAM} first generates a default model configuration with
#' \code{stockassessment}'s \code{setup.sam.data} and additional
#' configurations available in \code{conf} replace default configurations.
#' For details about possible configurations and
#' format see '\code{help("defcon", package = "stockassessment")} and
#' \url{https://github.com/fishfollower/SAM}.
#'
#' The function can handle input objects with multiple iterations (\code{iter}
#' dimension in \code{stk} and \code{idx}). If multiple iterations are provided
#' for \code{stk} but not for \code{idx}, \code{idx} will be inflated, and vice
#'  versa.
#' If the assessment fails for some iterations, the error messages are returned
#' for these iterations.
#' If argument \code{DoParallel} is set to \code{TRUE}, the individual
#' iterations are processed in parallel. This uses parallel computing provided
#' by the package \code{DoParallel}. The parallel workers need to be set up
#' and registered before calling the function. See \code{?DoParallel}.
#'
#' Argument \code{par_ini} allows the provision of initial parameter values for
#' SAM and can speed up the model. Either a single set of parameters can be
#' supplied and they are recycled if neccessary. Alternatively, a list of
#' initial parameters can be supplied, one for each iteration of the
#' stock/index. If the dimensions of initial values for numbers/fishing
#' mortality at age differ from the data, redundant years are automatically
#' removed and if years are missing, the values from the last provided year are
#' recycled.
#'
#' The default console output generated by SAM is not printed but saved. It is
#' stored as an attribute of the fit and can be accessed with
#' \code{attr(fit, "messages")}.
#'
#' Tagging data can be provided in the usual SAM format and should be stored
#' as an attribute of \code{attr(catch.n(stk), "recap")}. Additional tagging
#' configurations can be supplied as a list with the attribute
#' \code{attr(catch.n(stk), "recap_conf")}, e.g.
#' \code{list(map = list(logitRecapturePhi = factor(c(1, 1))))}, which are then
#' passed on to \code{sam.fit()}.
#'
#' Weights for the catch numbers can be supplied as an attribute of
#' \code{attr(catch.n(stk), "weight")} and should be formatted as \code{FLQuant}
#' objects.
#'
#' @section Warning:
#' This methods requires the \code{stockassessment} package and all its
#' dependencies to be installed. For details how to obtain
#' \code{stockassessment}, see \url{https://github.com/fishfollower/SAM/}.
#'
#' @param stk Object of class \linkS4class{FLStock} with stock and fishery data.
#' @param idx Object of class \linkS4class{FLIndices} or \linkS4class{FLIndex}
#'   object with survey index time series.
#' @param conf Optional configurations passed to SAM. Defaults to \code{NULL}.
#'   If provided, should be a list.
#' @param conf_full Use provided configuration object in full without ANY
#'   checking (see Details for more information).
#' @param  par_ini Optional starting parameters for SAM. See details for more
#'   information.
#' @param DoParallel Optional, defaults to \code{FALSE}. If set to \code{TRUE},
#'   will perform iterations of stock in parallel. See Details below for
#'   description.
#' @param NA_rm Remove trailing years with NAs, defaults to \code{TRUE}.
#' @param ... Additional arguments passed to \code{sam.fit()}, e.g.
#'   \code{newtonsteps}
#'
#' @return An object of class \code{sam} (for single iteration) or
#'   \code{sam_list} (list of \code{sam} objects for multiple iterations) with
#'   the model results.
#'
#' @examples
#' # fit SAM to North Sea cod
#' fit <- FLR_SAM(stk = cod4_stk, idx = cod4_idx)
#'
#' # use WGNSSK 2017 configuration
#' fit <- FLR_SAM(stk = cod4_stk, idx = cod4_idx, conf = cod4_conf_sam)
#'
#' # fit SAM to Irish Sea plaice
#' fit <- FLR_SAM(stk = ple7a_stk, idx = ple7a_idx, conf = ple7a_conf_sam)
#'
#' @export

setGeneric("FLR_SAM", function(stk, idx, conf = NULL, conf_full = FALSE,
                               par_ini = NULL, DoParallel = FALSE, NA_rm = TRUE,
                               ...) {
  standardGeneric("FLR_SAM")
})

### stk = FLStock, idx = FLIndices
#' @rdname FLR_SAM
setMethod(f = "FLR_SAM",
          signature = signature(stk = "FLStock", idx = "FLIndices"),
          definition = function(stk, idx, conf = NULL, conf_full = FALSE,
                                par_ini = NULL, DoParallel = FALSE,
                                NA_rm = TRUE, ...) {

  FLR_SAM_run(stk = stk, idx = idx, conf = conf, conf_full = conf_full,
              DoParallel = DoParallel, par_ini = par_ini,  NA_rm = NA_rm,
              ...)

})
### stk = FLStock, idx = FLIndex
#' @rdname FLR_SAM
setMethod(f = "FLR_SAM",
          signature = signature(stk = "FLStock", idx = "FLIndex"),
          definition = function(stk, idx, conf = NULL, conf_full = FALSE,
                                par_ini = NULL, DoParallel = FALSE,
                                NA_rm = TRUE, ...) {

  ### coerce FLIndex into FLIndices
  idx <- FLIndices(idx)

  ### run SPiCT
  FLR_SAM_run(stk = stk, idx = idx, conf = conf, conf_full = conf_full,
              DoParallel = DoParallel, par_ini = NULL,  NA_rm = NA_rm,
              ...)

})


### ------------------------------------------------------------------------ ###
### convert sam object into FLStock ####
### ------------------------------------------------------------------------ ###

### function for converting sam into FLStock
sam_to_FLStock <- function(object, ### sam object
                        uncertainty = FALSE, ### confidence intervals?
                        conf_level = 95, ### confidence interval level in %
                        catch_estimate = FALSE, ### use catch input or
                                                ### model estimates
                        correct_catch = FALSE, ### "correct" catch with
                                               ### estimated multiplier
                        ...) {

  ### check class type of input
  if (!isTRUE("sam" %in% class(object))) stop("object has to be class sam")

  ### calculate SD multiplier for uncertainty range
  if (isTRUE(uncertainty)) SD_mult <- stats::qnorm(0.5 + conf_level/200)

  ### get dimensions
  years <- object$data$years
  ages <- object$conf$minAge:object$conf$maxAge

  ### create FLQuant dummy
  qnt <- FLQuant(NA, dimnames = list(age = ages, year = years))
  ### create FLStock
  stk <- FLStock(qnt)

  ### stock numbers @ age
  stock.n <- exp(object$pl$logN)
  n_ages <- dim(stock.n)[1]
  n_yrs <- dim(stock.n)[2]
  stock.n(stk)[1:n_ages, 1:n_yrs] <- stock.n

  ### add range
  if (isTRUE(uncertainty)) {
    stock.n_sd <- qnt
    stock.n_sd[1:n_ages, 1:n_yrs] <- object$plsd$logN
    attr(stock.n(stk), "low") <- stock.n(stk) - stock.n_sd * SD_mult
    attr(stock.n(stk), "high") <- stock.n(stk) + stock.n_sd * SD_mult
  }

  ### harvest @ age
  harvest <- exp(object$pl$logF)
  ### duplicate linked ages
  harvest <- harvest[object$conf$keyLogFsta[1,] + 1, ]
  n_ages <- dim(harvest)[1]
  n_yrs <- dim(harvest)[2]
  harvest(stk)[1:n_ages, 1:n_yrs] <- harvest
  units(harvest(stk)) <- "f"

  ### add range
  if (isTRUE(uncertainty)) {
    harvest_sd <- qnt
    harvest_sd[1:n_ages,1:n_yrs]<-object$plsd$logF[object$conf$keyLogFsta[1,]+1,]
    attr(harvest(stk), "low") <- harvest(stk) - harvest_sd * SD_mult
    attr(harvest(stk), "high") <- harvest(stk) + harvest_sd * SD_mult
  }

  ### ---------------------------------------------------------------------- ###
  ### input data ####

  ### stock ####
  ### stock weight @ age
  stock.wt <- t(object$data$stockMeanWeight)
  n_ages <- dim(stock.wt)[1]
  n_yrs <- dim(stock.wt)[2]
  stock.wt(stk)[1:n_ages, 1:n_yrs] <- stock.wt
  ### stock biomass
  stock <- quantSums(qnt)
  stock[, 1:n_yrs] <- exp(object$sdrep$value[names(object$sdrep$value) == "logtsb"])
  stock(stk) <- stock

  ### add range
  if (isTRUE(uncertainty)) {
    stock_sd <- quantSums(qnt)
    stock_sd[, 1:n_yrs] <- object$sdrep$sd[names(object$sdrep$value) == "logtsb"]
    attr(stock(stk), "low") <- exp(log(stock(stk)) - stock_sd * SD_mult)
    attr(stock(stk), "high") <- exp(log(stock(stk)) + stock_sd * SD_mult)
  }


  ### catch ####
  ### catch weight @ age
  if (isTRUE(length(dim(object$data$catchMeanWeight)) == 3)) {
    catch.wt <- t(object$data$catchMeanWeight[,,, drop = TRUE])
  } else {
    catch.wt <- t(object$data$catchMeanWeight)
  }
  n_ages <- dim(catch.wt)[1]
  n_yrs <- dim(catch.wt)[2]
  catch.wt(stk)[1:n_ages, 1:n_yrs] <- catch.wt
  ### catch numbers @ age
  catch_fleets <- which(object$data$fleetTypes == 0)
  ### extract observations & estimates
  dat_catch <- cbind(object$data$aux, value = exp(object$data$logobs),
                     estimate = exp(object$rep$predObs))
  dat_catch <- dat_catch[dat_catch[, "fleet"] == catch_fleets, ]
  ### sum up catch numbers over all commercial fleets
  dat_catch <- stats::aggregate(cbind(value, estimate) ~ age + year, dat_catch,
                         FUN = sum)
  ### add missing ages/years
  dat_full <- expand.grid(age = unique(dat_catch$age),
                          year = unique(dat_catch$year))
  dat_catch <- merge(x = dat_catch, y = dat_full, all = TRUE)
  dat_catch <- dat_catch[order(dat_catch$year, dat_catch$age), ] ### sort

  ### insert catch: estimates or input values
  if (!isTRUE(catch_estimate)) {

    ### use input data
    catch.n(stk)[ac(unique(dat_catch$age)), ac(unique(dat_catch$year))] <-
      dat_catch$value

  } else {

    ### use value estimated by SAM
    catch.n(stk)[ac(unique(dat_catch$age)), ac(unique(dat_catch$year))] <-
      dat_catch$estimate

  }

  ### correct catch numbers if a catch multiplier was estimated
  if (isTRUE(correct_catch)) {

    ### warn if no catch multiplier available
    if (!(length(object$pl$logScale) > 1)) {
      warning("catch correction requested but no catch multiplier estimated")
    }

    ### get catch multiplier dimensions
    ages_mult <- object$conf$minAge:object$conf$maxAge
    yrs_mult <- object$conf$keyScaledYears
    ### create catch multiplier FLQuant
    catch_mult_data <- FLQuant(
      matrix(data = object$pl$logScale[(object$conf$keyParScaledYA + 1)],
             ncol = object$conf$noScaledYears,
             nrow = length(object$conf$minAge:object$conf$maxAge),
             byrow = TRUE),
      dimnames = list(year = object$conf$keyScaledYears,
                      age = object$conf$minAge:object$conf$maxAge))
    ### model values in log scale, exponentiate
    catch_mult_data <- exp(catch_mult_data)
    ### for simpler calculations, expand dimensions for stock
    catch_mult <- catch.n(stk) %=% 1
    catch_mult[ac(ages_mult), ac(yrs_mult)] <- catch_mult_data

    ### correct catch.n
    catch.n(stk) <- catch.n(stk) * catch_mult
    # ### split into landings and discards, based on landing fraction
    # ### done later
    # land_frac <- landings.n(stk) / catch.n(stk)
    # landings.n(stk) <- catch.n(stk) * land_frac
    # discards.n(stk) <- catch.n(stk) * (1 - land_frac)
    # ### update stock
    # catch(stk) <- computeCatch(stk)
    # landings(stk) <- computeLandings(stk)
    # discards(stk)<- computeDiscards(stk)

    # ### catch biomass
    # catch <- quantSums(qnt)
    # catch[, 1:n_yrs] <- exp(object$sdrep$value[names(object$sdrep$value) ==
    #                                              "logCatch"])
    # catch(stk) <- catch

  }

  ### total catch
  catch(stk) <- computeCatch(stk)

  ### catch range
  if (isTRUE(uncertainty) & isTRUE(catch_estimate)) {
    catch_sd <- quantSums(qnt)
    catch_sd[, 1:n_yrs] <- object$sdrep$sd[names(object$sdrep$value) == "logCatch"]
    attr(catch(stk), "low") <- exp(log(catch(stk)) - catch_sd * SD_mult)
    attr(catch(stk), "high") <- exp(log(catch(stk)) + catch_sd * SD_mult)
  }

  ### landings
  ### calculate with landings fraction of total catch
  if (isTRUE(length(dim(object$data$landFrac)) == 3)) {
    lfrac <- t(object$data$landFrac[,,, drop = TRUE])
  } else {
    lfrac <- t(object$data$landFrac)
  }
  n_ages <- dim(lfrac)[1]
  n_yrs <- dim(lfrac)[2]
  lfrac_qnt <- qnt ### FLQuant template
  lfrac_qnt[1:n_ages, 1:n_yrs] <- lfrac
  ### calculate landings numbers @ age
  landings.n(stk) <- catch.n(stk) * lfrac_qnt
  ### landings weights @ age
  if (isTRUE(length(dim(object$data$landMeanWeight)) == 3)) {
    landings.wt <- t(object$data$landMeanWeight[,,, drop = TRUE])
  } else {
    landings.wt <- t(object$data$landMeanWeight)
  }
  n_ages <- dim(landings.wt)[1]
  n_yrs <- dim(landings.wt)[2]
  landings.wt(stk)[1:n_ages, 1:n_yrs] <- landings.wt
  ### total landings
  landings(stk) <- computeLandings(stk)
  ### landings range
  if (isTRUE(uncertainty)) {
    ### landing fraction of total catch
    lfrac_qnt_sum <- landings(stk) / catch(stk)
    ### replace NAs with 1
    lfrac_qnt_sum[is.na(lfrac_qnt_sum)] <- 1
    ### calculate range
    ### assume that landing fraction is the same for estimate, low and high
    ### -> split catch into landings and discards
    attr(landings(stk), "low") <- attr(catch(stk), "low") * lfrac_qnt_sum
    attr(landings(stk), "high") <- attr(catch(stk), "high") * lfrac_qnt_sum
  }

  ### discards
  ### calculate discards number @ age
  discards.n(stk) <- catch.n(stk) * (1 - lfrac_qnt)
  ### discards weights @ age
  if (isTRUE(length(dim(object$data$disMeanWeight)) == 3)) {
    discards.wt <- t(object$data$disMeanWeight[,,, drop = TRUE])
  } else {
    discards.wt <- t(object$data$disMeanWeight)
  }
  n_ages <- dim(discards.wt)[1]
  n_yrs <- dim(discards.wt)[2]
  discards.wt(stk)[1:n_ages, 1:n_yrs] <- discards.wt
  ### total discards
  discards(stk) <- computeDiscards(stk)
  ### discard range
  if (isTRUE(uncertainty)) {
    ### calculate range
    ### assume that landing fraction is the same for estimate, low and high
    attr(discards(stk), "low") <- attr(catch(stk), "low") * (1 - lfrac_qnt_sum)
    attr(discards(stk), "high") <- attr(catch(stk), "high") * (1 - lfrac_qnt_sum)
  }

  ### natural mortality
  m <- t(object$data$natMor)
  n_ages <- dim(m)[1]
  n_yrs <- dim(m)[2]
  m(stk)[1:n_ages, 1:n_yrs] <- m

  ### maturity ogive
  mat <- t(object$data$propMat)
  n_ages <- dim(mat)[1]
  n_yrs <- dim(mat)[2]
  mat(stk)[1:n_ages, 1:n_yrs] <- mat

  ### proportion of F before spawning
  if (isTRUE(length(dim(object$data$propF)) == 3)) {
    harvest.spwn <- t(object$data$propF[,,, drop = TRUE])
  } else {
    harvest.spwn <- t(object$data$propF)
  }
  n_ages <- dim(harvest.spwn)[1]
  n_yrs <- dim(harvest.spwn)[2]
  harvest.spwn(stk)[1:n_ages, 1:n_yrs] <- harvest.spwn

  ### proportion of M before spawning
  m.spwn <- t(object$data$propM)
  n_ages <- dim(m.spwn)[1]
  n_yrs <- dim(m.spwn)[2]
  m.spwn(stk)[1:n_ages, 1:n_yrs] <- m.spwn

  ### set description
  desc(stk) <- "FLStock created from SAM model fit"

  ### set range
  ### plusgroup?
  range(stk)["plusgroup"] <- ifelse(isTRUE(object$conf$maxAgePlusGroup == 1),
                                    object$conf$maxAge, NA)
  ### fbar range
  range(stk)[c("minfbar", "maxfbar")] <- object$conf$fbarRange

  ### return FLStock
  return(stk)

}

### function for converting sam_list into FLStock
sam_list_to_FLStock <- function(object, uncertainty = FALSE, conf_level = 95,
                             ...) {

  ### check class type of input
  if (!isTRUE("sam_list" %in% class(object))) {
    stop("object has to be class sam_list")
  }

  ### if only non-sam classes, i.e. errors, stop here
  if (all(!sapply(object, is, "sam"))) {
    stop("input object does not contain any results")
  }

  ### convert all iterations into FLStock objects
  ### if error, return NA
  stk_iter <- lapply(object, function(i) {

    if (is(i, "sam")) {

      tmp <- sam_to_FLStock(object = i, uncertainty = uncertainty,
                         conf_level = conf_level, ...)

    } else {

      tmp <- NA

    }

  })

  ### positions (iters) with data
  pos_data <- which(sapply(stk_iter, is, "FLStock"))


  ### create template FLStock and expand iter dimension
  ### use first element which is FLStock
  ### create 1 iteration more
  stk <- propagate(stk_iter[[head(pos_data, 1)]],
                   iter = length(stk_iter) + 1, fill.iter = FALSE)
  ### delete first iteration
  stk <- FLCore::iter(stk, -1)

  ### insert values for all iterations
  for (i in which(sapply(stk_iter, is, "FLStock"))) {

    FLCore::iter(stk, i) <- stk_iter[[i]]

  }

  ### return FLStock including iterations
  return(stk)

}

### ------------------------------------------------------------------------ ###
### definition of methods for converting SAM results into FLStock ####
### ------------------------------------------------------------------------ ###

#showMethods("FLStock")
#showMethods(class = "sam")
### register "sam" class as formally defined class
setOldClass("sam")
setOldClass("sam_list")

#' Coerce SAM output into \code{FLStock} object
#'
#' This function takes the output from running the SAM stockassessment and
#' converts them into an \code{FLStock} object.
#'
#' If an \code{FLStock} is provided as \code{stk} argument, then only the
#' following slots are updated: \code{stock}, \code{stock.n} & \code{harvest},
#' and, if requested by \code{catch_estimate = TRUE}, \code{catch},
#' \code{catch.n}, \code{landings}, \code{landings.n}, \code{discards},
#' \code{discards.n} with model estimates. If the dimensions
#' (years, iterations) differ between the SAM results and the provided stock
#' template, the returned \code{FLStock} is expanded.
#'
#'
#' @param object Object of class \code{sam} with the results from a
#'   SAM stock assessment run. Alternatively, object of class \code{sam_list},
#'   i.e. a list of \code{sam} objects.
#' @param stk Optional. Object of class \linkS4class{FLStock}, to which the
#'   assessment results are added.
#' @param uncertainty If set to \code{TRUE}, the estimated uncertainty from
#' SAM will be added as attribute.
#' @param conf_level Confidence level used when uncertainty is returned.
#'   Defaults to 95 (percent).
#' @param catch_estimate If set to \code{TRUE}, the catch estimated by SAM will
#' be returned instead of the model input.
#' @param correct_catch If set to \code{TRUE}, returned catch are correct by
#'   catch multiplier estimated by SAM.
#'
#' @return An object of class \code{FLStock}.
#'
#' @examples
#' # fit SAM to North Sea cod
#' fit <- FLR_SAM(stk = cod4_stk, idx = cod4_idx, conf = cod4_conf_sam)
#'
#' # coerce the output into FLStock
#' stk <- SAM2FLStock(fit)
#'
#' # get catch estimates from model
#' stk <- SAM2FLStock(fit, catch_estimate = TRUE)
#'
#' @export

setGeneric("SAM2FLStock", function(object, stk, uncertainty = FALSE,
                                   conf_level = 95, catch_estimate = FALSE,
                                   correct_catch = FALSE) {
  standardGeneric("SAM2FLStock")
})

### object = sam, stk = missing
#' @rdname SAM2FLStock
setMethod(f = "SAM2FLStock",
          signature = signature(object = "sam", stk = "missing"),
          definition = function(object, stk, uncertainty = FALSE,
                                conf_level = 95, catch_estimate = FALSE,
                                correct_catch = FALSE) {

    sam_to_FLStock(object = object, stk = stk, uncertainty = uncertainty,
                   conf_level = conf_level, catch_estimate = catch_estimate,
                   correct_catch = correct_catch)

  })

### object = sam_list, stk = missing
#' @rdname SAM2FLStock
setMethod(f = "SAM2FLStock",
          signature = signature(object = "sam_list", stk = "missing"),
          definition = function(object, stk, uncertainty = FALSE,
                                conf_level = 95, catch_estimate = FALSE,
                                correct_catch = FALSE) {

    sam_list_to_FLStock(object = object, stk = stk, uncertainty = uncertainty,
                        conf_level = conf_level,
                        catch_estimate = catch_estimate,
                        correct_catch = correct_catch)

  })

### object = sam, stk = FLStock
#' @rdname SAM2FLStock
setMethod(f = "SAM2FLStock",
          signature = signature(object = "sam", stk = "FLStock"),
          definition = function(object, stk, uncertainty = FALSE,
                                conf_level = 95, catch_estimate = FALSE,
                                correct_catch = FALSE) {

  ### coerce SAM into FLStock
  stk_new <- sam_to_FLStock(object = object, stk = stk,
                            uncertainty = uncertainty,
                            conf_level = conf_level,
                            catch_estimate = catch_estimate,
                            correct_catch = correct_catch)

  ### adapt dimensions
  stks <- adapt_dims(stk, stk_new, fill.iter = TRUE)
  ### extract from list
  stk <- stks[[1]]
  stk_new <- stks[[2]]

  ### insert values
  stock(stk) <- stock(stk_new)
  stock.n(stk) <- stock.n(stk_new)
  harvest(stk) <- harvest(stk_new)

  ### update catch, if requested
  if (isTRUE(catch_estimate) | isTRUE(correct_catch)) {

    catch(stk) <- catch(stk_new)
    catch.n(stk) <- catch.n(stk_new)
    landings(stk) <- landings(stk_new)
    landings.n(stk) <- landings.n(stk_new)
    discards(stk) <- discards(stk_new)
    discards.n(stk) <- discards.n(stk_new)

  }

  ### return input stock, updated wit requested data
  return(stk)

})

### object = sam_list, stk = FLStock
#' @rdname SAM2FLStock
setMethod(f = "SAM2FLStock",
          signature = signature(object = "sam_list", stk = "FLStock"),
          definition = function(object, stk, uncertainty = FALSE,
                                conf_level = 95, catch_estimate = FALSE,
                                correct_catch = FALSE) {

  ### coerce SAM into FLStock
  stk_new <- sam_list_to_FLStock(object = object, stk = stk,
                                 uncertainty = uncertainty,
                                 conf_level = conf_level,
                                 catch_estimate = catch_estimate,
                                 correct_catch = correct_catch)

  ### adapt dimensions
  stks <- adapt_dims(stk, stk_new, fill.iter = TRUE)
  ### extract from list
  stk <- stks[[1]]
  stk_new <- stks[[2]]

  ### insert values
  stock(stk) <- stock(stk_new)
  stock.n(stk) <- stock.n(stk_new)
  harvest(stk) <- harvest(stk_new)

  ### update catch, if requested
  if (isTRUE(catch_estimate) | isTRUE(correct_catch)) {

    catch(stk) <- catch(stk_new)
    catch.n(stk) <- catch.n(stk_new)
    landings(stk) <- landings(stk_new)
    landings.n(stk) <- landings.n(stk_new)
    discards(stk) <- discards(stk_new)
    discards.n(stk) <- discards.n(stk_new)

  }

  ### return input stock, updated wit requested data
  return(stk)

})

### ------------------------------------------------------------------------ ###
### get SAM parameter values ####
### ------------------------------------------------------------------------ ###

### extract parameters of model fit
sam_getpar <- function(fit) {

  p <- fit$pl
  p$missing <- NULL
  attr(p, "what") <- NULL

  return(p)

}

#' Get parameter estimates from SAM model fit.
#'
#' This function extracts the parameter estimates from a SAM model fit. These
#' are useful e.g. as initial values in subsequent model fits and can improve
#' model convergence/computing time.
#'
#' @param fit A single SAM model fit of class \code{sam} or a list of fits of
#' class \code{sam_lst}.
#'
#' @return A list with the model parameters of the SAM or a list of them in case
#' of several supplied models
#'
#' @examples
#' ### fit SAM to North Sea cod
#' fit <- FLR_SAM(stk = cod4_stk, idx = cod4_idx, conf = cod4_conf_sam)
#'
#' ### extract parameters
#' pars <- getpars(fit)
#'
#' ### use them as starting values
#' fit2 <- FLR_SAM(stk = cod4_stk, idx = cod4_idx, conf = cod4_conf_sam,
#'                 par_ini = pars)
#' @export

setGeneric("getpars", function(fit) {
  standardGeneric("getpars")
})

### fit = sam
#' @rdname getpars
setMethod(f = "getpars",
          signature = signature(fit = "sam"),
          definition = function(fit) {

  sam_getpar(fit)

})
### fit = sam_list
#' @rdname getpars
setMethod(f = "getpars",
          signature = signature(fit = "sam_list"),
          definition = function(fit) {

  foreach(fit_i = fit, .errorhandling = "pass") %do% {

    sam_getpar(fit_i)

  }

})


### ------------------------------------------------------------------------ ###
### extract uncertainty from SAM object ####
### ------------------------------------------------------------------------ ###
### function for creating replicates based on estimation of uncertainty in SAM
### note: SAM works on a log scale and all reported parameters are also on a
###       scale, even standard deviations.
###       Therefore, the values returned from this function are exponentiated


#' Create replicates/iterations of SAM model fit based on variance-covariance
#' matrix
#'
#' This function use the uncertainty estimated by SAM to create
#' replicates/iterations of assessment results. The function uses the
#' variance-covariance matrix to quantify uncertainty.
#'
#' The returned objects are \code{FLQuant}s where the iteration dimension contains the
#' replicates. Each replicate is internally consistent, e.g. the fishing
#' mortality matches the stock numbers of the same replicate.
#'
#' The following metrics are returned:
#' \itemize{
#'   \item{\code{stock.n}} Stock numbers at age for all years,
#'                         \code{FLQuant}
#'   \item{\code{harvest}} Fishing mortalities at age for all years,
#'                         \code{FLQuant}
#'   \item{\code{survey_catchability}} Catchability at age for all years for all
#'                                     survey indices, list of \code{FLQuant}s
#'   \item{\code{catch_sd}} Standard deviation of the catch numbers at age,
#'                          time invariant, \code{FLQuant}
#'   \item{\code{survey_sd}} Standard deviation of all surveys at age, time
#'                           invariant, list of \code{FLQuants},
#'   \item{\code{survey_cov}} Covariance matrices of survey ages, one for each
#'     survey. Return object is a list of lists, the first level corresponds to
#'     the replicates, the second level to the surveys. If no covariance
#'     between ages is assumed in the SAM model, the diagonal in the covariance
#'     matrices is simply the square root of the standard deviation
#'     (in \code{survey_sd})
#'   \item{\code{proc_error}} Standard deviation of the stock numbers at age,
#'     time invariant, class \code{FLQuant}. This corresponds to the survival
#'     process error assumed/estimated, i.e. quantifies how much the actual
#'     stock numbers at age deviate from the deterministic catch equation.
#'   \item{\code{catch_n}} Estimates of catch numbers at age for all years. This
#'     differs from the assessment input values. If the SAM model fit contains
#'     catch multipliers, the values returned here are corrected for this.
#'     Class \code{FLQuant}.
#' }
#'
#'
#' @param fit A SAM model fit object of class \code{sam}.
#' @param n Number of replicates
#' @param print_screen If set to \code{TRUE}, print output of \code{TMB::sdreport} to screen.
#' @param idx_cov If set to \code{TRUE}, return covariance of survey index/indices.
#' @param catch_est If set to \code{TRUE}, return catch estimates from SAM.

#'
#' @return A list of FLQuants with the elements: stock.n, harvest,
#'   survey_catchability, catch_sd, survey_sd, survey_cov, proc_error, catch_n.
#'
#' @examples
#' ### fit SAM to North Sea cod
#' fit <- FLR_SAM(stk = cod4_stk, idx = cod4_idx, conf = cod4_conf_sam)
#'
#' ### create 2 replicates
#' reps <- SAM_uncertainty(fit, n = 2, idx_cov = TRUE, catch_est = TRUE)
#' @export

setGeneric("SAM_uncertainty", function(fit, n = 1000, print_screen = FALSE,
                                       idx_cov = FALSE, catch_est = FALSE) {
  standardGeneric("SAM_uncertainty")
})

### fit = sam
#' @rdname SAM_uncertainty
setMethod(f = "SAM_uncertainty",
          signature = signature(fit = "sam"),
          definition = function(fit, n = 1000, print_screen = FALSE,
                                idx_cov = FALSE, catch_est = FALSE) {

  ### check if required package is available
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop(paste("Package 'TMB' needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }

  ### index for fishing mortality ages
  idxF <- fit$conf$keyLogFsta[1, ] + 1# +dim(stk)[1]
  idxF <- idxF[idxF != 0] ### remove 0s

  ### index for F variances (usually some ages are bound)
  #idxNVar <- fit$conf$keyVarLogN

  ### get ages used for calculating fbar
  #bAges <- fit$conf$fbarRange
  #bAges <- do.call(':', as.list(bAges))

  ### index for stock numbers ages
  #idxN <- 1:ncol(natural.mortality)
  #idxN <- seq(min(fit$data$minAgePerFleet), max(fit$data$maxAgePerFleet))
  ### index for observation variances
  idxObs <- fit$conf$keyVarObs # starts at 0

  ##Resample estimated values to get N, F and q

  ### calculate standard deviations of model parameters
  . <- capture.output(sds <- TMB::sdreport(obj = fit$obj,
                                           par.fixed = fit$opt$par,
                                           getJointPrecision = TRUE))
  if (isTRUE(print_screen)) cat(paste0(., sep = "\n"))

  ### extract values for parameters
  est <- c(sds$par.fixed, sds$par.random)
  ### get covariance matrix of all model parameters
  cov <- solve(sds$jointPrecision)

  ### create random values based on estimation and covariance
  sim.states <- mvrnorm(n, est, cov) ### REPLACE with rmvnorm()
  ### matrix, columns are values, rows are requested samples
  table(colnames(sim.states))
  ### contains, among others, logF, logN...

  ### combine SAM estimate and random samples
  #dat <- rbind(est, sim.states)
  dat <- sim.states

  ### workaround if only one iteration/replicate
  if (!is.matrix(dat)) dat <- t(as.matrix(dat))

  ### ---------------------------------------------------------------------- ###
  ### stock ####
  ### ---------------------------------------------------------------------- ###

  ### stock characteristics
  min_age <- min(fit$data$minAgePerFleet[fit$data$fleetTypes == 0])
  max_age <- max(fit$data$maxAgePerFleet[fit$data$fleetTypes == 0])
  years <- fit$data$years
  ### FLQuant template for stock
  stk_template <- FLQuant(dimnames = list(age = min_age:max_age, year = years,
                                          iter = 1:n))

  ### numbers at age
  stock.n <- stk_template
  stock.n[] <- exp(t(dat[, colnames(dat) == "logN"]))

  ### F at age
  harvest <- stk_template
  ### insert values for estimated ages
  harvest[unique(idxF)] <- exp(t(dat[, colnames(dat) == "logF"]))
  ### duplicate ages, if some of them are bound
  harvest[] <- harvest[idxF]

  ### ---------------------------------------------------------------------- ###
  ### surveys ####
  ### ---------------------------------------------------------------------- ###

  ### survey specs
  idx_surveys <- which(fit$data$fleetTypes > 0) ### which observation are surveys
  ### age range of surveys
  survey_ages <- lapply(seq_along(idx_surveys), function(x) {
    seq(fit$data$minAgePerFleet[idx_surveys][x],
        fit$data$maxAgePerFleet[idx_surveys][x])
  })
  ### index for estimated parameters
  idx_LogFpar <- lapply(idx_surveys, function(x) {
    tmp <- fit$conf$keyLogFpar[x, ] + 1
    tmp <- tmp[tmp > 0]
  })

  sum(colnames(dat) == "logFpar") ### there are 9 parameters for cod
  ### 5 for Q1 (ages 1-5), 4 for Q3 (ages 1-4).
  survey_ages_idx <- split(seq(length(unlist(survey_ages))),
                           rep(seq(survey_ages), sapply(survey_ages, length)))

  ### get catchability at age (time-invariant) samples
  catchability <- lapply(seq_along(idx_surveys), function(x) {

    ### create FLQuant template
    tmp <- FLQuant(dimnames = list(age = survey_ages[[x]],
                                   year = "all", iter = 1:n))
    ### fill with catchability values
    tmp[] <-
      exp(t(dat[, colnames(dat) == "logFpar", drop = FALSE][, idx_LogFpar[[x]]]))

    return(tmp)

  })

  ### ---------------------------------------------------------------------- ###
  ### standard deviation - catch ####
  ### ---------------------------------------------------------------------- ###
  ### time-invariant

  ### template
  catch_sd <- FLQuant(dimnames = list(age = dimnames(stock.n)$age, year = "all",
                                      iter = 1:n))
  # sum(colnames(dat) == "logSdLogObs")

  ### index for catch sd (some ages are linked)
  catch_sd_idx <- idxObs[1, ][idxObs[1, ] > -1] + 1

  ### extract values
  catch_sd[] <-
    exp(t(dat[, colnames(dat) == "logSdLogObs", drop = FALSE][, catch_sd_idx]))
  ### logSdLogObs is the log of the SD of the log observations
  ### exponentiate values here to get SD (and not logSD)

  ### ---------------------------------------------------------------------- ###
  ### standard deviation - surveys ####
  ### ---------------------------------------------------------------------- ###
  ### time-invariant

  ### get catchability at age (time-invariant) samples
  survey_sd <- lapply(seq_along(idx_surveys), function(x) {

    ### create FLQuant template
    tmp <- FLQuant(dimnames = list(age = survey_ages[[x]],
                                   year = "all", iter = 1:n))
    ### index for sd (some ages are linked)
    idx_sd <- idxObs[idx_surveys[x], ]
    idx_sd <- idx_sd[idx_sd > -1] + 1

    ### fill with catchability values
    tmp[] <- exp(t(dat[, colnames(dat) == "logSdLogObs", drop = FALSE][, idx_sd]))

    return(tmp)

  })

  ### ---------------------------------------------------------------------- ###
  ### surveys - get covariance ####
  ### ---------------------------------------------------------------------- ###

  if (isTRUE(idx_cov)) {
    . <- capture.output(survey_cov <- foreach(iter_i = 1:n) %do% {

      fit$obj$fn(sim.states[iter_i, 1:length(sds$par.fixed)])
      cov <- fit$obj$report()$obsCov

      return(cov[idx_surveys])

    })

  } else {

    survey_cov <- NULL

  }

  ### ---------------------------------------------------------------------- ###
  ### stock.n uncertainty ####
  ### ---------------------------------------------------------------------- ###
  ### process error

  ### get index for ages
  idx_SdLogN <- fit$conf$keyVarLogN + 1

  ### FLQuant template
  SdLogN <- catch_sd %=% NA_integer_

  ### fill
  SdLogN[] <-
    exp(t(dat[, colnames(dat) == "logSdLogN", drop = FALSE][, idx_SdLogN]))

  ### ---------------------------------------------------------------------- ###
  ### get catch estimates ####
  ### ---------------------------------------------------------------------- ###

  if (isTRUE(catch_est)) {

    ### catch fleet index/indices
    catch_fleets <- which(fit$data$fleetTypes == 0)
    catch_desc <- fit$data$aux


    . <- capture.output(res_n <- foreach(iter_i = 1:n, .combine = c
    ) %do% {

      fit$obj$fn(sim.states[iter_i, 1:length(sds$par.fixed)])
      tmp <- cbind(fit$data$aux, est = fit$obj$report()$predObs)
      tmp <- tmp[tmp[, 2] == catch_fleets, ]
      tmp <- stats::aggregate(est ~ age + year, tmp,
                              FUN = sum)
      tmp_full <- expand.grid(age = unique(tmp$age),
                              year = unique(tmp$year))
      tmp <- merge(x = tmp, y = tmp_full, all = TRUE)
      tmp <- tmp[order(tmp$year, tmp$age), ]
      return(tmp$est)

    })
    catch_n <- FLQuant(NA,
                       dimnames = list(year = rownames(fit$data$catchMeanWeight),
                                       age = colnames(fit$data$catchMeanWeight),
                                       iter = seq(n)))
    catch_n[] <- exp(res_n)

    ### apply catch multiplier, if used
    if (fit$conf$noScaledYears > 0) {

      ### get catch multiplier dimensions
      ages_mult <- fit$conf$minAge:fit$conf$maxAge
      yrs_mult <- fit$conf$keyScaledYears
      ### get multiplier positions in created data set
      pos_mult <- which(colnames(dat) == "logScale")
      ### create catch multiplier FLQuant
      catch_mult_data <- FLQuant(NA,
                                 dimnames = list(year = fit$conf$keyScaledYears,
                                                 age = fit$conf$minAge:fit$conf$maxAge,
                                                 iter = seq(nrow(dat))))
      catch_mult_data[] <- rep(t(dat[, pos_mult]), each = length(ages_mult))
      ### model values in log scale, exponentiate
      catch_mult_data <- exp(catch_mult_data)
      ### for simpler calculations, use catch_n as template
      catch_mult <- catch_n %=% 1
      catch_mult[ac(ages_mult), ac(yrs_mult)] <- catch_mult_data

      ### correct catch.n
      catch_n <- catch_n * catch_mult

    }

  } else {

    catch_n <- NULL

  }

  return(list(stock.n = stock.n, harvest = harvest,
              survey_catchability = catchability, catch_sd = catch_sd,
              survey_sd = survey_sd, survey_cov = survey_cov,
              proc_error = SdLogN, catch_n = catch_n))

})
