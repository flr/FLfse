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

### load FLR objects and run SAM
### function not exported, use FLR_SAM instead
FLR_SAM_run <- function(stk, idx, conf = NULL,
                        force_list_output = FALSE,
                        DoParallel = FALSE ### compute iterations in parallel
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

  ### set parallel or serial processing of iterations
  `%do_tmp%` <- ifelse(isTRUE(DoParallel), `%dopar%`, `%do%`)

  ### go through iterations
  res_iter <- foreach(i = seq(dims(stk)$iter), .errorhandling = "pass",
                      .packages = c("FLCore", "stockassessment")) %do_tmp% {

    ### subset stock and index to current iter
    stk_i <- FLCore::iter(stk, i)
    idx_i <- lapply(idx, FLCore::iter, i)

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
        if (isTRUE(!is.numeric(dimnames(tmp)[[2]])) | ### non numeric ages
                   is.na(dimnames(tmp)[[2]]) | ### NA as age
                   dimnames(tmp)[[2]] == -1) {
          dimnames(tmp)[[2]] <- -1 ### recognized by SAM as SSB index
        }
      }

      ### add survey timing as attribute
      attr(tmp, "time") <- as.vector(range(x)[c("startf", "endf")])
      return(tmp)
    })

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
    if (!is.null(conf)) {

      ### find slots that can be used
      conf_names <- intersect(names(conf), names(conf_sam))

      ### go through slot names
      if (length(conf_names) > 0) {

        for (i in conf_names) {

          ### workaround for keyParScaledYA
          ### default 0x0 matrix, needs to extended
          if (i == "keyParScaledYA") {

            conf_sam[[i]] <- conf[[i]]

            ### otherwise enter only values
          } else {

            ### position where values supplied (i.e. not NA)
            pos <- NULL ### reset from previous iteration
            pos <- which(!is.na(conf[[i]]))
            ### insert values
            conf_sam[[i]][pos] <- conf[[i]][pos]

          }

        }

      }
    }

    ### define parameters for SAM
    par <- stockassessment::defpar(dat_sam, conf_sam)

    ### run SAM
    sam_msg <- capture.output(fit <- stockassessment::sam.fit(dat_sam, conf_sam,
                                                              par))

    ### save screen message(s) as attribute
    attr(x = fit, which = "messages") <- sam_msg

    ### return
    return(fit)

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
#' Trailing years without data in \code{stk} are removed.
#'
#' Survey indices are extracted from \code{idx}, using the \code{index}
#' slot(s).
#' SSB indices can be used. To define an index as SSB index, the first (age)
#' dimension of the \code{index} slot of \code{idx} has to be of length 1 and
#' the name of this dimension can be either missing (\code{NA}), non-numeric
#' or \code{-1}.
#'
#' Additional configurations can be passed as a list to SAM with the
#' \code{conf} argument. \code{FLR_SAM} first generates a default model
#' configuration with \code{stockassessment}'s \code{setup.sam.data}. If
#' additional configurations are passed to \code{FLR_SAM}, they will replace
#' the default configuration. For details about possible configurations and
#' format see '\code{help("defcon", package = "stockassessment")} and
#' \url{https://github.com/fishfollower/SAM}.
#'
#' The function can handle input objects with multiple iterations (\code{iter}
#' dimension in \code{stk} and \code{idx}). If multiple iterations are provided
#' for \code{stk} but not for \code{idx}, \code{idx} will be inflated, and vice
#'  versa.
#' If the assessment fails for some iterations, the error messages are return
#' for these iterations.
#' If argument \code{DoParallel} is set to \code{TRUE}, the individual
#' iterations are processed in parallel. This uses parallel computing provided
#' by the package \code{DoParallel}. The parallel workers need to be set up
#' and registered before calling the function. See \code{?DoParallel}.
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
#' @param DoParallel Optional, defaults to \code{FALSE}. If set to \code{TRUE},
#'   will perform iterations of stock in parallel. See Details below for
#'   description.
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

setGeneric("FLR_SAM", function(stk, idx, conf = NULL, DoParallel = FALSE) {
  standardGeneric("FLR_SAM")
})

### stk = FLStock, idx = FLIndices
#' @rdname FLR_SAM
setMethod(f = "FLR_SAM",
          signature = signature(stk = "FLStock", idx = "FLIndices"),
          definition = function(stk, idx, conf = NULL, DoParallel = FALSE) {

  FLR_SAM_run(stk = stk, idx = idx, conf = conf, DoParallel = DoParallel)

})
### stk = FLStock, idx = FLIndex
#' @rdname FLR_SAM
setMethod(f = "FLR_SAM",
          signature = signature(stk = "FLStock", idx = "FLIndex"),
          definition = function(stk, idx, conf = NULL, DoParallel = FALSE) {

  ### coerce FLIndex into FLIndices
  idx <- FLIndices(idx)

  ### run SPiCT
  FLR_SAM_run(stk = stk, idx = idx, conf = conf, DoParallel = DoParallel)

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
  catch.wt <- t(object$data$catchMeanWeight)
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

    catch.n(stk)[ac(unique(dat_catch$age)), ac(unique(dat_catch$year))] <-
      dat_catch$value

    } else {
    catch.n(stk)[ac(unique(dat_catch$age)), ac(unique(dat_catch$year))] <-
      dat_catch$estimate
    ### correct catch numbers if a catch multiplier was estimated
    if (length(object$pl$logScale) > 0) {
      ### catch multipliers
      catch_mult <- object$pl$logScale
      ### replicate multipliers as defined in configuration
      catch_mult <- catch_mult[(object$conf$keyParScaledYA + 1)]
      ### format into matrix
      catch_mult <- matrix(data = catch_mult, ncol = object$conf$noScaledYears,
                           nrow = length(object$conf$minAge:object$conf$maxAge),
                           byrow = TRUE)
      ### exponentiate
      catch_mult <- exp(catch_mult)
      ### coerce into FLQuant
      catch_qnt <- qnt ### dummy object
      catch_qnt[, ac(object$conf$keyScaledYears)] <- catch_mult
      ### replace NA with 1
      catch_qnt[is.na(catch_qnt)] <- 1
      ### multiply catch numbers
      catch.n(stk) <- catch.n(stk) * catch_qnt

      ### catch biomass
      catch <- quantSums(qnt)
      catch[, 1:n_yrs] <- exp(object$sdrep$value[names(object$sdrep$value) ==
                                                   "logCatch"])
      catch(stk) <- catch

    }

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
  lfrac <- t(object$data$landFrac)
  n_ages <- dim(lfrac)[1]
  n_yrs <- dim(lfrac)[2]
  lfrac_qnt <- qnt ### FLQuant template
  lfrac_qnt[1:n_ages, 1:n_yrs] <- lfrac
  ### calculate landings numbers @ age
  landings.n(stk) <- catch.n(stk) * lfrac_qnt
  ### landings weights @ age
  landings.wt <- t(object$data$landMeanWeight)
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
  discards.wt <- t(object$data$disMeanWeight)
  n_ages <- dim(discards.wt)[1]
  n_yrs <- dim(discards.wt)[2]
  discards.wt(stk)[1:n_ages, 1:n_yrs] <- discards.wt
  ### total discards
  discards(stk) <- computeDiscards(stk)
  ### discard range
  if (isTRUE(uncertainty)) {
    ### calculate range
    ### assume that landing fraction is the same for estimate, low and high
    attr(discards(stk), "low") <- attr(catch(stk), "low") * (1-lfrac_qnt_sum)
    attr(discards(stk), "high") <- attr(catch(stk), "high") * (1-lfrac_qnt_sum)
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
  harvest.spwn <- t(object$data$propF)
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
#' @param object Object of class \linkS4class{sam} with the results from a
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
                                   conf_level = 95, catch_estimate = FALSE) {
  standardGeneric("SAM2FLStock")
})

### object = sam, stk = missing
#' @rdname SAM2FLStock
setMethod(f = "SAM2FLStock",
          signature = signature(object = "sam", stk = "missing"),
          definition = function(object, stk, uncertainty = FALSE,
                                conf_level = 95, catch_estimate = FALSE) {

    sam_to_FLStock(object = object, stk = stk, uncertainty = uncertainty,
                   conf_level = conf_level, catch_estimate = catch_estimate)

  })

### object = sam_list, stk = missing
#' @rdname SAM2FLStock
setMethod(f = "SAM2FLStock",
          signature = signature(object = "sam_list", stk = "missing"),
          definition = function(object, stk, uncertainty = FALSE,
                                conf_level = 95, catch_estimate = FALSE) {

    sam_list_to_FLStock(object = object, stk = stk, uncertainty = uncertainty,
                        conf_level = conf_level, catch_estimate = catch_estimate)

  })

### object = sam, stk = FLStock
#' @rdname SAM2FLStock
setMethod(f = "SAM2FLStock",
          signature = signature(object = "sam", stk = "FLStock"),
          definition = function(object, stk, uncertainty = FALSE,
                                conf_level = 95, catch_estimate = FALSE) {

  ### coerce SAM into FLStock
  stk_new <- sam_to_FLStock(object = object, stk = stk,
                            uncertainty = uncertainty,
                            conf_level = conf_level,
                            catch_estimate = catch_estimate)

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
  if (isTRUE(catch_estimate)) {

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
                                conf_level = 95, catch_estimate = FALSE) {

  ### coerce SAM into FLStock
  stk_new <- sam_list_to_FLStock(object = object, stk = stk,
                                 uncertainty = uncertainty,
                                 conf_level = conf_level,
                                 catch_estimate = catch_estimate)

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
  if (isTRUE(catch_estimate)) {

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
