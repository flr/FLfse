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
FLR_SAM_run <- function(stk, idx, conf = NULL) {

  ### check if required package is available
  if (!requireNamespace("stockassessment", quietly = TRUE)) {
    stop(paste("Package 'stockassessment' needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }

  ### calculate landings fraction
  lf <- landings.n(stk) / catch.n(stk)
  ### assume only landings if NA
  lf[is.na(lf) | is.infinite(lf)] <- 1

  ### do some sanity checks ...
  ### fill empty landings/discards weights

  ### create SAM input objects
  dat_lst <- list(residual.fleet = catch.n(stk),
                  prop.mature = mat(stk),
                  stock.mean.weight = stock.wt(stk),
                  catch.mean.weight = catch.wt(stk),
                  dis.mean.weight = discards.wt(stk),
                  land.mean.weight = landings.wt(stk),
                  natural.mortality = m(stk),
                  prop.f = harvest.spwn(stk),
                  prop.m = m.spwn(stk),
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
  dat_idx <- lapply(idx, function(x){
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
  if (all(!is.na(range(stk)[c("minfbar", "maxfbar")]))) {
    conf_sam$fbarRange <- range(stk)[c("minfbar", "maxfbar")]
  }

  ### insert configuration, if supplied to function
  if (!is.null(conf)) {

    ### find slots that can be used
    conf_names <- intersect(names(conf), names(conf_sam))

    ### go through slot names
    if (length(conf_names) > 0) {

      for(i in conf_names) {

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
#' @section Warning:
#' This methods requires the \code{stockassessment} package and all its
#' dependencies to be installed. For details how to obtain
#' \code{stockassessment}, see \url{https://github.com/fishfollower/SAM/}.
#'
#' @param stk Object of class \linkS4class{FLStock} with stock and fishery data.
#' @param idx Object of class \linkS4class{FLIndices} or \linkS4class{FLIndex}
#'   object with survey index time series.
#' @param conf Optional configurations passed to SAM. Should be a list.
#'
#' @return An object of class \code{spictcls} with the model results.
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

setGeneric("FLR_SAM", function(stk, idx, conf = NULL) {
  standardGeneric("FLR_SAM")
})

### stk = FLStock, idx = FLIndices
#' @rdname FLR_SAM
setMethod(f = "FLR_SAM",
          signature = signature(stk = "FLStock", idx = "FLIndices"),
          definition = function(stk, idx, conf = NULL) {

  FLR_SAM_run(stk = stk, idx = idx, conf = conf)

})
### stk = FLStock, idx = FLIndex
#' @rdname FLR_SAM
setMethod(f = "FLR_SAM",
          signature = signature(stk = "FLStock", idx = "FLIndex"),
          definition = function(stk, idx, conf = NULL) {

  ### coerce FLIndex into FLIndices
  idx <- FLIndices(idx)

  ### run SPiCT
  FLR_SAM_run(stk = stk, idx = idx, conf = conf)

})


### ------------------------------------------------------------------------ ###
### convert sam object into FLStock ####
### ------------------------------------------------------------------------ ###

### function for converting sam into FLStock
sam2FLStock <- function(object, uncertainty = FALSE, conf_level = 95, ...) {

  ### check class type of input
  if (!isTRUE("sam" %in% class(object))) stop("object has to be class sam")

  ### calculate SD multiplier for uncertainty range
  if (isTRUE(uncertainty)) SD_mult <- stats::qnorm(0.5+conf_level/200)

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
  harvest <- harvest[object$conf$keyLogFsta[1,]+1, ]
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

  ### insert estimates
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
                         nrow = length(object$conf$minAge:object$conf$maxAge), byrow = TRUE)
    ### exponentiate
    catch_mult <- exp(catch_mult)
    ### coerce into FLQuant
    catch_qnt <- qnt ### dummy object
    catch_qnt[, ac(object$conf$keyScaledYears)] <- catch_mult
    ### replace NA with 1
    catch_qnt[is.na(catch_qnt)] <- 1
    ### multiply catch numbers
    catch.n(stk) <- catch.n(stk) * catch_qnt
  }
  ### save observations as attribute
  catch.n_obs <- qnt
  catch.n_obs[ac(unique(dat_catch$age)), ac(unique(dat_catch$year))] <-
    dat_catch$value
  attr(catch.n(stk), "observations") <- catch.n_obs
  ### catch biomass
  catch <- quantSums(qnt)
  catch[, 1:n_yrs] <- exp(object$sdrep$value[names(object$sdrep$value) == "logCatch"])
  catch(stk) <- catch
  ### total catch observations
  attr(catch(stk), "observations") <- quantSums(catch.wt(stk) *
                                                  attr(catch.n(stk), "observations"))

  ### catch range
  if (isTRUE(uncertainty)) {
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
  ### discardrange
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

# ### function for converting sam into FLStock
# sam2FLStock <- function(object, uncertainty = FALSE, ...) {
#
#   ### get dimensions
#   years <- object$data$years
#   ages <- object$conf$minAge:object$conf$maxAge
#
#   ### create FLQuant dummy
#   qnt <- FLQuant(NA, dimnames = list(age = ages, year = years))
#   ### create FLStock
#   stk <- FLStock(qnt)
#
#   ### fishing mortality
#   stk_lst <- list(harvest = faytable(object))
#
#
#   harvest <- matrix2FLQuant(faytable(object))
#   harvest(stk)[] <- harvest
#
#   #exp(object$pl$logF)
#
#   ### stock numbers @ age
#   stock.n <- exp(object$pl$logN)
#   n_ages <- dim(stock.n)[1]
#   n_yrs <- dim(stock.n)[2]
#   stock.n(stk)[1:n_ages, 1:n_yrs] <- stock.n
#
#   ### harvest @ age
#   harvest <- exp(object$pl$logF)
#   ### duplicate linked ages
#   harvest <- harvest[object$conf$keyLogFsta[1,]+1, ]
#   n_ages <- dim(harvest)[1]
#   n_yrs <- dim(harvest)[2]
#   harvest(stk)[1:n_ages, 1:n_yrs] <- harvest
#   units(harvest(stk)) <- "f"
#
#   ### ---------------------------------------------------------------------- ###
#   ### input data
#
#   ### stock
#   ### stock weight @ age
#   stock.wt <- t(object$data$stockMeanWeight)
#   n_ages <- dim(stock.wt)[1]
#   n_yrs <- dim(stock.wt)[2]
#   stock.wt(stk)[1:n_ages, 1:n_yrs] <- stock.wt
#   ### stock biomass
#   stock(stk) <- computeStock(stk)
#
#   ### catch
#   ### catch weight @ age
#   catch.wt <- t(object$data$catchMeanWeight)
#   n_ages <- dim(catch.wt)[1]
#   n_yrs <- dim(catch.wt)[2]
#   catch.wt(stk)[1:n_ages, 1:n_yrs] <- catch.wt
#   ### catch numbers @ age
#   catch_fleets <- which(object$data$fleetTypes == 0)
#   ### extract observations
#   dat_catch <- cbind(object$data$aux, value = exp(object$data$logobs))
#   dat_catch <- dat_catch[dat_catch[, "fleet"] == catch_fleets, ]
#   ### sum up catch numbers over all commercial fleets
#   dat_catch <- aggregate(value ~ age + year, dat_catch,  FUN = sum)
#   ### insert
#   catch.n(stk)[ac(unique(dat_catch$age)), ac(unique(dat_catch$year))] <- dat_catch$value
#   ### catch biomass
#   catch(stk) <- computeCatch(stk)
#
#   ### landings
#   ### calculate with landings fraction of total catch
#   lfrac <- t(object$data$landFrac)
#   n_ages <- dim(lfrac)[1]
#   n_yrs <- dim(lfrac)[2]
#   lfrac_qnt <- landings.n(stk) ### FLQuant template
#   lfrac_qnt[1:n_ages, 1:n_yrs] <- lfrac
#   ### calculate landings numbers @ age
#   landings.n(stk) <- catch.n(stk) * lfrac_qnt
#   ### landings weights @ age
#   landings.wt <- t(object$data$landMeanWeight)
#   n_ages <- dim(landings.wt)[1]
#   n_yrs <- dim(landings.wt)[2]
#   landings.wt(stk)[1:n_ages, 1:n_yrs] <- landings.wt
#   ### total landings
#   landings(stk) <- computeLandings(stk)
#
#   ### discards
#   ### calculate discards number @ age
#   discards.n(stk) <- catch.n(stk) * (1 - lfrac_qnt)
#   ### discards weights @ age
#   discards.wt <- t(object$data$disMeanWeight)
#   n_ages <- dim(discards.wt)[1]
#   n_yrs <- dim(discards.wt)[2]
#   discards.wt(stk)[1:n_ages, 1:n_yrs] <- discards.wt
#   ### total landings
#   discards(stk) <- computeDiscards(stk)
#
#   ### natural mortality
#   m <- t(object$data$natMor)
#   n_ages <- dim(m)[1]
#   n_yrs <- dim(m)[2]
#   m(stk)[1:n_ages, 1:n_yrs] <- m
#
#   ### maturity ogive
#   mat <- t(object$data$propMat)
#   n_ages <- dim(mat)[1]
#   n_yrs <- dim(mat)[2]
#   mat(stk)[1:n_ages, 1:n_yrs] <- mat
#
#   ### proportion of F before spawning
#   harvest.spwn <- t(object$data$propF)
#   n_ages <- dim(harvest.spwn)[1]
#   n_yrs <- dim(harvest.spwn)[2]
#   harvest.spwn(stk)[1:n_ages, 1:n_yrs] <- harvest.spwn
#
#   ### proportion of M before spawning
#   m.spwn <- t(object$data$propM)
#   n_ages <- dim(m.spwn)[1]
#   n_yrs <- dim(m.spwn)[2]
#   m.spwn(stk)[1:n_ages, 1:n_yrs] <- m.spwn
#
#   ### set description
#   desc(stk) <- "FLStock created from SAM model fit"
#
#   ### set range
#   ### plusgroup?
#   range(stk)["plusgroup"] <- ifelse(isTRUE(object$conf$maxAgePlusGroup == 1),
#                                     object$conf$maxAge, NA)
#   ### fbar range
#   range(stk)[c("minfbar", "maxfbar")] <- object$conf$fbarRange
#
#   ### iters
#   ### input data      1
#   ### model estimate  2
#   ### sd              3
#
#   ### input & model estimate
#   ### catch, catch.n, landings, landings.n, discards, discards.n
#
#   ### input only
#   ### m, mat, catch.wt, landings.wt, discards.wt, stock.wt, harvest.spwn,
#   ### m.spwn
#
#   ### what if the input stock has > 1 iteration?
#
#   ### estimate only
#   ### harvest, stock, stock.n
#
#   ### return FLStock
#   return(stk)
#
# }

#showMethods("FLStock")
#showMethods(class = "sam")
#setGeneric(name = "FLStock", def = function(object, ...) standardGeneric("FLStock"))
### register "sam" class a formally defined class
setOldClass("sam")
setMethod(f = "FLStock", signature = signature(object = "sam"),
          definition = sam2FLStock)

