## Function to parse a SAM forecast into in FLStock object

## add_FRL_forecast: flag to add a deterministic forecast in FLR using FLasher
## include base year: start the FLR forecast in the base year (often last assessment year), this will overwrite stock numbers, weights, etc. in the base year
## add_deterministic_SAM: add a SAM forecast without process error in recruitment, F and N processes
## geomean: use a geometric mean of historic recruitment values or an arithmetic mean (the latter is similar to a deterministic SAM forecast)
## keep_SAM_catch_weights: future catch weights may differ between SAM and FLasher::fwd
# SAM catchweights in projection years = average(catch weights in historic years over which to average)
# FLR catchweights in projection years = average(landings.sel) + average(landings.wt) + average(discards.sel) + average(discards.wt)SAMforecast2FLStock <- function(sam_forecast, add_FLR_forecast = T, include_base_year = F, add_deterministic_SAM = T, geomean = T, keep_SAM_catch_weights = T,...){

SAMforecast2FLStock <- function(sam_forecast, add_FLR_forecast = T, include_base_year = F, add_deterministic_SAM = T, geomean = T, keep_SAM_catch_weights = T,...){

  forecast_args <- attributes(sam_forecast)$forecast_arguments
  year.base     <- forecast_args$year.base
  proj.yrs      <- colnames(attr(sam_forecast,"shorttab"))
  ave.years     <- forecast_args$ave.years
  rec.years     <- forecast_args$rec.years
  nosim         <- forecast_args$nosim
  fit           <- attr(sam_forecast,"fit")
  max_yr        <- max(fit$data$years)

  if(!forecast_args$savesim){
    cat("running SAM forecast with simulation output \n")
    input <- forecast_args
    input <- input[!unlist(lapply(input, is.null))]
    input$savesim <- T
    input$fit <- attr(sam_forecast,"fit")
    input$estimate <- as.name(attr(forecast_args$estimate,"estimateLabel"))
    if(is.null(attr(forecast_args,"RNG_state")))rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed",attr(forecast_args,"RNG_state"), envir = .GlobalEnv)
    frc <- do.call("forecast",input)
    if(all.equal(attr(frc,"tab")[,1:24],attr(sam_forecast,"tab")[,1:24])){
      cat("DONE!")
      sam_forecast <- frc
    }
  }

  # create FLStock object from SAM fit
  Fkeys      <- fit$conf$keyLogFsta[1,] + 1
  FLstk      <- FLfse::SAM2FLStock(fit, catch_estimate = T)
  Flast      <- c(fbar(FLstk)[,ac(year.base)])

  FLstk_stf  <- trim(FLstk, year = range(FLstk)["minyear"]:(year.base-1))

  # create and fill forecast object with forecast assumptions on biology and fishery
  FLstk_stf <- stf(object = FLstk_stf, nyears = length(proj.yrs), wts.nyears = length(ave.years),
                   fbar.nyears = 1, f.rescale = TRUE, disc.nyears = length(ave.years))

  FLstk_stf@harvest[,proj.yrs]        <- FLstk@harvest[,ac(base_year)]

  FLstk_stf@mat[,proj.yrs]            <- yearMeans(FLstk@mat[,as.character(ave.years)])
  FLstk_stf@m[,proj.yrs]              <- yearMeans(FLstk@m[,as.character(ave.years)])
  FLstk_stf@harvest.spwn[,proj.yrs]   <- yearMeans(FLstk@harvest.spwn[,as.character(ave.years)])
  FLstk_stf@m.spwn[,proj.yrs]         <- yearMeans(FLstk@m.spwn[,as.character(ave.years)])
  FLstk_stf@stock.wt[,proj.yrs]       <- yearMeans(FLstk@stock.wt[,as.character(ave.years)])

  FLstk_stf@catch.wt[,proj.yrs]       <- yearMeans(FLstk@catch.wt[,as.character(ave.years)])
  FLstk_stf@discards.wt[,proj.yrs]    <- yearMeans(FLstk@discards.wt[,as.character(ave.years)])
  FLstk_stf@landings.wt[,proj.yrs]    <- yearMeans(FLstk@landings.wt[,as.character(ave.years)])

  FLstk_stf@discards.n[,proj.yrs]    <- yearMeans(FLstk@discards.n[,as.character(ave.years)]/FLstk@catch.n[,as.character(ave.years)])
  FLstk_stf@landings.n[,proj.yrs]    <- yearMeans(FLstk@landings.n[,as.character(ave.years)]/FLstk@catch.n[,as.character(ave.years)])

  if(keep_SAM_catch_weights){
    ## modify discard weights so the catch weights match with SAM
    FLstk_stf@discards.wt[,proj.yrs]   <- (FLstk_stf@catch.wt[,proj.yrs] - FLstk_stf@landings.n[,proj.yrs] * FLstk_stf@landings.wt[,proj.yrs])/FLstk_stf@discards.n[,proj.yrs]
  }

  FLstk_sam <- FLstk_stf

  FLstk_fwd <- NULL

  if(add_FLR_forecast){

    if(!include_base_year){
      FLstk_stf[,ac(base_year)]                <- FLstk[,ac(base_year)]
    }
    # define forecast scenario
    sc      <- do.call("rbind", forecast_args[c("catchval","fscale","catchval.exact","fval","nextssb","landval")])
    values  <- apply(sc,2,function(x)x[!is.na(x)])
    quants  <- apply(sc,2,function(x)rownames(sc)[!is.na(x)])
    quants[quants=="fscale"]         = "f"
    quants[quants=="catchval"]       = "catch"
    quants[quants=="catchval.exact"] = "catch"
    quants[quants=="fval"]           = "fbar"
    quants[quants=="nextssb"]        = "ssb"
    quants[quants=="landval"]        = "landings"

    # year.base, same fishing mortality
    quants[1] <- "fbar"
    values[1] <- Flast * values[1]

    relYear   <- rep(NA, length(proj.yrs))
    relYear[which(quants=="f")] <- as.numeric(proj.yrs[which(quants=="f")])-1
    ctrl <- fwdControl(
      data.frame(
        year = proj.yrs,
        value = values,
        quant = quants,
        relYear = relYear
      )[c(include_base_year,rep(T,length(proj.yrs)-1)),]
    )
    if(geomean){
      R <- c(rectable(fit)[ac(proj.yrs[1]),1],rep(exp(mean(log(rectable(fit)[as.character(rec.years),"Estimate"]))), length(proj.yrs)-1))[c(include_base_year,rep(T,length(proj.yrs)-1))]
      srPar <- FLPar(R,
                     dimnames = list(params="a", year = unlist(proj.yrs), iter = 1))
      srMod <- FLSR(model = "geomean", params = srPar)

      # projection
      FLstk_fwd <- FLasher::fwd(object = FLstk_stf, control = ctrl, sr = srMod)
    } else {
      R <- c(rectable(fit)[ac(proj.yrs[1]),1],rep(mean(rectable(fit)[as.character(rec.years),"Estimate"]), length(proj.yrs)-1))[c(include_base_year,rep(T,length(proj.yrs)-1))]
      srPar <- FLPar(R,
                     dimnames = list(params="a", year = unlist(proj.yrs), iter = 1))
      srMod <- FLSR(model = "geomean", params = srPar)

      # projection
      FLstk_fwd <- FLasher::fwd(object = FLstk_stf, control = ctrl, sr = srMod)
    }
  }

  FLstk_stf      <- propagate(FLstk_sam, iter = nosim, fill.iter = T)

  for(i in proj.yrs){
    sims_yr                  <- sam_forecast[proj.yrs == i]
    FLstk_stf@catch.n[,i]    <- sims_yr[[1]]$catchatage
    FLstk_stf@landings.n[,i] <- FLstk_stf@catch.n[,i] * yearMeans(FLstk@landings.n[,as.character(ave.years)]/FLstk@catch.n[,as.character(ave.years)])
    FLstk_stf@discards.n[,i] <- FLstk_stf@catch.n[,i] - FLstk_stf@landings.n[,i]

    FLstk_stf@stock.n[,i]    <- t(exp(sims_yr[[1]]$sim[,1:dim(FLstk_stf)[1]]))
    FLstk_stf@harvest[,i]    <- t(exp(sims_yr[[1]]$sim[,(dim(FLstk_stf)[1]+1):ncol(sims_yr[[1]]$sim)])[,Fkeys])
    FLstk_stf@catch[,i]      <- sims_yr[[1]]$catch
    FLstk_stf@landings[,i]   <- sims_yr[[1]]$land
    FLstk_stf@discards[,i]   <- FLstk_stf@catch[,i] - FLstk_stf@landings[,i]
  }

  if(add_deterministic_SAM){
    input <- forecast_args
    input <- input[!unlist(lapply(input, is.null))]
    input$deterministic <- T
    input$nosim <- 2
    input$savesim <- T
    input$fit <- attr(sam_forecast,"fit")
    input$estimate <- as.name(attr(forecast_args$estimate,"estimateLabel"))
    if(is.null(attr(forecast_args,"RNG_state")))rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed",attr(forecast_args,"RNG_state"), envir = .GlobalEnv)
    frc <- do.call("forecast",input)
    FLstk_sam      <- propagate(FLstk_sam, iter = 2, fill.iter = T)
    for(i in proj.yrs){
      sims_yr            <- frc[proj.yrs == i]
      FLstk_sam@catch.n[,i]    <- sims_yr[[1]]$catchatage
      FLstk_sam@landings.n[,i] <- FLstk_sam@catch.n[,i] * yearMeans(FLstk@landings.n[,as.character(ave.years)]/FLstk@catch.n[,as.character(ave.years)])
      FLstk_sam@discards.n[,i] <- FLstk_sam@catch.n[,i] - FLstk_sam@landings.n[,i]

      FLstk_sam@stock.n[,i]    <- t(exp(sims_yr[[1]]$sim[,1:dim(FLstk_sam)[1]]))
      FLstk_sam@harvest[,i]    <- t(exp(sims_yr[[1]]$sim[,(dim(FLstk_sam)[1]+1):ncol(sims_yr[[1]]$sim)])[,Fkeys])
      FLstk_sam@catch[,i]      <- sims_yr[[1]]$catch
      FLstk_sam@landings[,i]   <- sims_yr[[1]]$land
      FLstk_sam@discards[,i]   <- FLstk_sam@catch[,i] - FLstk_sam@landings[,i]
    }
  }
  if(is.null(FLstk_fwd)){
    if(add_deterministic_SAM){
      return(FLStocks(list("SAM" = FLstk_stf, "SAM_deterministic" = FLstk_sam)))
    } else {
      return(FLStocks(list("SAM" = FLstk_stf)))
    }
  } else {
    if(add_deterministic_SAM){
      return(FLStocks(list("SAM" = FLstk_stf, "FLR" = FLstk_fwd, "SAM_deterministic" = FLstk_sam)))
    } else {
      return(FLStocks(list("SAM" = FLstk_stf, "FLR" = FLstk_fwd)))
    }
  }
}
