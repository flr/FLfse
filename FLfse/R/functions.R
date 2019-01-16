### ------------------------------------------------------------------------ ###
### additional functions (internal) ####
### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### adapt dimensions in 2 FLR objects ####
### ------------------------------------------------------------------------ ###
### returned as list
adapt_dims <- function(obj1, obj2, fill.iter = FALSE) {

  ### save objects in list
  res <- list(obj1, obj2)

  ### do nothing if dimensions are identical
  if (identical(dims(obj1), dims(obj2))) return(res)

  ### min year
  min_year <- min(as.numeric(sapply(lapply(res, dims), "[", "minyear")))
  ### max year
  max_year <- max(as.numeric(sapply(lapply(res, dims), "[", "maxyear")))

  ### adapt year range
  res <- lapply(res, window, start = min_year)
  res <- lapply(res, window, end = max_year)

  ### number of iterations
  n_iter <- as.numeric(sapply(lapply(res, dims), "[", "iter"))

  ### check comparability of iterations and stop if not
  if (!identical(n_iter[1], n_iter[2])) {

    if (!any(n_iter == 1)) stop("incompatible iter dimension")

    ### adapt
    res[[which(n_iter == 1)]] <- propagate(res[[which(n_iter == 1)]],
                                           n_iter[which(n_iter != 1)],
                                    fill.iter = fill.iter)

  }

  return(res)

}
