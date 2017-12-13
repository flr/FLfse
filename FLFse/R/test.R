#
# # fit_fwd <- forecast(fit, fscale = c(1))
# # fit_fwd2 <- forecast(fit, fscale = c(2))
#
#
# library(stockassessment)
# library(FLCore)
# library(ggplotFL)
# library(spict)
#
# path <- "C:/Users/SF02/OneDrive - CEFAS/MA016N/"
# source(paste0(path, "scripts/functions.R"))
#
#
# ### ------------------------------------------------------------------------ ###
# ### NS cod from stockassessment.org
# ### ------------------------------------------------------------------------ ###
# ### WGNSSK 2017 assessment
#
# ### load stock, index & SAM configuration
# idx <- readRDS(paste0(path, "SAM/example_stocks/cod.27.47d_idx.rds"))
# stk <- readRDS(paste0(path, "SAM/example_stocks/cod.27.47d_stk.rds"))
# conf <- readRDS(paste0(path, "SAM/example_stocks/cod.27.47d_conf.rds"))
#
# ### fit SAM
# fit <- FLR_SAM(stk = stk, idx = idx, conf = conf)
#
# ### check convergence
# fit$fit$opt$convergence
#
# ### convert into FLStock
# stk_fit <- FLStock(fit$fit, uncertainty = TRUE)
# plot(stk_fit)
#
# ### ------------------------------------------------------------------------ ###
# ### Irish Sea plaice ple.27.7a from WGCSE 2017
# ### ------------------------------------------------------------------------ ###
# ### WGCSE 2017 assessment
#
# ### load stock, index & SAM configuration
# idx <- readRDS(paste0(path, "SAM/example_stocks/ple.27.7a_idx.rds"))
# stk <- readRDS(paste0(path, "SAM/example_stocks/ple.27.7a_stk.rds"))
# conf <- readRDS(paste0(path, "SAM/example_stocks/ple.27.7a_conf.rds"))
#
# ### fit SAM
# fit <- FLR_SAM(stk = stk, idx = idx, conf = conf)
#
# ### convert into FLStock
# stk_fit <- FLStock(fit$fit, uncertainty = TRUE)
# plot(stk_fit)
#
#
# ### ------------------------------------------------------------------------ ###
# ### SPiCT
# ### ------------------------------------------------------------------------ ###
#
#
# fit_spict <- FLR_SPiCT(stk = stk, idx = idx)
#
#
#
#fit_spict <- FLR_SPiCT(stk = cod4_stk, idx = cod4_idx)
