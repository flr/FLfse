NEWS
================

## Version 1.0.2

Changes to `SAM2FLStock`

-   Fixed bug when observations in the SAM model fit do not include
    years but instead an index for the years (1, 2, 3, …).

## Version 1.0.1

Changes to `SAM2FLStock`

-   This function can now handle multi-fleet SAM model fits. An example
    is given in the documentation for this method (see `?SAM2FLStock`).

-   Added `stock_only` argument. If set to `TRUE`, catch data (numbers,
    weights) are ignored and only stock data (numbers, SSB, etc.) are
    returned.

## Version 1.0

### Changes since version 0.x

Changes to `SAM2FLStock`

-   Fix when SAM model fit included age 0 but without catches for this
    age.

-   Fix setting of plusgroup.

-   New arguments `mat_est`, `stock.wt_est`, `catch.wt_est`, `m_est`
    (logical `TRUE`/`FALSE`) for returning estimates of biological data
    from SAM, if they are estimated by SAM.

-   New argument `spinoutyear` (logical `TRUE`/`FALSE`) to
    include/exclude estimates of biological data from SAM beyond last
    data year.

-   Generic methods for `SAM2FLStock` updated when argument `stk` is
    provided as a template. In the returned output, `stk` is used as
    template and all slots are updated, instead of only those
    corresponding to catch, fishing mortality and stock size.

Changes to `SAM_uncertainty`

-   Name of output list element with catch numbers including uncertainty
    changed from `catch_n` to `catch.n` so that it follows the same
    convention as `stock.n`.

General changes

-   Updated help files.
