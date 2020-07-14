#' North Sea cod
#'
#' Cod (\emph{Gadus morhua}) in Subarea 4, Division 7.d, and Subdivision 20
#' (North Sea, eastern English Channel, Skagerrak).
#' Assessment input data for North Sea cod as used by ICES WGNSSK 2018 and 2019.
#' Includes the stock, the indices and the SAM configuration. See the example(s)
#' below for how to run the assessment.
#'
#' @examples
#' # Replicate the 2018 assessment:
#' fit <- FLR_SAM(stk = cod4_stk,
#'                idx = cod4_idx,
#'                conf = cod4_conf_sam)
#'
#' # Replicate the 2019 assessment:
#' fit <- FLR_SAM(stk = cod4_stk_2019,
#'                idx = cod4_idx_2019,
#'                conf = cod4_conf_sam)
#'
#' @source \url{https://www.stockassessment.org}
"cod4_stk"

#' @rdname cod4_stk
"cod4_stk_2019"

#' @rdname cod4_stk
"cod4_stk_2019"

#' @rdname cod4_stk
"cod4_idx"

#' @rdname cod4_stk
"cod4_idx_2019"

#' @rdname cod4_stk
"cod4_conf_sam"

#' Irish Sea plaice
#'
#' Plaice (\emph{Pleuronectces platessa}) in Division 7.a (Irish Sea).
#' Assessment input data for Irish Sea plaice as used by ICES WGCSE 2017.
#' Includes the stock, the indices and the SAM configuration. See the example
#' below for how to run the assessment.
#'
#' @examples
#' # Replicate the 2017 assessment:
#' \dontrun{
#' fit <- FLR_SAM(stk = ple7a_stk,
#'                idx = ple7a_idx,
#'                conf = ple7a_conf_sam)
#' }
#'
"ple7a_stk"

#' @rdname ple7a_stk
"ple7a_idx"

#' @rdname ple7a_stk
"ple7a_conf_sam"


#' North Sea whiting
#'
#' Whiting (\emph{Merlangius merlangus}) in Subarea 4 and Division 7.d
#' (North Sea and eastern English Channel).
#' Assessment input data for North Sea whiting as used by ICES WGNSSK 2018.
#' Includes the stock, the indices and the SAM configuration. See the example
#' below for how to run the assessment.
#'
#' @source \url{https://www.stockassessment.org}
#'
#' #' @examples
#' # Replicate the 2017 assessment:
#' fit <- FLR_SAM(stk = whg4_stk,
#'                idx = whg4_idx,
#'                conf = whg4_conf_sam)
#'
"whg4_stk"

#' @rdname whg4_stk
"whg4_idx"

#' @rdname whg4_stk
"whg4_conf_sam"


#' North Sea haddock
#'
#' Haddock (\emph{Melanogrammus aeglefinus}) in Subarea 4, Division 6.a, and
#' Subdivision 20 (North Sea, West of Scotland, Skagerrak).
#' Assessment input data for North Sea haddock as used by ICES WGNSSK 2018.
#' Includes the stock, the indices and the SAM configuration. See the example
#' below for how to run the assessment.
#'
#' @source \url{https://www.stockassessment.org}
#'
#' @examples
#' # Replicate the 2018 assessment:
#' fit <- FLR_SAM(stk = had4_stk,
#'                idx = had4_idx,
#'                conf = had4_stk)
#'
"had4_stk"

#' @rdname had4_stk
"had4_idx"

#' @rdname had4_stk
"had4_conf_sam"

#' North Sea herring
#'
#' Herring (Clupea harengus) in Subarea 4 and divisions 3.a and 7.d,
#' autumn spawners (North Sea, Skagerrak and Kattegat, eastern English Channel).
#' Assessment input data for North Sea herring as used by ICES HAWG 2019.
#' Includes the stock, the indices and the SAM configuration.
#'
#' Please note that the herring SAM assessment requires a different version of
#' the stockassessment package. For details, see the readme at
#' https://github.com/shfischer/FLfse
#'
#' @source \url{https://github.com/ices-eg/wg_HAWG/tree/master/NSAS}
#'
#' #' @examples
#' \dontrun{
#' # NOTE: the herring SAM assessment requires a different version of
#' # the stockassessment package
#' fit <- FLR_SAM(stk = her4_stk,
#'                idx = her4_idx,
#'                conf = her4_conf_sam)
#' }
#'
"her4_stk"

#' @rdname her4_stk
"her4_idx"

#' @rdname her4_stk
"her4_conf_sam"


#' Northeast Atlantic mackerel 2019
#'
#' Mackerel  (\emph{Scomber  scombrus}) in subareas 1-8 and 14, and in
#' Division 9.a (the Northeast Atlantic and adjacent waters).
#' Assessment input data for Northeast Atlantic mackerel as used by ICES WGWIDE
#' 2019. Includes the stock, the indices and the SAM configuration. See the
#' example below for how to run the assessment.
#'
#' \code{mac.27.nea_stk_2019} includes tagging data, stored as an attribute in
#' \code{attr(catch.n(mac.27.nea_stk_2019), "recap")}. An additional tagging
#' configuration is stored in
#' \code{attr(catch.n(mac.27.nea_stk_2019), "recap_conf")}.
#' Also, weights for the catch numbers are supplied as an attribute in
#' \code{attr(catch.n(mac.27.nea_stk_2019), "weight")}. All additional data is
#' passed on automatically to the SAM assessment.
#'
#' @examples
#' # The 2019 WGWIDE Northeast Atlantic mackerel assessment can be reproduced
#' # with:
#' fit <- FLR_SAM(stk = mac.27.nea_stk_2019,
#'                idx = mac.27.nea_idx_2019,
#'                conf = mac.27.nea_conf_2019)
#'
#' @source \url{https://www.stockassessment.org}
"mac.27.nea_stk_2019"

#' @rdname mac.27.nea_stk_2019
"mac.27.nea_idx_2019"

#' @rdname mac.27.nea_stk_2019
"mac.27.nea_conf_2019"
