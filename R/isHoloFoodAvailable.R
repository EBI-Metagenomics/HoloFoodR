#' Check whether the HoloFood database is reachable
#'
#' @details
#' The HoloFood server is occasionally unresponsive. \code{isHoloFoodAvailable}
#' sends a lightweight request to the database and reports whether a valid
#' response was received. It never throws an error, so it can be used to guard
#' code that requires the database, e.g., in scripts, examples, tests or
#' vignettes.
#'
#' @param base.url \code{Character scalar} specifying the base URL of the
#' HoloFood API. (Default: \code{"https://www.holofooddata.org/api"})
#'
#' @param timeout \code{Numeric scalar} specifying the maximum time in seconds
#' to wait for a response. (Default: \code{30})
#'
#' @return \code{Logical scalar}: \code{TRUE} if the database responded
#' successfully, \code{FALSE} otherwise.
#'
#' @examples
#'
#' if( isHoloFoodAvailable() ){
#'     res <- doQuery("genome-catalogues")
#' }
#'
#' @export
#' @importFrom httr2 request req_url_query req_timeout req_error req_perform
#' @importFrom httr2 resp_status
isHoloFoodAvailable <- function(
        base.url = "https://www.holofooddata.org/api", timeout = 30){
    temp <- .check_input(base.url, list("character scalar"))
    temp <- .check_input(timeout, list("numeric scalar"))
    # Query a small endpoint. Any failure (no connection, timeout, HTTP error)
    # means the database is not available.
    available <- tryCatch({
        req <- request(paste0(base.url, "/genome-catalogues"))
        req <- req_url_query(req, format = "json")
        req <- req_timeout(req, timeout)
        req <- req_error(req, is_error = \(resp) FALSE)
        resp <- req_perform(req)
        resp_status(resp) == 200L
    }, error = function(e) FALSE)
    return(isTRUE(available))
}
