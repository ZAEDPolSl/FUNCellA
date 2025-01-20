#' Exemplary data of pathways.
#'
#' File containing 45 pathways form CellMarker database.
#'
#' @docType data
#' @name pathways
#' @usage data(pathways)
#' @format A list with 45 elements which represent pathways. Each element contain vector of gene names (characters) of each pathway.
#' Names represent CellMarker pathway ID for which information are given in pathways_info file.
#' @source \url{https://cran.r-project.org/web/packages/tmod/index.html}
"pathways"

#' Pathay description for exemplary data.
#'
#' File containing 45 pathways form CellMarker database.
#'
#' @docType data
#' @name pathways_info
#' @usage data(pathways_info)
#' @format A data frame with {45} rows and {3} variables:
#' \describe{
#'   \item{ID}{ID pathways ID which math to names of pathway data.}
#'   \item{Title}{Title pathways title}
#'   \item{DataBase}{DataBase information about database for obtaining the pathways.}
#' }
#' @source \url{https://cran.r-project.org/web/packages/tmod/index.html}
"pathways_info"
