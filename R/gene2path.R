#' Function to transform gene level information to pathway level by single-sample approaches.
#'
#' Function reduce genes to pathways.
#'
#' @param X matrix of data.
#' @param Y list of pathways
#' @param method methods for single-sample pathway enrichment analysis.
#'
#' @returns Function returns...
#'
#' @importFrom
#' @examples
#' \dontrun{
#'
#'
#' }
#'
#' @seealso \code{\link{gaussian_mixture_vector}}, \code{\link{EM_iter}}
#'
#' @export
gene2path<- function(X, Y, method = "ssGSEA"){
  # Check part
  if (!hasArg("X")){
    stop("No data.")}


  x<-1


}
