#' Function to transform gene level information to pathway level by single-sample approaches
#'
#' Function reduce genes to pathways.
#'
#' @param X matrix of data.
#' @param pathway list of pathways to analsysis
#' @param method methods for single-sample pathway enrichment analysis.
#'
#' @returns Function returns...
#'
#' @examples
#' \dontrun{
#'
#'
#' }
#'
#' @seealso \code{\link{gaussian_mixture_vector}}, \code{\link{EM_iter}}
#'
#' @export
gene2path<- function(X, pathway, method = "ssGSEA"){
  # Check part
  if (!hasArg("X")){
    stop("No data.")}
#
#   if (method=="ssGSEA"){
#     dfssGSEA <- gsva(as.matrix(data), pathway, method = "ssgsea", kcdf = "Gaussian", mx.diff = F)
#     dfssGSEA <- as.data.frame(dfssGSEA)

    switch(method,
      "ssGSEA"= df_path<-ssGSEA(X,pathway),
      "Mean" =  df_path<-Mean(X,pathway)
    )


return(df_path)
}
