#' Function to transform genes to pathway level by single-sample approaches
#'
#' Function reduce genes to pathways.
#'
#' @param X matrix of data
#' @param pathway list of pathways to analysis.
#' @param method methods for single-sample pathway enrichment analysis.
#' @param filt_cov numeric vector of length 1. Minimum % of pathway coverage after gene mapping. By default is 0 (no filtration), the 0.65 will indicate that only pathway with coverage of genes larger than 65% will be taken.
#' @param filt_min numeric vector of length 1. Minimum size of the resulting pathway after gene identifier mapping. By default, the minimum size is 15.
#' @param filt_max numeric vector of length 1. Maximum size of the resulting pathway after gene identifier mapping. By default, the minimum size is 500.
#'
#' @return Function returns...
#'
#' @examples
#' \dontrun{
#'
#'
#' }
#'
#' @import cli
#'
#' @seealso \code{\link{gaussian_mixture_vector}}, \code{\link{EM_iter}}
#'
#' @export
gene2path<- function(X, pathway, method = "ssGSEA",filt_cov=0,filt_min=15,filter_max=500){
  # parameters check part ----
  if (!hasArg("X")){
    stop("No data.")}

  if (!hasArg("pathway")){
    stop("No pathway list.")}

  if (filt_cov<0 | filt_cov>1){
    stop("filt_cov has to be value from 0 to 1")}

  method <- match.arg(method,choices=c('ssGSEA','Mean','BINA',"CERNO","ZSc",'JASMINE'))

  # check of filtration by coverage ----
  if (filt_cov!=0){
    pathway<-cover_check(X,pathway,filt_cov)
  } else {cli_alert_success("No filtration due to coverage")}

  # check of filtration by pathways size ----
  pathway<-min_max_filter(pathway,filt_min,filter_max)

  # run of gene to path transformation ----
  switch(method,
      "ssGSEA"= df_path<-ssGSEA(X,pathway),
      "Mean" =  df_path<-Mean_path(X,pathway)
    )


return(df_path)
}
