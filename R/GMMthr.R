#' Function for selecting thresholds ofor pathway activity
#'
#' Function reduce genes to pathways.
#'
#' @param gmms output from GMMdecomp function.
#'
#' @return Function returns...
#'
#' @import cli
#'
#' @export
GMMthr<-function(gmms){

  idprog <- cli_progress_bar("Calculating thresholds", total=nrow(X))
  for(i in 1:length(gmms)){

  }
  cli_progress_done(idprog)
  cli_alert_success('Thresholds calculated')


  return(ret)
}
