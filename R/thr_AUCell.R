#' Thresholding pathway activity scores using the AUCell approach
#'
#' Applies the AUCell package's methodology to identify thresholds that indicate
#' samples with relatively high pathway activity. This function determines, for each
#' pathway, which samples exhibit significant enrichment based on expression rankings.
#'
#' @param df_path A numeric matrix or data.frame of pathway activity scores,
#'   where rows correspond to pathways and columns to samples.
#'
#' @return A vector of threshold values.
#'   For detailed behavior, refer to the AUCell package documentation.
#'
#' @import cli AUCell
#'
#' @source \url{https://www.bioconductor.org/packages/release/bioc/html/AUCell.html}
#' @references Aibar, S., Bravo González-Blas, C., Moerman, T., Huynh-Thu, V.A.,
#'   Imrichová, H., Hulselmans, G., Rambow, F., Marine, J.C., Geurts, P., Aerts, J.,
#'   van den Oord, J., Kalender Atak, Z., Wouters, J., & Aerts, S. (2017).
#'   SCENIC: Single-Cell Regulatory Network Inference And Clustering. *Nature Methods*, 14, 1083–1086.
#'   \doi{10.1038/nmeth.4463}
#'
#' @import AUCell
#' @export
#'
thr_AUCell<-function(df_path){

  idprog <- cli_progress_bar("Calculating thresholds", total=nrow(df_path))
  ret <- c()
  for (i in 1:nrow(df_path)) {
    cli_progress_update(id = idprog)
    ret[i]<-AUCell:::.auc_assignmnetThreshold_v6(as.matrix(df_path[i,]),plotHist = F)$selected
  }

  cli_progress_done(idprog)
  cli_alert_success('Thresholds calculated')

  names(ret)<-rownames(df_path)
  return(ret)
}
