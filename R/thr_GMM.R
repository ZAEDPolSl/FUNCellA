#' Threshold selection for pathway activity using GMM decomposition
#'
#' This function selects thresholds indicating samples with relatively high pathway activity
#' based on Gaussian Mixture Model (GMM) decomposition. It processes GMM output from the
#' \code{dpGMM} package or the \code{GMMdecomp} function, extracting thresholds using different
#' strategies, including maximum threshold (Top1) and optional k-means-based refinement.
#'
#' @param gmms A list of GMM decomposition results, where each top-level list element corresponds
#'   to one pathway. Each element should contain:
#'   \describe{
#'     \item{\code{threshold}}{A numeric vector of threshold values derived from the GMM.}
#'     \item{\code{model}}{A matrix of GMM parameters (e.g., means, variances and alphas for each component).}
#'   }
#'
#' @return A named list of thresholding results per pathway. For each pathway, a sublist is returned with:
#'   \describe{
#'     \item{\code{Kmeans_thr}}{Threshold selected based on k-means clustering of GMM parameters.}
#'     \item{\code{Top1_thr}}{Top threshold (maximum of the available thresholds).}
#'     \item{\code{All_thr}}{All threshold values derived from the GMM.}
#'   }
#'
#' @import cli
#'
#' @export
thr_GMM<-function(gmms){

  idprog <- cli_progress_bar("Calculating thresholds", total=length(gmms))
  ret <- vector("list", length(gmms))

  for (i in 1:length(gmms)) {
    cli_progress_update(id = idprog)

    thr<-gmms[[i]]$threshold
    params<-gmms[[i]]$model

    # Top 1 threshold
    thrT<-max(thr,na.rm=T)

    # K-means threshold
    if (length(thr)>0){
      # Remove NA components
      idx <- !is.na(thr)
      thr <- thr[idx]

      if (sum(idx ) < nrow(params)) {
        params <- params[idx, , drop = FALSE]
      }

      if (nrow(params)==2){
        thrK<-thr
      } else if ((length(thr)+1)!=nrow(params)){
        thrK<-max(thr,na.rm=T)
      }else{
        thrK<-km_search(params,thr)
      }
    } else{
      thrK<--Inf
    }


    ret[[i]]<-list(Kmeans_thr=thrK,Top1_thr=thrT,All_thr=thr)

  }
  cli_progress_done(idprog)
  cli_alert_success('Thresholds calculated')

  names(ret)<-names(gmms)
  return(ret)
}
