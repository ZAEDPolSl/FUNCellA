#' Function for GMM decomposition of pathways
#'
#' Function reduce genes to pathways.
#'
#' @param X matrix of data enrichment scores rows pathways columns samples.
#' @param K maximum number of components for GMM decomposition (default: 10).
#' @param multiply logical, scaling values by 10 before fitting GMM (default: TRUE).
#' @param parallel logical, lunching parallel computing (default: FALSE).
#'
#' @return Function returns...
#'
#' @import dpGMM
#' @import cli
#'
#' @export
GMMdecomp <- function(X, K=10, multiply = TRUE, parallel = FALSE) {

  opt <- dpGMM::GMM_1D_opts
  opt$max_iter <- 1000
  opt$KS <- K
  opt$plot <- FALSE
  opt$quick_stop <- FALSE
  opt$SW <- 0.05
  opt$sigmas.dev <- 0

  cli_alert_info("Start calculating GMM decompositions")

  row_multiple <- function(row) {
    tmp <- as.numeric(row)
    if (multiply) {
      tmp <- tmp * 10
    }

    result <- runGMM(tmp, opts = opt)

    if (multiply) {
      result$model$mu <- result$model$mu / 10
      result$model$sigma <- result$model$sigma / 10
      result$threshold <- result$threshold / 10
      tmp <- tmp / 10
    }
    dist.plot <- generate_dist(tmp, result$model$alpha, result$model$mu, result$model$sigma, 1e4)
    result$fig <- plot_gmm_1D(tmp, dist.plot, Y = NULL, threshold = result$threshold, pal = "Blues")
    return(result)
  }

  idprog <- cli_progress_bar("Calculating GMM decompositions", total = nrow(X))
  results_list <- lapply(1:nrow(X), function(i) {
    cli_progress_update(id = idprog)
    row_multiple(X[i, ])
  })
  cli_progress_done(idprog)
  cli_alert_success('GMM calculated')

  names(results_list) <- rownames(X)
  return(results_list)
}
