#' Function for GMM decomposition of pathways
#'
#' Function reduce genes to pathways.
#'
#' @param X data.frame of data enrichment scores (rows: pathways, columns: samples).
#' @param K Maximum number of components for GMM decomposition (default: 10).
#' @param multiply Logical, whether to scale values by 10 before fitting GMM (default: TRUE).
#' @param IC Information criterion used to select the number of model components.
#' Possible methods are "AIC","AICc", "BIC" (default), "ICL-BIC" or "LR".
#' @param parallel Logical, whether to use parallel computing (default: FALSE). (Not implemented yet)
#'
#' @return Function returns...
#'
#' @import dpGMM
#' @import cli
#'
#' @export
GMMdecomp <- function(X, K=10, multiply = TRUE, IC="BIC",parallel = FALSE) {
  # ---- Parameter validation ----
  if (!is.data.frame(X)) {
    cli_abort(c("x" = "Input X must be a data.frame"))
  }

  if (!is.numeric(K) || K <= 0 || K %% 1 != 0) {
    cli_abort(c("x" = "K must be a positive integer."))
  }

  if (!is.logical(multiply)) {
    cli_abort(c("x" = "multiply must be a logical value (TRUE or FALSE)."))
  }

  IC_list <- c("AIC","AICc", "BIC", "ICL-BIC", "LR")
  if (!opts$IC %in% IC_list) {
    stop("Criterion not implemented. Please use AIC, AICc, BIC, ICL-BIC or LR")
  }

  if (!is.logical(parallel)) {
    cli_abort(c("x" = "parallel must be a logical value (TRUE or FALSE)."))
  }

  # ---- GMM options setup ----
  opt <- dpGMM::GMM_1D_opts
  opt$max_iter <- 1000
  opt$KS <- K
  opt$plot <- FALSE
  opt$quick_stop <- FALSE
  opt$SW <- 0.05
  opt$sigmas.dev <- 0
  opt$IC<-IC

  cli_alert_info("Start calculating GMM decompositions")

  # ---- Helper for row calculation ----
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

  # ---- GMM Calculation ----
  idprog <- cli_progress_bar("Calculating GMM decompositions", total = nrow(X))

  if (parallel) {
    if (.Platform$OS.type == "windows") {

    } else {
      future::plan(future::multicore)
    }
  } else{
    results_list <- lapply(1:nrow(X), function(i) {
      cli_progress_update(id = idprog)
      row_multiple(X[i, ])
    })
  }

  cli_progress_done(idprog)
  cli_alert_success('GMM calculated')

  names(results_list) <- rownames(X)
  return(results_list)
}
