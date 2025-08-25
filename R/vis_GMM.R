#' Visualize Gaussian Mixture Model (GMM) clustering results
#'
#' @description
#' This function generates visualizations of Gaussian Mixture Model (GMM) clustering
#' applied to pathway-level data. For each pathway (row in `df_path`), three plots
#' are returned: (1) points colored by GMM cluster assignment,
#' (2) points colored by pathway activity score (PAS),
#' and (3) the fitted GMM distribution with threshold lines.
#'
#' @param PR A numeric matrix or data.frame with two columns (first and second dimensions
#' of the reduced space, e.g. PCA, UMAP or t-SNE coordinates).
#' @param df_path A numeric matrix or data.frame of pathway activity scores,
#'   where rows correspond to pathways and columns to samples. Row names should correspond to pathway identifiers. Its is output of function gene2path().
#' @param res_gmm A list of GMM results, one per pathway,
#' each containing at least:
#'   \itemize{
#'     \item \code{cluster} – vector of cluster assignments per sample,
#'     \item \code{KS} – number of clusters,
#'     \item \code{fig} – a ggplot object representing the fitted GMM distribution.
#'   }
#' Output of GMMdecomp() function is valid.
#' @param thr_method Character string, thresholding method for cluster significance.
#' Must be one of:
#'   \itemize{
#'     \item `"KM"` – use GMM with k-means threshold,
#'     \item `"Top1"` – use GMM Top-1 threshold.
#'   }
#' Default is `"KM"`.
#'
#'
#' @return
#' A list of length equal to the number of pathways (`nrow(df_path)`).
#' Each element is itself a list with three ggplot objects:
#'   \itemize{
#'     \item \code{clusters} – scatterplot colored by cluster assignment.
#'     \item \code{PAS} – scatterplot colored by PAS values.
#'     \item \code{gmm} – fitted GMM distribution with threshold line.
#'   }
#'
#' @seealso [thr_GMM()] for computing thresholds.
#'
#' @examples
#' \dontrun{
#' # Example (assuming PR, df_path, and res_gmm already exist)
#' plots <- vis_GMM(PR, df_path, res_gmm, thr_method = "KM")
#'
#' # Access plots for first pathway
#' plots[[1]]$clusters
#' plots[[1]]$PAS
#' plots[[1]]$gmm
#' }
#'
#' @import cli
#' @import ggplot2
#'
#' @export
vis_GMM<-function(PR,df_path,res_gmm,thr_method="KM"){
  # ---- Arg chceck ----
  if (missing(PR) || missing(df_path) || missing(res_gmm)) {
    cli_abort("All arguments {.arg PR}, {.arg df_path}, {.arg res_gmm} must be provided.")
  }
  if (!is.data.frame(PR) && !is.matrix(PR)) {
    cli_abort("{.arg PR} must be a data.frame or matrix.")
  }
  if (!is.matrix(df_path) && !is.data.frame(df_path)) {
    cli_abort("{.arg df_path} must be a data.frame or matrix.")
  }
  if (!is.list(res_gmm)) {
    cli_abort("{.arg res_gmm} must be a list.")
  }

  # Dimension check
  if (nrow(df_path) != length(res_gmm)) {
    cli_abort("Number of rows in {.arg df_path} must equal the number of elements in {.arg res_gmm}.")
  }

  valid_methods <- c("KM", "Top1")
  if (!thr_method %in% valid_methods) {
    cli_abort("Invalid {.arg thr_method}. Must be one of: {toString(valid_methods)}.")
  }



  # ---- Main function body ----
  PR <- as.data.frame(PR)
  colnames(PR) <- c("V1", "V2")

  thrs <- thr_GMM(res_gmm)
  thr_col <- if (thr_method == "KM") "Kmeans_thr" else "Top1_thr"
  all_plots <- vector("list", nrow(df_path))

  cli_inform("Generating {nrow(df_path)} GMM visualizations...")
  for (i in seq_len(nrow(df_path))) {
    df<-PR
    tit <- rownames(df_path)[i]

    # Assign clusters
    df$col_factor <- factor(res_gmm[[i]]$cluster,levels = 1:res_gmm[[i]]$KS,
      labels = paste0("Cluster ", 1:res_gmm[[i]]$KS))

    # Path scores
    df$PAS <- as.numeric(df_path[i, ])

    # Significant vs non-significant clusters by thresholding
    clust_sig <- which(thrs[[i]]$All_thr == thrs[[i]][[thr_col]]) + 1
    col_sig <- colorRampPalette(c("#FFEB6B", "red4"))(max(res_gmm[[i]]$cluster) - clust_sig + 1)
    col_unsig <- colorRampPalette(rev(c("black", "grey75")))(clust_sig - 1)

    p1 <- ggplot(df, aes(x = V1, y = V2, color = col_factor)) +
      geom_point(size = 1) + theme_bw(base_size = 12) +
      labs(title = tit, color = "", x = "Dim. 1", y = "Dim. 2") +
      scale_color_manual(values = c(rev(col_unsig), col_sig)) +
      theme(legend.position = "bottom")

    p2 <- ggplot(df, aes(x = V1, y = V2, color = PAS)) +
      geom_point(size = 1) + theme_bw(base_size = 10) +
      labs(title = tit, color = "", x = "Dim. 1", y = "Dim. 2") +
      scale_color_viridis_c() +
      theme(legend.position = "bottom") +
      labs(color = "")

    z <- res_gmm[[i]]$fig
    n_layers <- length(z$layers)

    # adjust linewidth
    if (n_layers > 1) {z$layers[[n_layers - 1]]$aes_params$linewidth <- 1.2}

    z <- z +
      scale_color_manual(values = c(rev(col_unsig), col_sig, "grey90")) +
      labs(title = tit, color = "", x = "PAS value") +theme(legend.position = "none")

    # adjust vertical lines
    for (j in seq_along(z$layers)) {
      if ("GeomVline" %in% class(z$layers[[j]]$geom)) {
        z$layers[[j]]$aes_params$colour <- "grey65"}}

    z <- z + geom_vline(xintercept = thrs[[i]][[thr_col]],color = "deepskyblue3")
    all_plots[[i]] <- list(clusters = p1, PAS = p2, gmm = z)

    cli_progress_message("Processed {i}/{nrow(df_path)} pathways: {tit}")
  }

  cli_alert_success("Finished generating {nrow(df_path)} GMM visualizations.")
  names(all_plots)<-rownames(df_path)
  return(all_plots)

}
