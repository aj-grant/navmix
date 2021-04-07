#'Noise-augmented directional clustering
#'
#'Performs directional clustering by fitting a noise-augmented von Mises-Fisher mixture model
#'
#'@param x Matrix of values where rows represent observations and columns represent features.
#'@param K The number of clusters to fit.
#'@param select_K If TRUE (the default setting), the number of clusters will be chosen by BIC, with K the maximum number of clusters
#'considered. If FALSE, then a model with K clusters will be fit.
#'@param common_kapp If TRUE, then model will force the kappa parameter to be equal for all clusters, except the noise
#'cluster.
#'@param pj_ini The initial proportion of observations which belong in the noise cluster. Must be a number greater or
#'equal to 0 and strictly less than 1. The default value is 0.05. If set to 0, no observations will be placed in the
#'noise cluster.
#'@param no_ini The number of time the algorithm is run with different initialisations. Must be a number greater than
#'zero. The default value is 5.
#'@param tol The tolerance threshold for convergence of the EM algorithm. Must be a number greater than 0. The default
#'value is 1.0e-4.
#'@param max_iter The maximum number of iterations of the EM algorithm. Must be a number greater than 0. The default
#'value is 100.
#'@param plot Plots of the results will be produced if set to TRUE. Default is FALSE.
#'@param plot_heat Produces a heatmap of the results if plot is set to TRUE. The heatmap will also be returned as a
#'ggplot object.
#'@param plot_radial Produces (a) radial plot(s) of the results if plot is set to TRUE.
#'@param plot_radial_separate If set to FALSE (the default value), the fitted means of each cluster are plotted on the
#'same radial plot. If set to TRUE, they are plotted on separate radial plots.
#'@param radial_legend_pos Adjusts the position of the legend for a radial plots with all fitted means plotted together.
#'@param radial_separate_col Adjusts the format of the output of radial plots on separate plots.
#'@return Returned are the BIC values for each model fitted ($BIC), the final fitted model ($fit) and, if produced, the
#'heatmap as a ggplot object ($heatmap_plot). The fitted model has the following.
#'\item{mu}{A matrix where each column represents the mean of the fitted von Mises-Fisher distribution for each cluster.}
#'\item{kappa}{A row vector where each element represents the kappa parameter of the fitted von Mises-Fisher distribution
#'for each cluster.}
#'\item{g}{A matrix of probabilities for each observation belonging to each cluster. The value in the jth row and kth
#'column represents the probability that the jth observation belongs to the kth cluster.}
#'\item{z}{A vector of the cluster membership of each observation when allocated according to the cluster for which it
#'has the highest probability of membership (hard clustering).}
#'\item{bic}{The BIC for the fitted model.}
#'\item{l}{The value of the likelihood function at the estimated parameters.}

navmix = function(x, K = 10, select_K = TRUE, common_kappa = FALSE, pj_ini = 0.05, no_ini = 5, tol = 1.0e-4,
                  max_iter = 100, plot = FALSE, plot_heat = TRUE, plot_radial = TRUE, plot_radial_separate = FALSE,
                  radial_legend_pos = c(-2.5, 2.7), radial_separate_col = 2){
  if (is.numeric(pj_ini) == FALSE){stop('pj_ini must be a number greater or equal to 0 and strictly less than 1')}
    else if (pj_ini < 0 | pj_ini > 1) {stop('pj_ini must be a number greater or equal to 0 and strictly less than 1')}
  if (is.numeric(no_ini) == FALSE){stop('no_ini must be a number greater or equal to 1')}
    else if (no_ini < 1){stop('no_ini must be a number greater or equal to 1')}
  if (is.numeric(tol) == FALSE){stop('tol must be a number greater than 0')}
    else if (tol <= 0){stop('tol must be a number greater than 0')}
  if (is.numeric(max_iter) == FALSE){stop('max_iter must be a number greater than 0')}
    else if (max_iter <= 0){stop('max_iter must be a number greater than 0')}
  if (is.null(dim(x))){
    n = length(x)
    x = as.matrix(x, nrow = n)
    m = dim(x)[2]
  } else {
    n = dim(x)[1]
    m = dim(x)[2]
  }
  if (K > n){
    if (select_K == TRUE){
      K = n
      warning('Cannot fit more clusters than observations. Maximum K has been set to the number of observations.')
    } else {
      stop('Cannot fit more clusters than observations.')
    }
  }
  if (is.null(rownames(x))){
    snps_names = sprintf("snp %0d", seq(1:n))
  } else{
    snps_names = rownames(x)
  }
  if (is.null(colnames(x))){
    trait_names = sprintf("trait %0d", seq(1:m))
  } else{
    trait_names = colnames(x)
  }
  all_fit = vector(mode = "list", length = K)
  bic = vector(length = K)
  if (select_K == TRUE){
    for (j in 1:K){
      all_fit[[j]] = navmix_K(x, j, pj_ini = pj_ini, common_kappa = FALSE, no_ini = no_ini, tol = tol, max_iter = max_iter)
      bic[j] = all_fit[[j]]$bic
      if (j > 2){
        if (bic[j] < bic[(j-1)] & bic[(j-1)] < bic[(j-2)]){
          bic = bic[1:j]
          break
        }
      }
    }
  select_fit = which.max(bic)
  fit = all_fit[[select_fit]]
  } else {
    fit = navmix_K(x, K, pj_ini = pj_ini, common_kappa = FALSE, no_ini = no_ini, tol = tol, max_iter = max_iter)
    bic = fit$bic
  }
  rownames(fit$mu) = trait_names
  rownames(fit$g) = snps_names
  navmix_out = list(BIC = bic, fit = fit)
  if (plot == TRUE){
    x_prop = row_norm(x)
    rownames(x_prop) = rownames(x)
    colnames(x_prop) = colnames(x)
    noise_cl = ncol(fit$g)
    if(plot_heat == TRUE){
      heatmap_plot = hm_plot(x_prop[fit$z != noise_cl, ], fit$z[fit$z != noise_cl])
      navmix_out$heatmap_plot = heatmap_plot
      }
    if(plot_radial == TRUE){
      radial_plot = rad_plot(fit$mu[, -noise_cl], plot_radial_separate = plot_radial_separate, radial_legend_pos,
                             radial_separate_col)
      }
  }
  return(navmix_out)
}

hm_plot = function(B, z, reorder_traits = TRUE){
  heat_df = data.frame(B)
  names(heat_df) = colnames(B)
  if (reorder_traits == TRUE){
    d = dist(t(heat_df))
    hc = hclust(d)
    heat_colord = hc$order
  } else {
    heat_colord = seq(1, ncol(B))
    }
  heat_df = heat_df[, heat_colord]
  if (is.null(rownames(B))){heat_df$Variant = seq(1, length(z))} else {heat_df$variant = rownames(B)}
  heat_df$clust = z
  heat_df = pivot_longer(heat_df, 1:ncol(B))
  names(heat_df) = c("Variant", "clust", "Trait", "Value")
  heat_df$Trait = factor(heat_df$Trait, levels = colnames(B)[heat_colord])
  heat_df$Variant = factor(heat_df$Variant)
  heatmap_plot = ggplot(heat_df, aes(x = Variant, y = Trait, fill = Value)) +
    geom_tile() + facet_grid(cols = vars(clust), scales = "free_x", space = 'free') +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 7),
          axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7),
          legend.text = element_text(size = 7), legend.title = element_text(size = 7),
          strip.text = element_text(size = 7)) +
    scale_fill_distiller(palette="RdYlBu")
  print(heatmap_plot)
}

rad_plot = function(mu, plot_radial_separate = FALSE, plot_radial_par = NULL, radial_legend_pos = c(-2.5, 2.7),
                    radial_separate_col = 2){
  if (is.null(plot_radial_par)){
    current_par = par(cex.axis = 0.5, cex.lab = 0.5, cex.main = 0.75, oma = c(0, 0, 0, 0), mfrow = c(1, 1))
  } else {
      current_par = plot_radial_par
      }
  K = ncol(mu)
  if (K <= 8){line_col = brewer.pal(K, "Dark2")} else{line_col = hcl(seq(15, 375, length = (K + 1)))}
  if (plot_radial_separate == FALSE){
    rad_prop = radial.plot(t(mu), rp.type = "p", labels = rownames(mu), show.grid.labels = TRUE, line.col = line_col,
                           lwd = 1.5, radial.lim = c(-1, 1), label.prop = 1.3)
    legend(radial_legend_pos[1], legend_pos[2], seq(1, K), col=line_col, lty=1, cex = 0.5, title = "Cluster", bty = "n")
  } else {
    rlist = vector(mode = "list", length = K)
    par(mfrow = c(ceiling(K/separate_col), radial_separate_col))
    for (j in 1:K){
      rlist[[j]] = radial.plot(t(mu[, j]), rp.type = "p", labels = rownames(mu), show.grid.labels = TRUE,
                               line.col = line_col[j], lwd = 1.5, radial.lim = c(-1, 1), label.prop = 1.3)
      title(paste('Cluster', j, sep = ' '))
    }
  }
  par(current_par)
}

row_norm = function(x){
  x = matrix(x, nrow = nrow(x))
  x_norm = sapply(1:nrow(x), function(j){
    x[j, ] / sqrt(sum(x[j, ]^2))
  })
  t(x_norm)
}

col_norm = function(x){
  x = matrix(x, nrow = nrow(x))
  x_norm = sapply(1:ncol(x), function(j){
    x[, j] / sqrt(sum(x[, j]^2))
  })
}

C_vMF = function(kappa, d){
  nu = d / 2 - 1
  kappa^nu / ((2 * pi)^(nu+1) * besselI(kappa, nu))
}

navmix_K = function(x0, K, pj_ini = 0.05, common_kappa = FALSE, no_ini = 5, tol = 1.0e-4, max_iter = 100){
  n = nrow(x0)
  m = ncol(x0)
  if (pj_ini > 0){
    no_par = (m + 1) * K + m + 1 + abs(as.numeric(common_kappa) - 1) * (K - 1)
  } else {
    no_par = (m + 1) * K + abs(as.numeric(common_kappa) - 1) * (K - 1)
  }
  x = row_norm(x0)
  M = col_norm(matrix(colSums(x), ncol = 1))
  l = vector(length = no_ini)
  for (i in 1:no_ini){
    if (K > 1){
      clust_ini = skmeans(x0, K, m = 1)
      g = sapply(1:K, function(k){
        gk = rep(0, n)
        gk[which(clust_ini$cluster == k)] = 1
        gk
      })
    } else {
      g = rep(1, n)
    }
    g = cbind((1 - pj_ini) * g, rep(pj_ini, n))
    l[i] = -Inf
    e = tol + 1
    iter = 0
    while (e > tol){
      r = t(x) %*% g[, 1:K, drop = FALSE]
      mu = cbind(col_norm(r), M)
      if (common_kappa == TRUE){
        r1 = sum(sqrt(colSums(r^2))) / sum(g[, 1:K])
        kappa = rep((r1*m - r1^3) / (1 - r1^2), K)
      } else {
        kappa = sapply(1:K, function(k){
          r1 = min(sqrt(sum(r[, k]^2)) / sum(g[, k]), 1)
          min((r1*m - r1^3) / (1 - r1^2), 500)
        })
      }
      kappa[(K+1)] = 0.0001
      p = colSums(g) / n
      C = sapply(1:(K+1), function(k){C_vMF(kappa[k], m)}) * p
      G = exp(x %*% mu %*% diag(kappa)) %*% diag(C)
      g = G / rowSums(G)
      l1 = sum(log(rowSums(G)))
      e = abs(l1 - l[i])
      l[i] = l1
      iter = iter + 1
      if (is.na(e)){
        l[i] = -Inf
        break
      }
      if (iter >= max_iter){
        break
      }
    }
    bic = 2 * l[i] - no_par * log(n)
    if (i == 1){
      mu_out = mu
      kappa_out = kappa
      g_out = g
      z_out = sapply(1:n, function(j){which.max(g[j, ])})
      bic_out = bic
      l_out = l[i]
    } else if (l[i] > l[(i - 1)]){
      mu_out = mu
      kappa_out = kappa
      g_out = g
      z_out = sapply(1:n, function(j){which.max(g[j, ])})
      bic_out = bic
      l_out = l[i]
    }
  }
  kappa_out = matrix(kappa_out, nrow = 1)
  cluster_labels = c(sprintf("Cluster %0d", seq(1:(ncol(g)-1))), "Noise")
  colnames(mu_out) = colnames(kappa_out) = colnames(g_out) = cluster_labels
  return(list(mu = mu_out, kappa = kappa_out, g = g_out, z = z_out, bic = bic_out, l = l_out))
}
