#Functions used for reproducing simulations and applied example
sstat = function(Y, X, intercept = TRUE){
  n = length(Y)
  if (intercept == TRUE){
    xx = cbind(rep(1, n), X)
  }
  else {xx = X}
  mod = lm.fit(xx, Y)
  bhat= c(mod$coefficients[2])
  s = t(mod$residuals) %*% (mod$residuals) / (mod$df.residual)
  se = sqrt((c(s) * solve(t(xx) %*% xx))[2,2])
  return(list("bhat" = bhat, "se" = se))
}

mrest = function(mrob, log_scale = TRUE, y_lab_pos = NULL){
  est_ivw = mr_ivw(mrob)
  est_med = mr_median(mrob)
  est_conmix = mr_conmix(mrob)
  est_egger = mr_egger(mrob)
  presso_df = data.frame("betaY" = mrob@betaY, "betaX" = mrob@betaX, "betaYse" = mrob@betaYse, "betaXse" = mrob@betaXse)
  est_presso = mr_presso(BetaOutcome = "betaY", BetaExposure = "betaX", SdOutcome = "betaYse", SdExposure = "betaXse",
                         data = presso_df, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 12000,
                         SignifThreshold = 0.05)
  if (is.na(est_presso$'Main MR results'[2, 3])){
    Est = data.frame("Method" = factor(c("MR-IVW", "MR-Median", "MR-ConMix", "MR-PRESSO"), levels = c("MR-PRESSO", "MR-ConMix", "MR-Median", "MR-IVW")),
                     "Estimate" = c(est_ivw$Estimate, est_med$Estimate, est_conmix$Estimate, est_presso$'Main MR results'[1, 3]),
                     "CILower" = c(est_ivw$CILower, est_med$CILower, est_conmix$CILower, est_presso$'Main MR results'[1, 3] - qt(0.975, dim(presso_df)[1]-2) * est_presso$'Main MR results'[1, 4]),
                     "CIUpper" = c(est_ivw$CIUpper, est_med$CIUpper, est_conmix$CIUpper, est_presso$'Main MR results'[1, 3] + qt(0.975, dim(presso_df)[1]-2) * est_presso$'Main MR results'[1, 4])
    )
  } else {
    Est = data.frame("Method" = factor(c("MR-IVW", "MR-Median", "MR-ConMix", "MR-PRESSO"), levels = c("MR-PRESSO", "MR-ConMix", "MR-Median", "MR-IVW")),
                     "Estimate" = c(est_ivw$Estimate, est_med$Estimate, est_conmix$Estimate, est_presso$'Main MR results'[2, 3]),
                     "CILower" = c(est_ivw$CILower, est_med$CILower, est_conmix$CILower, est_presso$'Main MR results'[2, 3] - qt(0.975, dim(presso_df)[1]-2) * est_presso$'Main MR results'[2, 4]),
                     "CIUpper" = c(est_ivw$CIUpper, est_med$CIUpper, est_conmix$CIUpper, est_presso$'Main MR results'[2, 3] + qt(0.975, dim(presso_df)[1]-2) * est_presso$'Main MR results'[2, 4])
    )
  }
  if (log_scale == TRUE){
    if (is.null(y_lab_pos)) {y_lab_pos = max(Est$CIUpper) * 1.1}
    Plot = ggplot(Est, aes(x = Method, y = Estimate, ymin = CILower, ymax = CIUpper)) +
      geom_pointrange(size = 0.25) + coord_flip(clip = "off") + theme_classic() +
      geom_hline(yintercept = 0, lty = 2) +
      geom_text(aes(y = y_lab_pos,
                    label = paste(sprintf("%.2f", round(Estimate, 2)), " (", sprintf("%.2f", round(CILower, 2)), ", ", sprintf("%.2f", round(CIUpper, 2)), ")", sep = "")),
                hjust = 0, size = 2.5) +
      ylab("Log OR (95% CI)") +
      theme(plot.margin = unit(c(5.5, 60, 5.5, 5.5), "pt"), axis.title.y = element_blank(), axis.text = element_text(size = 6))
  } else  {
    Est_OR = Est
    Est_OR[, -1] = exp(Est_OR[, -1])
    if (is.null(y_lab_pos)) {y_lab_pos = max(Est_OR$CIUpper) * 1.1}
    Plot = ggplot(Est_OR, aes(x = Method, y = Estimate, ymin = CILower, ymax = CIUpper)) +
      geom_pointrange(size = 0.25) + coord_flip(clip = "off") + theme_classic() +
      geom_hline(yintercept = 1, lty = 2) +
      geom_text(aes(y = y_lab_pos,
                    label = paste(sprintf("%.2f", round(Estimate, 2)), " (", sprintf("%.2f", round(CILower, 2)), ", ", sprintf("%.2f", round(CIUpper, 2)), ")", sep = "")),
                hjust = 0, size = 2.5) +
      ylab("OR (95% CI)") +
      theme(plot.margin = unit(c(5.5, 60, 5.5, 5.5), "pt"), axis.title.y = element_blank(), axis.text = element_text(size = 6))
  }
  return(list("Est" = Est, "Plot" = Plot))
}
