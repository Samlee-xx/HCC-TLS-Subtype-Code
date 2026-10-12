###### Slide-bootstrap confidence intervals ######

slide_bootstrap_ci <- function(slide, y, score, n_boot, seed) {
  if (length(slide) != length(y) || length(y) != length(score)) {
    stop("Slide, outcome and score lengths differ")
  }
  if (anyNA(slide) || anyNA(y) || !all(is.finite(score))) {
    stop("Missing slide, outcome or score")
  }
  if (!setequal(unique(y), c(0, 1))) stop("Both reference classes are required")

  slide <- as.character(slide)
  slides <- unique(slide)
  rows.by.slide <- split(seq_along(slide), slide)
  bootstrap.auc <- rep(NA_real_, n_boot)

  set.seed(seed)
  for (i in seq_len(n_boot)) {
    sampled.slides <- sample(slides, length(slides), replace = TRUE)
    sampled.rows <- unlist(rows.by.slide[sampled.slides], use.names = FALSE)
    target <- y[sampled.rows]
    if (length(unique(target)) < 2L) next

    n.positive <- sum(target == 1L)
    n.negative <- sum(target == 0L)
    ranks <- rank(score[sampled.rows], ties.method = "average")
    positive.rank.sum <- sum(ranks[target == 1L])
    minimum.rank.sum <- n.positive * (n.positive + 1) / 2
    bootstrap.auc[i] <- (positive.rank.sum - minimum.rank.sum) / (n.positive * n.negative)
  }

  if (!any(is.finite(bootstrap.auc))) stop("No evaluable bootstrap samples")
  ci <- quantile(bootstrap.auc, c(0.025, 0.975), na.rm = TRUE, type = 7)
  unname(ci)
}


###### ROC figure ######

save_roc <- function(object, ci, title, path) {
  grDevices::pdf(path, width = 3.7, height = 3.7, useDingbats = FALSE)
  on.exit(grDevices::dev.off())

  false.positive.rate <- 1 - object$specificities
  sensitivity <- object$sensitivities
  plot(false.positive.rate, sensitivity,
       type = "l", col = "#2F78B7", lwd = 1.5,
       xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       xlab = "1 - Specificity", ylab = "Sensitivity", main = title)
  abline(0, 1, col = "grey75")

  auc.value <- as.numeric(pROC::auc(object))
  label <- sprintf("AUC: %.1f%%\n(%.1f%%-%.1f%%)",
                   100 * auc.value, 100 * ci[1], 100 * ci[2])
  text(0.98, 0.10, label, adj = c(1, 0), cex = 0.75)
}
