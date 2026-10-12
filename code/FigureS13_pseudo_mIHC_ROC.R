###### Figure S13E: pseudo-mIHC ROC curves ######

library(pROC)


###### input and output ######

script.file <- grep("^--file=", commandArgs(), value = TRUE)[1]
script.file <- sub("^--file=", "", script.file)
repo <- dirname(normalizePath(script.file))
source(file.path(repo, "tlsprofiler_roc.R"))

args <- commandArgs(trailingOnly = TRUE)
input.file <- file.path(repo, "result/FigureS13/pseudo_mIHC/tls_predictions.csv")
out.dir <- file.path(repo, "result/FigureS13/pseudo_mIHC/ROC")
if (length(args) >= 1) input.file <- args[1]
if (length(args) >= 2) out.dir <- args[2]
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

tls <- read.csv(input.file, check.names = FALSE, colClasses = c(slide = "character"))
required <- c("cohort", "slide", "tls_id", "target_TNC", "locked_lr_logit")
if (!all(required %in% names(tls))) stop("Missing TLS prediction columns")
if (anyNA(tls[, required])) stop("Missing TLS prediction values")
if (anyDuplicated(tls[, c("cohort", "tls_id")])) stop("Duplicate TLS identifiers")
if (!setequal(unique(tls$cohort), c("cohort2", "cohort3"))) stop("Expected Cohorts 2 and 3")


###### AUC and slide-bootstrap intervals ######

cohorts <- c("cohort2", "cohort3")
cohort.titles <- c("Cohort 2", "Cohort 3")
bootstrap.seeds <- c(20261024L, 20261124L)
auc.tables <- list()
coordinate.tables <- list()

for (i in seq_along(cohorts)) {
  cohort.tls <- tls[tls$cohort == cohorts[i], ]
  outcome <- cohort.tls$target_TNC
  score <- cohort.tls$locked_lr_logit
  if (!setequal(unique(outcome), c(0, 1))) stop("Both reference classes are required")
  if (!all(is.finite(score))) stop("Non-finite TLS scores")

  roc.curve <- roc(outcome, score, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- slide_bootstrap_ci(cohort.tls$slide, outcome, score, 2000L, bootstrap.seeds[i])

  auc.tables[[i]] <- data.frame(
    cohort = cohorts[i],
    n_tls = nrow(cohort.tls),
    n_slides = length(unique(cohort.tls$slide)),
    AUC = as.numeric(auc(roc.curve)),
    CI_low = ci[1],
    CI_high = ci[2]
  )
  coordinate.tables[[i]] <- data.frame(
    cohort = cohorts[i],
    threshold = roc.curve$thresholds,
    specificity = roc.curve$specificities,
    sensitivity = roc.curve$sensitivities
  )

  plot.file <- file.path(out.dir, paste0(cohorts[i], "_ROC.pdf"))
  save_roc(roc.curve, ci, cohort.titles[i], plot.file)
}

auc.summary <- do.call(rbind, auc.tables)
roc.coordinates <- do.call(rbind, coordinate.tables)
write.csv(auc.summary, file.path(out.dir, "roc_summary.csv"), row.names = FALSE)
write.csv(roc.coordinates, file.path(out.dir, "roc_coordinates.csv"), row.names = FALSE)
print(auc.summary)
