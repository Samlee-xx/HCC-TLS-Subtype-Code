###### Figure S13J: Cohort6 ROC curve ######

library(pROC)


###### input and output ######

script.file <- grep("^--file=", commandArgs(), value = TRUE)[1]
script.file <- sub("^--file=", "", script.file)
repo <- dirname(normalizePath(script.file))
source(file.path(repo, "tlsprofiler_roc.R"))

args <- commandArgs(trailingOnly = TRUE)
input.dir <- file.path(repo, "result/FigureS13/Cohort6")
if (length(args) >= 1) input.dir <- args[1]
out.dir <- file.path(input.dir, "ROC")
if (length(args) >= 2) out.dir <- args[2]
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

tls <- read.csv(file.path(input.dir, "tls_comparison.csv"),
                check.names = FALSE, colClasses = c(slide = "character"))
auc.summary <- read.csv(file.path(input.dir, "roc_summary.csv"), check.names = FALSE)

if (nrow(tls) != 184L || length(unique(tls$slide)) != 52L) stop("Invalid Cohort6 TLS set")
if (anyDuplicated(tls$tls_id)) stop("Duplicate TLS identifiers")
if (!setequal(unique(tls$reference_tls_subtype), c("BF-TLS", "TNC-TLS"))) {
  stop("Invalid reference TLS subtypes")
}


###### ROC and the slide-bootstrap interval from the comparison script ######

outcome <- as.integer(tls$reference_tls_subtype == "TNC-TLS")
score <- tls$automated_prob_TNC
if (anyNA(outcome) || !all(is.finite(score))) stop("Missing outcomes or non-finite scores")

roc.curve <- roc(outcome, score, levels = c(0, 1), direction = "<", quiet = TRUE)
auc.value <- as.numeric(auc(roc.curve))
if (nrow(auc.summary) != 1L) stop("Expected one Cohort6 ROC summary")
if (abs(auc.value - auc.summary$auc) > 1e-10) stop("AUC does not match the prediction table")

ci <- c(auc.summary$CI_low, auc.summary$CI_high)
if (anyNA(ci) || any(ci < 0 | ci > 1) || ci[1] > ci[2]) stop("Invalid confidence interval")

plot.file <- file.path(out.dir, "Cohort6_ROC.pdf")
save_roc(roc.curve, ci, "Cohort 6: method comparison", plot.file)

roc.coordinates <- data.frame(
  threshold = roc.curve$thresholds,
  specificity = roc.curve$specificities,
  sensitivity = roc.curve$sensitivities
)
write.csv(roc.coordinates, file.path(out.dir, "roc_coordinates.csv"), row.names = FALSE)
write.csv(auc.summary, file.path(out.dir, "roc_summary.csv"), row.names = FALSE)
print(auc.summary)
