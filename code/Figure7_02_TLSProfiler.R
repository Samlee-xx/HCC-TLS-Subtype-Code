###### Figure S13A-B: TLSProfiler training and cohort comparison ######

library(jsonlite)
library(pROC)


###### input and output ######

script.file <- grep("^--file=", commandArgs(), value = TRUE)[1]
script.file <- sub("^--file=", "", script.file)
repo <- dirname(normalizePath(script.file))
source(file.path(repo, "tlsprofiler_roc.R"))

args <- commandArgs(trailingOnly = TRUE)
data.dir <- file.path(repo, "data/09-TLSProfiler")
out.dir <- file.path(repo, "result/Figure7_TLSProfiler/02-final_model")
if (length(args) >= 1) data.dir <- args[1]
if (length(args) >= 2) out.dir <- args[2]
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

model.definition <- fromJSON(file.path(data.dir, "model/tlsprofiler_log1p_BTN.json"))
cutoff <- model.definition$probability_cutoff


###### B, T and neutrophil percentages ######

read.cohort <- function(number, expected.n) {
  filename <- paste0("cohort", number, "_tls_features.csv")
  path <- file.path(data.dir, "CODEX_reference", filename)
  dat <- read.csv(path, check.names = FALSE, colClasses = c(slide = "character"))

  required <- c("tls_id", "slide", "group", "B_pct", "T_pct", "Neutrophil_pct")
  if (!all(required %in% names(dat))) stop("Missing cohort input columns")
  if (nrow(dat) != expected.n || anyDuplicated(dat$tls_id)) stop("Invalid TLS cohort")
  if (anyNA(dat[, required])) stop("Missing cohort metadata or cell percentages")
  if (!all(dat$group %in% c("BF-TLS", "TNC-TLS"))) stop("Invalid TLS subtypes")

  percentages <- as.matrix(dat[, c("B_pct", "T_pct", "Neutrophil_pct")])
  if (!all(is.finite(percentages)) || any(percentages < 0 | percentages > 100)) {
    stop("Expected cell percentages on the 0-100 scale")
  }
  if (any(rowSums(percentages) > 100 + 1e-8)) stop("Cell percentages exceed 100")

  dat$log_B <- log1p(dat$B_pct)
  dat$log_T <- log1p(dat$T_pct)
  dat$log_Neutrophil <- log1p(dat$Neutrophil_pct)
  dat$y <- as.integer(dat$group == "TNC-TLS")
  dat
}

cohort2 <- read.cohort(2, 667L)
cohort3 <- read.cohort(3, 154L)
formula <- y ~ log_B + log_T + log_Neutrophil


###### Cohort2 leave-one-slide-out predictions ######

cohort2$probability_TNC <- NA_real_
for (slide in unique(cohort2$slide)) {
  test.rows <- cohort2$slide == slide
  training <- cohort2[!test.rows, ]
  held.out <- cohort2[test.rows, ]

  fit <- glm(formula, data = training, family = binomial())
  cohort2$probability_TNC[test.rows] <- predict(fit, newdata = held.out, type = "response")
}


###### Fit on Cohort2 and predict Cohort3 ######

tlsprofiler <- glm(formula, data = cohort2, family = binomial())
cohort3$probability_TNC <- predict(tlsprofiler, newdata = cohort3, type = "response")

coefficient.names <- c("B_pct", "T_pct", "Neutrophil_pct")
reference.coefficients <- unlist(model.definition$coefficients[coefficient.names])
reference.coefficients <- c(model.definition$intercept, reference.coefficients)
coefficient.error <- unname(coef(tlsprofiler)) - reference.coefficients
if (max(abs(coefficient.error)) > 1e-8) stop("Fitted coefficients differ from the model definition")

cohort2$predicted_TLS_subtype <- ifelse(cohort2$probability_TNC >= cutoff, "TNC-TLS", "BF-TLS")
cohort3$predicted_TLS_subtype <- ifelse(cohort3$probability_TNC >= cutoff, "TNC-TLS", "BF-TLS")


###### ROC curves and slide-bootstrap intervals ######

roc2 <- roc(cohort2$y, cohort2$probability_TNC, levels = c(0, 1), direction = "<", quiet = TRUE)
roc3 <- roc(cohort3$y, cohort3$probability_TNC, levels = c(0, 1), direction = "<", quiet = TRUE)
ci2 <- slide_bootstrap_ci(cohort2$slide, cohort2$y, cohort2$probability_TNC, 5000L, 20260907L)
ci3 <- slide_bootstrap_ci(cohort3$slide, cohort3$y, cohort3$probability_TNC, 5000L, 20260908L)

auc.summary <- data.frame(
  cohort = c("Cohort 2", "Cohort 3"),
  evaluation = c("slide-LOSO", "cross-cohort"),
  n_TLS = c(nrow(cohort2), nrow(cohort3)),
  n_slides = c(length(unique(cohort2$slide)), length(unique(cohort3$slide))),
  AUC = c(as.numeric(auc(roc2)), as.numeric(auc(roc3))),
  CI_low = c(ci2[1], ci3[1]),
  CI_high = c(ci2[2], ci3[2]),
  probability_cutoff = cutoff
)
save_roc(roc2, ci2, "Cohort 2: slide-LOSO", file.path(out.dir, "cohort2_ROC.pdf"))
save_roc(roc3, ci3, "Cohort 3", file.path(out.dir, "cohort3_ROC.pdf"))


###### save model and prediction tables ######

coefficient.table <- data.frame(
  term = names(coef(tlsprofiler)),
  coefficient = unname(coef(tlsprofiler))
)
write.csv(coefficient.table, file.path(out.dir, "model_coefficients.csv"), row.names = FALSE)
write.csv(auc.summary, file.path(out.dir, "roc_summary.csv"), row.names = FALSE)
write.csv(cohort2, file.path(out.dir, "cohort2_loso_predictions.csv"), row.names = FALSE)
write.csv(cohort3, file.path(out.dir, "cohort3_predictions.csv"), row.names = FALSE)
saveRDS(list(model = tlsprofiler, definition = model.definition),
        file.path(out.dir, "TLSProfiler_log1p_BTN_model.rds"))
print(auc.summary)
