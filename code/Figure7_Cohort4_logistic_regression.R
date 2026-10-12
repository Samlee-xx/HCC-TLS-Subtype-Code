###### Cohort 4 patient-level logistic regression ######

library(dplyr)
library(ggplot2)
library(broom)

script.file <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1])
repo <- dirname(normalizePath(script.file))
args <- commandArgs(trailingOnly = TRUE)
data.dir <- file.path(repo, 'data/11-CODEX_Cohort4')
output.dir <- file.path(repo, 'result/Figure7_Cohort4/02-logistic_regression')
tls.file <- file.path(repo, 'result/Figure7_Cohort4/01-TLS_clustering/Cohort4_TLS_subtype_assignment.csv')
if (length(args) >= 1) data.dir <- args[1]
if (length(args) >= 2) output.dir <- args[2]
if (length(args) >= 3) tls.file <- args[3]
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(file.path(data.dir, 'patient_metadata.csv'), na.strings = c('', 'NA'))
tls <- read.csv(tls.file)
if (nrow(dat) != 103 || anyDuplicated(dat$patient_id)) stop('Expected 103 unique patients')
if (nrow(tls) != 55 || anyDuplicated(tls$tls_id)) stop('Expected 55 unique TLSs')
if (anyNA(tls$patient_id) || !all(tls$patient_id %in% dat$patient_id)) stop('Unmatched TLS patient')
if (anyNA(tls$TLS_subtype) || !all(tls$TLS_subtype %in% c('BF-TLS', 'TNC-TLS'))) stop('Invalid TLS subtype')

###### patient TLS burden and subtype ######

patient.tls <- tls %>%
  group_by(patient_id) %>%
  summarise(
    tls_count = n(),
    bf_tls_count = sum(TLS_subtype == 'BF-TLS'),
    tnc_tls_count = sum(TLS_subtype == 'TNC-TLS'),
    .groups = 'drop'
  )
reference <- dat[, c('patient_id', 'tls_count', 'bf_tls_count', 'tnc_tls_count', 'tls_subtype')]
dat <- dat[, !names(dat) %in% c('tls_count', 'bf_tls_count', 'tnc_tls_count', 'tls_subtype')]
dat <- left_join(dat, patient.tls, by = 'patient_id')
for (column in c('tls_count', 'bf_tls_count', 'tnc_tls_count')) {
  dat[[column]][is.na(dat[[column]])] <- 0L
}
dat$tls_subtype <- ifelse(dat$bf_tls_count > 0, 'BF-TLS',
                          ifelse(dat$tnc_tls_count > 0, 'TNC-TLS', 'Non-TLS'))
reference <- reference[match(dat$patient_id, reference$patient_id), ]
for (column in c('tls_count', 'bf_tls_count', 'tnc_tls_count', 'tls_subtype')) {
  if (!identical(as.character(dat[[column]]), as.character(reference[[column]]))) {
    stop(paste('Patient metadata disagrees with TLS assignments:', column))
  }
}

dat$tls_subtype <- factor(dat$tls_subtype, levels = c('Non-TLS', 'TNC-TLS', 'BF-TLS'))
dat$sex2 <- factor(dat$sex2, levels = c('Male', 'Female'))
dat$treatment2 <- factor(dat$treatment2, levels = c('Single ICI', 'Combination'))
dat$bclc_c <- factor(dat$bclc_c, levels = c('A/B', 'C'))
dat$vascular2 <- factor(dat$vascular2, levels = c('No', 'Yes'))
dat$multifocal_primary <- factor(dat$multifocal_primary, levels = c('Single/non-multiple', 'Multiple'))
dat$age_10 <- dat$age / 10
dat$tumor_diameter_10mm <- dat$tumor_diameter_mm / 10
dat$additional_tls_count <- pmax(dat$tls_count - 1, 0)
# All infiltration fractions use the total number of cells as the denominator.
dat$total_immune_10 <- dat$total_immune_ratio * 10
dat$lymphoid_total_10 <- dat$lymphoid_total_ratio * 10
dat$myeloid_total_10 <- dat$myeloid_total_ratio * 10
write.csv(dat, file.path(output.dir, 'Cohort4_103_patient_analysis_table.csv'), row.names = FALSE)
write.csv(as.data.frame(table(dat$tls_subtype, dat$response_group)),
          file.path(output.dir, 'Cohort4_patient_response_counts.csv'), row.names = FALSE)

###### univariable models and variable-level likelihood-ratio tests ######

formulas <- list(
  'TLS subtype' = responder ~ tls_subtype,
  'TLS burden' = responder ~ tls_count,
  'Age' = responder ~ age_10,
  'Sex' = responder ~ sex2,
  'Treatment' = responder ~ treatment2,
  'BCLC stage' = responder ~ bclc_c,
  'Vascular invasion' = responder ~ vascular2,
  'Tumor diameter' = responder ~ tumor_diameter_10mm,
  'Multifocal tumor' = responder ~ multifocal_primary,
  'Total immune infiltration' = responder ~ total_immune_10,
  'Lymphoid infiltration' = responder ~ lymphoid_total_10,
  'Myeloid infiltration' = responder ~ myeloid_total_10
)
univariable <- list()
for (name in names(formulas)) {
  formula <- formulas[[name]]
  complete <- dat[complete.cases(dat[, all.vars(formula)]), ]
  fit <- glm(formula, data = complete, family = binomial())
  null <- glm(responder ~ 1, data = complete, family = binomial())
  lrt.p <- anova(null, fit, test = 'LRT')[2, 'Pr(>Chi)']
  univariable[[name]] <- tidy(fit) %>%
    filter(term != '(Intercept)') %>%
    mutate(model = name, n_model = nobs(fit), variable_LRT_p = lrt.p,
           OR = exp(estimate),
           CI_low = exp(estimate - qnorm(0.975) * std.error),
           CI_high = exp(estimate + qnorm(0.975) * std.error))
}
univariable <- bind_rows(univariable)
write.csv(univariable, file.path(output.dir, 'Cohort4_univariable_logistic_regression.csv'), row.names = FALSE)
variable.tests <- unique(univariable[, c('model', 'n_model', 'variable_LRT_p')])
variable.tests$included_at_p_lt_0_10 <- variable.tests$variable_LRT_p < 0.10
write.csv(variable.tests, file.path(output.dir, 'Cohort4_variable_LRT.csv'), row.names = FALSE)

###### multivariable model ######

selected.models <- variable.tests$model[variable.tests$included_at_p_lt_0_10]
if (length(selected.models) == 0) stop('No variable met the LRT P < 0.10 criterion')
selected.terms <- vapply(formulas[selected.models], function(formula) all.vars(formula)[2], character(1))
# Use the number of TLSs beyond the first when adjusting for TLS subtype.
selected.terms[selected.terms == 'tls_count'] <- 'additional_tls_count'
model.formula <- reformulate(unname(selected.terms), response = 'responder')
complete <- dat[complete.cases(dat[, all.vars(model.formula)]), ]
fit.multi <- glm(model.formula, data = complete, family = binomial())
if (!fit.multi$converged) stop('Multivariable model did not converge')
multivariable <- tidy(fit.multi) %>%
  filter(term != '(Intercept)') %>%
  mutate(OR = exp(estimate),
         CI_low = exp(estimate - qnorm(0.975) * std.error),
         CI_high = exp(estimate + qnorm(0.975) * std.error),
         n_model = nobs(fit.multi))
write.csv(multivariable, file.path(output.dir, 'Cohort4_multivariable_logistic_regression.csv'), row.names = FALSE)
write.csv(data.frame(n_patient = nobs(fit.multi), n_responder = sum(complete$responder == 1),
                     n_non_responder = sum(complete$responder == 0)),
          file.path(output.dir, 'Cohort4_model_sample_size.csv'), row.names = FALSE)

###### forest plot ######

term.label <- c(
  'tls_subtypeBF-TLS' = 'BF-TLS vs Non-TLS',
  'tls_subtypeTNC-TLS' = 'TNC-TLS vs Non-TLS',
  'additional_tls_count' = 'TLS burden, per additional TLS',
  'vascular2Yes' = 'Vascular invasion: Yes vs No',
  'tumor_diameter_10mm' = 'Tumor diameter, per 10 mm',
  'multifocal_primaryMultiple' = 'Multifocal tumor: Yes vs No',
  'lymphoid_total_10' = 'Lymphoid / total cells, per 10%'
)
multivariable$label <- factor(term.label[multivariable$term], levels = rev(unname(term.label)))
forest <- ggplot(multivariable, aes(x = OR, y = label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = 'grey60') +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = 'y', width = 0.15, colour = 'grey40') +
  geom_point(size = 2, colour = '#377EB8') +
  scale_x_log10() +
  labs(x = 'Covariate-adjusted odds ratio for response', y = NULL) +
  theme_classic()
ggsave(file.path(output.dir, 'Cohort4_multivariable_forest.pdf'), forest, width = 7.5, height = 4.5)
print(table(dat$tls_subtype, dat$response_group))
print(multivariable[, c('term', 'OR', 'CI_low', 'CI_high', 'p.value', 'n_model')])
