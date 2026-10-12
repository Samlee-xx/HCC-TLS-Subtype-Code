###### HookNet-TLS analysis in Cohort 7 ######

options(stringsAsFactors = FALSE)


###### directories ######

command.args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep('^--file=', command.args, value = TRUE)

if (length(file.arg) == 1) {
  script.dir <- dirname(normalizePath(
    sub('^--file=', '', file.arg),
    winslash = '/',
    mustWork = TRUE
  ))
} else {
  script.dir <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
}

hook.file <- file.path(
  script.dir,
  'data',
  '07-HookNet_TLS',
  'Cohort7_HookNet_TLS_patient_classification.csv'
)

output.dir <- file.path(
  script.dir,
  'result',
  'Figure7_TLSProfiler',
  '03-HookNet_TLS'
)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)


###### analysis settings ######

analysis.parameters <- data.frame(
  parameter = c(
    'TLS confidence value',
    'TLS confidence count',
    'GC confidence value',
    'GC confidence count'
  ),
  value = c(192L, 62L, 245L, 60L)
)

write.csv(
  analysis.parameters,
  file.path(output.dir, 'HookNet_TLS_analysis_parameters.csv'),
  row.names = FALSE
)


###### read the patient-level table ######

hook <- read.csv(
  hook.file,
  check.names = FALSE,
  fileEncoding = 'UTF-8-BOM'
)

required.columns <- c(
  'sample_id',
  'response_group',
  'hooknet_tls_status',
  'tls_count',
  'gc_positive_tls_count',
  'gc_count'
)

stopifnot(all(required.columns %in% colnames(hook)))
stopifnot(nrow(hook) == 120)
stopifnot(!anyDuplicated(hook$sample_id))
stopifnot(!anyNA(hook[, required.columns]))
stopifnot(all(hook$response_group %in% c('R', 'NR')))
stopifnot(all(hook$hooknet_tls_status %in% c(
  'non-TLS', 'GC-negative', 'GC-positive'
)))

hook$hooknet_tls_status <- factor(
  hook$hooknet_tls_status,
  levels = c('non-TLS', 'GC-negative', 'GC-positive')
)

hook$response_group <- factor(
  hook$response_group,
  levels = c('R', 'NR')
)

write.csv(
  hook,
  file.path(output.dir, 'Cohort7_HookNet_TLS_patient_classification.csv'),
  row.names = FALSE
)


###### observed counts and observed-to-expected ratios ######

count.table <- table(
  HookNet_group = hook$hooknet_tls_status,
  Response_group = hook$response_group
)
count.table

count.result <- data.frame(
  HookNet_group = rownames(count.table),
  R = as.integer(count.table[, 'R']),
  NR = as.integer(count.table[, 'NR']),
  row.names = NULL
)

expected.table <- outer(
  rowSums(count.table),
  colSums(count.table)
) / sum(count.table)

observed.expected <- data.frame(
  HookNet_group = rep(rownames(count.table), each = ncol(count.table)),
  Response_group = rep(colnames(count.table), times = nrow(count.table)),
  observed_n = as.integer(t(count.table)),
  expected_n = as.numeric(t(expected.table))
)

observed.expected$observed_expected_ratio <- with(
  observed.expected,
  observed_n / expected_n
)

write.csv(
  count.result,
  file.path(output.dir, 'HookNet_TLS_response_counts.csv'),
  row.names = FALSE
)

write.csv(
  observed.expected,
  file.path(output.dir, 'HookNet_TLS_observed_expected.csv'),
  row.names = FALSE
)


###### Fisher exact tests ######

compare.with.non.tls <- function(status, count.table) {

  test.table <- matrix(
    c(
      count.table[status, 'R'],
      count.table['non-TLS', 'R'],
      count.table[status, 'NR'],
      count.table['non-TLS', 'NR']
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Response = c('R', 'NR'),
      HookNet_group = c(status, 'non-TLS')
    )
  )

  fisher.result <- fisher.test(test.table, alternative = 'two.sided')

  odds.ratio <- (
    count.table[status, 'R'] * count.table['non-TLS', 'NR']
  ) / (
    count.table[status, 'NR'] * count.table['non-TLS', 'R']
  )

  data.frame(
    comparison = paste(status, 'versus non-TLS'),
    status_R = as.integer(count.table[status, 'R']),
    status_NR = as.integer(count.table[status, 'NR']),
    non_TLS_R = as.integer(count.table['non-TLS', 'R']),
    non_TLS_NR = as.integer(count.table['non-TLS', 'NR']),
    odds_ratio_for_response = as.numeric(odds.ratio),
    fisher_P_two_sided = fisher.result$p.value,
    row.names = NULL
  )
}

test.result <- do.call(
  rbind,
  lapply(
    c('GC-negative', 'GC-positive'),
    compare.with.non.tls,
    count.table = count.table
  )
)

test.result

write.csv(
  test.result,
  file.path(output.dir, 'HookNet_TLS_Fisher_test.csv'),
  row.names = FALSE
)


###### summary ######

print(analysis.parameters)
print(count.result)
print(observed.expected)
print(test.result)
