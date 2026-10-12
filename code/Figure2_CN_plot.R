# CN plots for the final 667 TLSs. Run Figure2_CN.py first.
# Figure 2G: Wilcoxon test and BH correction across nine CNs.

script.arg = grep('^--file=', commandArgs(), value = TRUE)
if (length(script.arg) != 1) stop('Run this file with Rscript.')
script.dir = dirname(normalizePath(sub('^--file=', '', script.arg)))
args = commandArgs(trailingOnly = TRUE)
data.dir = if (length(args) >= 1) args[1] else file.path(script.dir, 'data')
result.dir = if (length(args) >= 2) args[2] else file.path(script.dir, 'result/Figure2_CN')
extra.library = Sys.getenv('CN_R_LIBRARY', '')
if (nzchar(extra.library)) .libPaths(c(extra.library, .libPaths()))

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(pheatmap)
library(scales)

dir.create(file.path(result.dir, "02-figures/01-CN-celltype-enrichment"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result.dir, "02-figures/02-CN-layer-ROE"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result.dir, "02-figures/03-CN-alluvial"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result.dir, "02-figures/04-CN-boxplot"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(result.dir, "02-figures/05-CN-layer-proportion"), recursive = TRUE, showWarnings = FALSE)

cn.order = paste0('CN', 1:9)
celltype.order = c(
  'B', 'CD4T', 'CD8T', 'T-B', 'Cholangiocyte', 'DC', 'Endothelial',
  'Fibroblast', 'HEV', 'Lymphatic', 'Macrophage', 'Neutrophil', 'Unassigned'
)
tls.group.order = c('BF-TLS', 'TNC-TLS')
layer.order = c('Edge', 'Transition', 'Core')


# ---------- 01 CN cell-type enrichment heatmap ----------
cn.mean = read.csv(
  file.path(data.dir, "02-CODEX_Cohort2/CN", "cn_celltype_mean_k9_final.csv"),
  check.names = FALSE
)

colnames(cn.mean)[1] = 'CN'
cn.mean$CN = gsub('CN_', 'CN', cn.mean$CN)
cn.mean = cn.mean[match(cn.order, cn.mean$CN), ]

final.plot = cn.mean[, celltype.order]
rownames(final.plot) = cn.mean$CN
final.plot = as.matrix(final.plot)

final.plot.z = scale(final.plot)

write.csv(
  final.plot,
  file.path(result.dir, "01-tables/CN_celltype_mean_matrix_final.csv")
)

write.csv(
  final.plot.z,
  file.path(result.dir, "01-tables/CN_celltype_zscore_matrix_final.csv")
)

colors <- colorRampPalette(c("#51a16a", "white", "red"))(99)
breaks = seq(-2, 2, length.out = 100)

pdf(
  file = file.path(result.dir, "02-figures/01-CN-celltype-enrichment/CN_celltype_enrichment_heatmap.pdf"),
  height = 3.2,
  width = 5.2
)

a = pheatmap(
  final.plot.z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 15,
  cellheight = 15,
  color = colors,
  breaks = breaks
)

dev.off()

final.plot.z.range1 = final.plot.z
final.plot.z.range1[final.plot.z.range1 > 1] = 1
final.plot.z.range1[final.plot.z.range1 < -1] = -1

pdf(
  file = file.path(result.dir, "02-figures/01-CN-celltype-enrichment/CN_celltype_enrichment_heatmap_range_minus1_to_1.pdf"),
  height = 3.2,
  width = 5.2
)

a = pheatmap(
  final.plot.z.range1,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 15,
  cellheight = 15,
  color = colors,
  breaks = seq(-1, 1, length.out = 100)
)

dev.off()


# ---------- 02 CN enrichment in TLS layers ----------
total.meta = read.csv(
  file.path(result.dir, "01-tables/cn_tls_layer_cell_table_final.csv"),
  check.names = FALSE
)

total.meta = total.meta %>%
  filter(!is.na(tls_subtype), !is.na(CN), !is.na(Boundary_Class)) %>%
  filter(CN %in% cn.order, Boundary_Class %in% layer.order)

total.meta$CN = factor(total.meta$CN, levels = cn.order)
total.meta$tls_subtype = factor(total.meta$tls_subtype, levels = tls.group.order)
total.meta$Boundary_Class = factor(total.meta$Boundary_Class, levels = layer.order)

roe.source = data.frame()
roe.matrix = NULL
p.matrix = NULL

for (tls.group in tls.group.order) {
  total.meta.sub = subset(total.meta, tls_subtype == tls.group)
  class.boundary = as.matrix(table(total.meta.sub$CN, total.meta.sub$Boundary_Class))
  class.boundary = class.boundary[cn.order, layer.order]

  roe = class.boundary
  p.value = class.boundary

  for (i in 1:dim(class.boundary)[1]) {
    for (j in 1:dim(class.boundary)[2]) {
      onecn_onelayer = class.boundary[i, j]
      onelayer_othercn = sum(class.boundary[i, ]) - onecn_onelayer
      onecn_otherlayer = sum(class.boundary[, j]) - onecn_onelayer
      othercn_otherlayer = sum(class.boundary) - onecn_onelayer - onelayer_othercn - onecn_otherlayer

      x = matrix(
        c(onecn_onelayer, onecn_otherlayer, onelayer_othercn, othercn_otherlayer),
        nrow = 2,
        ncol = 2
      )

      chisq.result = suppressWarnings(chisq.test(x))
      roe[i, j] = onecn_onelayer / chisq.result$expected[1, 1]
      p.value[i, j] = chisq.result$p.value

      roe.source = rbind(
        roe.source,
        data.frame(
          tls_subtype = tls.group,
          CN = rownames(class.boundary)[i],
          Boundary_Class = colnames(class.boundary)[j],
          n_cell = onecn_onelayer,
          expected = chisq.result$expected[1, 1],
          roe = onecn_onelayer / chisq.result$expected[1, 1],
          p_value = chisq.result$p.value
        )
      )
    }
  }

  colnames(roe) = paste(tls.group, colnames(roe), sep = '_')
  colnames(p.value) = paste(tls.group, colnames(p.value), sep = '_')
  roe.matrix = cbind(roe.matrix, roe)
  p.matrix = cbind(p.matrix, p.value)
}

roe.matrix = as.matrix(roe.matrix)
storage.mode(roe.matrix) = 'numeric'
p.matrix = as.matrix(p.matrix)
storage.mode(p.matrix) = 'numeric'

write.csv(
  roe.source,
  file.path(result.dir, "01-tables/CN_layer_ROE_source_data_final.csv"),
  row.names = FALSE
)

write.csv(
  roe.matrix,
  file.path(result.dir, "01-tables/CN_layer_ROE_matrix_final.csv")
)

write.csv(
  p.matrix,
  file.path(result.dir, "01-tables/CN_layer_ROE_pvalue_matrix_final.csv")
)

roe.plot = roe.matrix
roe.plot[roe.plot > 2] = 2
roe.plot[roe.plot < 0] = 0

pdf(
  file = file.path(result.dir, "02-figures/02-CN-layer-ROE/CN_layer_ROE_heatmap.pdf"),
  height = 3.2,
  width = 4.3
)

a = pheatmap(
  as.matrix(roe.plot),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = c(3),
  cellwidth = 15,
  cellheight = 15,
  color = colorRampPalette(c("gray90", "white", "red"))(100),
  breaks = seq(0, 2, length.out = 101),
  labels_col = c('Edge', 'Transition', 'Core', 'Edge', 'Transition', 'Core'),
  main = 'BF-TLS              TNC-TLS'
)

dev.off()


# ---------- 03 CN alluvial plot ----------
cn.prop = read.csv(
  file.path(result.dir, "01-tables/cn_tls_proportion_long_final.csv"),
  check.names = FALSE
)

colnames(cn.prop)[1] = 'tls_id'
cn.prop$CN = factor(cn.prop$CN, levels = cn.order)
cn.prop$tls_subtype = factor(cn.prop$tls_subtype, levels = tls.group.order)

ct.tls.mean = cn.prop %>%
  group_by(tls_subtype, CN) %>%
  summarise(percentage = mean(percentage), .groups = 'drop')

write.csv(
  ct.tls.mean,
  file.path(result.dir, "01-tables/CN_proportion_alluvial_source_data_final.csv"),
  row.names = FALSE
)

colorlist <- c(
  "CN1" = "#ef8a8c",
  "CN2" = "#2c7fb8",
  "CN3" = "#359a35",
  "CN4" = "#e31a1c",
  "CN5" = "#b2d58a",
  "CN6" = "#a6cee3",
  "CN7" = "#f6bd6b",
  "CN8" = "#f47c13",
  "CN9" = "#6a3d9a"
)

p = ggplot(ct.tls.mean, aes(x = tls_subtype, y = percentage, fill = CN,
                            stratum = CN, alluvium = CN)) +
  geom_stratum(width = 0.5, color = "white") +
  geom_alluvium(alpha = 0.5, width = 0.5, curve_type = "linear") +
  scale_fill_manual(values = colorlist) +
  ylab('Percentage') +
  xlab(NULL) +
  theme_classic()

ggsave(
  file.path(result.dir, "02-figures/03-CN-alluvial/CN_proportion_alluvial_by_TLS_subtype.pdf"),
  p,
  width = 4.6,
  height = 3.4
)


# ---------- 04 CN proportion boxplot ----------
prop_df = cn.prop %>%
  filter(!is.na(tls_id), !is.na(tls_subtype), !is.na(CN))

pvalue.table = prop_df %>%
  group_by(CN) %>%
  summarise(
    p_value = wilcox.test(percentage[tls_subtype == 'BF-TLS'], percentage[tls_subtype == 'TNC-TLS'], exact = FALSE, correct = TRUE, alternative = 'two.sided')$p.value,
    BF_TLS_n = sum(tls_subtype == 'BF-TLS'),
    TNC_TLS_n = sum(tls_subtype == 'TNC-TLS'),
    BF_TLS_median = median(percentage[tls_subtype == 'BF-TLS']),
    TNC_TLS_median = median(percentage[tls_subtype == 'TNC-TLS']),
    BF_TLS_mean = mean(percentage[tls_subtype == 'BF-TLS']),
    TNC_TLS_mean = mean(percentage[tls_subtype == 'TNC-TLS']),
    .groups = 'drop'
  )

pvalue.table$BH_FDR = p.adjust(pvalue.table$p_value, method = 'BH')
pvalue.table$test = 'Two-sided Wilcoxon rank-sum; tie and continuity correction'
pvalue.table$correction_family = 'All 9 CN comparisons'
pvalue.table$plot_label = paste0('BH FDR = ', format.pval(pvalue.table$BH_FDR, digits = 2, eps = 1e-20))
set.seed(20260618)

write.csv(
  prop_df,
  file.path(result.dir, "01-tables/CN_proportion_boxplot_source_data_final.csv"),
  row.names = FALSE
)

write.csv(
  pvalue.table,
  file.path(result.dir, "01-tables/CN_proportion_boxplot_pvalue_final.csv"),
  row.names = FALSE
)

p = ggplot(prop_df, aes(x = tls_subtype, y = percentage / 100)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = tls_subtype), width = 0.15, size = 0.8, alpha = 0.5, show.legend = FALSE) +
  facet_wrap(~ CN, scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.18))) +
  scale_color_manual(values = c('BF-TLS' = '#c095e4', 'TNC-TLS' = '#72b043')) +
  labs(x = NULL, y = "Proportion within TLS", title = "CN proportion across TLS subtypes") +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = NA),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) +
  geom_text(
    data = pvalue.table, aes(x = 1.5, y = Inf, label = plot_label),
    inherit.aes = FALSE, vjust = 1.5, size = 3
  )


ggsave(
  file.path(result.dir, "02-figures/04-CN-boxplot/CN_proportion_boxplot_by_TLS_subtype.pdf"),
  p,
  width = 6,
  height = 7
)


# ---------- 05 CN percentage within each TLS layer ----------
cn.layer.prop = read.csv(
  file.path(result.dir, "01-tables/CN_layer_proportion_by_TLS_final.csv"),
  check.names = FALSE
)

cn.layer.prop = cn.layer.prop %>%
  filter(!is.na(percentage)) %>%
  filter(CN %in% cn.order, tls_subtype %in% tls.group.order) %>%
  filter(delaunay_layer %in% c('Periphery', 'Transition', 'Core'))

cn.layer.prop$CN = factor(cn.layer.prop$CN, levels = cn.order)
cn.layer.prop$tls_subtype = factor(cn.layer.prop$tls_subtype, levels = tls.group.order)
cn.layer.prop$delaunay_layer = factor(
  cn.layer.prop$delaunay_layer,
  levels = c('Periphery', 'Transition', 'Core')
)

cn.layer.summary = cn.layer.prop %>%
  group_by(tls_subtype, delaunay_layer, CN) %>%
  summarise(
    mean_percentage = mean(percentage),
    sd_percentage = sd(percentage),
    n_TLS = n(),
    sem_percentage = sd_percentage / sqrt(n_TLS),
    .groups = 'drop'
  )

write.csv(
  cn.layer.summary,
  file.path(result.dir, "01-tables/CN_layer_proportion_mean_final.csv"),
  row.names = FALSE
)

cn.layer.pvalue = data.frame()

for (one.cn in cn.order) {
  for (one.layer in c('Periphery', 'Transition', 'Core')) {
    one.table = subset(cn.layer.prop, CN == one.cn & delaunay_layer == one.layer)
    bf.value = one.table$percentage[one.table$tls_subtype == 'BF-TLS']
    tnc.value = one.table$percentage[one.table$tls_subtype == 'TNC-TLS']

    p.t = t.test(bf.value, tnc.value)$p.value
    p.w = wilcox.test(bf.value, tnc.value)$p.value

    cn.layer.pvalue = rbind(
      cn.layer.pvalue,
      data.frame(
        CN = one.cn,
        delaunay_layer = one.layer,
        BF_TLS_mean = mean(bf.value),
        TNC_TLS_mean = mean(tnc.value),
        difference_TNC_minus_BF = mean(tnc.value) - mean(bf.value),
        BF_TLS_n = length(bf.value),
        TNC_TLS_n = length(tnc.value),
        p_value_ttest = p.t,
        p_value_wilcoxon = p.w
      )
    )
  }
}

cn.layer.pvalue$fdr_ttest = p.adjust(cn.layer.pvalue$p_value_ttest, method = 'BH')

write.csv(
  cn.layer.pvalue,
  file.path(result.dir, "01-tables/CN_layer_proportion_difference_pvalue_final.csv"),
  row.names = FALSE
)

mean.mat = matrix(
  NA,
  nrow = length(cn.order),
  ncol = length(tls.group.order) * 3
)

rownames(mean.mat) = cn.order
colnames(mean.mat) = c(
  'BF-TLS_Edge', 'BF-TLS_Transition', 'BF-TLS_Core',
  'TNC-TLS_Edge', 'TNC-TLS_Transition', 'TNC-TLS_Core'
)

for (i in 1:nrow(cn.layer.summary)) {
  one.cn = as.character(cn.layer.summary$CN[i])
  one.group = as.character(cn.layer.summary$tls_subtype[i])
  one.layer = as.character(cn.layer.summary$delaunay_layer[i])

  if (one.layer == 'Periphery') {
    one.layer.label = 'Edge'
  } else {
    one.layer.label = one.layer
  }

  mean.mat[one.cn, paste(one.group, one.layer.label, sep = '_')] =
    cn.layer.summary$mean_percentage[i]
}

write.csv(
  mean.mat,
  file.path(result.dir, "01-tables/CN_layer_proportion_mean_matrix_final.csv")
)

pdf(
  file = file.path(result.dir, "02-figures/05-CN-layer-proportion/CN_layer_proportion_mean_heatmap.pdf"),
  height = 3.2,
  width = 4.3
)

a = pheatmap(
  mean.mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = c(3),
  cellwidth = 15,
  cellheight = 15,
  color = colorRampPalette(c("gray90", "white", "red"))(100),
  breaks = seq(0, 25, length.out = 101),
  labels_col = c('Edge', 'Transition', 'Core', 'Edge', 'Transition', 'Core'),
  main = 'BF-TLS              TNC-TLS'
)

dev.off()

print('Saved CN cell-type enrichment heatmap.')
print('Saved CN layer ROE heatmap.')
print('Saved CN alluvial plot.')
print('Saved CN proportion boxplot.')
print('Saved CN layer proportion heatmap.')

writeLines(capture.output(sessionInfo()), file.path(result.dir, 'R_sessionInfo.txt'))
