###### Figure 2L: subtype probabilities across serial sections ######

script.file <- grep("^--file=", commandArgs(), value = TRUE)[1]
script.file <- sub("^--file=", "", script.file)
repo <- dirname(normalizePath(script.file))
args <- commandArgs(trailingOnly = TRUE)
result.dir <- file.path(repo, "result/Figure2_serial")
if (length(args) >= 1) result.dir <- args[1]

tls.names <- paste0("TLS", 1:7)
data <- read.csv(file.path(result.dir, "serial_TLS_predictions_TLS.csv"))
sections <- seq(1, 49, 3)
probability <- matrix(NA_real_, nrow = 7, ncol = 17)

for (i in seq_along(tls.names)) {
  one.tls <- data[data$TLS == tls.names[i], ]
  one.tls <- one.tls[order(one.tls$section_number), ]
  if (nrow(one.tls) != 17 || !identical(as.integer(one.tls$section_number), as.integer(sections))) {
    stop("Incomplete TLS-section grid")
  }
  probability[i, ] <- one.tls$prob_TNC_logistic
}
if (anyNA(probability) || any(probability < 0 | probability > 1)) stop("Invalid TLS probabilities")

palette <- colorRampPalette(
  c("#2166AC", "#67A9CF", "#F2F2F2", "#EF8A62", "#B2182B"), space = "Lab"
)(1001)
tile.colors <- palette[round(probability * 1000) + 1]


###### heatmap and probability legend ######

draw.heatmap <- function() {
  par(fig = c(0, 0.82, 0, 1), mar = c(3.2, 4.2, 3.6, 0.5), family = "sans", xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0.5, 17.5), ylim = c(7.5, 0.5))
  x <- rep(1:17, each = 7)
  y <- rep(1:7, times = 17)
  rect(x - 0.48, y - 0.46, x + 0.48, y + 0.46, col = tile.colors, border = "white", lwd = 0.6)
  axis(1, at = 1:17, labels = sections, tck = -0.04, cex.axis = 0.88, mgp = c(2, 0.35, 0))
  axis(2, at = 1:7, labels = paste("TLS", 1:7), las = 1, tick = FALSE, cex.axis = 0.88)
  mtext("Sections", side = 1, line = 1.9, adj = 0)
  mtext("Stability analysis of TLS subtypes", side = 3, line = 2.25, cex = 1.05)
  legend("top", inset = c(0, -0.22), legend = c("BF-TLS", "TNC-TLS"),
         fill = c("#2166AC", "#B2182B"), border = NA, horiz = TRUE, bty = "n", xpd = NA, cex = 0.9)

  par(fig = c(0.90, 0.925, 0.23, 0.72), mar = c(0, 0, 0, 0), new = TRUE, xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  edges <- seq(0, 1, length.out = 1002)
  rect(0, edges[-1002], 1, edges[-1], col = palette, border = NA)
  axis(4, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), las = 1, tick = FALSE, cex.axis = 0.8, mgp = c(0, 0.1, 0))
  mtext("TNC-TLS", side = 3, line = 0.3, cex = 0.8)
  mtext("BF-TLS", side = 1, line = 0.3, cex = 0.8)
  mtext("Predicted probability", side = 2, line = 0.4, cex = 0.85)
}

pdf(file.path(result.dir, "Figure2L_serial_TLS_probability.pdf"), width = 8.1, height = 3.25, useDingbats = FALSE)
draw.heatmap()
dev.off()
png(file.path(result.dir, "Figure2L_serial_TLS_probability.png"), width = 2430, height = 975, res = 300)
draw.heatmap()
dev.off()
