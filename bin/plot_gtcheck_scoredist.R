#!/usr/bin/env Rscript

## R-script to plot score distributions obtained with bcftools gtcheck

#library(svglite)

standardize.vector <- function(x)
{
  z <- (x - mean(x))/sd(x)
}

standardize.table <- function(dat)
{
  dz <- dat
  n_col = dim(dat)[2]
  for (j in c(1:n_col)) {
    dz[,j] <- standardize.vector(dat[,j])
  }
  dz
}

plot.densities <- function(
  z.scores,
  xlim = c(-4,4), ylim = c(0,0.5),
  colour = "blue",
  axis_label_type = 1,
  file = "bcftools_gtcheck_score_dist.pdf"
  )
{
  do_plot <- length(file) > 1
  # h1 <- hist(x1, breaks = 40, xlim = xlim, ylim = ylim, plot = F)
  x = sort(z.scores, decreasing = F)
  ds <- density(x)

  #svglite(file = file, width = 8, height = 8)
  if (do_plot) {
    cat("Opening file ", file, " for plotting ...\n")
    pdf(file = file, width = 8, height = 8)
  }
  # first plot densities on top of each other
  plot(ds, ann = F, xlim = xlim, ylim = ylim, col = cols[1])
  lines(x=c(x[1],x[1]), y=ylim, col = cols[1], lty = 1, lwd = 1)
  lines(x=c(x[2],x[2]), y=ylim, col = cols[1], lty = 2, lwd = 1)

  if (axis_label_type == 1 | axis_label_type == 2) {
    title( main = "discrepancy score distribution")
  }
  if (axis_label_type == 1 | axis_label_type == 3 | axis_label_type == 5) {
    title(ylab = 'density')
  }
  if (axis_label_type == 5 | axis_label_type == 6) {
    title(xlab = "z-score", sub = "bcftools gtmatch score")
  }

  ##title(
  ##  xlab = "z-score", ylab = "density",
  ##  main = "discrepancy score distribution",
  ##  sub = "bcftools gtmatch score"
  ##  )
  if (do_plot) {
    dev.off()
  }
}

## main
args = commandArgs(trailingOnly=TRUE)

score.file <- args[1]
#outfil.prfx <- args[2]
oufn <- args[2]

cols <- c("blue", "darkgreen")

dat <- read.csv(score.file, row.names = 1)
n_donors = dim(dat)[2]

z.scores <- standardize.table(dat)

pdf(file = oufn, width = 8.25, height = 11.75, paper = "a4", pointsize = 9)
par(mfrow = c(n_donors, 2))

for (i in c(1:n_donors)) {
  if (i == 1) {
    typ = 1
  } else if (i == n_donors) {
    typ = 5
  } else {
    typ = 3
  }
  donor.name = names(dat)[i]
  x <- as.vector(z.scores[,i])
  #plot.densities(x, colour = cols[1], xlim = range(x), ylim = c(0,0.6), file = paste(outfil.prfx, donor.name, "all.pdf", sep = "_"))
  #plot.densities(x, colour = cols[1], xlim = c(-4,4), ylim = c(0,0.6), file = paste(outfil.prfx, donor.name, "centre.pdf", sep = "_"))
  plot.densities(x, colour = cols[1], xlim = range(x), ylim = c(0,0.6), axis_label_type = typ, file ="")
  plot.densities(x, colour = cols[1], xlim = c(-4,4), ylim = c(0,0.6), axis_label_type = typ + 1, file ="")
}

dev.off()
