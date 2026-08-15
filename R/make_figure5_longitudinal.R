#!/usr/bin/env Rscript
# Figure 5: longitudinal change over the ~4-year interval.
# Reads the numerical results produced by R/reproduce_manuscript.R.
# Significance markers use the Bonferroni-adjusted p-values reported there.

regions <- c("AUD","SM","VV")
lc <- read.csv(file.path("outputs","longitudinal_change.csv"), stringsAsFactors=FALSE)
tc <- read.csv(file.path("outputs","tissue_fraction_change.csv"), stringsAsFactors=FALSE)

get_met <- function(met, measure, value) {
  d <- lc[lc$metabolite == met & lc$measure == measure, ]
  d <- d[match(regions, d$region), ]
  d[[value]]
}
get_tissue <- function(tissue, value) {
  d <- tc[tc$tissue == tissue, ]
  d <- d[match(regions, d$region), ]
  d[[value]]
}

gaba_unc   <- get_met("GABA","UNC","percent_change")
gaba_unc_p <- get_met("GABA","UNC","p_bonferroni")
gaba_cr    <- get_met("GABA","CR","percent_change")
gaba_cr_p  <- get_met("GABA","CR","p_bonferroni")
glx_unc    <- get_met("Glx","UNC","percent_change")
glx_unc_p  <- get_met("Glx","UNC","p_bonferroni")
glx_cr     <- get_met("Glx","CR","percent_change")
glx_cr_p   <- get_met("Glx","CR","p_bonferroni")

# Tissue values are fractions in the numerical output; convert to percentage points.
GM    <- 100 * get_tissue("GM_FRA","mean_change")
GM_p  <- get_tissue("GM_FRA","p_bonferroni")
WM    <- 100 * get_tissue("WM_FRA","mean_change")
WM_p  <- get_tissue("WM_FRA","p_bonferroni")
CSF   <- 100 * get_tissue("CSF_FRA","mean_change")
CSF_p <- get_tissue("CSF_FRA","p_bonferroni")

col_unc <- "#1b7a74"; col_cr <- "#e0a458"
tcol <- c(GM="#5b8a72", WM="#b8b0a1", CSF="#6a8ec9")
stars <- function(p) ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*","")))
ablab <- function(l) mtext(l, side=3, line=0.4, adj=0, font=2, cex=1.05)

pdf(file.path("outputs","longitudinal_change.pdf"), width=11, height=3.9, pointsize=11)
layout(matrix(1:3, nrow=1), widths=c(1,1,1.05))
par(mar=c(4.0,4.4,2.6,1.0), mgp=c(2.6,0.7,0), las=1, family="sans")

grp <- function(M, P, cols, names, ylim, ylab, title, tag, leg) {
  bp <- barplot(M, beside=TRUE, names.arg=names, col=cols, border=NA,
                ylim=ylim, ylab=ylab); abline(h=0, col="grey40")
  barplot(M, beside=TRUE, col=cols, border=NA, add=TRUE, axes=FALSE)
  off <- diff(ylim)*0.04
  for (i in seq_along(M)) {
    v <- M[i]; s <- stars(P[i]); if (s=="") next
    text(bp[i], v + ifelse(v>=0, off, -off), s, cex=0.9, col="grey25")
  }
  legend("topleft", leg, fill=cols, border=NA, bty="n", cex=0.82)
  ablab(tag); mtext(title, side=3, line=0.5, cex=0.95)
}

grp(rbind(gaba_unc,gaba_cr), rbind(gaba_unc_p,gaba_cr_p), c(col_unc,col_cr), regions,
    c(-22,6), "Change over interval (%)", "GABA+", "A",
    c("Water-referenced","Creatine-referenced"))
grp(rbind(glx_unc,glx_cr), rbind(glx_unc_p,glx_cr_p), c(col_unc,col_cr), regions,
    c(-12,16), "Change over interval (%)", "Glx", "B",
    c("Water-referenced","Creatine-referenced"))
grp(rbind(GM,WM,CSF), rbind(GM_p,WM_p,CSF_p), tcol, regions,
    c(-3,2.2), "Change (percentage points)", "Tissue fractions", "C",
    c("Gray matter","White matter","CSF"))

invisible(dev.off()); cat("wrote outputs/longitudinal_change.pdf\n")
