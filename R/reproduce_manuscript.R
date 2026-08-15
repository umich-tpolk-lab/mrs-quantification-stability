#!/usr/bin/env Rscript
# Reproduction code for:
# "How quantification choices shape MRS measures of GABA and Glx:
#  Longitudinal stability and change in older adults"
#
# Run from the repository root:
#   Rscript R/reproduce_manuscript.R
#
# This script uses only the de-identified, manuscript-specific files in data/.
# It reproduces the primary numerical analyses reported in the manuscript and
# writes machine-readable results to outputs/.

options(stringsAsFactors = FALSE)

DATA_FILE <- file.path("data", "mrs_measurements.csv")
QC_FILE <- file.path("data", "gannet_qc_metrics.csv")
PARTICIPANTS_FILE <- file.path("data", "participants.csv")
OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

REGIONS <- c("AUD","SM","VV")
MEASURES <- c("UNC","CSF","TC","ATC","CR")
METABOLITES <- c(GABA="G", Glx="GLX")
B_BOOT <- 10000
SEED <- 42

# ----------------------- Data and QC -------------------------------------------
load_public_data <- function() {
  ms <- read.csv(DATA_FILE, check.names = FALSE)
  qc <- read.csv(QC_FILE, check.names = FALSE)
  participants <- read.csv(PARTICIPANTS_FILE, check.names = FALSE)

  stopifnot(length(unique(ms$participant_id)) == 56)
  stopifnot(nrow(ms) == 112)
  stopifnot(nrow(qc) == 668)
  stopifnot(all(ms$session %in% c(1,2)))
  list(ms=ms, qc=qc, participants=participants)
}

robust_outlier <- function(x, k=8) {
  med <- median(x, na.rm=TRUE)
  md <- mad(x, center=med, na.rm=TRUE)  # R default scale = 1.4826
  if (!is.finite(md) || md == 0) return(rep(FALSE, length(x)))
  is.finite(x) & abs(x-med)/md > k
}

apply_qc <- function(ms, qc, linewidth_max=15) {
  qc$hr <- paste0(qc$hemisphere, qc$region)
  lw_bad <- is.finite(qc$lw_cr_hz) & qc$lw_cr_hz > linewidth_max
  ratio_bad_g <- robust_outlier(qc$gaba_cr, 8)
  ratio_bad_x <- robust_outlier(qc$glx_cr, 8)

  # Fit error is reported but is not used as an exclusion gate.
  qc$exclude_G <- lw_bad | ratio_bad_g
  qc$exclude_GLX <- lw_bad | ratio_bad_x

  for (met in c("G","GLX")) {
    bad <- qc[qc[[paste0("exclude_", met)]], c("participant_id","session","hr")]
    if (nrow(bad)) {
      for (i in seq_len(nrow(bad))) {
        rows <- ms$participant_id == bad$participant_id[i] & ms$session == bad$session[i]
        for (me in MEASURES) {
          col <- paste(met, me, bad$hr[i], sep="_")
          if (col %in% names(ms)) ms[rows, col] <- NA_real_
        }
      }
    }
  }

  # Bilateral metabolite values are the mean of QC-passing hemispheres.
  for (met in c("G","GLX")) for (me in MEASURES) for (rg in REGIONS) {
    lcol <- paste(met, me, paste0("L",rg), sep="_")
    rcol <- paste(met, me, paste0("R",rg), sep="_")
    bcol <- paste(met, me, rg, sep="_")
    ms[[bcol]] <- rowMeans(cbind(ms[[lcol]], ms[[rcol]]), na.rm=TRUE)
    ms[[bcol]][is.nan(ms[[bcol]])] <- NA_real_
  }

  # Bilateral tissue fractions use both hemispheres (as in the source mastersheet).
  for (fr in c("GM_FRA","WM_FRA","CSF_FRA")) for (rg in REGIONS) {
    lcol <- paste0(fr, "_L", rg)
    rcol <- paste0(fr, "_R", rg)
    bcol <- paste0(fr, "_", rg)
    ms[[bcol]] <- rowMeans(cbind(ms[[lcol]], ms[[rcol]]), na.rm=TRUE)
    ms[[bcol]][is.nan(ms[[bcol]])] <- NA_real_
  }

  attr(ms, "qc_counts") <- c(
    available=nrow(qc),
    excluded_GABA=sum(qc$exclude_G),
    excluded_Glx=sum(qc$exclude_GLX),
    linewidth=sum(lw_bad),
    ratio_GABA=sum(ratio_bad_g),
    ratio_Glx=sum(ratio_bad_x)
  )
  ms
}

# ----------------------- ICC helpers -------------------------------------------
icc_c1 <- function(t1, t2) {
  ok <- is.finite(t1) & is.finite(t2)
  t1 <- t1[ok]; t2 <- t2[ok]
  n <- length(t1)
  if (n < 3) return(c(icc=NA_real_, n=n))
  X <- cbind(t1,t2); k <- 2
  grand <- mean(X)
  row_m <- rowMeans(X); col_m <- colMeans(X)
  MSR <- k * sum((row_m-grand)^2)/(n-1)
  SST <- sum((X-grand)^2)
  SSE <- SST - k*sum((row_m-grand)^2) - n*sum((col_m-grand)^2)
  MSE <- SSE / ((n-1)*(k-1))
  c(icc=(MSR-MSE)/(MSR+(k-1)*MSE), n=n)
}

pair_of <- function(ms, col, ids=NULL) {
  if (is.null(ids)) ids <- sort(unique(ms$participant_id))
  d <- ms[ms$participant_id %in% ids, c("participant_id","session",col)]
  names(d) <- c("s","w","v")
  d$v[is.finite(d$v) & d$v <= 0] <- NA_real_
  w <- reshape(d, idvar="s", timevar="w", direction="wide")
  idx <- match(ids, w$s)
  cbind(w[["v.1"]][idx], w[["v.2"]][idx])
}

bootstrap_icc <- function(M, B=B_BOOT) {
  Mc <- M[complete.cases(M),,drop=FALSE]
  n <- nrow(Mc)
  obs <- icc_c1(Mc[,1], Mc[,2])[["icc"]]
  bv <- replicate(B, {
    ix <- sample.int(n, n, replace=TRUE)
    icc_c1(Mc[ix,1], Mc[ix,2])[["icc"]]
  })
  ci <- quantile(bv, c(.025,.975), na.rm=TRUE)
  c(icc=obs, lo=ci[[1]], hi=ci[[2]], n=n)
}

bootstrap_cross_roi_mean <- function(ms, met, me, B=B_BOOT, ids=NULL) {
  if (is.null(ids)) ids <- sort(unique(ms$participant_id))
  Ms <- setNames(lapply(REGIONS, function(rg)
    pair_of(ms, paste(met, me, rg, sep="_"), ids)), REGIONS)

  obs_each <- sapply(Ms, function(M) icc_c1(M[,1], M[,2])[["icc"]])
  obs <- mean(obs_each, na.rm=TRUE)

  bv <- replicate(B, {
    ix <- sample.int(length(ids), length(ids), replace=TRUE)
    vals <- sapply(Ms, function(M) icc_c1(M[ix,1], M[ix,2])[["icc"]])
    mean(vals, na.rm=TRUE)
  })
  ci <- quantile(bv, c(.025,.975), na.rm=TRUE)
  c(icc=obs, lo=ci[[1]], hi=ci[[2]], n=length(ids))
}

make_table1 <- function(ms) {
  set.seed(SEED)
  rows <- list(); z <- 1
  for (metname in names(METABOLITES)) {
    met <- METABOLITES[[metname]]
    for (me in MEASURES) {
      for (rg in REGIONS) {
        r <- bootstrap_icc(pair_of(ms, paste(met,me,rg,sep="_")))
        rows[[z]] <- data.frame(metabolite=metname, region=rg, measure=me,
                                icc=r["icc"], ci_low=r["lo"], ci_high=r["hi"], n=r["n"])
        z <- z + 1
      }
      r <- bootstrap_cross_roi_mean(ms, met, me)
      rows[[z]] <- data.frame(metabolite=metname, region="Cross-ROI mean", measure=me,
                              icc=r["icc"], ci_low=r["lo"], ci_high=r["hi"], n=r["n"])
      z <- z + 1
    }
  }
  do.call(rbind, rows)
}

# ----------------------- Paired ICC comparisons -------------------------------
boot_compare <- function(Ma, Mb, B=B_BOOT) {
  n <- nrow(Ma)  # aligned to all 56 participants
  obs <- icc_c1(Ma[,1],Ma[,2])[["icc"]] - icc_c1(Mb[,1],Mb[,2])[["icc"]]
  bd <- replicate(B, {
    ix <- sample.int(n,n,replace=TRUE)
    icc_c1(Ma[ix,1],Ma[ix,2])[["icc"]] - icc_c1(Mb[ix,1],Mb[ix,2])[["icc"]]
  })
  ci <- quantile(bd, c(.025,.975), na.rm=TRUE)
  list(diff=obs, lo=ci[[1]], hi=ci[[2]], distribution=bd)
}

TABLE2_PAIRS <- list(
  c("UNC","CSF"), c("CSF","TC"), c("TC","ATC"),
  c("UNC","ATC"), c("UNC","CR"), c("CR","ATC")
)

make_table2 <- function(ms) {
  set.seed(SEED)
  rows <- list(); z <- 1
  for (metname in names(METABOLITES)) {
    met <- METABOLITES[[metname]]
    for (pr in TABLE2_PAIRS) for (rg in REGIONS) {
      a <- pair_of(ms, paste(met,pr[1],rg,sep="_"))
      b <- pair_of(ms, paste(met,pr[2],rg,sep="_"))
      r <- boot_compare(a,b)
      rows[[z]] <- data.frame(metabolite=metname, region=rg,
                              comparison=paste(pr[1],pr[2],sep=" - "),
                              difference=r$diff, ci_low=r$lo, ci_high=r$hi,
                              significant=!(r$lo <= 0 && r$hi >= 0))
      z <- z + 1
    }
  }
  do.call(rbind, rows)
}

# ----------------------- Variance decomposition -------------------------------
variance_components <- function(v1,v2) {
  ok <- is.finite(v1) & is.finite(v2); v1 <- v1[ok]; v2 <- v2[ok]
  n <- length(v1)
  if (n < 3) return(c(cv_between=NA, cv_within=NA, icc=NA, n=n))
  X <- cbind(v1,v2); k <- 2; grand <- mean(X)
  rm <- rowMeans(X); cm <- colMeans(X)
  MSR <- k*sum((rm-grand)^2)/(n-1)
  SSE <- sum((X-grand)^2) - k*sum((rm-grand)^2) - n*sum((cm-grand)^2)
  MSE <- SSE/(n-1)
  varB <- max((MSR-MSE)/k,0); varW <- MSE
  c(cv_between=100*sqrt(varB)/grand,
    cv_within=100*sqrt(varW)/grand,
    icc=(MSR-MSE)/(MSR+MSE), n=n)
}

make_variance_components <- function(ms) {
  rows <- list(); z <- 1
  for (metname in names(METABOLITES)) {
    met <- METABOLITES[[metname]]
    for (me in MEASURES) {
      per_roi <- sapply(REGIONS, function(rg) {
        M <- pair_of(ms, paste(met,me,rg,sep="_"))
        variance_components(M[,1],M[,2])
      })
      rows[[z]] <- data.frame(
        metabolite=metname, measure=me,
        cv_between=mean(per_roi["cv_between",],na.rm=TRUE),
        cv_within=mean(per_roi["cv_within",],na.rm=TRUE),
        icc=mean(per_roi["icc",],na.rm=TRUE)
      )
      z <- z + 1
    }
  }
  do.call(rbind, rows)
}

# ----------------------- Tissue-explained variance -----------------------------
make_tissue_variance <- function(ms) {
  w1 <- ms[ms$session == 1,]
  rows <- list(); z <- 1
  for (metname in names(METABOLITES)) {
    met <- METABOLITES[[metname]]
    for (me in c("UNC","ATC")) for (rg in REGIONS) {
      ycol <- paste(met,me,rg,sep="_")
      d <- w1[,c(ycol,paste0("GM_FRA_",rg),paste0("WM_FRA_",rg),paste0("CSF_FRA_",rg))]
      names(d) <- c("y","GM","WM","CSF")
      d$y[is.finite(d$y) & d$y <= 0] <- NA_real_
      d <- d[complete.cases(d),]
      d$tGM <- d$GM/(d$GM+d$WM)
      # Fractions close to the closure constraint; R's lm drops an aliased term.
      r_full <- summary(lm(y ~ GM + WM + CSF, data=d))$r.squared
      r_csf  <- summary(lm(y ~ CSF, data=d))$r.squared
      r_tgm  <- summary(lm(y ~ tGM, data=d))$r.squared
      rows[[z]] <- data.frame(metabolite=metname, measure=me, region=rg,
                              r2_tissue=r_full, r2_csf=r_csf, r2_tgm=r_tgm)
      z <- z + 1
    }
  }
  do.call(rbind, rows)
}

make_residualized_icc <- function(ms) {
  rows <- list(); z <- 1
  for (metname in names(METABOLITES)) {
    met <- METABOLITES[[metname]]
    for (rg in REGIONS) {
      resid_wave <- function(wv) {
        d <- ms[ms$session == wv,
                c("participant_id",paste(met,"UNC",rg,sep="_"),
                  paste0("GM_FRA_",rg),paste0("WM_FRA_",rg),paste0("CSF_FRA_",rg))]
        names(d) <- c("s","y","GM","WM","CSF")
        d <- d[complete.cases(d),]
        data.frame(s=d$s, r=resid(lm(y ~ GM + WM + CSF, data=d)))
      }
      r1 <- resid_wave(1); r2 <- resid_wave(2)
      rr <- merge(r1,r2,by="s")
      icc_resid <- icc_c1(rr$r.x,rr$r.y)[["icc"]]
      U <- pair_of(ms,paste(met,"UNC",rg,sep="_"))
      A <- pair_of(ms,paste(met,"ATC",rg,sep="_"))
      rows[[z]] <- data.frame(metabolite=metname, region=rg,
                              icc_uncorrected=icc_c1(U[,1],U[,2])[["icc"]],
                              icc_tissue_residualized=icc_resid,
                              icc_alpha_corrected=icc_c1(A[,1],A[,2])[["icc"]])
      z <- z + 1
    }
  }
  do.call(rbind, rows)
}

# ----------------------- Longitudinal change ----------------------------------
make_longitudinal_change <- function(ms) {
  rows <- list(); z <- 1
  for (metname in names(METABOLITES)) {
    met <- METABOLITES[[metname]]
    for (me in c("UNC","CR")) for (rg in REGIONS) {
      M <- pair_of(ms,paste(met,me,rg,sep="_"))
      M <- M[complete.cases(M),,drop=FALSE]
      d <- M[,2]-M[,1]
      tt <- t.test(d)
      rows[[z]] <- data.frame(
        metabolite=metname, measure=me, region=rg, n=nrow(M),
        mean_change=mean(d), percent_change=100*mean(d)/mean(M[,1]),
        t=unname(tt$statistic), p_raw=tt$p.value
      )
      z <- z + 1
    }
  }
  out <- do.call(rbind,rows)
  out$p_bonferroni <- p.adjust(out$p_raw, method="bonferroni")
  out
}

make_tissue_change <- function(ms) {
  rows <- list(); z <- 1
  for (fr in c("GM_FRA","WM_FRA","CSF_FRA")) for (rg in REGIONS) {
    M <- pair_of(ms,paste(fr,rg,sep="_"))
    M <- M[complete.cases(M),,drop=FALSE]
    d <- M[,2]-M[,1]
    tt <- t.test(d)
    rows[[z]] <- data.frame(
      tissue=fr, region=rg, n=nrow(M),
      mean_change=d |> mean(), t=unname(tt$statistic), p_raw=tt$p.value
    )
    z <- z + 1
  }
  out <- do.call(rbind,rows)
  out$p_bonferroni <- p.adjust(out$p_raw,method="bonferroni")
  out
}

# ----------------------- Sensitivity analyses ---------------------------------
qc_sweep <- function(ms_raw,qc,thresholds=c(12,15,16,17,18,20)) {
  rows <- list(); z <- 1
  for (thr in thresholds) {
    m <- apply_qc(ms_raw,qc,linewidth_max=thr)
    cnt <- attr(m,"qc_counts")
    for (metname in names(METABOLITES)) {
      met <- METABOLITES[[metname]]
      vals <- sapply(c("UNC","CR","ATC"), function(me)
        mean(sapply(REGIONS,function(rg) {
          M <- pair_of(m,paste(met,me,rg,sep="_"))
          icc_c1(M[,1],M[,2])[["icc"]]
        }),na.rm=TRUE))
      rows[[z]] <- data.frame(
        linewidth_threshold=thr, metabolite=metname,
        excluded=if (met=="G") cnt[["excluded_GABA"]] else cnt[["excluded_Glx"]],
        ICC_UNC=vals[["UNC"]], ICC_CR=vals[["CR"]], ICC_ATC=vals[["ATC"]]
      )
      z <- z + 1
    }
  }
  do.call(rbind,rows)
}

short_interval_sensitivity <- function(ms, participants, cutoff=4.02) {
  ids <- participants$participant_id[participants$interval_years <= cutoff]
  rows <- list(); z <- 1
  for (metname in names(METABOLITES)) {
    met <- METABOLITES[[metname]]
    for (me in MEASURES) {
      vals <- sapply(REGIONS,function(rg) {
        M <- pair_of(ms,paste(met,me,rg,sep="_"),ids=sort(ids))
        icc_c1(M[,1],M[,2])[["icc"]]
      })
      rows[[z]] <- data.frame(metabolite=metname, measure=me,
                              n_participants=length(ids),
                              cutoff_years=cutoff,
                              mean_interval_years=mean(participants$interval_years[participants$participant_id %in% ids]),
                              cross_roi_mean_icc=mean(vals,na.rm=TRUE))
      z <- z + 1
    }
  }
  do.call(rbind,rows)
}

# ----------------------- Meta-analysis ----------------------------------------
# The manuscript reports random-effects meta-analysis (REML, Knapp-Hartung)
# and a multivariate sensitivity analysis using the bootstrap covariance among
# the six metabolite-by-region ICC differences.
meta_analysis_contrast <- function(ms, a, b, B=B_BOOT) {
  if (!requireNamespace("metafor", quietly=TRUE)) {
    warning("Package 'metafor' is not installed; meta-analysis outputs skipped.")
    return(NULL)
  }
  ids <- sort(unique(ms$participant_id))
  cells <- expand.grid(metabolite=names(METABOLITES), region=REGIONS,
                       stringsAsFactors=FALSE)
  yi <- numeric(nrow(cells))
  for (j in seq_len(nrow(cells))) {
    met <- METABOLITES[[cells$metabolite[j]]]
    Ma <- pair_of(ms,paste(met,a,cells$region[j],sep="_"),ids)
    Mb <- pair_of(ms,paste(met,b,cells$region[j],sep="_"),ids)
    yi[j] <- icc_c1(Ma[,1],Ma[,2])[["icc"]] - icc_c1(Mb[,1],Mb[,2])[["icc"]]
  }

  # Precompute the aligned matrices once; only row indices change in the bootstrap.
  Ma_list <- vector("list", nrow(cells))
  Mb_list <- vector("list", nrow(cells))
  for (j in seq_len(nrow(cells))) {
    met <- METABOLITES[[cells$metabolite[j]]]
    Ma_list[[j]] <- pair_of(ms,paste(met,a,cells$region[j],sep="_"),ids)
    Mb_list[[j]] <- pair_of(ms,paste(met,b,cells$region[j],sep="_"),ids)
  }

  set.seed(SEED)
  boot <- matrix(NA_real_, nrow=B, ncol=nrow(cells))
  for (i in seq_len(B)) {
    ix <- sample.int(length(ids),length(ids),replace=TRUE)
    for (j in seq_len(nrow(cells))) {
      Ma <- Ma_list[[j]]; Mb <- Mb_list[[j]]
      boot[i,j] <- icc_c1(Ma[ix,1],Ma[ix,2])[["icc"]] -
                   icc_c1(Mb[ix,1],Mb[ix,2])[["icc"]]
    }
  }
  V <- cov(boot,use="pairwise.complete.obs")
  vi <- diag(V)

  fit <- metafor::rma.uni(yi=yi, vi=vi, method="REML", test="knha")
  mod <- metafor::rma.uni(yi=yi, vi=vi, mods=~factor(cells$metabolite),
                          method="REML", test="knha")
  cells$cell <- paste(cells$metabolite,cells$region,sep="_")
  mv <- metafor::rma.mv(yi=yi, V=V, random=~1|cell, data=cells, method="REML", test="t")

  data.frame(
    comparison=paste(a,b,sep=" - "),
    model=c("univariate_REML_KH","multivariate_bootstrap_V"),
    estimate=c(as.numeric(coef(fit)[1]),as.numeric(coef(mv)[1])),
    ci_low=c(fit$ci.lb,mv$ci.lb), ci_high=c(fit$ci.ub,mv$ci.ub),
    p=c(fit$pval,mv$pval),
    metabolite_moderator_p=c(mod$QMp,NA_real_)
  )
}

# ----------------------- Run ---------------------------------------------------
d <- load_public_data()
ms_raw <- d$ms
ms <- apply_qc(ms_raw,d$qc,15)
qc_counts <- attr(ms,"qc_counts")
write.csv(data.frame(metric=names(qc_counts),value=as.numeric(qc_counts)),
          file.path(OUT_DIR,"qc_counts.csv"),row.names=FALSE)

table1 <- make_table1(ms)
write.csv(table1,file.path(OUT_DIR,"table1_icc.csv"),row.names=FALSE)

table2 <- make_table2(ms)
write.csv(table2,file.path(OUT_DIR,"table2_pairwise_icc.csv"),row.names=FALSE)

vc <- make_variance_components(ms)
write.csv(vc,file.path(OUT_DIR,"variance_components.csv"),row.names=FALSE)

tv <- make_tissue_variance(ms)
write.csv(tv,file.path(OUT_DIR,"tissue_explained_variance.csv"),row.names=FALSE)

ri <- make_residualized_icc(ms)
write.csv(ri,file.path(OUT_DIR,"tissue_residualized_icc.csv"),row.names=FALSE)

lc <- make_longitudinal_change(ms)
write.csv(lc,file.path(OUT_DIR,"longitudinal_change.csv"),row.names=FALSE)

tc <- make_tissue_change(ms)
write.csv(tc,file.path(OUT_DIR,"tissue_fraction_change.csv"),row.names=FALSE)

qs <- qc_sweep(ms_raw,d$qc)
write.csv(qs,file.path(OUT_DIR,"supplement_qc_sweep.csv"),row.names=FALSE)

si <- short_interval_sensitivity(ms,d$participants)
write.csv(si,file.path(OUT_DIR,"short_interval_sensitivity.csv"),row.names=FALSE)

# Meta-analysis is the only analysis requiring a non-base package.
meta_rows <- lapply(list(c("UNC","CSF"),c("UNC","ATC"),c("CSF","TC"),c("TC","ATC"),c("UNC","CR")),
                    function(x) meta_analysis_contrast(ms,x[1],x[2]))
meta_rows <- Filter(Negate(is.null),meta_rows)
if (length(meta_rows))
  write.csv(do.call(rbind,meta_rows),file.path(OUT_DIR,"meta_analysis.csv"),row.names=FALSE)

cat("Reproduction complete. Results written to outputs/.\n")
cat(sprintf("QC: %d spectra available; %d GABA+ and %d Glx excluded at 15 Hz + 8-MAD rule.\n",
            qc_counts[["available"]], qc_counts[["excluded_GABA"]], qc_counts[["excluded_Glx"]]))
