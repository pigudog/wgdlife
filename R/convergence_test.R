## ============================================================================================
## Does cancer type explain LESS of the somatic event spectrum AFTER WGD than BEFORE it?
## Paired within-sample design on the WGD+ timing layer.
##   estimator 1  R2  = SS_between(cancer type) / SS_total on rarefied relative-abundance
##                      spectra.  With d = Euclidean this is the PERMANOVA R2 (Anderson 2001),
##                      computed in O(N*P) rather than O(N^2).
##   estimator 2  NMI = event-level I(cancer type ; event) / H(cancer type)
##   estimator 3  leave-one-out nearest-centroid macro recall (classification-based)
## Inference: (a) bootstrap over samples WITHIN cancer type (composition + n_g fixed) for a CI
##            on delta = post - pre; (b) within-sample permutation of the pre/post label as the
##            null (preserves each sample's event multiset AND its n_pre / n_post exactly).
## ============================================================================================
suppressPackageStartupMessages({library(data.table); library(parallel)})
RNGkind("L'Ecuyer-CMRG"); set.seed(20260924)   # fixed seed of the published runs
NCORE <- 8L
D3   <- file.path(Sys.getenv("WGDLIFE_DATA", "data"), "")   # directory of the two timed-event tables (see R/README.md)
LONG <- fread(paste0(D3, "fig3_convergence_events_long.tsv"))
UNI  <- fread(paste0(D3, "fig3_convergence_event_universe.tsv"))

r2_stats <- function(X, g, ng, G) {
  N <- nrow(X); mu <- colMeans(X)
  gmu <- rowsum(X, g, reorder = TRUE) / as.numeric(ng)
  sst <- sum(X * X) - N * sum(mu * mu)
  ssb <- sum(as.numeric(ng) * rowSums(sweep(gmu, 2, mu)^2))
  c(r2 = ssb / sst, r2adj = 1 - ((sst - ssb) / (N - G)) / (sst / (N - 1)))
}
nmi_stats <- function(M, g) {
  C <- rowsum(M, g, reorder = TRUE); P <- C / sum(C)
  pi_ <- rowSums(P); pj <- colSums(P); E <- outer(pi_, pj); nz <- P > 0
  I <- sum(P[nz] * log(P[nz] / E[nz])); Ht <- -sum(pi_[pi_ > 0] * log(pi_[pi_ > 0]))
  c(nmi = I / Ht)
}
lonc <- function(X, g, ng) {
  gs <- rowsum(X, g, reorder = TRUE); G <- nrow(gs); N <- nrow(X); Pn <- ncol(X)
  gmu <- gs / as.numeric(ng); x2 <- rowSums(X * X); D <- matrix(0, N, G)
  for (k in seq_len(G)) {
    own <- g == rownames(gs)[k]; ck <- matrix(gmu[k, ], N, Pn, byrow = TRUE)
    ck[own, ] <- (matrix(gs[k, ], sum(own), Pn, byrow = TRUE) - X[own, , drop = FALSE]) / (ng[k] - 1)
    D[, k] <- x2 - 2 * rowSums(X * ck) + rowSums(ck * ck)
  }
  rec <- tapply(rownames(gs)[max.col(-D, ties.method = "first")] == g, g, mean)
  list(bacc = mean(rec), rec = rec)
}

## K: integer depth, "min" = per-sample min(n_pre, n_post), "all" = no rarefaction
one_draw <- function(dt, K, nev, N, permute_tier = FALSE) {
  d <- copy(dt)
  if (permute_tier) d[, tier := sample(tier), by = sid]
  if (!identical(K, "all")) {
    d[, r := runif(.N)]; setorder(d, sid, tier, r); d[, k := seq_len(.N), by = .(sid, tier)]
    if (identical(K, "min")) { d[, kk := min(.N), by = .(sid, tier)]; d[, kk := min(kk), by = sid]
                               d <- d[k <= kk] } else d <- d[k <= K]
  }
  lapply(c(pre = "pre", post = "post"), function(tt) {
    s <- d[tier == tt]; M <- matrix(0, N, nev); M[cbind(s$sid, s$eid)] <- 1; M })
}
run_stats <- function(dr, meta, idx = NULL, bacc = FALSE) {
  g <- meta$g; ng <- meta$ng; G <- meta$G
  if (!is.null(idx)) g <- g[idx]
  out <- lapply(dr, function(M) {
    if (!is.null(idx)) M <- M[idx, , drop = FALSE]
    rs <- rowSums(M); X <- M / pmax(rs, 1)
    v <- c(r2_stats(X, g, ng, G), nmi_stats(M, g))
    if (bacc) v <- c(v, bacc = lonc(X, g, ng)$bacc)
    v })
  c(pre = out$pre, post = out$post, setNames(out$post - out$pre, paste0("d_", names(out$pre))))
}

prep <- function(platform_, currency_, minev, K, kinds = NULL, min_type = 100L, min_each = 3L) {
  L <- LONG[platform == platform_ & currency == currency_]
  if (!is.null(kinds)) L <- L[kind %in% kinds]
  keep <- UNI[platform == platform_ & currency == currency_ & pre >= minev & post >= minev, event]
  L <- L[event %in% keep]
  amb <- L[, uniqueN(cancer_type), by = sample][V1 > 1, sample]   # patient-collapsed: drop
  L <- L[!sample %in% amb]                                        # the few ambiguous-type patients
  cnt <- dcast(L, sample ~ tier, value.var = "event", fun.aggregate = length)
  L <- L[sample %in% cnt[pre >= min_each & post >= min_each, sample]]
  nt <- L[, uniqueN(sample), by = cancer_type]
  L <- L[cancer_type %in% nt[V1 >= min_type, cancer_type]]
  L[, sid := as.integer(factor(sample))]
  ev <- sort(unique(L$event)); L[, eid := match(event, ev)]
  sm <- unique(L[, .(sid, cancer_type)])[order(sid)]
  g <- sm$cancer_type; ng <- as.integer(table(g))
  stopifnot(nrow(sm) == uniqueN(L$sample))
  list(dt = L[, .(sid, eid, tier)], ev = ev, n_amb_dropped = length(amb),
       meta = list(g = g, ng = ng, G = length(ng), N = nrow(sm)),
       ntab = L[, uniqueN(sample), by = cancer_type][order(-V1)])
}

analyse <- function(tag, platform_, currency_ = "strict", minev = 20L, K = 3L, kinds = NULL,
                    min_each = 3L, n_point = 200L, n_boot = 1000L, n_null = 1000L) {
  P <- prep(platform_, currency_, minev, K, kinds, min_each = min_each)
  dt <- P$dt; meta <- P$meta; nev <- length(P$ev); N <- meta$N
  Kl <- if (is.character(K)) K else as.integer(K)
  cat(sprintf("[%s] N=%d P=%d G=%d rows=%d (dropped %d ambiguous-type patients)\n",
              tag, N, nev, meta$G, nrow(dt), P$n_amb_dropped))
  if (N < 200 || meta$G < 3) { cat("  SKIPPED: underpowered\n"); return(NULL) }
  pt <- simplify2array(mclapply(seq_len(n_point), function(i)
      run_stats(one_draw(dt, Kl, nev, N), meta, bacc = TRUE), mc.cores = NCORE))
  point <- rowMeans(pt)
  ibg <- split(seq_len(N), meta$g)
  bt <- simplify2array(mclapply(seq_len(n_boot), function(i) {
      idx <- unlist(lapply(ibg, function(v) sample(v, length(v), replace = TRUE)), use.names = FALSE)
      run_stats(one_draw(dt, Kl, nev, N), meta, idx = idx) }, mc.cores = NCORE))
  nl <- simplify2array(mclapply(seq_len(n_null), function(i)
      run_stats(one_draw(dt, Kl, nev, N, permute_tier = TRUE), meta), mc.cores = NCORE))
  mk <- function(st, b, nn) data.table(tag = tag, platform = platform_, currency = currency_,
      minev = minev, K = as.character(K), kinds = if (is.null(kinds)) "all" else paste(kinds, collapse = "+"),
      N = N, n_events_universe = nev, G = meta$G, stat = st,
      pre = point[[paste0("pre.", st)]], post = point[[paste0("post.", st)]],
      delta = point[[paste0("d_", st)]],
      rel_drop = 1 - point[[paste0("post.", st)]] / point[[paste0("pre.", st)]],
      boot_lo = if (is.null(b)) NA_real_ else unname(quantile(b, .025)),
      boot_hi = if (is.null(b)) NA_real_ else unname(quantile(b, .975)),
      boot_frac_delta_ge0 = if (is.null(b)) NA_real_ else mean(b >= 0),
      null_mean = if (is.null(nn)) NA_real_ else mean(nn),
      null_sd = if (is.null(nn)) NA_real_ else sd(nn),
      null_p_two = if (is.null(nn)) NA_real_ else
        (1 + sum(abs(nn - mean(nn)) >= abs(point[[paste0("d_", st)]] - mean(nn)))) / (1 + length(nn)),
      z_vs_null = if (is.null(nn)) NA_real_ else (point[[paste0("d_", st)]] - mean(nn)) / sd(nn))
  res <- rbindlist(c(lapply(c("r2", "r2adj", "nmi"), function(st) mk(st, bt[paste0("d_", st), ],
                                                                    nl[paste0("d_", st), ])),
                     list(mk("bacc", NULL, NULL))))
  list(res = res, prep = P, point = point, boot = bt, null = nl)
}
