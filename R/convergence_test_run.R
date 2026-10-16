args <- commandArgs(trailingOnly = FALSE)
here <- sub("^--file=", "", args[grepl("^--file=", args)])
source(file.path(if (length(here)) dirname(here[1]) else "R", "convergence_test.R"))
OUT <- D3
runs <- list()
add <- function(...) { a <- analyse(...); if (!is.null(a)) runs[[length(runs) + 1]] <<- a; invisible(a) }

## ---------------------------------------------------------------- MSK: PRIMARY
PRI <- add("msk_primary", "msk", "strict", 20L, 3L, n_point = 200L, n_boot = 2000L, n_null = 2000L)

## ---------------------------------------------------------------- MSK: confounder controls
NB <- 500L
add("msk_K5",           "msk", "strict", 20L, 5L,  min_each = 5L, n_point = 100L, n_boot = NB, n_null = NB)
add("msk_Kmin_paired",  "msk", "strict", 20L, "min",              n_point = 100L, n_boot = NB, n_null = NB)
add("msk_Kall_noraref", "msk", "strict", 20L, "all",              n_point = 20L,  n_boot = NB, n_null = NB)
add("msk_univ50",       "msk", "strict", 50L, 3L,                 n_point = 100L, n_boot = NB, n_null = NB)
add("msk_univ10",       "msk", "strict", 10L, 3L,                 n_point = 100L, n_boot = NB, n_null = NB)
add("msk_currency_all", "msk", "all",    20L, 3L,                 n_point = 100L, n_boot = NB, n_null = NB)
add("msk_armloss_only", "msk", "strict", 20L, 3L, kinds = "arm_loss", n_point = 100L, n_boot = NB, n_null = NB)

## ---------------------------------------------------------------- TCGA replication
add("tcga_primary",     "tcga", "strict", 20L, 3L, n_point = 200L, n_boot = 2000L, n_null = 2000L)
add("tcga_Kall_noraref","tcga", "strict", 20L, "all",              n_point = 20L,  n_boot = NB, n_null = NB)
add("tcga_armloss_only","tcga", "strict", 20L, 3L, kinds = "arm_loss", n_point = 100L, n_boot = NB, n_null = NB)

main <- rbindlist(lapply(runs, `[[`, "res"))
fwrite(main, paste0(OUT, "fig3_convergence_main.tsv"), sep = "\t")
print(main[stat == "r2", .(tag, N, P = n_events_universe, G, pre = round(pre, 5), post = round(post, 5),
                           delta = round(delta, 5), rel_drop = round(rel_drop, 3),
                           lo = round(boot_lo, 5), hi = round(boot_hi, 5),
                           null_p = null_p_two, z = round(z_vs_null, 1))])

## ---------------------------------------------------------------- primary: per-type + per-event
P <- PRI$prep; meta <- P$meta; nev <- length(P$ev); N <- meta$N
set.seed(99)
acc_rec <- list(); acc_disp <- list(); acc_ev <- list()
for (i in 1:200) {
  dr <- one_draw(P$dt, 3L, nev, N)
  for (tt in c("pre", "post")) {
    M <- dr[[tt]]; X <- M / pmax(rowSums(M), 1)
    lc <- lonc(X, meta$g, meta$ng)
    acc_rec[[length(acc_rec) + 1]] <- data.table(tier = tt, cancer_type = names(lc$rec), recall = as.numeric(lc$rec))
    gmu <- rowsum(X, meta$g, reorder = TRUE) / as.numeric(meta$ng)
    mu  <- colMeans(X)
    d_own <- sqrt(rowSums((X - gmu[match(meta$g, rownames(gmu)), , drop = FALSE])^2))
    acc_disp[[length(acc_disp) + 1]] <- data.table(tier = tt, cancer_type = meta$g, d_within = d_own)[
      , .(within_disp = mean(d_within)), by = .(tier, cancer_type)][
      , between_disp := sqrt(rowSums((gmu[match(cancer_type, rownames(gmu)), , drop = FALSE] -
                                      matrix(mu, .N, length(mu), byrow = TRUE))^2))][]
    ## per-event: chi-square contribution of the type x event table
    C <- rowsum(M, meta$g, reorder = TRUE); E <- outer(rowSums(C), colSums(C)) / sum(C)
    acc_ev[[length(acc_ev) + 1]] <- data.table(tier = tt, event = P$ev,
                                               chi2 = colSums((C - E)^2 / pmax(E, 1e-9)), n = colSums(C))
  }
}
byt <- rbindlist(acc_rec)[, .(recall = mean(recall)), by = .(tier, cancer_type)]
byd <- rbindlist(acc_disp)[, .(within_disp = mean(within_disp), between_disp = mean(between_disp)),
                           by = .(tier, cancer_type)]
byt <- merge(byt, byd, by = c("tier", "cancer_type"))
byt <- merge(byt, P$ntab[, .(cancer_type, n_samples = V1)], by = "cancer_type")
byt[, ratio_between_within := between_disp / within_disp]
bywide <- dcast(byt, cancer_type + n_samples ~ tier,
                value.var = c("recall", "within_disp", "between_disp", "ratio_between_within"))
fwrite(bywide[order(-n_samples)], paste0(OUT, "fig3_convergence_bytype.tsv"), sep = "\t")
print(bywide[order(-n_samples)])

bye <- rbindlist(acc_ev)[, .(chi2 = mean(chi2), n = mean(n)), by = .(tier, event)]
bye <- dcast(bye, event ~ tier, value.var = c("chi2", "n"))
fwrite(bye[order(-chi2_pre)], paste0(OUT, "fig3_convergence_byevent.tsv"), sep = "\t")

## null / bootstrap draws of the primary, for the figure
fwrite(rbind(data.table(kind = "bootstrap", d_r2 = PRI$boot["d_r2", ], d_nmi = PRI$boot["d_nmi", ]),
             data.table(kind = "null_within_sample_perm", d_r2 = PRI$null["d_r2", ], d_nmi = PRI$null["d_nmi", ])),
       paste0(OUT, "fig3_convergence_primary_draws.tsv"), sep = "\t")
cat("DONE\n")
