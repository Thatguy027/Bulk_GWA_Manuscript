## Shared scoring core for the panel optimiser -------------------------------
##
## Sourced by DIAG_pool_optimizer.R and DIAG_pool_optimizer_size_sweep.R so the
## two cannot drift apart. The caller sets the tuning constants (KSTAR, MAC_MIN,
## R2_MAX, Q) and calls pool_opt_attach() once, which binds the large objects
## into the global environment rather than copying them per call.
##
## Names are deliberately unabbreviated and upper case. The obvious short ones
## collide with base: `cm` is grDevices::cm, whose binding is locked, so a
## superassignment to it fails at run time rather than at parse time.
## ---------------------------------------------------------------------------

pool_opt_attach <- function(P) {
  ## Carrier pairs sorted by strain, with per-strain offsets. Selecting the
  ## pairs belonging to a panel is then |S| contiguous reads instead of a
  ## logical pass over all 10.9 M rows, which is the whole cost of the loop.
  data.table::setorder(P$car, s)
  e <- globalenv()
  assign("N_UNIV",  length(P$strains), e)
  assign("N_RARE",  P$n_rare,          e)
  assign("CAR_M",   P$car$m,           e)
  assign("CAR_S",   P$car$s,           e)
  assign("CAR_DIV", P$car$in_div,      e)
  off <- c(0L, cumsum(tabulate(P$car$s, nbins = length(P$strains))))
  assign("CAR_START", off[-length(off)] + 1L, e)
  assign("CAR_LEN",   diff(off),              e)
  assign("GENO", P$gp, e)
  assign("GRAM", P$G,  e)
  assign("N_MARK", ncol(P$gp), e)
  invisible(TRUE)
}

## Privateness: one pass over the carrier list. A marker is private to the panel
## when exactly one of its carriers is in it. Returns counts for EVERY strain in
## the universe; index by S to get the panel's own.
priv_counts <- function(S) {
  k  <- sequence(CAR_LEN[S], CAR_START[S])
  ms <- CAR_M[k]; ss <- CAR_S[k]; ds <- CAR_DIV[k]
  cnt <- tabulate(ms, nbins = N_RARE)
  p <- cnt[ms] == 1L
  list(nondiv = tabulate(ss[p & !ds], nbins = N_UNIV),
       div    = tabulate(ss[p &  ds], nbins = N_UNIV))
}

## Mapping: the panel kinship is a double-centred submatrix of the one
## precomputed cross-product. Centring each marker within the panel is exactly
## post-multiplication by H = I - J/n, so H GRAM[S,S] H is the marker-centred
## kinship and no marker is ever re-centred. Its top eigenvectors are already
## orthogonal to the intercept, so markers project onto them uncentred.
map_score <- function(S) {
  n  <- length(S)
  Gs <- GRAM[S, S, drop = FALSE]
  rm_ <- rowMeans(Gs)
  K  <- Gs - rm_ - rep(rm_, each = n) + mean(Gs)
  ev <- eigen(K, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  Qm <- ev$vectors[, seq_len(Q), drop = FALSE]
  A  <- GENO[S, , drop = FALSE]
  s1 <- colSums(A); s2 <- colSums(A * A)
  tot  <- s2 - s1 * s1 / n                    # centred sum of squares
  expl <- colSums(crossprod(Qm, A)^2)
  r2  <- ifelse(tot > 1e-9, expl / tot, 1)
  mac <- pmin(s1, n - s1)
  list(n_mappable = sum(mac >= MAC_MIN & r2 <= R2_MAX),
       n_mac      = sum(mac >= MAC_MIN),
       pc1_share  = lam[1] / sum(lam),
       eff_dim    = sum(lam)^2 / sum(lam^2))
}

score <- function(S) {
  pc <- priv_counts(S); ms <- map_score(S)
  nd <- pc$nondiv[S]
  c(ms, list(min_nondiv = min(nd), med_nondiv = median(nd),
             n_below = sum(nd < KSTAR), med_div = median(pc$div[S])))
}
feasible <- function(S) min(priv_counts(S)$nondiv[S]) >= KSTAR

## Swap-based simulated annealing. Proposals that break the identifiability
## floor are rejected outright rather than penalised, so the search never leaves
## the feasible set once inside it.
pool_opt_anneal <- function(S, iters, T0 = 40, report = NULL) {
  cur <- score(S); best <- cur; bestS <- S; trace <- numeric(iters)
  n <- length(S)
  for (it in seq_len(iters)) {
    Temp <- T0 * (1 - it / iters) + 1e-6
    out <- S[sample.int(n, 1)]
    ins <- sample(setdiff(seq_len(N_UNIV), S), 1)
    T_ <- c(setdiff(S, out), ins)
    st <- score(T_)
    if (st$min_nondiv >= KSTAR &&
        (st$n_mappable > cur$n_mappable ||
         runif(1) < exp((st$n_mappable - cur$n_mappable) / Temp))) {
      S <- T_; cur <- st
      if (cur$n_mappable > best$n_mappable) { best <- cur; bestS <- S }
    }
    trace[it] <- best$n_mappable
    if (!is.null(report) && it %% report == 0)
      cat(format(Sys.time(), "[%H:%M:%S] "),
          sprintf("  iter %5d  best %d mappable\n", it, best$n_mappable), sep = "")
  }
  list(S = bestS, best = best, trace = trace)
}

## Farthest-point seed: a maximally dispersed panel, which is the right starting
## point because dispersion is what both halves of the problem want.
pool_opt_seed <- function(n, fixed = integer(0)) {
  d2 <- outer(diag(GRAM), diag(GRAM), "+") - 2 * GRAM
  S <- integer(n)
  if (length(fixed)) {
    S[seq_along(fixed)] <- fixed
    dmin <- apply(d2[fixed, , drop = FALSE], 2, min)
    start <- length(fixed) + 1L
  } else {
    S[1] <- which.max(rowSums(d2)); dmin <- d2[S[1], ]; start <- 2L
  }
  if (start <= n) for (i in start:n) {
    dmin[S[S > 0]] <- -Inf
    S[i] <- which.max(dmin); dmin <- pmin(dmin, d2[S[i], ])
  }
  S
}

## Repair an infeasible panel: swap out the strain that most violates the floor
## for the candidate that best relieves it. Privateness is NOT monotone in the
## panel -- adding a strain can destroy another's -- so this is re-checked each
## pass rather than computed once.
pool_opt_repair <- function(S, uni_nondiv, passes = 40L, ncand = 25L) {
  for (pass in seq_len(passes)) {
    pc <- priv_counts(S)$nondiv
    bad <- S[pc[S] < KSTAR]
    if (!length(bad)) break
    out <- bad[which.min(pc[bad])]
    cand <- setdiff(seq_len(N_UNIV), S)
    cand <- cand[order(uni_nondiv[cand], decreasing = TRUE)][seq_len(min(ncand, length(cand)))]
    bestmin <- -1; bestS <- NULL
    for (cc in cand) {
      T_ <- c(setdiff(S, out), cc)
      mn <- min(priv_counts(T_)$nondiv[T_])
      if (mn > bestmin) { bestmin <- mn; bestS <- T_ }
    }
    S <- bestS
    if (bestmin >= KSTAR) break
  }
  S
}
