#!/usr/bin/env Rscript
# generate_fixtures.R
#
# Writes the JSON reference fixtures for the Layer 1 tests (SPEC.md §5): R's own
# distribution functions evaluated on fixed grids. The Vitest suites in tests/layer1/
# check that the TypeScript routines in src/stats/ reproduce these values.
#
# Usage, from the repository root:
#   Rscript reference/generate_fixtures.R
# Requires the jsonlite package.
#
# The fixtures in tests/fixtures/ are committed, so CI does not need R.
# Never edit a fixture by hand (SPEC.md §10, rule 2): change this script and rerun it.

suppressPackageStartupMessages(library(jsonlite))

# Any warning not deliberately caught below is a bug in this script: stop on it.
options(warn = 2)

out_dir <- file.path("tests", "fixtures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

to_json <- function(x) as.character(toJSON(x, digits = I(17), auto_unbox = TRUE))

# Writes one fixture file: a small header plus one JSON object per case and line.
# Numbers are written with 17 significant digits, so every double survives the
# round trip exactly; the script checks this before writing.
write_fixture <- function(name, description, cases) {
  stopifnot(is.data.frame(cases), nrow(cases) > 0)
  stopifnot(all(vapply(cases, function(col) all(is.finite(col)), logical(1))))
  rows <- vapply(seq_len(nrow(cases)), function(i) to_json(as.list(cases[i, , drop = FALSE])), "")
  json <- paste0(
    "{\n",
    '  "description": ', to_json(description), ",\n",
    '  "generated_by": "reference/generate_fixtures.R",\n',
    '  "r_version": ', to_json(R.version.string), ",\n",
    '  "n_cases": ', nrow(cases), ",\n",
    '  "cases": [\n    ', paste(rows, collapse = ",\n    "), "\n  ]\n",
    "}\n"
  )
  back <- fromJSON(json)$cases
  for (col in names(cases)) {
    stopifnot(identical(as.numeric(back[[col]]), as.numeric(cases[[col]])))
  }
  path <- file.path(out_dir, paste0(name, ".json"))
  writeLines(json, path, sep = "")
  message(sprintf("wrote %-34s %5d cases", path, nrow(cases)))
}

# Evaluates expr and returns 1 if it raised a warning, else 0. The warning is muffled.
warned <- function(expr) {
  flag <- 0
  withCallingHandlers(expr, warning = function(w) {
    flag <<- 1
    invokeRestart("muffleWarning")
  })
  flag
}

# ---------------------------------------------------------------------------
# Normal distribution: pnorm, qnorm
# ---------------------------------------------------------------------------

x <- sort(unique(c(
  seq(-8, 8, by = 0.125),
  c(-1, 1) * 1e-10, c(-1, 1) * 1e-5, c(-1, 1) * 1e-3,
  c(-1, 1) * 0.67448975, c(-1, 1) * 0.6745, # branch boundary of Cody's algorithm
  c(-1, 1) * sqrt(32), c(-1, 1) * 5.66, # branch boundary
  as.vector(outer(c(-1, 1), c(9, 10, 12, 15, 20, 25, 30, 35, 37)))
)))
write_fixture(
  "pnorm",
  "pnorm(x) and pnorm(x, lower.tail = FALSE)",
  data.frame(x = x, lower = pnorm(x), upper = pnorm(x, lower.tail = FALSE))
)

p_tail <- c(10^-(12:3), 2.5e-3, 5e-3, 0.01, 0.025, 0.05, 0.074, 0.075, 0.076,
            exp(-25) * c(0.999, 1.001)) # exp(-25): branch boundary r = 5 of AS 241
p <- sort(unique(c(p_tail, seq(0.1, 0.9, by = 0.05), 0.925, 1 - p_tail)))
write_fixture(
  "qnorm",
  "qnorm(p) and qnorm(p, lower.tail = FALSE), p in [1e-12, 1 - 1e-12]",
  data.frame(p = p, lower = qnorm(p), upper = qnorm(p, lower.tail = FALSE))
)

# ---------------------------------------------------------------------------
# Gamma and beta functions: lgamma, lbeta, pbeta
# ---------------------------------------------------------------------------

gx <- c(1e-8, 1e-3, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 1, 1.001, 1.01, 1.5,
        1.9, 1.99, 2, 2.01, 2.5, 3, 3.7, 5, 7.5, 10, 10.5, 14.5, 15, 15.5, 20, 33.3, 50,
        99.5, 100, 250.25, 500, 500.5, 1000, 2500.5, 5001, 1e4, 1e5, 1e6)
write_fixture("lgamma", "lgamma(x)", data.frame(x = gx, value = lgamma(gx)))

ab <- expand.grid(a = c(0.5, 1, 2.5, 7, 10, 15.5, 50, 500, 2500.5),
                  b = c(0.5, 1, 2.5, 7, 10, 15.5, 50, 500, 2500.5))
write_fixture("lbeta", "lbeta(a, b)", data.frame(a = ab$a, b = ab$b, value = lbeta(ab$a, ab$b)))

bg <- expand.grid(x = c(1e-6, 1e-3, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999),
                  a = c(0.5, 1, 2.5, 10, 50, 500),
                  b = c(0.5, 1, 2.5, 10, 50, 500))
write_fixture(
  "pbeta",
  "pbeta(x, a, b) and pbeta(x, a, b, lower.tail = FALSE)",
  data.frame(x = bg$x, a = bg$a, b = bg$b,
             lower = pbeta(bg$x, bg$a, bg$b),
             upper = pbeta(bg$x, bg$a, bg$b, lower.tail = FALSE))
)

# ---------------------------------------------------------------------------
# Student t distribution: pt, dt, qt
# ---------------------------------------------------------------------------

# The df values required by SPEC.md §5, plus non-integer df as produced by Welch's test.
t_df <- c(1, 2, 5, 10, 30, 100, 1000, 1.5, 7.3, 47.8)

tg <- expand.grid(
  q = c(-100, -40, -20, -10, -6, -5, -4, -3, -2.5, -2, -1.5, -1, -0.5, -0.1, 0,
        0.1, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 10, 20, 40, 100),
  df = t_df
)
write_fixture(
  "pt",
  "pt(q, df), pt(q, df, lower.tail = FALSE) and dt(q, df)",
  data.frame(q = tg$q, df = tg$df,
             lower = pt(tg$q, tg$df),
             upper = pt(tg$q, tg$df, lower.tail = FALSE),
             density = dt(tg$q, tg$df))
)

p_t <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.45)
qg <- expand.grid(p = sort(c(p_t, 0.5, 1 - p_t)), df = t_df)
write_fixture(
  "qt",
  "qt(p, df) and qt(p, df, lower.tail = FALSE), p in [1e-8, 1 - 1e-8]",
  data.frame(p = qg$p, df = qg$df,
             lower = qt(qg$p, qg$df),
             upper = qt(qg$p, qg$df, lower.tail = FALSE))
)

# ---------------------------------------------------------------------------
# Noncentral t distribution: pt(q, df, ncp)
# ---------------------------------------------------------------------------

# R warns "full precision may not have been achieved in 'pnt{final}'" when the lower
# tail it computes is within 1e-10 of 1, i.e. the complementary tail has lost
# *relative* precision. The absolute accuracy, which is what Layer 1 tests (abs 1e-8),
# is unaffected, so those cases are kept and flagged with r_warning = 1.
ng <- expand.grid(
  q = c(-10, -5, -3, -2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 10, 15, 20, 40),
  ncp = c(0, 0.5, 1, 2, 5, 10, -1, -5),
  df = c(1, 2, 5, 10, 30, 100, 1000)
)
nct <- data.frame(q = ng$q, df = ng$df, ncp = ng$ncp)
nct$lower <- NA_real_
nct$upper <- NA_real_
nct$r_warning <- 0
for (i in seq_len(nrow(nct))) {
  w1 <- warned(nct$lower[i] <- pt(nct$q[i], nct$df[i], nct$ncp[i]))
  w2 <- warned(nct$upper[i] <- pt(nct$q[i], nct$df[i], nct$ncp[i], lower.tail = FALSE))
  nct$r_warning[i] <- max(w1, w2)
}
write_fixture(
  "pnt",
  "pt(q, df, ncp) and pt(q, df, ncp, lower.tail = FALSE); r_warning = 1 if R warned",
  nct
)

# ---------------------------------------------------------------------------
# Discrete distributions: dbinom(log = TRUE), dhyper(log = TRUE)
# ---------------------------------------------------------------------------

binom <- do.call(rbind, lapply(c(1, 2, 5, 10, 20, 50, 100, 500, 1000, 5000), function(n) {
  do.call(rbind, lapply(c(0.001, 0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98, 0.999), function(p) {
    sd <- sqrt(n * p * (1 - p))
    xs <- round(c(0, 1, 2, n * p + c(-3, -1, 0, 1, 3) * sd, n / 3, n / 2, n - 2, n - 1, n))
    xs <- sort(unique(xs[xs >= 0 & xs <= n]))
    data.frame(x = xs, n = n, p = p)
  }))
}))
binom$value <- dbinom(binom$x, binom$n, binom$p, log = TRUE)
write_fixture("dbinom_log", "dbinom(x, n, p, log = TRUE), n up to 5000", binom)

# Urn with m white and n black balls, k draws; x = white balls drawn. Only supports
# with at least two points (a one-point support has log-pmf exactly 0).
hyper <- do.call(rbind, lapply(c(1, 5, 20, 100, 1000, 5000), function(m) {
  do.call(rbind, lapply(c(1, 7, 20, 150, 1000, 5000), function(n) {
    N <- m + n
    ks <- round(c(1, 2, N * 0.05, N * 0.3, N * 0.5, N * 0.9, N - 1))
    ks <- sort(unique(ks[ks >= 1 & ks <= N - 1]))
    do.call(rbind, lapply(ks, function(k) {
      lo <- max(0, k - n)
      hi <- min(k, m)
      if (lo == hi) return(NULL)
      mode <- floor((k + 1) * (m + 1) / (N + 2))
      xs <- round(c(lo, lo + 1, (lo + hi) / 2, mode, hi - 1, hi))
      xs <- sort(unique(xs[xs >= lo & xs <= hi]))
      data.frame(x = xs, m = m, n = n, k = k)
    }))
  }))
}))
hyper$value <- dhyper(hyper$x, hyper$m, hyper$n, hyper$k, log = TRUE)
write_fixture("dhyper_log", "dhyper(x, m, n, k, log = TRUE), m and n up to 5000", hyper)
