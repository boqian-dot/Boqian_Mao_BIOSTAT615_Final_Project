/**
 * Log probability mass functions of the binomial and hypergeometric distributions.
 *
 * C. Loader (2000), "Fast and accurate computation of binomial probabilities"
 * (manuscript; the method behind R's dbinom()). The binomial log-pmf is written as
 *
 *   log f(x; n, p) = stirlerr(n) − stirlerr(x) − stirlerr(n − x)
 *                    − bd0(x, np) − bd0(n − x, nq) − ½ log(2π x (n − x) / n),
 *
 * where bd0(x, M) = x log(x/M) + M − x is evaluated by a series when x ≈ M. Unlike
 * log C(n, x) + x log p + (n − x) log q, this has no cancellation between large
 * terms, so it stays accurate for large n, and working on the log scale means it
 * never underflows.
 *
 * The hypergeometric pmf is a ratio of three binomial pmfs with a common p (Loader 2000):
 *   P(X = x) = f(x; m, p) · f(k − x; n, p) / f(k; m + n, p),  p = k / (m + n).
 */

import { DBL_MIN, LN_2PI } from './constants'
import { stirlerr } from './gamma'

/** Deviance term bd0(x, M) = x log(x/M) + M − x, accurate also when x ≈ M. */
export function bd0(x: number, M: number): number {
  if (Math.abs(x - M) < 0.1 * (x + M)) {
    // With v = (x − M)/(x + M): x log(x/M) = 2x(v + v³/3 + v⁵/5 + …), and M − x = −v(x + M).
    let v = (x - M) / (x + M)
    let s = (x - M) * v
    if (Math.abs(s) < DBL_MIN) return s
    let ej = 2 * x * v
    v *= v
    for (let j = 1; j < 1000; j++) {
      ej *= v
      const s1 = s + ej / (2 * j + 1)
      if (s1 === s) return s1
      s = s1
    }
  }
  return x * Math.log(x / M) + M - x
}

/**
 * Binomial log-pmf with success probability p and q = 1 − p passed separately
 * (so a caller that knows q exactly can avoid the rounding of 1 − p).
 * Assumes integer 0 ≤ x ≤ n is checked by the caller, or returns −∞ outside.
 */
export function dbinomRawLog(x: number, n: number, p: number, q: number): number {
  if (p === 0) return x === 0 ? 0 : -Infinity
  if (q === 0) return x === n ? 0 : -Infinity
  if (x < 0 || x > n) return -Infinity
  if (x === 0) return n === 0 ? 0 : p < 0.1 ? n * Math.log1p(-p) : n * Math.log(q)
  if (x === n) return q < 0.1 ? n * Math.log1p(-q) : n * Math.log(p)
  const lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q)
  const lf = LN_2PI + Math.log(x) + Math.log1p(-x / n)
  return lc - 0.5 * lf
}

/** log P(X = x) for X ~ Binomial(n, p). R: dbinom(x, n, p, log = TRUE). */
export function dbinomLog(x: number, n: number, p: number): number {
  if (!Number.isInteger(n) || n < 0 || Number.isNaN(p) || p < 0 || p > 1 || Number.isNaN(x)) return NaN
  if (!Number.isInteger(x) || x < 0 || x > n) return -Infinity
  return dbinomRawLog(x, n, p, 1 - p)
}

/**
 * log P(X = x) for X ~ Hypergeometric: the number of white balls in k draws without
 * replacement from an urn with m white and n black balls. R: dhyper(x, m, n, k, log = TRUE).
 */
export function dhyperLog(x: number, m: number, n: number, k: number): number {
  if (![m, n, k].every((v) => Number.isInteger(v) && v >= 0) || k > m + n || Number.isNaN(x)) return NaN
  if (!Number.isInteger(x) || x < 0 || x > k || x > m || k - x > n) return -Infinity
  if (k === 0) return 0
  const total = m + n
  const p = k / total
  const q = (total - k) / total
  return dbinomRawLog(x, m, p, q) + dbinomRawLog(k - x, n, p, q) - dbinomRawLog(k, total, p, q)
}
