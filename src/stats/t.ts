/**
 * Student t distribution: CDF, density, and quantile function.
 *
 * pt — through the incomplete beta function: P(|T| > |t|) = I_{df/(df+t²)}(df/2, 1/2)
 *   (Abramowitz & Stegun 1964, eq. 26.7.1 and 26.5.27). For df > 4e5 the normal
 *   approximation of A&S eq. 26.7.8 is used, as in R.
 *
 * qt — starting value from G. W. Hill (1970), "Algorithm 396: Student's
 *   t-quantiles", Communications of the ACM 13, 619–620, then refined by Newton
 *   steps with a second-order (Taylor) correction until the relative change is
 *   below 1e-14. The closed forms for df = 1 (Cauchy) and df = 2 are used directly.
 */

import { pbetaBoth } from './beta'
import { lbeta } from './gamma'
import { pnorm, qnorm } from './normal'

/** Density of the t distribution with df degrees of freedom. R: dt(). */
export function dt(x: number, df: number): number {
  if (Number.isNaN(x) || Number.isNaN(df) || df <= 0) return NaN
  if (!Number.isFinite(x)) return 0
  // f(x) = (1 + x²/df)^(−(df+1)/2) / (√df · B(1/2, df/2))
  return Math.exp(-((df + 1) / 2) * Math.log1p((x * x) / df) - 0.5 * Math.log(df) - lbeta(0.5, df / 2))
}

/** t distribution CDF, P(T ≤ q), or P(T > q) when lowerTail is false. R: pt(). */
export function pt(q: number, df: number, lowerTail = true): number {
  if (Number.isNaN(q) || Number.isNaN(df) || df <= 0) return NaN
  if (!Number.isFinite(q)) return q < 0 === lowerTail ? 0 : 1
  if (!Number.isFinite(df)) return pnorm(q, lowerTail)
  if (df > 4e5) {
    const s = 1 / (4 * df)
    return pnorm((q * (1 - s)) / Math.sqrt(1 + q * q * 2 * s), lowerTail)
  }
  const q2 = q * q
  // Two-sided tail P(|T| > |q|); both arguments of the beta are formed without cancellation.
  const [twoSided] = pbetaBoth(df / (df + q2), q2 / (df + q2), df / 2, 0.5)
  const farTail = twoSided / 2 // P(T > |q|)
  // The requested tail is the small one when it lies beyond q (q > 0 upper, q < 0 lower).
  const requestedIsFar = q > 0 !== lowerTail
  return requestedIsFar ? farTail : 1 - farTail
}

/** Hill's (1970) approximation to the t > 0 with P(T > t) = alpha, for df ≥ 1. */
function hillStart(alpha: number, df: number): number {
  const P = 2 * alpha // two-sided probability
  const a = 1 / (df - 0.5)
  const b = 48 / (a * a)
  let c = ((20700 * a) / b - 98) * a - 16
  c = c * a + 96.36
  const d = ((94.5 / (b + c) - 3) / b + 1) * Math.sqrt((a * Math.PI) / 2) * df
  let y = Math.pow(d * P, 2 / df)
  if (y > 0.05 + a) {
    // Asymptotic inverse expansion about the normal quantile.
    const x = qnorm(0.5 * P)
    y = x * x
    if (df < 5) c += 0.3 * (df - 4.5) * (x + 0.6)
    c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c
    y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x
    y = Math.expm1(a * y * y)
  } else {
    y =
      ((1 / (((df + 6) / (df * y) - 0.089 * d - 0.822) * (df + 2) * 3) + 0.5 / (df + 4)) * y - 1) *
        ((df + 1) / (df + 2)) +
      1 / y
  }
  return Math.sqrt(df * y)
}

/** The t > 0 with P(T > t) = alpha, for 0 < alpha < 0.5. */
function upperQuantile(alpha: number, df: number): number {
  if (df === 1) return 1 / Math.tan(Math.PI * alpha)
  if (df === 2) return (1 - 2 * alpha) / Math.sqrt(2 * alpha * (1 - alpha))

  let t = df >= 1 ? hillStart(alpha, df) : 1 / Math.tan(Math.PI * alpha)
  if (!(t > 0) || !Number.isFinite(t)) t = -qnorm(alpha)

  // G(t) = P(T > t) is convex and decreasing for t > 0, so Newton's method converges
  // monotonically once an iterate lies below the root.
  const maxIter = df >= 1 ? 50 : 200
  for (let it = 0; it < maxIter; it++) {
    const f = dt(t, df)
    if (!(f > 0)) break
    const x = (pt(t, df, false) - alpha) / f
    if (!Number.isFinite(x)) break
    // Second-order Taylor step, using f'(t)/f(t) = −t(df + 1)/(df + t²).
    const step = x * (1 + (x * t * (df + 1)) / (2 * (t * t + df)))
    const next = t + step
    t = next > 0 ? next : t / 2
    if (Math.abs(step) <= 1e-14 * t) break
  }
  return t
}

/** t distribution quantile: the q with P(T ≤ q) = p, or P(T > q) = p when lowerTail is false. R: qt(). */
export function qt(p: number, df: number, lowerTail = true): number {
  if (Number.isNaN(p) || Number.isNaN(df) || p < 0 || p > 1 || df <= 0) return NaN
  if (!Number.isFinite(df)) return qnorm(p, lowerTail)
  // Reduce to the smaller tail: |q| solves P(T > |q|) = alpha, with alpha ≤ 0.5.
  const alpha = p <= 0.5 ? p : 1 - p
  const negative = p <= 0.5 === lowerTail
  if (alpha === 0.5) return 0
  if (alpha === 0) return negative ? -Infinity : Infinity
  const t = upperQuantile(alpha, df)
  return negative ? -t : t
}
