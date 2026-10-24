/**
 * Noncentral t distribution CDF.
 *
 * R. V. Lenth (1989), "Algorithm AS 243: Cumulative distribution function of the
 * non-central t distribution", Applied Statistics 38, 185–189. It sums the twin
 * series of J. Guenther (1978), J. Statist. Comput. Simul. 6, 199–206:
 *
 *   P(T ≤ t; δ) = Φ(−δ) + ½ Σ_j [ p_j I_x(j + ½, df/2) + q_j I_x(j + 1, df/2) ],  x = t²/(t² + df),
 *
 * with Poisson weights p_j = e^(−δ²/2)(δ²/2)^j / j! and q_j = δ e^(−δ²/2)(δ²/2)^j / (√2 Γ(j + 3/2)).
 * The incomplete beta terms are updated by recurrence, and the sum stops when
 * Lenth's error bound falls below 1e-12. Negative t uses
 * P(T ≤ t; δ) = 1 − P(T ≤ −t; −δ).
 *
 * As in R's pt(q, df, ncp), the normal approximation of Abramowitz & Stegun
 * eq. 26.7.10 is used when df > 4e5 or δ² > 2·log(2)·1021 (where e^(−δ²/2) would
 * underflow).
 */

import { pbetaBoth } from './beta'
import { SQRT_2_OVER_PI } from './constants'
import { lbeta } from './gamma'
import { pnorm } from './normal'
import { pt } from './t'

const MAX_ITER = 1000
const ERR_MAX = 1e-12

/**
 * Noncentral t CDF with df degrees of freedom and noncentrality ncp:
 * P(T ≤ q), or P(T > q) when lowerTail is false. R: pt(q, df, ncp).
 */
export function pnt(q: number, df: number, ncp: number, lowerTail = true): number {
  if (Number.isNaN(q) || Number.isNaN(df) || Number.isNaN(ncp) || df <= 0) return NaN
  if (ncp === 0) return pt(q, df, lowerTail)
  if (!Number.isFinite(q)) return q < 0 === lowerTail ? 0 : 1

  // Work with tt ≥ 0; for q < 0 use P(T ≤ q; δ) = 1 − P(T ≤ −q; −δ).
  const negate = q < 0
  const tt = negate ? -q : q
  const del = negate ? -ncp : ncp
  // The series below gives P(T ≤ tt; del); `wantLower` says whether that is the tail requested.
  const wantLower = lowerTail !== negate

  if (df > 4e5 || del * del > 2 * Math.LN2 * 1021) {
    const s = 1 / (4 * df)
    return pnorm((tt * (1 - s) - del) / Math.sqrt(1 + tt * tt * 2 * s), wantLower)
  }

  let tnc = 0
  if (tt > 0) {
    const tt2 = tt * tt
    const x = tt2 / (tt2 + df)
    const y = df / (tt2 + df) // 1 − x, without cancellation
    const logY = -Math.log1p(tt2 / df)
    const lambda = del * del
    let p = 0.5 * Math.exp(-0.5 * lambda)
    let qw = SQRT_2_OVER_PI * p * del
    // s = ½ − Σ p_j tracks the Poisson mass not yet used (for the error bound).
    let s = 0.5 - p
    if (s < 1e-7) s = -0.5 * Math.expm1(-0.5 * lambda)
    let a = 0.5
    const b = 0.5 * df
    const rxb = Math.exp(b * logY) // (1 − x)^b
    let xodd = pbetaBoth(x, y, a, b)[0] // I_x(1/2, b)
    let godd = 2 * rxb * Math.exp(a * Math.log(x) - lbeta(0.5, b))
    // I_x(1, b) = 1 − (1 − x)^b, computed with expm1 for accuracy at small x.
    let xeven = -Math.expm1(b * logY)
    let geven = b * x * rxb
    tnc = p * xodd + qw * xeven

    for (let it = 1; it <= MAX_ITER; it++) {
      a += 1
      xodd -= godd
      xeven -= geven
      godd *= (x * (a + b - 1)) / a
      geven *= (x * (a + b - 0.5)) / (a + 0.5)
      p *= lambda / (2 * it)
      qw *= lambda / (2 * it + 1)
      tnc += p * xodd + qw * xeven
      s -= p
      // s < 0 means rounding error has used up the Poisson mass: stop.
      if (s < -1e-10) break
      if (s <= 0 && it > 1) break
      const errbd = 2 * s * (xodd - godd)
      if (Math.abs(errbd) < ERR_MAX) break
    }
  }
  tnc += pnorm(-del)
  const lower = Math.min(tnc, 1)
  return wantLower ? lower : 1 - lower
}
