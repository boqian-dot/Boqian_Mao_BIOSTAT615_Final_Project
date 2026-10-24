/**
 * Regularized incomplete beta function I_x(a, b) = P(X ≤ x) for X ~ Beta(a, b).
 *
 * Continued fraction for I_x(a, b) (Abramowitz & Stegun 1964, eq. 26.5.8), evaluated
 * with the modified Lentz algorithm (W. J. Lentz 1976, Applied Optics 15, 668–671;
 * I. J. Thompson & A. R. Barnett 1986, J. Comput. Phys. 64, 490–509). The fraction
 * converges quickly for x < (a + 1)/(a + b + 2); otherwise the symmetry
 * I_x(a, b) = 1 − I_{1−x}(b, a) is used. Both tails are returned, so a small tail
 * probability is never formed as 1 − (a number close to 1).
 */

import { lbeta } from './gamma'

const MAX_ITER = 10000
const CF_EPS = 1e-15
const TINY = 1e-300

/** Continued fraction part of I_x(a, b) (A&S 26.5.8), by the modified Lentz method. */
function betaContinuedFraction(x: number, a: number, b: number): number {
  const qab = a + b
  const qap = a + 1
  const qam = a - 1
  let c = 1
  let d = 1 - (qab * x) / qap
  if (Math.abs(d) < TINY) d = TINY
  d = 1 / d
  let h = d
  for (let m = 1; m <= MAX_ITER; m++) {
    const m2 = 2 * m
    // Even term: d_{2m} = m(b − m)x / ((a + 2m − 1)(a + 2m))
    let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2))
    d = 1 + aa * d
    if (Math.abs(d) < TINY) d = TINY
    c = 1 + aa / c
    if (Math.abs(c) < TINY) c = TINY
    d = 1 / d
    h *= d * c
    // Odd term: d_{2m+1} = −(a + m)(a + b + m)x / ((a + 2m)(a + 2m + 1))
    aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2))
    d = 1 + aa * d
    if (Math.abs(d) < TINY) d = TINY
    c = 1 + aa / c
    if (Math.abs(c) < TINY) c = TINY
    d = 1 / d
    const del = d * c
    h *= del
    if (Math.abs(del - 1) < CF_EPS) break
  }
  return h
}

/**
 * Both tails of the Beta(a, b) distribution at x: [I_x(a, b), 1 − I_x(a, b)].
 * The caller passes y = 1 − x separately, so that it can be computed without
 * cancellation (e.g. y = df / (df + t²) in the t distribution).
 */
export function pbetaBoth(x: number, y: number, a: number, b: number): [number, number] {
  if (Number.isNaN(x) || Number.isNaN(y) || !(a > 0) || !(b > 0)) return [NaN, NaN]
  if (x <= 0) return [0, 1]
  if (y <= 0) return [1, 0]
  // log of x^a (1 − x)^b / B(a, b), the factor in front of the continued fraction
  const logFront = a * Math.log(x) + b * Math.log(y) - lbeta(a, b)
  if (x < (a + 1) / (a + b + 2)) {
    const lower = (Math.exp(logFront) * betaContinuedFraction(x, a, b)) / a
    return [lower, 1 - lower]
  }
  const upper = (Math.exp(logFront) * betaContinuedFraction(y, b, a)) / b
  return [1 - upper, upper]
}

/** Beta(a, b) CDF, P(X ≤ x), or P(X > x) when lowerTail is false. R: pbeta(). */
export function pbeta(x: number, a: number, b: number, lowerTail = true): number {
  const [lower, upper] = pbetaBoth(x, 1 - x, a, b)
  return lowerTail ? lower : upper
}
