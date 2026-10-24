/**
 * Log-gamma, Stirling's error term, and log-beta.
 *
 * lgamma — Lanczos (1964), "A precision approximation of the gamma function",
 *   SIAM J. Numer. Anal. Ser. B 1, 86–96, with g = 7 and the 9 coefficients of
 *   P. Godfrey (2001), "A note on the computation of the convergent Lanczos complex
 *   gamma approximation". The reflection formula Γ(x)Γ(1 − x) = π / sin(πx) is used
 *   for x < 0.5.
 *
 * stirlerr — the error of Stirling's formula, as used by C. Loader (2000), "Fast and
 *   accurate computation of binomial probabilities": tabulated for n = 0.5, 1, …, 15,
 *   and the asymptotic series 1/(12n) − 1/(360n³) + … for n > 15.
 *
 * lbeta — log B(a, b). When an argument is ≥ 10, the large lgamma terms are cancelled
 *   analytically with Stirling's formula, so the result keeps full accuracy (the same
 *   decomposition as SLATEC DLBETA and DiDonato & Morris 1992, ACM TOMS 18, 360–373).
 */

import { LN_SQRT_2PI } from './constants'

const LANCZOS_G = 7
const LANCZOS_COEF = [
  0.99999999999980993227684700473478, 676.520368121885098567009190444019,
  -1259.13921672240287047156078755283, 771.3234287776530788486528258894,
  -176.61502916214059906584551354, 12.507343278686904814458936853,
  -0.13857109526572011689554707, 9.984369578019570859563e-6, 1.50563273514931155834e-7,
]

/** log Γ(x) for x > 0. R: lgamma(). */
export function lgamma(x: number): number {
  if (Number.isNaN(x) || x <= 0) return NaN
  if (x === 1 || x === 2) return 0
  if (x === Infinity) return Infinity
  if (x < 0.5) {
    // Reflection: log Γ(x) = log(π / sin(πx)) − log Γ(1 − x)
    return Math.log(Math.PI / Math.sin(Math.PI * x)) - lgamma(1 - x)
  }
  const z = x - 1
  let a = LANCZOS_COEF[0]
  for (let i = 1; i < LANCZOS_COEF.length; i++) a += LANCZOS_COEF[i] / (z + i)
  const t = z + LANCZOS_G + 0.5
  return LN_SQRT_2PI + (z + 0.5) * Math.log(t) - t + Math.log(a)
}

/**
 * stirlerr(n) = log Γ(n + 1) − [(n + ½) log n − n + log √(2π)] for n = 0.5, 1, …, 15.
 * Index 2n. Computed to 25 significant digits from exact factorials (and
 * Γ(k + ½) = (2k)! √π / (4ᵏ k!)) with 60-digit decimal arithmetic. Entry 0 is unused.
 */
const STIRLERR_HALVES = [
  0, 1.5342640972002734529138394e-1, 8.1061466795327258219670264e-2,
  5.4814121051917653896138702e-2, 4.1340695955409294093822081e-2,
  3.316287351993628748511051e-2, 2.7677925684998339148789293e-2,
  2.3746163656297495971330279e-2, 2.0790672103765093111522772e-2,
  1.8488450532673185230779357e-2, 1.6644691189821192163194865e-2,
  1.5134973221917378873513837e-2, 1.3876128823070747998745727e-2,
  1.2810465242920226924250655e-2, 1.1896709945891770095055724e-2,
  1.1104559758206917326630755e-2, 1.0411265261972096497478567e-2,
  9.7994161261588032983903734e-3, 9.2554621827127329177286366e-3,
  8.7687001341393854629550473e-3, 8.3305634333628712564693187e-3,
  7.9341145643140205472495625e-3, 7.5736754879518407949720242e-3,
  7.2445543013203831795461966e-3, 6.9428401072095298656641527e-3,
  6.6652470327076824423561809e-3, 6.4089941880042070684396311e-3,
  6.1717122630394576475346048e-3, 5.951370112758847735624416e-3,
  5.7462165130101156820261025e-3, 5.55473355196280137103869e-3,
]

// Coefficients of the asymptotic series for stirlerr(n).
const S0 = 1 / 12
const S1 = 1 / 360
const S2 = 1 / 1260
const S3 = 1 / 1680
const S4 = 1 / 1188

/**
 * Error of Stirling's approximation to n!:
 * stirlerr(n) = log Γ(n + 1) − (n + ½) log n + n − log √(2π), for n > 0.
 */
export function stirlerr(n: number): number {
  if (n <= 15) {
    const twice = n + n
    if (twice === Math.trunc(twice)) return STIRLERR_HALVES[twice]
    return lgamma(n + 1) - (n + 0.5) * Math.log(n) + n - LN_SQRT_2PI
  }
  // Number of series terms chosen so that the truncation error is below 1e-16 (Loader 2000).
  const nn = n * n
  if (n > 500) return (S0 - S1 / nn) / n
  if (n > 80) return (S0 - (S1 - S2 / nn) / nn) / n
  if (n > 35) return (S0 - (S1 - (S2 - S3 / nn) / nn) / nn) / n
  return (S0 - (S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n
}

/** log B(a, b) = log Γ(a) + log Γ(b) − log Γ(a + b), for a, b > 0. R: lbeta(). */
export function lbeta(a: number, b: number): number {
  if (Number.isNaN(a) || Number.isNaN(b) || a <= 0 || b <= 0) return NaN
  const p = Math.min(a, b)
  const q = Math.max(a, b)
  if (q === Infinity) return -Infinity
  // Writing log Γ(x) = (x − ½) log x − x + log √(2π) + stirlerr(x) and cancelling terms:
  if (p >= 10) {
    const corr = stirlerr(p) + stirlerr(q) - stirlerr(p + q)
    return -0.5 * Math.log(q) + LN_SQRT_2PI + corr + (p - 0.5) * Math.log(p / (p + q)) + q * Math.log1p(-p / (p + q))
  }
  if (q >= 10) {
    const corr = stirlerr(q) - stirlerr(p + q)
    return lgamma(p) + corr + p - p * Math.log(p + q) + (q - 0.5) * Math.log1p(-p / (p + q))
  }
  return lgamma(p) + lgamma(q) - lgamma(p + q)
}
