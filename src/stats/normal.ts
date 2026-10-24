/**
 * Standard normal distribution: CDF and quantile function.
 *
 * pnorm — W. J. Cody (1969), "Rational Chebyshev approximations for the error
 *   function", Mathematics of Computation 23, 631–637. The coefficients are those
 *   of Cody's routine ANORM in the SPECFUN package (netlib, public domain), which
 *   R's pnorm() also uses. Accuracy is close to double precision on the whole line,
 *   and both tails are computed directly, so small tail probabilities keep their
 *   relative accuracy.
 *
 * qnorm — M. J. Wichura (1988), "Algorithm AS 241: The percentage points of the
 *   normal distribution", Applied Statistics 37, 477–484 (routine PPND16, relative
 *   accuracy about 1e-16 for p ≥ 1e-300).
 */

import { INV_SQRT_2PI, SQRT_32 } from './constants'

// Cody's ANORM coefficients: |x| ≤ 0.674 (A, B), 0.674 < |x| ≤ sqrt(32) (C, D), |x| > sqrt(32) (P, Q).
const A = [
  2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582, 18154.981253343561249,
  0.065682337918207449113,
]
const B = [47.20258190468824187, 976.09855173777669322, 10260.932208618978205, 45507.789335026729956]
const C = [
  0.39894151208813466764, 8.8831497943883759412, 93.506656132177855979, 597.27027639480026226,
  2494.5375852903726711, 6848.1904505362823326, 11602.651437647350124, 9842.7148383839780218,
  1.0765576773720192317e-8,
]
const D = [
  22.266688044328115691, 235.38790178262499861, 1519.377599407554805, 6485.558298266760755,
  18615.571640885098091, 34900.952721145977266, 38912.003286093271411, 19685.429676859990727,
]
const P = [
  0.21589853405795699, 0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466,
  2.9112874951168792e-5, 0.02307344176494017303,
]
const Q = [
  1.28426009614491121, 0.468238212480865118, 0.0659881378689285515, 0.00378239633202758244,
  7.29751555083966205e-5,
]

/**
 * exp(-y²/2), computed as exp(-a²/2) · exp(-(y - a)(y + a)/2) where a is y rounded
 * toward zero to a multiple of 1/16. This avoids the rounding error of forming y²
 * directly when y is large (Cody 1969).
 */
function expHalfSquare(y: number): number {
  const a = Math.trunc(y * 16) / 16
  const del = (y - a) * (y + a)
  return Math.exp(-a * a * 0.5) * Math.exp(-del * 0.5)
}

/**
 * Both tails of the standard normal distribution at x: [P(Z ≤ x), P(Z > x)].
 * Equivalent to R's pnorm(x) and pnorm(x, lower.tail = FALSE).
 */
export function pnormBoth(x: number): [number, number] {
  if (Number.isNaN(x)) return [NaN, NaN]
  const y = Math.abs(x)

  if (y <= 0.67448975) {
    // |x| ≤ qnorm(3/4): rational approximation to erf(x / sqrt(2)). Cody's original
    // boundary is 0.66291; R uses qnorm(3/4), and so do we, so that the branches
    // match the reference (the fixtures include points on both sides).
    let xnum = 0
    let xden = 0
    if (y > Number.EPSILON / 2) {
      const xsq = x * x
      xnum = A[4] * xsq
      xden = xsq
      for (let i = 0; i < 3; i++) {
        xnum = (xnum + A[i]) * xsq
        xden = (xden + B[i]) * xsq
      }
    }
    const temp = (x * (xnum + A[3])) / (xden + B[3])
    return [0.5 + temp, 0.5 - temp]
  }

  // Otherwise compute the smaller tail, Φ(-|x|), and get the other one by subtraction.
  let tail: number
  if (y <= SQRT_32) {
    let xnum = C[8] * y
    let xden = y
    for (let i = 0; i < 7; i++) {
      xnum = (xnum + C[i]) * y
      xden = (xden + D[i]) * y
    }
    tail = expHalfSquare(y) * ((xnum + C[7]) / (xden + D[7]))
  } else {
    const xsq = 1 / (x * x)
    let xnum = P[5] * xsq
    let xden = xsq
    for (let i = 0; i < 4; i++) {
      xnum = (xnum + P[i]) * xsq
      xden = (xden + Q[i]) * xsq
    }
    const r = (xsq * (xnum + P[4])) / (xden + Q[4])
    tail = expHalfSquare(y) * ((INV_SQRT_2PI - r) / y)
  }
  return x > 0 ? [1 - tail, tail] : [tail, 1 - tail]
}

/** Standard normal CDF, P(Z ≤ x), or P(Z > x) when lowerTail is false. R: pnorm(). */
export function pnorm(x: number, lowerTail = true): number {
  const [lower, upper] = pnormBoth(x)
  return lowerTail ? lower : upper
}

// Wichura's AS 241 coefficients, in increasing powers of r (paper notation a0..a7, b0..b7, …).
// Central region |p - 0.5| ≤ 0.425:
const AS241_A = [
  3.387132872796366608, 133.14166789178437745, 1971.5909503065514427, 13731.693765509461125,
  45921.953931549871457, 67265.770927008700853, 33430.575583588128105, 2509.0809287301226727,
]
const AS241_B = [
  1, 42.313330701600911252, 687.1870074920579083, 5394.1960214247511077, 21213.794301586595867,
  39307.89580009271061, 28729.085735721942674, 5226.495278852545925,
]
// Intermediate tail, r = sqrt(-log(min(p, 1 - p))) ≤ 5:
const AS241_C = [
  1.42343711074968357734, 4.6303378461565452959, 5.7694972214606914055, 3.64784832476320460504,
  1.27045825245236838258, 0.24178072517745061177, 0.0227238449892691845833, 7.7454501427834140764e-4,
]
const AS241_D = [
  1, 2.05319162663775882187, 1.6763848301838038494, 0.68976733498510000455, 0.14810397642748007459,
  0.0151986665636164571966, 5.475938084995344946e-4, 1.05075007164441684324e-9,
]
// Far tail, r > 5:
const AS241_E = [
  6.6579046435011037772, 5.4637849111641143699, 1.7848265399172913358, 0.29656057182850489123,
  0.026532189526576123093, 0.0012426609473880784386, 2.71155556874348757815e-5,
  2.01033439929228813265e-7,
]
const AS241_F = [
  1, 0.59983220655588793769, 0.13692988092273580531, 0.0148753612908506148525,
  7.868691311456132591e-4, 1.8463183175100546818e-5, 1.4215117583164458887e-7,
  2.04426310338993978564e-15,
]

/** Evaluates c[0] + c[1]·x + … + c[n]·xⁿ by Horner's rule. */
function horner(c: readonly number[], x: number): number {
  let s = c[c.length - 1]
  for (let i = c.length - 2; i >= 0; i--) s = s * x + c[i]
  return s
}

/**
 * Standard normal quantile: the z with P(Z ≤ z) = p, or P(Z > z) = p when lowerTail
 * is false. R: qnorm(). Wichura's AS 241 (PPND16).
 */
export function qnorm(p: number, lowerTail = true): number {
  if (Number.isNaN(p) || p < 0 || p > 1) return NaN
  if (p === 0) return lowerTail ? -Infinity : Infinity
  if (p === 1) return lowerTail ? Infinity : -Infinity

  const pLower = lowerTail ? p : 1 - p
  const pUpper = lowerTail ? 1 - p : p
  const q = pLower - 0.5

  if (Math.abs(q) <= 0.425) {
    const r = 0.180625 - q * q
    return (q * horner(AS241_A, r)) / horner(AS241_B, r)
  }

  // r = sqrt(-log(smaller tail probability)), taken from whichever tail was supplied.
  const r = Math.sqrt(-Math.log(q > 0 ? pUpper : pLower))
  const val =
    r <= 5
      ? horner(AS241_C, r - 1.6) / horner(AS241_D, r - 1.6)
      : horner(AS241_E, r - 5) / horner(AS241_F, r - 5)
  return q < 0 ? -val : val
}
