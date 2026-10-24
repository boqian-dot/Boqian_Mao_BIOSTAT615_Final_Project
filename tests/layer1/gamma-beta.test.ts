// Layer 1: log-gamma, log-beta and the incomplete beta function vs. R.
// SPEC §5 lists no tolerance for these building blocks; the ones below were set
// before the first run and are recorded in docs/dev-log.md:
//   lgamma, lbeta: |error| ≤ 1e-12 · max(1, |value|)  (relative, or absolute near the zeros at 1 and 2)
//   pbeta:         abs 1e-10, the tolerance of pt(), which is computed from it

import { describe, expect, it } from 'vitest'
import { lbeta, lgamma, pbeta } from '../../src/stats'
import { absError, loadFixture, mismatches } from '../helpers'

const LOG_TOL = 1e-12
const PBETA_ABS = 1e-10

const mixedError = (a: number, e: number) => Math.abs(a - e) / Math.max(1, Math.abs(e))

interface LgammaCase {
  x: number
  value: number
}
interface LbetaCase {
  a: number
  b: number
  value: number
}
interface PbetaCase {
  x: number
  a: number
  b: number
  lower: number
  upper: number
}

describe('lgamma vs R lgamma()', () => {
  it('matches on x in [1e-8, 1e6]', () => {
    const bad = mismatches(loadFixture<LgammaCase>('lgamma').cases, {
      label: (c) => `x=${c.x}`,
      actual: (c) => lgamma(c.x),
      expected: (c) => c.value,
      error: mixedError,
      tolerance: () => LOG_TOL,
    })
    expect(bad).toEqual([])
  })

  it('is exactly 0 at 1 and 2, and NaN for x ≤ 0', () => {
    expect(lgamma(1)).toBe(0)
    expect(lgamma(2)).toBe(0)
    expect(lgamma(0)).toBeNaN()
    expect(lgamma(-1.5)).toBeNaN()
  })
})

describe('lbeta vs R lbeta()', () => {
  it('matches on a, b in [0.5, 2500.5]', () => {
    const bad = mismatches(loadFixture<LbetaCase>('lbeta').cases, {
      label: (c) => `a=${c.a}, b=${c.b}`,
      actual: (c) => lbeta(c.a, c.b),
      expected: (c) => c.value,
      error: mixedError,
      tolerance: () => LOG_TOL,
    })
    expect(bad).toEqual([])
  })
})

describe('pbeta vs R pbeta()', () => {
  const cases = loadFixture<PbetaCase>('pbeta').cases
  it.each(['lower', 'upper'] as const)('%s tail', (tail) => {
    const bad = mismatches(cases, {
      label: (c) => `x=${c.x}, a=${c.a}, b=${c.b}`,
      actual: (c) => pbeta(c.x, c.a, c.b, tail === 'lower'),
      expected: (c) => c[tail],
      error: absError,
      tolerance: () => PBETA_ABS,
    })
    expect(bad).toEqual([])
  })
})
