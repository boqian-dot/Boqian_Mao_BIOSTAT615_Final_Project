// Layer 1 (SPEC §5): normal CDF and quantile vs. R's pnorm() and qnorm().
// Tolerance: abs 1e-12 in the centre, rel 1e-9 in the tails. A case is in the
// tail when its probability (the CDF value, or min(p, 1 - p) for the quantile)
// is below 1e-3.

import { describe, expect, it } from 'vitest'
import { pnorm, qnorm } from '../../src/stats'
import { absError, loadFixture, mismatches, relError } from '../helpers'

const TAIL = 1e-3
const ABS_CENTRAL = 1e-12
const REL_TAIL = 1e-9

interface PnormCase {
  x: number
  lower: number
  upper: number
}
interface QnormCase {
  p: number
  lower: number
  upper: number
}

const pnormCases = loadFixture<PnormCase>('pnorm').cases
const qnormCases = loadFixture<QnormCase>('qnorm').cases
const tails = ['lower', 'upper'] as const

describe('pnorm vs R pnorm()', () => {
  it.each(tails)('%s tail', (tail) => {
    const bad = mismatches(pnormCases, {
      label: (c) => `x=${c.x}`,
      actual: (c) => pnorm(c.x, tail === 'lower'),
      expected: (c) => c[tail],
      error: (a, e) => (e < TAIL ? relError(a, e) : absError(a, e)),
      tolerance: (c) => (c[tail] < TAIL ? REL_TAIL : ABS_CENTRAL),
    })
    expect(bad).toEqual([])
  })

  it('fixture covers the probabilities 1e-12 … 1 - 1e-12 in both tails', () => {
    expect(Math.min(...pnormCases.map((c) => c.lower))).toBeLessThan(1e-12)
    expect(Math.min(...pnormCases.map((c) => c.upper))).toBeLessThan(1e-12)
  })
})

describe('qnorm vs R qnorm()', () => {
  it.each(tails)('%s tail', (tail) => {
    const inTail = (c: QnormCase) => Math.min(c.p, 1 - c.p) < TAIL
    const bad = mismatches(qnormCases, {
      label: (c) => `p=${c.p}`,
      actual: (c) => qnorm(c.p, tail === 'lower'),
      expected: (c) => c[tail],
      error: (a, e, c) => (inTail(c) ? relError(a, e) : absError(a, e)),
      tolerance: (c) => (inTail(c) ? REL_TAIL : ABS_CENTRAL),
    })
    expect(bad).toEqual([])
  })

  it('fixture covers p in [1e-12, 1 - 1e-12]', () => {
    const ps = qnormCases.map((c) => c.p)
    expect(Math.min(...ps)).toBeLessThanOrEqual(1e-12)
    expect(Math.max(...ps)).toBeGreaterThanOrEqual(1 - 1e-12)
  })

  it('returns ±Infinity at 0 and 1 and NaN outside [0, 1]', () => {
    expect(qnorm(0)).toBe(-Infinity)
    expect(qnorm(1)).toBe(Infinity)
    expect(qnorm(-0.1)).toBeNaN()
    expect(qnorm(1.1)).toBeNaN()
  })
})
