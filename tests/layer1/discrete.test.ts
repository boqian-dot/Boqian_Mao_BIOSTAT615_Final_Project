// Layer 1 (SPEC §5): log binomial and hypergeometric pmfs vs. R's
// dbinom(log = TRUE) and dhyper(log = TRUE), n up to 5000. Tolerance: rel 1e-10.

import { describe, expect, it } from 'vitest'
import { dbinomLog, dhyperLog } from '../../src/stats'
import { loadFixture, mismatches, relError } from '../helpers'

const REL = 1e-10

interface BinomCase {
  x: number
  n: number
  p: number
  value: number
}
interface HyperCase {
  x: number
  m: number
  n: number
  k: number
  value: number
}

describe('dbinomLog vs R dbinom(log = TRUE)', () => {
  const cases = loadFixture<BinomCase>('dbinom_log').cases

  it('fixture reaches n = 5000', () => {
    expect(Math.max(...cases.map((c) => c.n))).toBeGreaterThanOrEqual(5000)
  })

  it('matches on n ≤ 5000, p in [0.001, 0.999]', () => {
    const bad = mismatches(cases, {
      label: (c) => `x=${c.x}, n=${c.n}, p=${c.p}`,
      actual: (c) => dbinomLog(c.x, c.n, c.p),
      expected: (c) => c.value,
      error: relError,
      tolerance: () => REL,
    })
    expect(bad).toEqual([])
  })

  it('handles the support boundary and degenerate p', () => {
    expect(dbinomLog(-1, 10, 0.3)).toBe(-Infinity)
    expect(dbinomLog(11, 10, 0.3)).toBe(-Infinity)
    expect(dbinomLog(0, 10, 0)).toBe(0)
    expect(dbinomLog(10, 10, 1)).toBe(0)
    expect(dbinomLog(3, 10, 1.5)).toBeNaN()
  })
})

describe('dhyperLog vs R dhyper(log = TRUE)', () => {
  const cases = loadFixture<HyperCase>('dhyper_log').cases

  it('fixture reaches m = n = 5000', () => {
    expect(Math.max(...cases.map((c) => c.m))).toBeGreaterThanOrEqual(5000)
    expect(Math.max(...cases.map((c) => c.n))).toBeGreaterThanOrEqual(5000)
  })

  it('matches on m, n ≤ 5000', () => {
    const bad = mismatches(cases, {
      label: (c) => `x=${c.x}, m=${c.m}, n=${c.n}, k=${c.k}`,
      actual: (c) => dhyperLog(c.x, c.m, c.n, c.k),
      expected: (c) => c.value,
      error: relError,
      tolerance: () => REL,
    })
    expect(bad).toEqual([])
  })

  it('returns −∞ outside the support', () => {
    // m = 5 white, n = 3 black, k = 6 draws: support is 3 … 5
    expect(dhyperLog(2, 5, 3, 6)).toBe(-Infinity)
    expect(dhyperLog(6, 5, 3, 6)).toBe(-Infinity)
    expect(dhyperLog(0, 5, 3, 0)).toBe(0)
  })
})
