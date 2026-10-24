// Layer 1 (SPEC §5): Student t CDF and quantile vs. R's pt() and qt().
// Tolerances: pt abs 1e-10; qt rel 1e-9 for p in [1e-8, 1 - 1e-8].
// dt (the density, used by qt's Newton steps) has no SPEC tolerance; rel 1e-11 was
// set before the first run and is recorded in docs/dev-log.md.

import { describe, expect, it } from 'vitest'
import { dt, pt, qt } from '../../src/stats'
import { absError, loadFixture, mismatches, relError } from '../helpers'

const PT_ABS = 1e-10
const QT_REL = 1e-9
const DT_REL = 1e-11
const SPEC_DF = [1, 2, 5, 10, 30, 100, 1000]

interface PtCase {
  q: number
  df: number
  lower: number
  upper: number
  density: number
}
interface QtCase {
  p: number
  df: number
  lower: number
  upper: number
}

const ptCases = loadFixture<PtCase>('pt').cases
const qtCases = loadFixture<QtCase>('qt').cases
const dfValues = [...new Set(ptCases.map((c) => c.df))]
const tails = ['lower', 'upper'] as const

describe('pt vs R pt()', () => {
  it('fixture includes every df required by SPEC §5', () => {
    expect(dfValues).toEqual(expect.arrayContaining(SPEC_DF))
  })

  describe.each(dfValues)('df = %s', (df) => {
    it.each(tails)('%s tail', (tail) => {
      const bad = mismatches(
        ptCases.filter((c) => c.df === df),
        {
          label: (c) => `q=${c.q}, df=${c.df}`,
          actual: (c) => pt(c.q, c.df, tail === 'lower'),
          expected: (c) => c[tail],
          error: absError,
          tolerance: () => PT_ABS,
        },
      )
      expect(bad).toEqual([])
    })
  })
})

describe('dt vs R dt()', () => {
  it('matches the density', () => {
    const bad = mismatches(ptCases, {
      label: (c) => `x=${c.q}, df=${c.df}`,
      actual: (c) => dt(c.q, c.df),
      expected: (c) => c.density,
      error: relError,
      tolerance: () => DT_REL,
    })
    expect(bad).toEqual([])
  })
})

describe('qt vs R qt()', () => {
  it('fixture covers p in [1e-8, 1 - 1e-8] and every SPEC df', () => {
    const ps = qtCases.map((c) => c.p)
    expect(Math.min(...ps)).toBeLessThanOrEqual(1e-8)
    expect(Math.max(...ps)).toBeGreaterThanOrEqual(1 - 1e-8)
    expect([...new Set(qtCases.map((c) => c.df))]).toEqual(expect.arrayContaining(SPEC_DF))
  })

  describe.each(dfValues)('df = %s', (df) => {
    it.each(tails)('%s tail', (tail) => {
      const bad = mismatches(
        qtCases.filter((c) => c.df === df),
        {
          label: (c) => `p=${c.p}, df=${c.df}`,
          actual: (c) => qt(c.p, c.df, tail === 'lower'),
          expected: (c) => c[tail],
          error: relError,
          tolerance: () => QT_REL,
        },
      )
      expect(bad).toEqual([])
    })
  })
})
