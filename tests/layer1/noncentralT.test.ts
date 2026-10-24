// Layer 1 (SPEC §5): noncentral t CDF vs. R's pt(q, df, ncp). Tolerance: abs 1e-8.
//
// Cases with r_warning = 1 are those where R warned "full precision may not have
// been achieved in 'pnt{final}'". R raises it when the lower tail it computed is
// within 1e-10 of 1, i.e. the complementary tail has lost relative precision; the
// absolute accuracy tested here is not affected, so these cases are tested too.

import { describe, expect, it } from 'vitest'
import { pnt } from '../../src/stats'
import { absError, loadFixture, mismatches } from '../helpers'

const PNT_ABS = 1e-8
const SPEC_DF = [1, 2, 5, 10, 30, 100, 1000]
const SPEC_NCP = [0, 0.5, 1, 2, 5, 10]

interface PntCase {
  q: number
  df: number
  ncp: number
  lower: number
  upper: number
  r_warning: number
}

const cases = loadFixture<PntCase>('pnt').cases
const dfValues = [...new Set(cases.map((c) => c.df))]

describe('pnt vs R pt(q, df, ncp)', () => {
  it('fixture includes every df and ncp required by SPEC §5', () => {
    expect(dfValues).toEqual(expect.arrayContaining(SPEC_DF))
    expect([...new Set(cases.map((c) => c.ncp))]).toEqual(expect.arrayContaining(SPEC_NCP))
  })

  describe.each(dfValues)('df = %s', (df) => {
    it.each(['lower', 'upper'] as const)('%s tail', (tail) => {
      const bad = mismatches(
        cases.filter((c) => c.df === df),
        {
          label: (c) => `q=${c.q}, df=${c.df}, ncp=${c.ncp}`,
          actual: (c) => pnt(c.q, c.df, c.ncp, tail === 'lower'),
          expected: (c) => c[tail],
          error: absError,
          tolerance: () => PNT_ABS,
        },
      )
      expect(bad).toEqual([])
    })
  })
})
