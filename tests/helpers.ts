import { readFileSync } from 'node:fs'

/** A JSON fixture written by reference/generate_fixtures.R. */
export interface Fixture<Case> {
  description: string
  generated_by: string
  r_version: string
  n_cases: number
  cases: Case[]
}

export function loadFixture<Case>(name: string): Fixture<Case> {
  const url = new URL(`./fixtures/${name}.json`, import.meta.url)
  const fixture = JSON.parse(readFileSync(url, 'utf8')) as Fixture<Case>
  if (fixture.cases.length !== fixture.n_cases) {
    throw new Error(`fixture ${name}: expected ${fixture.n_cases} cases, found ${fixture.cases.length}`)
  }
  return fixture
}

export function absError(actual: number, expected: number): number {
  return Math.abs(actual - expected)
}

/** Relative error; when the reference is exactly 0, only an exact 0 passes. */
export function relError(actual: number, expected: number): number {
  if (expected === 0) return actual === 0 ? 0 : Infinity
  return Math.abs(actual - expected) / Math.abs(expected)
}

export interface Mismatch {
  input: string
  actual: number
  expected: number
  error: number
  tolerance: number
}

export interface Check<Case> {
  /** Describes the input, e.g. "q=2, df=5"; shown when the check fails. */
  label: (c: Case) => string
  actual: (c: Case) => number
  expected: (c: Case) => number
  /** Error measure and tolerance for this case (they may depend on the case). */
  error: (actual: number, expected: number, c: Case) => number
  tolerance: (c: Case) => number
}

/**
 * Returns the cases whose error exceeds the tolerance (or is NaN), worst first,
 * at most `limit` of them. A passing check returns [], so a test can assert
 * `expect(mismatches(...)).toEqual([])` and a failure lists the offending inputs.
 */
export function mismatches<Case>(cases: readonly Case[], check: Check<Case>, limit = 10): Mismatch[] {
  const bad: Mismatch[] = []
  for (const c of cases) {
    const actual = check.actual(c)
    const expected = check.expected(c)
    const error = check.error(actual, expected, c)
    const tolerance = check.tolerance(c)
    if (!(error <= tolerance)) bad.push({ input: check.label(c), actual, expected, error, tolerance })
  }
  const severity = (m: Mismatch) => (Number.isNaN(m.error) ? Infinity : m.error / m.tolerance)
  bad.sort((a, b) => severity(b) - severity(a))
  return bad.slice(0, limit)
}
