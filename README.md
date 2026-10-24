# TrialSizer

Real-time sample size and power for two-arm clinical trials. BIOSTAT 615 (Fall 2026)
final project, Boqian Mao. The full specification is in [SPEC.md](SPEC.md).

> TrialSizer is a planning aid. It does not replace a statistician, and the final sample
> size for a real trial should be confirmed with one.

**Status:** week 1 of 7. The project scaffold and the verified numerical core
(`src/stats/`) are in place; the app itself is a placeholder page. Live site (once
GitHub Pages is enabled): https://boqian-dot.github.io/Boqian_Mao_BIOSTAT615_First_Project/

## Local development

Requires Node.js 24 (see `.nvmrc`).

```sh
npm ci              # install the exact dependency versions from package-lock.json
npm run dev         # start the development server
npm test            # run all Vitest suites
npm run lint        # ESLint
npm run typecheck   # TypeScript, strict mode
npm run build       # production build into dist/
```

## How correctness is checked

The TypeScript routines are tested against R, which is the trusted reference
(SPEC §5). `reference/generate_fixtures.R` evaluates R's own functions on fixed grids
and writes JSON files to `tests/fixtures/`; the Vitest suites in `tests/` compare the
TypeScript results with them. The fixtures are committed, so CI does not need R. To
regenerate them (requires R with the `jsonlite` package):

```sh
npm run fixtures    # = Rscript reference/generate_fixtures.R
```

Fixtures are never edited by hand.

## Credits

The numerical routines in `src/stats/` are implementations of published algorithms:

- Normal CDF: W. J. Cody (1969), *Math. Comp.* 23, 631–637; coefficients from Cody's
  SPECFUN routine ANORM (netlib, public domain).
- Normal quantile: M. J. Wichura (1988), Algorithm AS 241, *Applied Statistics* 37, 477–484.
- Log-gamma: C. Lanczos (1964), *SIAM J. Numer. Anal. B* 1, 86–96, with P. Godfrey's
  (2001) coefficients for g = 7.
- Binomial and hypergeometric probabilities: C. Loader (2000), "Fast and accurate
  computation of binomial probabilities".
- Incomplete beta: continued fraction of Abramowitz & Stegun (1964) eq. 26.5.8, evaluated
  by the modified Lentz method (Lentz 1976; Thompson & Barnett 1986).
- t quantile: G. W. Hill (1970), Algorithm 396, *Comm. ACM* 13, 619–620.
- Noncentral t CDF: R. V. Lenth (1989), Algorithm AS 243, *Applied Statistics* 38, 185–189.

The reference values come from R (R Core Team), whose implementations of the same
algorithms served as the comparison standard.
