# Development log

One entry per working session with the coding agent (SPEC §10, rule 4). This log is the
raw material for `AI_USAGE.md`.

---

## 2026-10-24 · Session 1 · Week 1: scaffold and Layer 1 numerics

**Goal.** The week-1 milestone of SPEC §8: a Vite + React + TypeScript (strict) scaffold
with ESLint, Vitest and a GitHub Actions workflow that deploys a placeholder page to
GitHub Pages; the `src/stats/` numerical routines; `reference/generate_fixtures.R`;
Layer 1 tests passing; start this log.

**Agent.** Claude Code (Claude Opus 5.5) in VS Code.

**Key prompt (verbatim).**

> 按照计划开始做第一步

("Following the plan, start on the first step.") It was sent with SPEC §8–9 selected in
the editor; the agent took "the first step" to be the week-1 row of the milestone table.

### What the agent changed

| Commit | Content |
|--------|---------|
| `fd87323` | SPEC.md |
| `f94e169` | Scaffold: `package.json`, TypeScript configs (`strict: true`), `vite.config.ts` (GitHub Pages base path, Vitest), `eslint.config.js` (forbids React, UI imports and DOM globals in `src/stats/`), placeholder page with the disclaimer |
| `f95c8a9` | `reference/generate_fixtures.R` and 10 fixture files in `tests/fixtures/` (4,020 cases) |
| `1def0a4` | `src/stats/normal.ts`: `pnorm` (Cody), `qnorm` (AS 241); `tsconfig.stats.json`, which type-checks `src/stats/` without the DOM library; test helpers |
| `b3b6fe9` | `src/stats/gamma.ts` (`lgamma` Lanczos, `stirlerr`, `lbeta`), `src/stats/beta.ts` (incomplete beta) |
| `fdcba22` | `src/stats/t.ts`: `pt`, `dt`, `qt` (Hill's start + Newton) |
| `a059a30` | `src/stats/noncentralT.ts`: `pnt` (AS 243) |
| `c825af7` | `src/stats/discrete.ts`: `dbinomLog`, `dhyperLog` (Loader) |
| `32a3391` | Placeholder page shows z₀.₉₇₅ computed by the in-house `qnorm` |
| `e5e623b` | `.github/workflows/ci.yml` (lint, type-check, test, build; deploy to Pages only from `main` and only if all pass) and `fixtures.yml` (manual fixture regeneration with R) |
| next commit | README (local commands, how fixtures are made, algorithm credits) and this log |

### Decisions made in this session (not fixed by the SPEC)

1. **Tolerances for building blocks that SPEC §5 does not list.** Set before the first
   test run: `lgamma` and `lbeta` |error| ≤ 1e-12 · max(1, |value|) (relative, but
   absolute near the zeros of log Γ at 1 and 2); `pbeta` abs 1e-10, the tolerance of
   `pt`, which is computed from it; `dt` rel 1e-11.
2. **"Central" vs. "tails" for the normal tolerance.** A case is in the tail when its
   probability (the CDF value, or min(p, 1 − p) for `qnorm`) is below 1e-3. Tails are
   tested at rel 1e-9, everything else at abs 1e-12.
3. **Fixture grids go beyond the SPEC minimum.** Normal tails out to |x| = 37 (≈ 1e-300);
   non-integer df (1.5, 7.3, 47.8) for `pt`/`qt`, because Welch df are not integers;
   negative noncentrality (−1, −5) for `pnt`; upper tails for every CDF.
4. **Binomial and hypergeometric probabilities use Loader's saddle-point method,** not
   differences of `lgamma` values. It avoids cancellation at large n and is the method R
   uses. Lanczos `lgamma` remains the log-gamma routine (used by `stirlerr` and `lbeta`).
5. **Implemented from the published papers, not translated from R's C source.** R is
   GPL; the project is planned to be MIT-licensed. R is used only as the reference.
6. **Small, numerically equivalent changes to the published algorithms** (each is noted
   in the code): AS 243's I_x(1, b) = 1 − (1 − x)^b is computed with `expm1`; Loader's
   x = 0 and x = n cases use `log1p`; `lbeta` cancels the large log Γ terms analytically
   with Stirling's series; Cody's first-branch boundary is R's 0.67448975 instead of 0.66291.
7. **Scope of "numerical routines" in week 1:** the routines tested by Layer 1. The PRNG,
   random variates and the integer n search are part of the simulation engine and come in
   week 2 with their Layer 3 tests.
8. **Git.** The project folder sat inside a git repository rooted at the home directory.
   A separate repository was initialised in the project folder, with `origin` set to the
   GitHub repository. Nothing has been pushed yet. The proposal files
   (`project_proposal.html`, `project_proposal_template.pdf`) are left untracked.

### What went wrong and how it was corrected

- **Node.js was not installed.** The agent installed the official Node 24.22.1 LTS binary
  (SHA-256 checked against nodejs.org) into `~/.local/node`, without editing the shell profile.
- **create-vite 9 now ships oxlint instead of ESLint** and does not set `strict`
  explicitly. The SPEC requires ESLint and strict mode, so the configs were written by hand.
- **TypeScript 7.0 is the newest release, but typescript-eslint supports only < 6.1.**
  TypeScript is pinned to ~6.0.3.
- **Bug in the R script.** `c(-1, 1) * c(9, 10, …)` recycled the vectors instead of giving
  ± each value; R printed a warning. Fixed with `outer()`. The script now runs with
  `options(warn = 2)`, so any unexpected warning stops it.
- **R's `pt(q, df, ncp)` warned "full precision may not have been achieved in
  'pnt{final}'" on 191 of 1,120 grid points.** R raises this when the lower tail it
  computed is within 1e-10 of 1, that is, when the complementary tail has lost *relative*
  precision. Absolute accuracy, which is what Layer 1 tests (abs 1e-8), is not affected.
  The points were kept and flagged with `r_warning = 1` in the fixture.
- **ESLint's `no-loss-of-precision` rule flagged 125 coefficient literals,** which are
  published with up to 20 significant digits. The rule is turned off for `src/stats/` only,
  with a comment; the coefficients stay verbatim so they can be checked against the papers.
- **All 76 tests passed on the first run.** A first-run pass could also mean the tests
  are too weak, so the agent ran the mutation check below.

### How the result was verified

- `npm run lint`, `npm run typecheck`, `npm test` (5 files, 76 tests) and `npm run build`
  pass. Each commit was rechecked in a clean git worktree.
- Largest error against R over all fixture cases:

  | Routine | Tolerance | Largest error |
  |---------|-----------|---------------|
  | `pnorm`, centre / tails | abs 1e-12 / rel 1e-9 | 6.9e-18 / 1.7e-16 |
  | `qnorm`, centre / tails | abs 1e-12 / rel 1e-9 | 8.9e-16 / 0 (identical) |
  | `lgamma` | 1e-12 · max(1, \|value\|) | 1.8e-15 |
  | `lbeta` | 1e-12 · max(1, \|value\|) | 2.0e-15 |
  | `pbeta` | abs 1e-10 | 3.5e-14 (x = 0.5, a = b = 500) |
  | `pt` | abs 1e-10 | 7.7e-15 |
  | `dt` | rel 1e-11 | 5.3e-14 |
  | `qt` | rel 1e-9 | 8.6e-14 |
  | `pnt` | abs 1e-8 | 1.7e-13 |
  | `dbinomLog` | rel 1e-10 | 4.8e-16 |
  | `dhyperLog` | rel 1e-10 | 2.2e-13 |

- **Mutation check:** four deliberate bugs, each introduced alone and then reverted. The
  suite caught every one:

  | Mutation | Tests failed |
  |----------|--------------|
  | AS 241 central constant 0.180625 → 0.180626 | 2 |
  | AS 243 error bound 1e-12 → 1e-4 | 14 |
  | Loader `bd0` series term 2j + 1 → 2j + 2 | 2 |
  | `qt` Newton refinement switched off (Hill's start only) | 12 |

- **DOM guard:** a probe file in `src/stats/` that used `document` was rejected by both
  ESLint and `tsc -p tsconfig.stats.json`.
- **Built site:** served with `vite preview`; `index.html`, the JS bundle, the CSS and
  the favicon all load under `/Boqian_Mao_BIOSTAT615_First_Project/`.

### Open items

- Push to GitHub and set Settings → Pages → Source to "GitHub Actions"; then confirm
  that CI passes and the placeholder page deploys. **Not verified yet.**
- Week 2: continuous endpoint (formulas, exact power, integer n search), PRNG and
  variates, the simulation engine in a Web Worker, Layers 2–3 for continuous.
