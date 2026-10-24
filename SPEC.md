# SPEC: TrialSizer — Real-Time Sample Size and Power for Two-Arm Clinical Trials

This document is the build specification for **TrialSizer**. It is addressed to the agentic coding tool that will build the project. It was produced from a structured interview with the author; every design decision and its rationale is recorded in the decision log at the end. When this spec is ambiguous, ask the author, then record the answer in the decision log.

## 0. Project facts

- **Author:** Boqian Mao
- **Course:** BIOSTAT 615 (Fall 2026), individual final project
- **Track:** A, interactive web app
- **Repository:** https://github.com/boqian-dot/Boqian_Mao_BIOSTAT615_First_Project (public)
- **Milestones:** proposal Mon Oct 26, 2026 · checkpoint (deployed, working MVP) Mon Nov 23, 2026 · final submission (repository + video) Thu Dec 17, 2026, 3:30 PM ET

## 1. Goal and users

**Goal.** Help clinical researchers, and other people with limited statistics training, quickly arrive at a reasonable, defensible sample size for a two-arm clinical trial. They should also understand what that number means and how fragile it is.

**Users.** Clinical researchers, residents, research coordinators, and students planning a trial or writing a grant. They know their clinical numbers (e.g., "control-group blood pressure is about 140; we hope to lower it by 10") but often not the SD, the right test, or what "80% power" means.

**What the app does.** The user describes the trial in clinical terms. The app returns the required sample size, computed three ways: a textbook formula, an exact computation, and a Monte Carlo simulation of the trial. It shows power, p-values, and confidence intervals from simulated trials in real time as inputs change, and it shows how the answer changes if the assumptions are wrong.

**Statistical-computing core.** The project implements in-house:
- the distribution functions (including the noncentral t)
- exact power for discrete tests by full enumeration, using an efficient margin-grouped algorithm
- a seeded simulation engine

All three are verified against R. The app makes visible where the textbook normal-approximation formulas fail, e.g., binary endpoints at small n or with rare events.

**Disclaimer (shown in the UI).** TrialSizer is a planning aid. It does not replace a statistician, and the final sample size for a real trial should be confirmed with one.

## 2. User workflow

The app follows the order in which researchers normally plan a trial: fix α and the target power, state the expected effect and variability, then solve for n. It has two modes; **Mode 1 is the default landing view.**

**Mode 1: Find the sample size (default).** The user enters:

- significance level α (default 0.05, two-sided; one-sided is available in the Advanced panel)
- target power (default 0.80)
- expected treatment effect, in clinical units (see "Effect inputs")
- variability (continuous endpoint only), with an SD helper for users who do not know the SD
- allocation ratio k (see §3, Notation)
- expected dropout rate

The app returns:
- the required number of analyzable participants per arm and in total
- the number to enroll after dropout inflation
- the three-way comparison (formula / exact / simulation)

**Mode 2: Power for a given sample size.** This mode is for studies whose size is fixed by budget or recruitment.
- The user enters the **total number of participants who can be enrolled**, N, plus the other design inputs.
- The app removes the expected dropout and splits the analyzable participants by the allocation ratio: n_C = ⌊N(1 − d) / (1 + k)⌋ and n_T = ⌊N(1 − d) · k / (1 + k)⌋.
- It returns the power, computed three ways.

### Effect inputs

Effects are entered **in clinical (raw) units**, because the target users think in clinical numbers, not in standardized effect sizes.

- **Continuous endpoint:**
  - Inputs: the expected difference Δ (treatment − control) and the SD, e.g., "lower by 10 mmHg, SD 20 mmHg". The control-group mean is an optional input used only for labels and plots.
  - The app shows the implied Cohen's d = Δ / SD **for reference only**; Cohen's d is never the primary input.
  - The sign of Δ does not affect two-sided power, but it matters for one-sided tests and for the wording of the justification paragraph.
- **Binary endpoint:** the event (or response) rate in the control group and in the treatment group, e.g., 30% vs. 20%. The app shows the absolute difference, relative risk, and odds ratio for reference.

### SD helper (continuous endpoint)

An "I don't know the SD" button opens a small panel that derives the SD from numbers that usually appear in papers or pilot data:

| User has | Formula | Note |
|----------|---------|------|
| SE of a group mean, with group size n | SD = SE × √n | |
| 95% CI for a group mean, with group size n | SD = √n × (upper − lower) / (2 × t₀.₉₇₅,ₙ₋₁) | the t quantile handles small pilot studies correctly |
| Typical range of values (min to max) | SD ≈ range / 4 | rough; the UI labels it as a rough guess |

The derived SD fills the SD field, and the user can still edit it.

## 3. Statistical method

**Design:** a two-arm, parallel-group superiority trial with a single primary endpoint.

| Endpoint | Effect measure | Default test | Advanced option | CI for the effect |
|----------|----------------|--------------|-----------------|-------------------|
| Continuous | difference in means Δ = μ_T − μ_C | Student's t test (common SD) | Welch's t test ("the two arms have different SDs") | t-based CI matching the chosen test |
| Binary | risk difference Δ = p_T − p_C | Pearson chi-square test, no continuity correction (≡ pooled z test) | Fisher's exact test | Newcombe hybrid score interval (Wilson-based) |

- **Defaults:** two-sided α = 0.05, target power 0.80.
- **Advanced options:** test choice, unequal SDs, one-sided α, number of simulation replicates, and RNG seed. They sit in a collapsed **Advanced** panel, so a novice only ever sees the defaults.
- **Binary CI:** the Newcombe interval is used instead of the Wald interval because Wald undercovers badly at small n.

### Notation

- Allocation ratio **k = n_T / n_C**. k = 1 means 1:1; k = 2 means two treated participants per control.
- z_q is the standard normal q-quantile, and 1 − β is the target power.
- Formulas below are written for a two-sided test, with z₁₋α/₂. A one-sided test uses z₁₋α instead.

### Closed-form sample-size formulas (used as the starting point)

**Continuous, common SD σ:**

n_C = (1 + 1/k) · σ² · (z₁₋α/₂ + z₁₋β)² / Δ²,  n_T = k · n_C

**Continuous, unequal SDs σ_C and σ_T (Welch):**

n_C = (σ_C² + σ_T² / k) · (z₁₋α/₂ + z₁₋β)² / Δ²,  n_T = k · n_C

**Binary (pooled z / chi-square),** with p̄ = (p_C + k·p_T) / (1 + k) and q = 1 − p:

n_C = [ z₁₋α/₂ · √( p̄·q̄·(1 + 1/k) ) + z₁₋β · √( p_C·q_C + p_T·q_T / k ) ]² / Δ²,  n_T = k · n_C

Sample sizes are always rounded **up** to whole numbers per arm, with n_T = ⌈k · n_C⌉.

### Exact power (no normal approximation)

- **Student's t test.** Use the noncentral t distribution with df = n_C + n_T − 2 and noncentrality δ = Δ / (σ·√(1/n_C + 1/n_T)). Power is P(T > t₁₋α/₂,df) + P(T < −t₁₋α/₂,df), counting both tails.
  - The reported n is the smallest integer n_C (with n_T = ⌈k·n_C⌉) that reaches the target power.
  - The search starts from the closed-form n.
- **Welch's t test.** Use the noncentral t approximation with
  - noncentrality δ = Δ / √(σ_C²/n_C + σ_T²/n_T)
  - Welch–Satterthwaite df ν = (σ_C²/n_C + σ_T²/n_T)² / [ (σ_C²/n_C)²/(n_C − 1) + (σ_T²/n_T)²/(n_T − 1) ]

  This is the reported value. The simulation is the reference that checks it.
- **Binary tests, by full enumeration.**
  - For fixed (n_C, n_T), each outcome (x_C, x_T) has probability Binom(x_C; n_C, p_C) · Binom(x_T; n_T, p_T), computed in log space.
  - The exact power of a test is the sum of these probabilities over the outcomes the test rejects. This gives the **true** power of both the chi-square test and Fisher's exact test, and the true type I error when p_T = p_C.
  - The chi-square test costs O(1) per table.
  - Fisher's test: group the tables by margin total m = x_C + x_T. Within a margin, the conditional distribution is hypergeometric and unimodal, so a two-pointer pass over tail cumulative sums gives the two-sided p-values of every table in that margin in O(support).
  - The whole computation is therefore O(n_C · n_T).
- **Sawtooth power.** Because outcomes are discrete, exact power for binary tests is **not monotone** in n. See the reporting rule below.

### Role of the three computations

Every result is computed three ways and shown side by side. **The exact computation is the authority.**

| Computation | Role | Shown as |
|-------------|------|----------|
| Closed-form formula | fast starting point for the search; what most textbooks and calculators report | "Formula: n = …, power = …" |
| Exact computation (noncentral t; full enumeration for binary) | **determines the reported sample size** | the headline number |
| Monte Carlo simulation | independent check of the exact result, and the raw material for the trial-level views (estimates, p-values, CIs of individual simulated trials) | "Simulation: power = … ± MC SE" |

**Why three computations.** Closed-form formulas rely on approximations: normality, and the large-sample normal approximation for proportions. The exact computation shows where those approximations fail. The simulation independently confirms the exact result and lets the user *see* individual trials. It is not a second copy of the formula.

**Simulation settings:**

- **Replicates:** default R = 10,000 per design, changeable in the Advanced panel.
- **Monte Carlo SE:** SE = √(p̂(1 − p̂) / R), always shown next to the simulated power (≈ 0.004 at power 0.80 with R = 10,000).
- **Seeded PRNG:** the seed is shown and editable, so a result can be reproduced exactly.
- **Patient-level data:** simulated by default, i.e., normal outcomes for continuous endpoints and binomial counts for binary endpoints.
  - **Performance fallback** if benchmarks show continuous simulation is too slow at large n: simulate the sufficient statistics directly instead (x̄ ~ N(μ, σ²/n) and s² ~ σ²·χ²ₙ₋₁/(n − 1), using Marsaglia–Tsang gamma variates).
- **Stored per replicate:** estimate, p-value, and CI, for the visualizations.

**Discrepancy messages:**

- **Formula vs. exact.** If the formula n differs from the exact n, or the formula power at the reported n differs from the exact power by more than 0.02, show a plain-language note, e.g., "The textbook formula underestimates the sample size for this design by 6 per arm."
- **Simulation vs. exact.** If the simulated power differs from the exact power by more than 3 MC SE, show a warning. This should not happen in normal operation, so it signals a bug. The test suite checks the same condition (§5).

**Reporting rule for sawtooth power (binary exact computation).** Report two numbers:

1. the smallest n whose exact power reaches the target, and
2. the smallest n from which exact power **stays** at or above the target for every larger n in the search window. The window runs up to max(2 × formula n, formula n + 50) per arm.

Number 2 is the headline ("recommended") value; number 1 is shown as a note. For continuous endpoints, power is monotone in n, so the two coincide.

## 4. Interface

- **Layout:** a single-page app in English. On desktop, the **input panel** is on the left and the **results area** on the right. On phones they stack, inputs first.
- **Disclaimer:** the §1 disclaimer is always visible (e.g., in the footer).

### Input panel

1. **Endpoint selector:** Continuous | Binary.
2. **Mode toggle:** "Find sample size" (default) | "Power for a given N".
3. **Design inputs (§2):**
   - α
   - target power (Mode 1) or total enrollable N (Mode 2)
   - effect inputs in clinical units
   - SD (continuous only), with the "I don't know the SD" button
   - allocation ratio: presets 1:1, 2:1, 3:1, plus a custom k
   - dropout rate
4. **Outcome name and unit** (e.g., "systolic blood pressure", "mmHg"), used in the justification paragraph and on axis labels.
5. **Advanced panel (collapsed):**
   - test choice (Student / Welch; chi-square / Fisher)
   - separate SDs per arm (Welch)
   - one-sided α
   - number of simulation replicates R
   - RNG seed

Every input has a one-line plain-language hint, e.g., "SD: how much individual patients typically vary around the group average." Invalid inputs (e.g., p_T = p_C, SD ≤ 0, rates outside (0, 1), Δ = 0 in Mode 1) produce an inline message. They never cause a crash or a NaN.

### Results area: four views (all in the MVP)

1. **Result card**
   - **Headline.** Mode 1: the recommended number of analyzable participants per arm and in total, plus the number to enroll after dropout inflation (n_enrolled = ⌈n / (1 − dropout)⌉ per arm). Mode 2: the power.
   - **Comparison row:** formula / exact / simulation (± MC SE), with any discrepancy or sawtooth notes (§3).
   - **Justification paragraph:** a copyable sample-size justification in English, ready to paste into a protocol or grant application. Example: *"A total of 128 participants (64 per arm, 1:1 allocation) provides 80% power to detect a difference of 10 mmHg in mean systolic blood pressure, assuming a common standard deviation of 20 mmHg, using a two-sided Student's t test at α = 0.05. Allowing for 10% dropout, 144 participants (72 per arm) will be enrolled."* The paragraph adapts to the endpoint, test, mode, and inputs, and has a "Copy" button.
2. **Power curve**
   - Plots power against n per arm (control arm), from about ¼ to 2× the recommended (or entered) n.
   - Formula power and exact power are two lines. A horizontal line marks the target power and a vertical marker marks the recommended (or entered) n.
   - The simulated power at that n is a point with ±2 MC SE error bars.
   - For binary exact tests, the sawtooth must be visible, never smoothed away.
   - Hovering shows the values.
3. **Sensitivity analysis:** what happens if the assumptions are wrong. A table with these rows:
   - SD × {0.8, 1.0, 1.2, 1.5} (continuous), or the control rate shifted by ±5 and ±10 percentage points with the absolute difference held fixed (binary; rates clipped to (0, 1))
   - effect × {0.8, 1.0, 1.2}
   - dropout {planned, planned + 10 percentage points}

   Each row shows (a) the n that would then be required and (b) the power the planned n would actually have. Rows where the planned n falls below the target power are highlighted. These rows use the same computation method as the headline.
4. **Simulated trials:** the trial-level view that explains what "80% power" means.
   - A CI plot of 100 simulated trials, taken from the R replicates: one horizontal CI per trial, colored significant / not significant, with reference lines at 0 and at the true Δ.
   - A caption, e.g., "In 81 of these 100 simulated trials, the result was statistically significant."
   - A histogram of p-values over all R replicates, with a line at α.
   - A "Draw another 100" button.

### "Real-time" behavior

- **Closed-form results** update on every input change, with no perceptible delay.
- **Exact computation and simulation** run in a Web Worker, off the main thread.
  - Inputs are debounced (~200 ms).
  - A new input cancels any stale job.
  - Each view shows a small "updating…" indicator until its fresh result arrives. Sensitivity rows may fill in progressively.
- **Performance target:** for typical designs (≤ 500 per arm), all four views refresh within about 1 s of the last input change on a mid-range laptop.

## 5. Correctness checks

The code is written by an AI agent, and the author is not fluent in TypeScript. So **correctness must be established by automated tests, not by reading code.** There are four layers, and CI runs all of them.

### Layer 1: numerical routines vs. R

`reference/generate_fixtures.R` evaluates R's own functions on a grid and writes JSON fixtures to `tests/fixtures/`. The fixtures are committed, so CI does not need R. Vitest compares the TypeScript routines against them.

| Routine | R reference | Grid (at least) | Tolerance |
|---------|-------------|-----------------|-----------|
| normal CDF / quantile | `pnorm`, `qnorm` | p ∈ [1e-12, 1 − 1e-12] | abs 1e-12 (central), rel 1e-9 (tails) |
| t CDF | `pt` | df ∈ {1, 2, 5, 10, 30, 100, 1000} | abs 1e-10 |
| t quantile | `qt` | same df, p ∈ [1e-8, 1 − 1e-8] | rel 1e-9 |
| noncentral t CDF | `pt(q, df, ncp)` | df as above × ncp ∈ {0, 0.5, 1, 2, 5, 10} | abs 1e-8 |
| log binomial / hypergeometric pmf | `dbinom(log = TRUE)`, `dhyper(log = TRUE)` | n up to 5,000 | rel 1e-10 |

### Layer 2: power, sample size, tests, and CIs vs. R

- **Continuous, Student, k = 1.** Power and n must match `power.t.test(..., strict = TRUE)`; `strict = TRUE` makes R count both tails, as this app does. The reported n must equal `ceiling(power.t.test(power = …, strict = TRUE)$n)` exactly.
- **Continuous, any k, and Welch.** The R script computes the reference power directly with `pt(…, ncp = …)`, using the same df and noncentrality. `power.t.test()` handles only equal n.
- **Binary, closed form, k = 1.** Must match `power.prop.test()`, which uses the same pooled formula.
- **Binary, exact power, by brute force in R.**
  - For every table (x_C, x_T), call `fisher.test()` and `chisq.test(correct = FALSE)` (with warnings suppressed).
  - Sum `dbinom(x_C, n_C, p_C) * dbinom(x_T, n_T, p_T)` over the rejected tables.
  - This is slow but obviously correct, and it is the reference for the fast margin-grouped TypeScript algorithm.
  - Exact power must agree to abs 1e-8. Reported sample sizes must be identical, including both sawtooth numbers.
- **Fisher two-sided p-value definition.** Follow R's rule: sum the probabilities of all tables whose hypergeometric probability is ≤ the observed probability × (1 + 1e-7). Without this rule, p-values for tied tables will not match R.
- **Fixture designs.** At least 40, spanning:
  - k ∈ {1, 2, 3} and α ∈ {0.05, 0.01}
  - one- and two-sided tests
  - n from 5 to 100 per arm (brute force in R is slow)
  - rare events (p ∈ [0.02, 0.10])
  - at least three designs where the sawtooth makes the two reported n differ
- **CIs.**
  - t-based CIs must match `t.test(…, var.equal = TRUE / FALSE)$conf.int` on fixed small datasets.
  - The Newcombe CI must match an R computation built from the Wilson intervals of `prop.test(x, n, correct = FALSE)$conf.int`, combined by Newcombe's (1998) hybrid score method.
- **SD helper and dropout arithmetic** (both modes): unit tests with hand-computed values.
- **Justification paragraph:** snapshot tests for every endpoint × mode × test combination.

### Layer 3: statistical properties of the simulator

These tests are stochastic but use **fixed seeds**, so CI is deterministic.

- **Type I error.** With Δ = 0 and normal data, the Student-t rejection rate is within α ± 3 MC SE. For the binary tests, the simulated rejection rate is within 3 MC SE of the **exact** type I error from enumeration, which is not α at small n.
- **Power.** On the fixture designs, simulated power is within 3 MC SE of the exact power.
- **CI coverage.** t-CI coverage on normal data is within 0.95 ± 3 MC SE. Newcombe coverage at moderate n (≥ 30 per arm, rates in [0.1, 0.9]) lies in [0.93, 0.97].
- **Random variates.**
  - uniform PRNG: mean and variance checks
  - normal variates: a Kolmogorov–Smirnov test against Φ on 10⁵ draws
  - binomial variates (and gamma variates, if used): mean, variance, and a goodness-of-fit test against R's `dbinom` / `pgamma`
- **Reproducibility.** The same seed and inputs give bit-identical results.

### Layer 4: published worked examples

- **Required sanity check:** Δ = 10, σ = 20, α = 0.05 two-sided, power 0.80 → **64 per arm**.
- **Published examples:** at least three worked examples from a standard reference, e.g., Chow, Shao, Wang & Lokhnygina, *Sample Size Calculations in Clinical Research* (3rd ed.), with page numbers cited in the test.
  - At least one continuous-endpoint example and one binary-endpoint example.
  - Where the book uses a different formula variant (e.g., a continuity correction), the test documents the variant and checks the app's corresponding computation.

### Continuous integration

- A GitHub Actions workflow runs lint, type-check, and every Vitest suite on each push and pull request.
- Deployment to GitHub Pages runs only when they pass.
- A second, manually triggered workflow (`workflow_dispatch`) may regenerate the fixtures with R.

## 6. Tech stack and hosting

- **Language:** TypeScript, strict mode.
- **Framework:** React + Vite. A static single-page app with no backend.
- **Charts:** Recharts. Only the CI plot may switch to a hand-written SVG component, if Recharts cannot draw it cleanly.
- **Numerics:** implemented in-house in `src/stats/`, with no general-purpose stats dependency. Common JS stats libraries lack the noncentral t, and writing these routines is part of the statistical-computing substance. Required routines:
  - normal CDF and quantile (e.g., Wichura's AS 241 for the quantile)
  - log-gamma (Lanczos), and log-space binomial and hypergeometric probabilities, so large n does not underflow
  - Student t CDF and quantile
  - noncentral t CDF (e.g., Lenth's AS 243)
  - a seeded PRNG (e.g., xoshiro128** or PCG32)
  - normal variates (Box–Muller or the polar method)
  - binomial variates (inversion for small n·p; otherwise a standard exact algorithm such as BTPE)
  - gamma variates (Marsaglia–Tsang), only if the sufficient-statistics fallback is used
  - integer search and root finding for the sample size

  Each routine cites its algorithm in a comment and is tested against R (§5).
- **Concurrency:** one Web Worker for the exact computation and the simulation (§4).
- **Stretch only:** KaTeX for rendering formulas on the Methods page.
- **Hosting:** GitHub Pages, deployed by GitHub Actions on every push to `main` once checks pass. Zero cost, no server, no cold start.
- **Local development:** `npm ci`, `npm run dev`, `npm test`. Document these in the README.

### Repository layout

```
/                       README.md, LICENSE (MIT unless the author decides otherwise), AI_USAGE.md, SPEC.md
src/stats/              distributions, PRNG, formulas, exact power, simulation (pure functions, no UI)
src/worker/             Web Worker entry point and typed messages
src/components/         input panel, SD helper, the four result views
src/text/               justification-paragraph templates
tests/                  Vitest suites and JSON reference fixtures
reference/              R scripts that generate the fixtures
docs/dev-log.md         development log (see §10)
.github/workflows/      CI and deployment to GitHub Pages
```

`src/stats/` must not import React or touch the DOM, so it can be tested and reasoned about in isolation.

## 7. Scope

### MVP: deployed and working by the Nov 23 checkpoint

- Continuous and binary endpoints, with the default and advanced tests (§3).
- Both modes: find the sample size, and power for a given N (§2).
- Effect inputs in clinical units, the SD helper, the allocation ratio, and dropout handling (§2).
- Formula, exact, and simulation results side by side, with discrepancy and sawtooth messages (§3).
- All four result views (§4).
- Test layers 1–3 passing in CI, with deployment gated on them (§5).
- Live on GitHub Pages.

**Minimum bar if behind schedule.** The checkpoint must not fall below this: both endpoints in Mode 1, with the result card and the power curve, test layers 1–2 passing, and the app deployed. Anything missing beyond that is reported honestly in the checkpoint report.

### Stretch goals, in priority order

Start these only after the MVP passes layers 1–3 for both endpoints.

1. **Shareable link and Methods page.**
   - All inputs, including the seed, are encoded in the URL hash, which works on GitHub Pages. A "Copy link" button copies the link. Opening a link restores the exact design and results.
   - A Methods page lists every formula (in KaTeX), the algorithms used, references, limitations, and the disclaimer.
2. **Non-inferiority and equivalence designs**, reusing the exact and simulation engines.
   - **Non-inferiority:** H₀: Δ ≤ −M vs. H₁: Δ > −M, with a user-specified margin M, a direction ("higher is better" / "lower is better"), and one-sided α (default 0.025). The expected true Δ defaults to 0.
   - **Equivalence:** two one-sided tests (TOST) with margins ±M.
   - **Tests:** continuous designs use shifted t tests. Binary non-inferiority uses the Farrington–Manning score test, whose exact power comes from the same enumeration engine.
   - **Verification:** where no closed-form exact power exists (e.g., continuous TOST), simulation is the reference. Also validate against published worked examples.

### Out of scope

Not built, not even as stretch goals:

- time-to-event endpoints
- interim / group-sequential analyses
- more than two arms
- cluster-randomized and crossover designs
- covariate adjustment
- Bayesian designs
- user accounts and any backend

## 8. Milestones

| Week | Dates (2026) | Deliverable |
|------|--------------|-------------|
| 1 | Oct 27 – Nov 2 | Scaffold: Vite + React + TS strict, ESLint, Vitest, and a GitHub Actions workflow that deploys a placeholder page to GitHub Pages. `src/stats/` numerical routines and `reference/generate_fixtures.R`. **Layer 1 tests passing.** Start `docs/dev-log.md`. |
| 2 | Nov 3 – 9 | Continuous endpoint: formulas, exact power (noncentral t), integer n search, seeded simulation engine in the Web Worker. Layers 2–3 for continuous. Minimal UI: input panel and result card. |
| 3 | Nov 10 – 16 | Binary endpoint: formula, margin-grouped exact enumeration (chi-square, Fisher), sawtooth rule, Newcombe CI, simulation. Layers 2–3 for binary. Power curve view. **Benchmark exact enumeration at large n** (risk R2). |
| 4 | Nov 17 – 23 | Sensitivity view, simulated-trials view, justification paragraph, SD helper, Mode 2. **Checkpoint (Nov 23): deployed MVP + report.** Decide whether to write the agentic-coding tutorial (§10). |
| 5 | Nov 24 – 30 | Address checkpoint feedback. Usability pass with at least one non-statistician. Layer 4 (published examples). |
| 6 | Dec 1 – 7 | Stretch 1 (shareable link + Methods page), then stretch 2 (non-inferiority / equivalence). |
| 7 | Dec 8 – 17 | README (live link, screenshots, how to run, how the fixtures are made, credits), LICENSE, final `AI_USAGE.md`, video. **Final submission: Thu Dec 17, 3:30 PM ET.** |

## 9. Risks and mitigations

| # | Risk | Mitigation |
|---|------|------------|
| R1 | The in-house noncentral t CDF is inaccurate for large noncentrality or extreme df. | Use a published algorithm (AS 243) and Layer 1 tests against `pt(ncp = )`. If accuracy cannot be met in some region, restrict inputs there (power ≈ 1 there anyway) and document it. |
| R2 | Exact binary enumeration costs O(n_C · n_T) per design. The sawtooth search and the sensitivity rows multiply that, so large n may be too slow for "real time". | Use exact enumeration only when n_C · n_T ≤ 10⁶ (≈ 1,000 per arm at 1:1). The sawtooth search may need a lower cap; set it from the week-3 benchmark. Above the cap, report the formula n checked by simulation, clearly labelled; the normal approximation is accurate there. For sensitivity rows, a single exact power evaluation at the planned n is cheap; only the required-n searches are expensive, and they may fill in progressively. |
| R3 | The two sawtooth numbers confuse novice users. | The headline shows a single recommended n. The other number appears as a plain-language note, and the power curve shows why. |
| R4 | Garbage in, garbage out: users guess the effect or SD badly. | Sensitivity view, SD helper, input hints, and the visible disclaimer (§1). |
| R5 | The author cannot verify AI-written TypeScript by reading it. | R reference fixtures and property tests (§5). `src/stats/` is kept pure and commented with algorithm citations. |
| R6 | Fisher two-sided p-values disagree with R on ties. | Use R's relative-tolerance rule (§5, Layer 2). |
| R7 | Scope creep before the MVP is verified. | Stretch work starts only after layers 1–3 pass for both endpoints. |
| R8 | Recharts cannot draw the CI plot cleanly. | Switch that one view to a hand-written SVG component. |

## 10. Working rules for the coding agent

1. **Build in milestone order (§8).** Write tests together with each `src/stats/` routine, never after the UI.
2. **Never weaken a test to make it pass.** Do not loosen a tolerance, edit a fixture by hand, or skip a test without the author's approval. Record every such decision and its reason in `docs/dev-log.md`.
3. **Make small, focused commits** with descriptive messages. Run lint, type-check, and tests before every commit.
4. **Keep the development log.** Append one entry to `docs/dev-log.md` for each working session with the agent. Each entry records:
   - date and goal
   - key prompts (verbatim where they shaped the result)
   - what the agent changed (files or commits)
   - what went wrong and how it was corrected
   - how the result was verified

   The log is the raw material for `AI_USAGE.md`. That document must describe the real workflow honestly, including failures and course corrections.
5. **The agentic-coding tutorial is optional.** Whether to write one (a step-by-step document or GitHub Pages site that lets a reader reproduce the development process) is decided at the Nov 23 checkpoint, based on progress. The dev log is kept so that this option stays open.
6. **Credit third-party work.** Credit any code, algorithm, or data you build on in the README, and respect its license.
7. **Ask when the spec is ambiguous, and record the answer.** Ask the author rather than guessing, then add the answer to the decision log below.

## Decision log

| # | Question | Decision | Rationale |
|---|----------|----------|-----------|
| 1 | Core workflow: what does the user enter, and what does the app return? | Solve for n by default (Mode 1), switchable to "power for a given N" (Mode 2). | Matches how researchers actually plan: α and power are fixed first, then the effect and SD, then n. Directly serves the goal of a reasonable sample size. Mode 2 covers fixed-budget studies. |
| 2 | Which endpoint types are in the MVP? | Continuous and binary. (Survival was initially listed as a stretch goal; see decision 9.) | Covers most two-arm trials. Binary endpoints give exact computation and simulation a real role, because formula approximations fail at small n and for rare events. Continuous-only risks being "a wrapper around `power.t.test`"; survival is too much work before Nov 23. |
| 3 | How do users enter the effect and variability? | Clinical units: difference + SD for continuous, two event rates for binary. An SD helper derives the SD from an SE, a 95% CI, or a range. Cohen's d is shown for reference only. | Target users know clinical numbers but often not the SD. Standardized "small / medium / large" presets lose clinical meaning and are widely criticized for trial planning. |
| 4 | Which tests, and which defaults? | Two-sided α = 0.05, power 0.80. Continuous: Student's t, with Welch in the Advanced panel. Binary: chi-square without correction, with Fisher's exact test in the Advanced panel. The Advanced panel is collapsed by default. Newcombe CI for binary endpoints. | Novices get sensible defaults without facing choices they do not understand. Offering Fisher's exact test keeps the approximate-vs.-exact comparison available. Newcombe has better small-n coverage than Wald. |
| 5 | What are the roles of the formula, the exact computation, and the simulation? | All three are shown side by side. The exact computation sets the reported n. The simulation (R = 10,000, seeded, with MC SE) checks it independently and feeds the trial-level views. Discrepancies are flagged. Sawtooth rule: the headline is the smallest n from which power stays ≥ the target. | Gives the simulation a real job and makes it visible where the textbook formula fails, which is the statistical-computing core of the project. The exact result has no Monte Carlo noise, so the reported n is stable. |
| 6 | Which outputs are in the MVP? | All four: the result card with the copyable justification paragraph, the power curve, the sensitivity analysis, and the simulated-trials view. | Each answers a different user question: "How many?", "How does n trade off against power?", "What if my guesses are wrong?", "What does 80% power mean?" |
| 7 | Which framework and hosting? | TypeScript + React + Vite, static on GitHub Pages. Simulation in a Web Worker. Numerics written in-house and tested against R. | Only fully in-browser JS meets the "real-time" requirement with no server, no cold start, and zero cost. Writing the distributions is statistical-computing work that R reference values can verify. Trade-off accepted: the author is less fluent in JS than in R. |
| 8 | How is correctness checked? | Four automated layers: numerics vs. R; power, n, tests, and CIs vs. R (including brute-force exact power in R); statistical properties of the simulator (seeded); published worked examples. All run in GitHub Actions and gate deployment. | The author cannot verify AI-written TypeScript line by line, and R is the trusted reference the author knows. Property tests catch simulator bugs that point comparisons miss. |
| 9 | Which stretch goals? | (1) Shareable link + Methods page, then (2) non-inferiority / equivalence. Survival and interim analyses are out of scope; this supersedes the survival note in decision 2. | Highest value per unit of effort for the target users, and non-inferiority / equivalence reuses the existing engines. Survival and group-sequential designs are too large for the ~3.5 weeks after the checkpoint. |
| 10 | What is the project called? | **TrialSizer**: "Real-Time Sample Size and Power for Two-Arm Clinical Trials". | Short and memorable, and it says what the tool does. |
| 11 | How is AI usage documented? | A structured `docs/dev-log.md` entry per session, from day 1; `AI_USAGE.md` is compiled from it. Whether to write the optional agentic-coding tutorial is decided at the Nov 23 checkpoint. | `AI_USAGE.md` is worth 15 points and must be honest and specific, which is hard to reconstruct after the fact. The log keeps the tutorial option open without committing a week to it now. |
