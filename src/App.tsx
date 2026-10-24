import { qnorm } from './stats'

// Placeholder page for the week-1 deployment. The real input panel and result
// views replace this from week 2 on (SPEC §8).
export default function App() {
  return (
    <div className="page">
      <header>
        <h1>TrialSizer</h1>
        <p className="subtitle">Real-time sample size and power for two-arm clinical trials</p>
      </header>

      <main>
        <p>
          TrialSizer is under construction. It will compute the sample size for a two-arm trial
          three ways (a textbook formula, an exact computation, and a Monte Carlo simulation) and
          show how the answer changes when the assumptions are wrong.
        </p>
        <p className="build-check">
          Numerical core check: z<sub>0.975</sub> = {qnorm(0.975).toFixed(6)}
        </p>
      </main>

      <footer>
        TrialSizer is a planning aid. It does not replace a statistician, and the final sample size
        for a real trial should be confirmed with one.
      </footer>
    </div>
  )
}
