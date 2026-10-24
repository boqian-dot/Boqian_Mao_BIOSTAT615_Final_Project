/**
 * TrialSizer numerical core. Pure functions only: no React, no DOM (SPEC §6).
 * Function names follow the R functions each routine is tested against.
 */

export { pbeta, pbetaBoth } from './beta'
export { lbeta, lgamma, stirlerr } from './gamma'
export { pnt } from './noncentralT'
export { pnorm, pnormBoth, qnorm } from './normal'
export { dt, pt, qt } from './t'
