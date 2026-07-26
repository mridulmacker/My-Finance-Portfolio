# Numerical Methods for Derivative Pricing — A Comparative Study

A self-contained quantitative-finance project that prices one reference contract
(6-month ATM equity option; S₀ = K = 35, σ = 45%, r = 3%, q = 4%) with **seven
numerical techniques**, then uses the structure of the cross-model agreements and
disagreements to draw quantitative conclusions about each method.

## Contents
- `report.pdf` — 9-page research note with all results, figures, and analysis
- `pricing_library.py` — single dependency-light module (NumPy/SciPy) that
  reproduces every number in the report; run `python pricing_library.py`

## Headline results (all seeded and reproducible)

| Quantity | Method | Value |
|---|---|---|
| European call | Black–Scholes (analytic) | 4.2621 |
| | CRR two-point average, n = 500 | 4.26211 (error 2×10⁻⁵) |
| | Implicit FD, 2000×1000 | 4.2615 |
| | MC antithetic, M = 200k | 4.266 ± 0.014 |
| American call | CRR n = 5,000 (benchmark) | 4.2879 |
| | Projected implicit FD | 4.2873 |
| | Longstaff–Schwartz sandwich | [4.278, 4.315] |
| American floating lookback | 3-step path tree | 5.980 |
| Down-and-out call (continuous) | closed form / bridge MC | 2.348 / 2.345 ± 0.015 |
| Down-and-out call (daily monitored) | MC / BGK closed form | 2.633 ± 0.015 / 2.647 |

## What makes this more than a pricing exercise

1. **Binomial oscillation exploited, not endured.** The CRR error oscillates with
   the parity of n (an O(1/n) envelope); averaging consecutive step counts cancels
   it, cutting the error by ~100× at no extra cost.
2. **Barrier "MC error" diagnosed as monitoring bias.** The 12% gap between a
   daily-monitored simulation and the continuous closed form is bias, not noise —
   confirmed two independent ways (Broadie–Glasserman–Kou barrier shift matches the
   discrete price; a Brownian-bridge estimator recovers the continuous one).
3. **Level vs. dynamics separated with a controlled experiment.** With common
   random numbers and the initial rate held at the baseline, stochastic interest
   rates move the price by < 0.1 cents while a rate-*level* change of 2 points moves
   it by 15 cents — quantifying which refinement actually matters at this horizon.
4. **CEV calibrated before compared.** The elasticity α is shown to be a *skew*
   parameter, not a level parameter: with ATM local vol matched across α, the ATM
   price is nearly invariant while implied vols span 48.1% → 43.3% across strikes
   28–42 (the equity skew).
5. **Both Longstaff–Schwartz biases exhibited.** The in-sample estimator (foresight
   bias, high) and out-of-sample estimator (policy suboptimality, low) bracket the
   lattice benchmark at every polynomial order — and a linear basis is shown to
   collapse the American premium entirely.
6. **Every scheme validated before use.** The FD code is certified against the
   European closed form (error 6×10⁻⁴) before its American output is trusted; every
   Monte Carlo estimate carries a standard error.

## Method inventory
Black–Scholes closed forms + Greeks · CRR binomial (European & American) ·
non-recombining path tree for a path-dependent American lookback · Monte Carlo with
antithetic variates and moment matching · correlated stochastic-short-rate MC with
pathwise discounting · barrier pricing (closed form, BGK correction, Brownian
bridge) · calibrated CEV + implied-vol inversion · Heston-CEV stochastic volatility
· implicit finite differences (Thomas solver, American projection) ·
Longstaff–Schwartz (two-pass).

*Originally developed from a Derivative Securities coursework project; fully
rebuilt, extended, and independently verified.*
