"""
Numerical Methods for Derivative Pricing — companion library.

Prices a common reference contract (S0 = K = 35, sigma = 45%, r = 3%, q = 4%,
T = 0.5y) with seven techniques, cross-validating each against a closed form
or a converged benchmark:

    black_scholes         analytic benchmark + Greeks + put-call parity
    crr_tree              CRR binomial tree (European / American) + two-point averaging
    lookback_tree         American floating-strike lookback on a 3-step path tree
    mc_european           exact-terminal Monte Carlo: plain / antithetic / moment matching
    mc_stochastic_rates   correlated GBM short rate, pathwise stochastic discounting
    barrier_pricing       down-and-out call: closed form, discrete MC, BGK shift,
                          Brownian-bridge correction
    cev_pricing           CEV with ATM-vol calibration + implied-vol skew extraction
    heston_cev            CEV diffusion with Heston variance (full-truncation Euler)
    fd_implicit           implicit finite differences, Thomas solver, American projection
    lsm                   Longstaff-Schwartz with in-sample (high-biased) and
                          out-of-sample (low-biased) estimators

Design principles:
  * every Monte Carlo result carries a standard error;
  * comparisons across model variants use common random numbers (shared seed);
  * every scheme is validated on a case with a known answer before being
    trusted on the case without one.

Requires: numpy, scipy.  Run `python pricing_library.py` to reproduce the
headline numbers of the report.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded
from scipy.special import comb
from scipy.stats import norm

# ----------------------------------------------------------------------------
# Baseline contract
# ----------------------------------------------------------------------------
S0, K, SIGMA, R, Q, T = 35.0, 35.0, 0.45, 0.03, 0.04, 0.5


# ----------------------------------------------------------------------------
# 1. Black-Scholes benchmark and Greeks
# ----------------------------------------------------------------------------
def black_scholes(S0=S0, K=K, sigma=SIGMA, r=R, q=Q, T=T) -> dict:
    """Analytic call/put values and Greeks under BSM with dividend yield q."""
    sqT = np.sqrt(T)
    d1 = (np.log(S0 / K) + (r - q + 0.5 * sigma**2) * T) / (sigma * sqT)
    d2 = d1 - sigma * sqT
    disc_q, disc_r = np.exp(-q * T), np.exp(-r * T)
    call = S0 * disc_q * norm.cdf(d1) - K * disc_r * norm.cdf(d2)
    put = K * disc_r * norm.cdf(-d2) - S0 * disc_q * norm.cdf(-d1)
    nd1 = norm.pdf(d1)
    return {
        "call": call,
        "put": put,
        "d1": d1,
        "d2": d2,
        "delta_call": disc_q * norm.cdf(d1),
        "delta_put": -disc_q * norm.cdf(-d1),
        "gamma": disc_q * nd1 / (S0 * sigma * sqT),
        "vega": S0 * disc_q * nd1 * sqT,  # per unit of vol
        "theta_call": (-S0 * disc_q * nd1 * sigma / (2 * sqT)
                       + q * S0 * disc_q * norm.cdf(d1)
                       - r * K * disc_r * norm.cdf(d2)),  # per year
        "rho_call": K * T * disc_r * norm.cdf(d2),
        "parity_gap": (call - put) - (S0 * disc_q - K * disc_r),  # should be ~0
    }


# ----------------------------------------------------------------------------
# 2. CRR binomial tree
# ----------------------------------------------------------------------------
def crr_tree(n: int, american: bool = False, **kw) -> float:
    """CRR call value. European uses the closed binomial sum; American uses
    backward induction with an early-exercise check at every node."""
    S0_, K_, sigma, r, q, T_ = (kw.get(k, v) for k, v in
                                zip("S0 K sigma r q T".split(),
                                    (S0, K, SIGMA, R, Q, T)))
    dt = T_ / n
    u = np.exp(sigma * np.sqrt(dt))
    d = 1.0 / u
    p = (np.exp((r - q) * dt) - d) / (u - d)
    if not american:
        j = np.arange(n + 1)
        ST = S0_ * u**j * d ** (n - j)
        w = comb(n, j) * p**j * (1 - p) ** (n - j)
        return float(np.exp(-r * T_) * np.sum(w * np.maximum(ST - K_, 0.0)))
    disc = np.exp(-r * dt)
    j = np.arange(n + 1)
    S = S0_ * u**j * d ** (n - j)
    V = np.maximum(S - K_, 0.0)
    for i in range(n - 1, -1, -1):
        S = S[1: i + 2] / u
        V = disc * (p * V[1: i + 2] + (1 - p) * V[: i + 1])
        V = np.maximum(V, S - K_)  # early-exercise projection
    return float(V[0])


def crr_averaged(n: int) -> float:
    """Two-point average 0.5*[C(n)+C(n+1)]: cancels the odd-even oscillation,
    accelerating convergence by roughly two orders of magnitude."""
    return 0.5 * (crr_tree(n) + crr_tree(n + 1))


# ----------------------------------------------------------------------------
# 3. American floating-strike lookback (3-step path tree)
# ----------------------------------------------------------------------------
def lookback_tree(n: int = 3) -> float:
    """American floating-strike lookback call, exercise value S_t - min(S).
    The running minimum breaks recombination, so the state is (S, Smin) and
    the tree is traversed path-wise (2^n paths — fine for small n)."""
    dt = T / n
    u = np.exp(SIGMA * np.sqrt(dt))
    d = 1.0 / u
    p = (np.exp((R - Q) * dt) - d) / (u - d)
    disc = np.exp(-R * dt)

    def value(t, S, Smin):
        exercise = S - Smin
        if t == n:
            return max(exercise, 0.0)
        cont = disc * (p * value(t + 1, S * u, min(Smin, S * u))
                       + (1 - p) * value(t + 1, S * d, min(Smin, S * d)))
        return max(cont, exercise)

    return value(0, S0, S0)


# ----------------------------------------------------------------------------
# 4a. Plain Monte Carlo with variance-reduction variants
# ----------------------------------------------------------------------------
def mc_european(M: int, mode: str = "plain", seed: int = 0) -> tuple[float, float]:
    """European call by exact terminal sampling (no path discretisation).
    mode: 'plain' | 'avt' (antithetic) | 'mm' (moment matching).
    Returns (estimate, standard_error). For AVT the SE is computed on the
    pair averages, which is the correct (correlation-aware) formula."""
    g = np.random.default_rng(seed)
    if mode == "avt":
        Z = g.standard_normal(M // 2)
        Z = np.concatenate([Z, -Z])
    else:
        Z = g.standard_normal(M)
        if mode == "mm":
            Z = (Z - Z.mean()) / Z.std()
    ST = S0 * np.exp((R - Q - 0.5 * SIGMA**2) * T + SIGMA * np.sqrt(T) * Z)
    pay = np.exp(-R * T) * np.maximum(ST - K, 0.0)
    if mode == "avt":
        pairs = 0.5 * (pay[: M // 2] + pay[M // 2:])
        return pairs.mean(), pairs.std(ddof=1) / np.sqrt(M // 2)
    return pay.mean(), pay.std(ddof=1) / np.sqrt(M)


# ----------------------------------------------------------------------------
# 4b. Stochastic short rate, correlated with the stock
# ----------------------------------------------------------------------------
def mc_stochastic_rates(M: int, r0: float = R, mu_r: float = 0.02,
                        sigma_r: float = 0.03, rho: float = 0.1,
                        n_steps: int = 126, seed: int = 0) -> tuple[float, float]:
    """European call with a GBM short rate correlated with the stock; each
    payoff is discounted along its own rate path.  Keep r0 equal to the
    deterministic baseline to isolate the effect of rate *randomness* from a
    rate *level* change (the controlled-experiment point of the report).
    Setting sigma_r = 0 gives the matched deterministic-rate control."""
    g = np.random.default_rng(seed)
    dt = T / n_steps
    S = np.full(M, S0)
    rt = np.full(M, r0)
    integral = np.zeros(M)
    for _ in range(n_steps):
        z1, z2 = g.standard_normal(M), g.standard_normal(M)
        dWr = np.sqrt(dt) * z1
        dWS = rho * dWr + np.sqrt(1 - rho**2) * np.sqrt(dt) * z2
        integral += rt * dt
        S = S * np.exp((rt - Q - 0.5 * SIGMA**2) * dt + SIGMA * dWS)
        rt = rt * np.exp((mu_r - 0.5 * sigma_r**2) * dt + sigma_r * dWr)
    pay = np.exp(-integral) * np.maximum(S - K, 0.0)
    return pay.mean(), pay.std(ddof=1) / np.sqrt(M)


# ----------------------------------------------------------------------------
# 4c. Down-and-out barrier: closed form, discrete MC, BGK, Brownian bridge
# ----------------------------------------------------------------------------
def barrier_closed_form(B: float = 32.0, S0_=S0, K_=K, sigma=SIGMA,
                        r=R, q=Q, T_=T) -> float:
    """Continuous-monitoring down-and-out call, B <= K (Hull Ch. 26).
    Note y1 uses B^2/(S0*K) — a frequently misquoted argument."""
    lam = (r - q + 0.5 * sigma**2) / sigma**2
    sqT = np.sqrt(T_)
    x1 = np.log(S0_ / K_) / (sigma * sqT) + lam * sigma * sqT
    y1 = np.log(B**2 / (S0_ * K_)) / (sigma * sqT) + lam * sigma * sqT
    return (S0_ * np.exp(-q * T_) * norm.cdf(x1)
            - K_ * np.exp(-r * T_) * norm.cdf(x1 - sigma * sqT)
            - S0_ * np.exp(-q * T_) * (B / S0_) ** (2 * lam) * norm.cdf(y1)
            + K_ * np.exp(-r * T_) * (B / S0_) ** (2 * lam - 2)
            * norm.cdf(y1 - sigma * sqT))


def barrier_bgk_shift(B: float, n_monitor: int, sigma=SIGMA) -> float:
    """Broadie-Glasserman-Kou (1997): a discretely monitored down barrier at B
    prices like a *continuous* barrier at B*exp(-0.5826*sigma*sqrt(dt))."""
    return B * np.exp(-0.5826 * sigma * np.sqrt(T / n_monitor))


def barrier_mc(M: int, n: int, B: float = 32.0, bridge: bool = False,
               seed: int = 0) -> tuple[float, float]:
    """Down-and-out call MC.  bridge=False: knock-out checked only at the n
    grid dates (discrete monitoring — biased high vs. the continuous form).
    bridge=True: paths are additionally killed with the Brownian-bridge
    crossing probability between grid points, restoring continuous
    monitoring inside a discrete simulation."""
    g = np.random.default_rng(seed)
    dt = T / n
    logB = np.log(B)
    logS = np.full(M, np.log(S0))
    alive = np.ones(M, bool)
    for _ in range(n):
        Z = g.standard_normal(M)
        logS_new = logS + (R - Q - 0.5 * SIGMA**2) * dt + SIGMA * np.sqrt(dt) * Z
        if bridge:
            p_cross = np.exp(-2 * (logS - logB) * (logS_new - logB)
                             / (SIGMA**2 * dt))
            p_cross = np.where((logS > logB) & (logS_new > logB),
                               np.clip(p_cross, 0, 1), 1.0)
            alive &= g.random(M) > p_cross
        else:
            alive &= logS_new > logB
        logS = logS_new
    pay = np.exp(-R * T) * np.where(alive, np.maximum(np.exp(logS) - K, 0.0), 0.0)
    return pay.mean(), pay.std(ddof=1) / np.sqrt(M)


# ----------------------------------------------------------------------------
# 5. CEV — calibrate the volatility LEVEL before comparing elasticities
# ----------------------------------------------------------------------------
def cev_pricing(M: int, n: int, alpha: float, calibrate: bool = True,
                barrier: float | None = None, seed: int = 0) -> tuple[float, float]:
    """CEV call via Euler.  With calibrate=True, sigma_cev = SIGMA*S0^(1-alpha)
    so that every alpha shares the same 45% ATM instantaneous volatility —
    without this, changing alpha changes the vol level (0.45*35^-0.5 = 7.6%
    at alpha = 0.5) and the comparison is meaningless.  Use a common seed
    across alphas so *differences* are estimated with high precision."""
    sigma_cev = SIGMA * S0 ** (1 - alpha) if calibrate else SIGMA
    g = np.random.default_rng(seed)
    dt = T / n
    S = np.full(M, S0)
    alive = np.ones(M, bool)
    for _ in range(n):
        Z = g.standard_normal(M)
        S = np.maximum(S + (R - Q) * S * dt
                       + sigma_cev * np.maximum(S, 1e-8) ** alpha
                       * np.sqrt(dt) * Z, 1e-8)
        if barrier is not None:
            alive &= S > barrier
    pay = np.maximum(S - K, 0.0)
    if barrier is not None:
        pay = np.where(alive, pay, 0.0)
    pay = np.exp(-R * T) * pay
    return pay.mean(), pay.std(ddof=1) / np.sqrt(M)


def cev_implied_skew(alpha: float = 0.5, strikes=(28, 30, 32, 34, 35, 36, 38, 40, 42),
                     M: int = 200_000, n: int = 252, seed: int = 55) -> dict:
    """Price calibrated-CEV calls across strikes and invert Black-Scholes to
    exhibit the implied-volatility skew that the elasticity generates."""
    g = np.random.default_rng(seed)
    dt = T / n
    sigma_cev = SIGMA * S0 ** (1 - alpha)
    S = np.full(M, S0)
    for _ in range(n):
        Z = g.standard_normal(M)
        S = np.maximum(S + (R - Q) * S * dt
                       + sigma_cev * np.sqrt(np.maximum(S, 1e-8))
                       * np.sqrt(dt) * Z, 1e-8)

    def bs_call(Kx, s):
        d1 = (np.log(S0 / Kx) + (R - Q + 0.5 * s**2) * T) / (s * np.sqrt(T))
        return (S0 * np.exp(-Q * T) * norm.cdf(d1)
                - Kx * np.exp(-R * T) * norm.cdf(d1 - s * np.sqrt(T)))

    def implied_vol(price, Kx, lo=0.01, hi=2.0):
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if bs_call(Kx, mid) < price else (lo, mid)
        return 0.5 * (lo + hi)

    out = {}
    for Kx in strikes:
        price = np.exp(-R * T) * np.maximum(S - Kx, 0.0).mean()
        out[Kx] = implied_vol(price, Kx)
    return out


# ----------------------------------------------------------------------------
# 5b. CEV diffusion with Heston stochastic variance
# ----------------------------------------------------------------------------
def heston_cev(M: int, n: int, alpha: float, kappa: float = 3.0,
               theta: float = SIGMA**2, v0: float = SIGMA**2, xi: float = 0.6,
               rho: float = -0.3, seed: int = 0) -> tuple[float, float]:
    """CEV local-vol shape driven by a Heston variance process (full-truncation
    Euler).  Centre the variance on the baseline (theta = v0 = sigma^2) so the
    comparison to Black-Scholes isolates the stochastic-volatility effect."""
    g = np.random.default_rng(seed)
    dt = T / n
    S = np.full(M, S0)
    v = np.full(M, v0)
    scale = S0 ** (1 - alpha)
    for _ in range(n):
        z1, z2 = g.standard_normal(M), g.standard_normal(M)
        dW1 = np.sqrt(dt) * z1
        dW2 = rho * dW1 + np.sqrt(1 - rho**2) * np.sqrt(dt) * z2
        v_pos = np.maximum(v, 0.0)
        S = np.maximum(S + (R - Q) * S * dt
                       + np.sqrt(v_pos) * scale
                       * np.maximum(S, 1e-8) ** alpha * dW1, 1e-8)
        v = np.maximum(v + kappa * (theta - v_pos) * dt
                       + xi * np.sqrt(v_pos) * dW2, 0.0)
    pay = np.exp(-R * T) * np.maximum(S - K, 0.0)
    return pay.mean(), pay.std(ddof=1) / np.sqrt(M)


# ----------------------------------------------------------------------------
# 6. Implicit finite differences with American projection
# ----------------------------------------------------------------------------
def fd_implicit(M: int, N: int, american: bool, Smax_mult: float = 4.0) -> float:
    """Fully implicit scheme for the BSM PDE (Hull coefficients), solved with
    a banded tridiagonal solver each step; the American constraint is enforced
    by projection V <- max(V, S-K) after each solve.  Validate the European
    version against the closed form before trusting the American number."""
    Smax = Smax_mult * K
    dt = T / N
    S = np.linspace(0.0, Smax, M + 1)
    V = np.maximum(S - K, 0.0)
    j = np.arange(1, M)
    a = 0.5 * dt * ((R - Q) * j - SIGMA**2 * j**2)      # sub-diagonal coeff
    b = 1 + dt * (SIGMA**2 * j**2 + R)                  # diagonal
    c = -0.5 * dt * ((R - Q) * j + SIGMA**2 * j**2)     # super-diagonal
    ab = np.zeros((3, M - 1))
    ab[0, 1:] = c[:-1]
    ab[1, :] = b
    ab[2, :-1] = a[1:]
    for step in range(N):
        tau = (step + 1) * dt  # time-to-maturity at this level
        rhs = V[1:M].copy()
        upper_bc = Smax * np.exp(-Q * tau) - K * np.exp(-R * tau)
        rhs[-1] -= c[-1] * upper_bc          # lower BC is 0, no adjustment
        V[1:M] = solve_banded((1, 1), ab, rhs)
        V[0], V[M] = 0.0, upper_bc
        if american:
            V = np.maximum(V, S - K)
    return float(np.interp(S0, S, V))


# ----------------------------------------------------------------------------
# 7. Longstaff-Schwartz with both bias directions
# ----------------------------------------------------------------------------
def _simulate_gbm_paths(M: int, n: int, seed: int) -> np.ndarray:
    g = np.random.default_rng(seed)
    Z = g.standard_normal((n, M))
    incr = (R - Q - 0.5 * SIGMA**2) * (T / n) + SIGMA * np.sqrt(T / n) * Z
    logS = np.log(S0) + np.cumsum(incr, axis=0)
    return np.exp(np.vstack([np.full(M, np.log(S0)), logS]))


def lsm(M: int, n: int, power: int, seed: int = 0,
        two_pass: bool = False) -> tuple[float, float]:
    """American call by Longstaff-Schwartz regression on (S/K)^k, k=0..power.
    two_pass=False: classic same-path estimator (biased HIGH: the exercise
    rule exploits in-sample regression noise — 'foresight bias').
    two_pass=True: fit the policy on one path set, price it on an independent
    set (biased LOW: any estimated policy is suboptimal).  Reporting both
    brackets the true value (Glasserman 2004)."""
    disc = np.exp(-R * T / n)

    def fit(S):
        cash = np.maximum(S[-1] - K, 0.0)
        coefs = [None] * n
        for t in range(n - 1, 0, -1):
            cash *= disc
            itm = S[t] > K
            if itm.sum() < power + 2:
                continue
            A = np.vander(S[t, itm] / K, power + 1, increasing=True)
            coef, *_ = np.linalg.lstsq(A, cash[itm], rcond=None)
            coefs[t] = coef
            ex = S[t, itm] - K
            stop = ex > A @ coef
            cash[np.where(itm)[0][stop]] = ex[stop]
        return coefs, (cash * disc)

    if not two_pass:
        _, pv = fit(_simulate_gbm_paths(M, n, seed))
        return pv.mean(), pv.std(ddof=1) / np.sqrt(M)

    coefs, _ = fit(_simulate_gbm_paths(M, n, seed))
    S = _simulate_gbm_paths(M, n, seed + 1)          # independent paths
    alive = np.ones(M, bool)
    pv = np.zeros(M)
    for t in range(1, n):
        if coefs[t] is None:
            continue
        itm = alive & (S[t] > K)
        if not itm.any():
            continue
        A = np.vander(S[t, itm] / K, power + 1, increasing=True)
        stop = (S[t, itm] - K) > A @ coefs[t]
        idx = np.where(itm)[0][stop]
        pv[idx] = (S[t, idx] - K) * np.exp(-R * t * T / n)
        alive[idx] = False
    pv[alive] = np.maximum(S[-1, alive] - K, 0.0) * np.exp(-R * T)
    return pv.mean(), pv.std(ddof=1) / np.sqrt(M)


# ----------------------------------------------------------------------------
# Reproduce the report's headline numbers
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    bs = black_scholes()
    print(f"BS call {bs['call']:.4f} | put {bs['put']:.4f} | "
          f"parity gap {bs['parity_gap']:.2e}")

    print(f"CRR n=1000: {crr_tree(1000):.4f} | "
          f"two-point avg n=500: {crr_averaged(500):.5f}")

    print(f"American lookback (3-step): {lookback_tree():.4f}")

    v, se = mc_european(200_000, "avt", seed=1200200)
    print(f"MC antithetic 200k: {v:.4f} ({se:.4f})")

    base = mc_stochastic_rates(100_000, r0=0.03, seed=11)
    ctrl = mc_stochastic_rates(100_000, r0=0.03, sigma_r=0.0, rho=0.0, seed=11)
    print(f"Stoch-rate {base[0]:.4f} vs deterministic control {ctrl[0]:.4f} "
          f"(CRN difference {base[0]-ctrl[0]:+.4f})")

    cf = barrier_closed_form()
    naive = barrier_mc(200_000, 126, seed=21)
    bridge = barrier_mc(200_000, 126, bridge=True, seed=21)
    bgk_cf = barrier_closed_form(B=barrier_bgk_shift(32.0, 126))
    print(f"DO barrier: CF {cf:.4f} | daily MC {naive[0]:.4f} ({naive[1]:.4f}) "
          f"| BGK CF {bgk_cf:.4f} | bridge MC {bridge[0]:.4f} ({bridge[1]:.4f})")

    for a in (0.5, 0.9, 1.0):
        v, se = cev_pricing(100_000, 252, a, seed=31)
        print(f"CEV alpha={a} (calibrated): {v:.4f} ({se:.4f})")

    eu = fd_implicit(2000, 1000, american=False)
    am = fd_implicit(2000, 1000, american=True)
    print(f"FD 2000x1000: European {eu:.4f} | American {am:.4f} | "
          f"CRR-5000 American benchmark {crr_tree(5000, american=True):.4f}")

    hi = lsm(200_000, 126, 3, seed=71)
    lo = lsm(150_000, 126, 3, seed=81, two_pass=True)
    print(f"LSM P=3 sandwich: [{lo[0]:.4f}, {hi[0]:.4f}] "
          f"(s.e. {lo[1]:.4f}/{hi[1]:.4f})")
