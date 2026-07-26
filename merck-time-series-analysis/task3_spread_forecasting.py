"""Task 3 - Forecasting the MRK daily high-low log price spread.

The spread S_t = log(H_t) - log(L_t) is a range-based volatility proxy
(Parkinson, 1980). Three models are estimated on the log spread:

  * Log-AR(1)      - benchmark, coefficients fixed at the initial window
  * Log-HAR        - Corsi (2009), re-estimated by OLS in every window
  * Log-ARMA(p,q)  - BIC-selected on the initial window, coefficients fixed

One-day-ahead forecasts are back-transformed with a log-normal bias
correction exp(y_hat + sigma^2/2) and evaluated on the level scale with
MSPE and QLIKE loss, Mincer-Zarnowitz regressions with a Wald test of
(c, beta) = (0, 1), and pairwise Diebold-Mariano tests.

The fixed-coefficient treatment of AR(1)/ARMA replicates the original
MATLAB design (a documented computational-constraint limitation); pass
--refit-arma to re-estimate them every `--refit-every` days instead.

Usage:
    python src/task3_spread_forecasting.py --data data/mrk_ohlc.csv
"""

from __future__ import annotations

import argparse
import itertools
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.tsa.arima.model import ARIMA

from utils import diebold_mariano, load_prices

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# HAR
# ---------------------------------------------------------------------------

def har_design(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Daily/weekly/monthly HAR regressors for a log-spread window."""
    n = len(y)
    y_d = y[21:n - 1]
    w = pd.Series(y).rolling(5).mean().to_numpy()[21:n - 1]
    m = pd.Series(y).rolling(22).mean().to_numpy()[21:n - 1]
    X = np.column_stack([np.ones_like(y_d), y_d, w, m])
    target = y[22:]
    return X, target


def har_rolling_forecasts(y_log: np.ndarray, is_len: int) -> np.ndarray:
    """One-step-ahead HAR forecasts, OLS re-estimated in each window."""
    T = len(y_log)
    out = np.empty(T - is_len)
    for t in range(T - is_len):
        win = y_log[t:t + is_len]
        X, target = har_design(win)
        beta, *_ = np.linalg.lstsq(X, target, rcond=None)
        resid = target - X @ beta
        sigma2 = resid @ resid / (len(resid) - 4)
        x_next = np.array([1.0, win[-1], win[-5:].mean(), win[-22:].mean()])
        out[t] = np.exp(x_next @ beta + sigma2 / 2)
    return out


# ---------------------------------------------------------------------------
# AR(1) / ARMA with fixed or refitted coefficients
# ---------------------------------------------------------------------------

def arima_rolling_forecasts(y_log: np.ndarray, is_len: int, order: tuple,
                            refit_every: int | None = None) -> np.ndarray:
    """One-step forecasts with rolled conditioning data.

    refit_every=None reproduces the original design: fit once on the initial
    window, hold coefficients and the bias-correction variance fixed, and
    roll only the conditioning data via `ARIMAResults.apply`.
    """
    T = len(y_log)
    n_fc = T - is_len
    out = np.empty(n_fc)
    fitted = ARIMA(y_log[:is_len], order=order).fit()
    sigma2 = float(fitted.params[-1])          # innovation variance
    for t in range(n_fc):
        if refit_every and t > 0 and t % refit_every == 0:
            fitted = ARIMA(y_log[t:t + is_len], order=order).fit()
            sigma2 = float(fitted.params[-1])
        window = y_log[t:t + is_len]
        applied = fitted.apply(window)
        out[t] = np.exp(float(np.asarray(applied.forecast(1)).ravel()[0]) + sigma2 / 2)
    return out


def select_arma_order(y: np.ndarray, max_p: int = 5, max_q: int = 5) -> tuple:
    best, best_bic = (1, 0, 0), np.inf
    for p, q in itertools.product(range(max_p + 1), range(max_q + 1)):
        if p == 0 and q == 0:
            continue
        try:
            bic = ARIMA(y, order=(p, 0, q)).fit().bic
        except Exception:
            continue
        if bic < best_bic:
            best, best_bic = (p, 0, q), bic
    return best


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def losses(actual: np.ndarray, fc: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mspe = (actual - fc) ** 2
    ratio = actual / fc
    qlike = ratio - np.log(ratio) - 1.0
    return mspe, qlike


def mincer_zarnowitz(actual: np.ndarray, fc: np.ndarray) -> dict:
    X = np.column_stack([np.ones_like(fc), fc])
    b, *_ = np.linalg.lstsq(X, actual, rcond=None)
    resid = actual - X @ b
    sigma2 = resid @ resid / (len(actual) - 2)
    V = sigma2 * np.linalg.inv(X.T @ X)
    diff = b - np.array([0.0, 1.0])
    wald = float(diff @ np.linalg.inv(V) @ diff)
    p = float(1 - stats.chi2.cdf(wald, 2))
    ss_res = resid @ resid
    ss_tot = np.sum((actual - actual.mean()) ** 2)
    return {"intercept": float(b[0]), "slope": float(b[1]),
            "wald": wald, "p": p, "r2": float(1 - ss_res / ss_tot)}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(data_path: str, is_frac: float, refit_arma: bool, refit_every: int,
        outdir: str = "outputs") -> None:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    df = load_prices(data_path)
    spread = np.log(df["High"].to_numpy()) - np.log(df["Low"].to_numpy())
    spread = spread[spread > 0]
    y_log = np.log(spread)
    T = len(spread)
    is_len = int(np.floor(is_frac * T))
    actual = spread[is_len:]
    print(f"Observations: {T}  initial window: {is_len}  "
          f"out-of-sample forecasts: {T - is_len}")

    # preliminary plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].plot(spread, lw=0.4)
    ax[0].set_title("MRK Daily High-Low Log Spread")
    ax[1].hist(spread, bins=80)
    ax[1].set_title("Spread Distribution")
    fig.tight_layout()
    fig.savefig(out / "t3_spread.png", dpi=150)

    # order selection on the initial window
    order = select_arma_order(y_log[:is_len])
    print(f"BIC-selected ARMA order on initial window: "
          f"({order[0]},{order[2]})")

    # forecasts
    print("Forecasting: Log-HAR (re-estimated each window) ...")
    fc_har = har_rolling_forecasts(y_log, is_len)
    print("Forecasting: Log-AR(1) ...")
    fc_ar = arima_rolling_forecasts(
        y_log, is_len, (1, 0, 0),
        refit_every if refit_arma else None)
    print(f"Forecasting: Log-ARMA({order[0]},{order[2]}) ...")
    fc_arma = arima_rolling_forecasts(
        y_log, is_len, order,
        refit_every if refit_arma else None)

    models = {"Log-HAR": fc_har,
              f"Log-ARMA({order[0]},{order[2]})": fc_arma,
              "Log-AR(1)": fc_ar}

    # losses + MZ
    rows = []
    loss_store = {}
    for name, fc in models.items():
        mspe, qlike = losses(actual, fc)
        loss_store[name] = (mspe, qlike)
        mz = mincer_zarnowitz(actual, fc)
        rows.append({"model": name, "MSPE(1e-4)": mspe.mean() * 1e4,
                     "QLIKE": qlike.mean(),
                     "MZ_intercept(1e-3)": mz["intercept"] * 1e3,
                     "MZ_slope": mz["slope"], "Wald": mz["wald"],
                     "Wald_p": mz["p"], "MZ_R2": mz["r2"]})
    tab = pd.DataFrame(rows)
    print("\nFORECAST EVALUATION")
    print(tab.round(4).to_string(index=False))

    # DM tests
    names = list(models)
    print("\nDIEBOLD-MARIANO TESTS")
    for loss_name, idx in (("MSPE", 0), ("QLIKE", 1)):
        print(f"--- {loss_name} loss ---")
        for a, b in itertools.combinations(names, 2):
            stat, p = diebold_mariano(loss_store[a][idx], loss_store[b][idx])
            better = b if stat > 0 else a
            verdict = f"{better} better" if p < 0.05 else "no significant difference"
            print(f"{a:18s} vs {b:18s}  DM={stat:8.4f}  p={p:.4g}  ({verdict})")

    # forecast plot
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(actual, lw=0.4, label="Realised spread", alpha=0.6)
    ax.plot(fc_har, lw=0.6, label="Log-HAR forecast")
    ax.legend()
    ax.set_title("One-Day-Ahead Spread Forecasts (Log-HAR)")
    fig.tight_layout()
    fig.savefig(out / "t3_forecasts.png", dpi=150)
    print("\nFigures written to", out.resolve())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/mrk_ohlc.csv")
    ap.add_argument("--is-frac", type=float, default=0.30,
                    help="initial estimation window as a fraction of the sample")
    ap.add_argument("--refit-arma", action="store_true",
                    help="re-estimate AR(1)/ARMA coefficients periodically "
                         "(removes the frozen-coefficient limitation of the "
                         "original design)")
    ap.add_argument("--refit-every", type=int, default=20)
    ap.add_argument("--outdir", default="outputs")
    args = ap.parse_args()
    run(args.data, args.is_frac, args.refit_arma, args.refit_every, args.outdir)
