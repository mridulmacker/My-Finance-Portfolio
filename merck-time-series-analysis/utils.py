"""Shared utilities for the MRK time-series analysis.

Performance summaries, signal post-processing, and HAC-based inference
(Newey-West mean-differential test and the Ledoit-Wolf (2008) Sharpe test).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

TRADING_DAYS = 252


# ---------------------------------------------------------------------------
# Performance
# ---------------------------------------------------------------------------

def summarise_returns(ret: np.ndarray, trades: np.ndarray | None = None,
                      rf_daily: float = 0.0) -> dict:
    """Cumulative/annualised return, annualised vol, Sharpe, trade count."""
    ret = np.asarray(ret, dtype=float)
    cum = float(np.prod(1.0 + ret) - 1.0)
    ann_ret = float(np.mean(ret) * TRADING_DAYS)
    ann_vol = float(np.std(ret, ddof=1) * np.sqrt(TRADING_DAYS))
    sharpe = np.nan
    if ann_vol > 0:
        sharpe = float((np.mean(ret) - rf_daily) / np.std(ret, ddof=1)
                       * np.sqrt(TRADING_DAYS))
    n_trades = float(np.sum(trades)) if trades is not None else 0.0
    return {"cum_return": cum, "ann_return": ann_ret, "ann_vol": ann_vol,
            "sharpe": sharpe, "n_trades": n_trades,
            "max_dd": max_drawdown(ret)}


def max_drawdown(ret: np.ndarray) -> float:
    wealth = np.cumprod(1.0 + np.asarray(ret, dtype=float))
    peak = np.maximum.accumulate(wealth)
    return float(np.max((peak - wealth) / peak))


# ---------------------------------------------------------------------------
# Signal post-processing
# ---------------------------------------------------------------------------

def trade_units(signal: np.ndarray) -> np.ndarray:
    """Absolute position change per day (first entry counts as a trade)."""
    signal = np.asarray(signal, dtype=float)
    out = np.zeros_like(signal)
    if signal.size:
        out[0] = abs(signal[0])
        out[1:] = np.abs(np.diff(signal))
    return out


def apply_min_holding_period(signal: np.ndarray, min_hold: int) -> np.ndarray:
    """Suppress sign flips until a position has been held `min_hold` days."""
    sig = np.asarray(signal, dtype=float).copy()
    if sig.size <= 1 or min_hold <= 1:
        return sig
    hold = 1
    for i in range(1, len(sig)):
        if np.sign(sig[i]) != np.sign(sig[i - 1]):
            if hold >= min_hold:
                hold = 1
            else:
                sig[i] = sig[i - 1]
                hold += 1
        else:
            hold += 1
    return sig


def make_signals_raw(fc: np.ndarray, thr1: float,
                     thr2: float | None = None) -> np.ndarray:
    """Threshold rule with optional doubled position beyond thr2."""
    fc = np.asarray(fc, dtype=float)
    sig = np.zeros_like(fc)
    sig[fc > thr1] = 1
    sig[fc < -thr1] = -1
    if thr2 is not None and thr2 > thr1:
        sig[fc > thr2] = 2
        sig[fc < -thr2] = -2
    return sig


# ---------------------------------------------------------------------------
# HAC inference
# ---------------------------------------------------------------------------

def _hac_variance(x: np.ndarray) -> float:
    """Newey-West long-run variance of the mean with Bartlett kernel,
    bandwidth floor(T^(1/3)) — matching the original MATLAB implementation."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    xbar = x.mean()
    max_lag = int(np.floor(n ** (1 / 3)))
    gamma0 = np.var(x, ddof=1)
    acc = 0.0
    for j in range(1, max_lag + 1):
        w = 1.0 - j / (max_lag + 1.0)
        acc += 2.0 * w * np.sum((x[j:] - xbar) * (x[:-j] - xbar)) / n
    return (gamma0 + acc) / n


def hac_mean_diff_test(d: np.ndarray) -> tuple[float, float]:
    """H0: E[d] = 0 with HAC standard errors. Returns (stat, p)."""
    var = _hac_variance(d)
    if var <= 0:
        return np.nan, np.nan
    stat = np.mean(d) / np.sqrt(var)
    return float(stat), float(2 * (1 - stats.norm.cdf(abs(stat))))


def ledoit_wolf_sharpe_test(r1: np.ndarray, r2: np.ndarray) -> tuple[float, float]:
    """Ledoit-Wolf (2008) HAC-corrected test of equal Sharpe ratios."""
    r1, r2 = np.asarray(r1, float), np.asarray(r2, float)
    m1, m2 = r1.mean(), r2.mean()
    s1, s2 = r1.std(ddof=1), r2.std(ddof=1)
    if s1 == 0 or s2 == 0:
        return np.nan, np.nan
    sr1, sr2 = m1 / s1, m2 / s2
    psi = ((r1 - m1) / s1 - (sr1 / 2) * (((r1 - m1) ** 2) / s1 ** 2 - 1)
           - (r2 - m2) / s2 + (sr2 / 2) * (((r2 - m2) ** 2) / s2 ** 2 - 1))
    var = _hac_variance(psi)
    if var <= 0:
        return np.nan, np.nan
    stat = (sr1 - sr2) / np.sqrt(var)
    return float(stat), float(2 * (1 - stats.norm.cdf(abs(stat))))


def diebold_mariano(loss1: np.ndarray, loss2: np.ndarray) -> tuple[float, float]:
    """DM test on a loss differential.

    All forecasts in this project are one-step-ahead, so under the null the
    differential is serially uncorrelated and the zero-lag variance suffices
    (Diebold and Mariano, 1995).
    """
    d = np.asarray(loss1, float) - np.asarray(loss2, float)
    var_d = np.var(d, ddof=0) / len(d)
    if var_d <= 0:
        return np.nan, np.nan
    stat = d.mean() / np.sqrt(var_d)
    return float(stat), float(2 * (1 - stats.norm.cdf(abs(stat))))


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------

def load_prices(path: str) -> pd.DataFrame:
    """Load a daily price CSV with a Date column; sorts and drops NA rows."""
    df = pd.read_csv(path, parse_dates=["Date"])
    df = df.dropna().sort_values("Date").reset_index(drop=True)
    return df
