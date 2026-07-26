"""Task 1 - Statistical properties of MRK daily log returns.

Stationarity (ADF), normality (Jarque-Bera, Q-Q, histogram), and the
white-noise hypothesis (ACF/PACF, Ljung-Box, periodogram / Welch PSD).

Usage:
    python src/task1_return_properties.py --data data/mrk_close.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal as sp_signal
from scipy import stats
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.stats.stattools import jarque_bera
from statsmodels.tsa.stattools import adfuller

from utils import load_prices


def run(data_path: str, outdir: str = "outputs") -> None:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    df = load_prices(data_path)
    log_p = np.log(df["Close"].to_numpy())
    ret = np.diff(log_p)               # decimal log returns
    ret_pct = 100 * ret

    # ------------------------------------------------------------------ plots
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].plot(log_p, lw=0.7)
    ax[0].set_title("MRK Log Prices")
    ax[1].plot(ret_pct, lw=0.4)
    ax[1].set_title("MRK Daily Log Returns (%)")
    fig.tight_layout()
    fig.savefig(out / "t1_prices_returns.png", dpi=150)

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    plot_acf(ret, lags=100, ax=ax[0], title="Sample ACF of Log Returns")
    plot_pacf(ret, lags=100, ax=ax[1], method="ywm",
              title="Sample PACF of Log Returns")
    fig.tight_layout()
    fig.savefig(out / "t1_acf_pacf.png", dpi=150)

    fig, ax = plt.subplots(figsize=(6, 6))
    stats.probplot(ret, dist="norm", plot=ax)
    ax.set_title("Q-Q Plot of MRK Log Returns")
    fig.tight_layout()
    fig.savefig(out / "t1_qq_plot.png", dpi=150)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(ret, bins=100, density=True, alpha=0.75)
    x = np.linspace(ret.min(), ret.max(), 400)
    ax.plot(x, stats.norm.pdf(x, ret.mean(), ret.std(ddof=1)), "k-", lw=1.5)
    ax.set_title("Histogram of Daily Log Returns with Fitted Normal")
    ax.set_xlabel("Daily Log Return")
    ax.set_ylabel("Density")
    fig.tight_layout()
    fig.savefig(out / "t1_histogram.png", dpi=150)

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    f_p, pxx_p = sp_signal.periodogram(ret)
    f_w, pxx_w = sp_signal.welch(ret, nperseg=min(1024, len(ret) // 4))
    ax[0].semilogy(f_p[1:], pxx_p[1:], lw=0.4)
    ax[0].set_title("Periodogram PSD Estimate")
    ax[1].semilogy(f_w[1:], pxx_w[1:], lw=0.8)
    ax[1].set_title("Welch PSD Estimate")
    for a in ax:
        a.set_xlabel("Normalised Frequency")
    fig.tight_layout()
    fig.savefig(out / "t1_psd.png", dpi=150)

    # ------------------------------------------------------------------ tests
    print("=" * 60)
    print("DESCRIPTIVE STATISTICS (log returns)")
    print("=" * 60)
    print(f"n = {len(ret)}  mean = {ret.mean():.6f}  std = {ret.std(ddof=1):.6f}")
    print(f"skewness = {stats.skew(ret):.4f}  "
          f"kurtosis (Pearson) = {stats.kurtosis(ret, fisher=False):.2f}")

    print("\n" + "=" * 60)
    print("AUGMENTED DICKEY-FULLER TESTS (lags by BIC)")
    print("=" * 60)
    for name, series, reg in (("Log Prices", log_p, "ct"),
                              ("Log Returns", ret, "c")):
        stat, pval, usedlag, nobs, crit, _ = adfuller(series, maxlag=20,
                                                      regression=reg,
                                                      autolag="BIC")
        concl = "Stationary" if pval < 0.05 else "Non-Stationary"
        print(f"{name:12s}  stat = {stat:9.4f}  p = {pval:.4f}  "
              f"lag = {usedlag:2d}  -> {concl}")

    print("\n" + "=" * 60)
    print("JARQUE-BERA NORMALITY TEST")
    print("=" * 60)
    jb, jb_p, skw, kurt = jarque_bera(ret)
    print(f"JB = {jb:.4f}  p = {jb_p:.4g}  "
          f"-> {'Reject' if jb_p < 0.05 else 'Fail to reject'} normality at 5%")

    print("\n" + "=" * 60)
    print("LJUNG-BOX Q-TEST, LAGS 1-20 (returns)")
    print("=" * 60)
    lb = acorr_ljungbox(ret, lags=range(1, 21), return_df=True)
    lb["reject_5pct"] = lb["lb_pvalue"] < 0.05
    print(lb.round(4).to_string())

    print("\n" + "=" * 60)
    print("LJUNG-BOX Q-TEST, LAGS 1-20 (squared returns - ARCH check)")
    print("=" * 60)
    lb2 = acorr_ljungbox(ret ** 2, lags=[10, 20], return_df=True)
    print(lb2.round(4).to_string())
    print("\nFigures written to", out.resolve())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/mrk_close.csv")
    ap.add_argument("--outdir", default="outputs")
    args = ap.parse_args()
    run(args.data, args.outdir)
