"""Task 2 - Rolling-window ARMA / ARMA-GARCH trading strategies on MRK.

Pipeline (mirrors the original MATLAB implementation):
  1. 60/20/20 train / development / test split of daily log returns.
  2. ARMA(p, q) order selection on the training sample, p, q in {0, 1, 2},
     by AIC / BIC / HQIC voting.
  3. GARCH(1,1) diagnostics on ARMA residuals (Ljung-Box on squared resids).
  4. Development-sample grid search over threshold rules (raw / vol-scaled,
     optional doubled positions), 10-day minimum holding period, 10 bps
     transaction costs.
  5. One-shot test-sample evaluation vs a dual moving-average rule
     (Taylor, 2007) and buy-and-hold, with Newey-West HAC mean tests,
     Ledoit-Wolf Sharpe tests, and a transaction-cost sensitivity sweep.

Note on estimation: the original MATLAB estimated ARMA-GARCH jointly.
Here the mean and variance equations are estimated in two steps
(statsmodels ARIMA, then `arch` GARCH(1,1) on residuals) - an asymptotically
valid QMLE approach whose point estimates can differ slightly from joint MLE.

Usage:
    python src/task2_trading_strategies.py --data data/mrk_close.csv
    # quick smoke test on a short series:
    python src/task2_trading_strategies.py --data data/mrk_close.csv --fast
"""

from __future__ import annotations

import argparse
import itertools
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from arch import arch_model
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.tsa.arima.model import ARIMA

from utils import (apply_min_holding_period, hac_mean_diff_test,
                   ledoit_wolf_sharpe_test, load_prices, make_signals_raw,
                   max_drawdown, summarise_returns, trade_units)

warnings.filterwarnings("ignore")

RF_ANNUAL = 0.04
RF_DAILY = RF_ANNUAL / 252
TC_DEFAULT = 0.001          # 10 bps per unit of position change
MIN_HOLD = 10
REFIT_EVERY = 20


# ---------------------------------------------------------------------------
# Model selection
# ---------------------------------------------------------------------------

def ic_table(train: np.ndarray, max_p: int = 2, max_q: int = 2) -> pd.DataFrame:
    rows = []
    for p, q in itertools.product(range(max_p + 1), range(max_q + 1)):
        if p == 0 and q == 0:
            continue
        try:
            res = ARIMA(train, order=(p, 0, q)).fit()
            rows.append({"p": p, "q": q, "AIC": res.aic, "BIC": res.bic,
                         "HQIC": res.hqic, "params": p + q + 2})
        except Exception:
            continue
    tab = pd.DataFrame(rows)
    for c in ("AIC", "BIC", "HQIC"):
        tab[f"pass_{c}"] = tab[c] <= tab[c].min() + 2.0
    tab["votes"] = tab[[f"pass_{c}" for c in ("AIC", "BIC", "HQIC")]].sum(axis=1)
    return tab.sort_values("AIC").reset_index(drop=True)


def garch_diagnostics(train: np.ndarray, p: int, q: int) -> dict:
    """Ljung-Box on squared ARMA residuals + two-step GARCH(1,1) fit."""
    arma = ARIMA(train, order=(p, 0, q)).fit()
    resid = arma.resid
    lb = acorr_ljungbox(resid ** 2, lags=[10], return_df=True)
    g = arch_model(resid * 100, vol="GARCH", p=1, q=1, mean="Zero").fit(disp="off")
    omega = g.params["omega"] / 100 ** 2
    alpha, beta = g.params["alpha[1]"], g.params["beta[1]"]
    sigma_target = float(np.mean(g.conditional_volatility) / 100)
    return {"lb_stat": float(lb["lb_stat"].iloc[0]),
            "lb_p": float(lb["lb_pvalue"].iloc[0]),
            "omega": float(omega), "alpha": float(alpha), "beta": float(beta),
            "sigma_target": sigma_target}


# ---------------------------------------------------------------------------
# Rolling one-step-ahead forecasts
# ---------------------------------------------------------------------------

def rolling_forecasts(log_ret: np.ndarray, start: int, end: int, win: int,
                      p: int, q: int, refit_every: int = REFIT_EVERY,
                      use_garch: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Mean forecast (and conditional vol if `use_garch`) for t in [start, end).

    The model is re-fitted every `refit_every` days on the trailing `win`
    observations; between refits, fixed coefficients are applied to the
    updated conditioning data via `ARIMAResults.apply`.
    """
    n = end - start
    fc = np.zeros(n)
    vol = np.zeros(n)
    fitted = None
    g_params = None
    last_vol = float(np.std(log_ret[start - win:start], ddof=1))

    for i in range(n):
        t = start + i
        y = log_ret[t - win:t]
        if fitted is None or i % refit_every == 0:
            fitted = ARIMA(y, order=(p, 0, q)).fit()
            if use_garch:
                g = arch_model(fitted.resid * 100, vol="GARCH", p=1, q=1,
                               mean="Zero").fit(disp="off")
                g_params = (g.params["omega"], g.params["alpha[1]"],
                            g.params["beta[1]"])
            applied = fitted
        else:
            applied = fitted.apply(y)   # fixed params, rolled conditioning data
        fc[i] = float(np.asarray(applied.forecast(1)).ravel()[0])

        if use_garch and g_params is not None:
            resid = applied.resid * 100
            omega, alpha, beta = g_params
            # filter conditional variance forward through the window
            h = np.var(resid, ddof=1)
            for e in resid:
                h = omega + alpha * e ** 2 + beta * h
            vol[i] = np.sqrt(max(h, 1e-12)) / 100
            last_vol = vol[i]
        else:
            vol[i] = max(float(np.std(y, ddof=1)), 1e-12)
    return fc, vol


# ---------------------------------------------------------------------------
# Strategy evaluation
# ---------------------------------------------------------------------------

def evaluate_rule(fc, cond_vol, simple_ret, mode, thr, dd_mult, tc):
    thr2 = dd_mult * thr if dd_mult > 1 else None
    base = fc / np.maximum(cond_vol, 1e-12) if mode == "volscaled" else fc
    sig = apply_min_holding_period(make_signals_raw(base, thr, thr2), MIN_HOLD)
    trades = trade_units(sig)
    net = sig * simple_ret - tc * trades
    s = summarise_returns(net, trades, RF_DAILY)
    s.update({"mode": mode, "thr": thr, "dd_mult": dd_mult,
              "flat_pct": 100 * np.mean(sig == 0)})
    return s, sig, net


def dev_grid_search(fc, cond_vol, dev_ret, tc, raw_grid, vol_grid,
                    dd_mults=(1, 2, 3), min_trades=25):
    rows = []
    for mode, grid in (("raw", raw_grid), ("volscaled", vol_grid)):
        for thr, dd in itertools.product(grid, dd_mults):
            s, sig, _ = evaluate_rule(fc, cond_vol, dev_ret, mode, thr, dd, tc)
            d = np.sign(sig)
            trivial = (np.all(d == 0) or s["n_trades"] < min_trades
                       or np.all(d >= 0) or np.all(d <= 0))
            rows.append({**s, "trivial": trivial})
    tab = pd.DataFrame(rows).sort_values("sharpe", ascending=False)
    return tab.reset_index(drop=True)


def ma_signal(price: np.ndarray, S: int, L: int, B: float) -> np.ndarray:
    s = pd.Series(price)
    short = s.rolling(S).mean()
    long_ = s.rolling(L).mean()
    r = (short - long_) / long_
    sig = np.where(r > B, 1.0, np.where(r < -B, -1.0, 0.0))
    sig[: L - 1] = 0.0
    return sig


def ma_benchmark(price: np.ndarray, tc: float,
                 S_list=(1, 5, 10), L_list=(50, 100), B_list=(0.0, 0.01)):
    """Taylor (2007) three-state MA rule, 80/20 split, best net OOS Sharpe."""
    ret = np.diff(np.log(price))
    oos_start = int(np.floor(0.80 * len(price)))
    best = None
    for S, L, B in itertools.product(S_list, L_list, B_list):
        if S >= L:
            continue
        sig_full = ma_signal(price, S, L, B)[:-1]      # position for next day
        oos_sig = sig_full[oos_start - 1:]
        oos_ret = ret[oos_start - 1:]
        gross = oos_sig * oos_ret
        change = np.concatenate([[0.0], np.abs(np.diff(oos_sig))])
        net = gross - tc * change
        vol = net.std(ddof=1)
        sharpe = net.mean() / vol if vol > 0 else np.nan
        row = {"S": S, "L": L, "B": B, "sharpe_daily": sharpe,
               "n_trades": int((change > 0).sum())}
        if best is None or (np.isfinite(sharpe)
                            and sharpe > best["sharpe_daily"]):
            best = row
    return best


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(data_path: str, tc: float, fast: bool, outdir: str = "outputs") -> None:
    Path(outdir).mkdir(parents=True, exist_ok=True)
    df = load_prices(data_path)
    price = df["Close"].to_numpy()
    log_ret = np.diff(np.log(price))
    simple_ret = np.exp(log_ret) - 1.0
    N = len(log_ret)

    T_train = int(np.floor(0.60 * N))
    T_dev = int(np.floor(0.20 * N))
    dev_s, dev_e = T_train, T_train + T_dev
    test_s, test_e = dev_e, N
    win = T_train
    refit = REFIT_EVERY if not fast else max(5, REFIT_EVERY)

    print(f"Split - train: {T_train}  dev: {T_dev}  test: {N - dev_e}")

    # 1. order selection ----------------------------------------------------
    train = log_ret[:T_train]
    tab = ic_table(train)
    print("\nARMA order selection (training sample):")
    print(tab.round(2).to_string(index=False))
    best_p, best_q = int(tab.iloc[0]["p"]), int(tab.iloc[0]["q"])
    print(f"Selected mean equation: ARMA({best_p},{best_q})")

    # 2. GARCH diagnostics --------------------------------------------------
    diag = garch_diagnostics(train, best_p, best_q)
    print(f"\nLjung-Box Q(10) on squared ARMA residuals: "
          f"stat={diag['lb_stat']:.2f}  p={diag['lb_p']:.4g}")
    print(f"GARCH(1,1): omega={diag['omega']:.2e}  alpha={diag['alpha']:.4f}  "
          f"beta={diag['beta']:.4f}  persistence={diag['alpha']+diag['beta']:.4f}")

    # 3. development-sample forecasts + grid search -------------------------
    raw_grid = [1e-4, 2e-4, 5e-4, 1e-3] if fast else \
        [1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7.5e-4, 1e-3, 1.5e-3, 2e-3, 3e-3, 5e-3]
    vol_grid = [0.005, 0.02, 0.05] if fast else \
        [0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15, 0.20]

    print("\nGenerating development-sample rolling forecasts ...")
    fc_a_dev, vol_a_dev = rolling_forecasts(log_ret, dev_s, dev_e, win,
                                            best_p, best_q, refit, False)
    fc_b_dev, vol_b_dev = rolling_forecasts(log_ret, dev_s, dev_e, win,
                                            best_p, best_q, refit, True)
    dev_ret = simple_ret[dev_s:dev_e]

    grid_a = dev_grid_search(fc_a_dev, vol_a_dev, dev_ret, tc, raw_grid, vol_grid)
    grid_b = dev_grid_search(fc_b_dev, vol_b_dev, dev_ret, tc, raw_grid, vol_grid)
    best_a = grid_a[~grid_a["trivial"]].iloc[0]
    best_b = grid_b[~grid_b["trivial"]].iloc[0]
    bh_dev = summarise_returns(dev_ret, None, RF_DAILY)
    print(f"\nBest dev config A: mode={best_a['mode']} thr={best_a['thr']:g} "
          f"dd={best_a['dd_mult']:.0f}  Sharpe={best_a['sharpe']:.4f}")
    print(f"Best dev config B: mode={best_b['mode']} thr={best_b['thr']:g} "
          f"dd={best_b['dd_mult']:.0f}  Sharpe={best_b['sharpe']:.4f}")
    print(f"Dev buy-and-hold Sharpe: {bh_dev['sharpe']:.4f}")

    # 4. MA benchmark -------------------------------------------------------
    ma_best = ma_benchmark(price, tc)
    print(f"\nMA benchmark (80/20 split): S={ma_best['S']} L={ma_best['L']} "
          f"B={ma_best['B']}  daily Sharpe={ma_best['sharpe_daily']:.4f}")

    # 5. test-sample evaluation --------------------------------------------
    print("\nGenerating test-sample rolling forecasts ...")
    fc_a, vol_a = rolling_forecasts(log_ret, test_s, test_e, win,
                                    best_p, best_q, refit, False)
    fc_b, vol_b = rolling_forecasts(log_ret, test_s, test_e, win,
                                    best_p, best_q, refit, True)
    test_ret = simple_ret[test_s:test_e]

    sA, sigA, netA = evaluate_rule(fc_a, vol_a, test_ret, best_a["mode"],
                                   best_a["thr"], int(best_a["dd_mult"]), tc)
    sB, sigB, netB = evaluate_rule(fc_b, vol_b, test_ret, best_b["mode"],
                                   best_b["thr"], int(best_b["dd_mult"]), tc)

    sig_ma_full = ma_signal(price, ma_best["S"], ma_best["L"], ma_best["B"])
    sig_ma = sig_ma_full[test_s:test_e]
    tr_ma = trade_units(sig_ma)
    net_ma = sig_ma * test_ret - tc * tr_ma
    sMA = summarise_returns(net_ma, tr_ma, RF_DAILY)
    sBH = summarise_returns(test_ret, None, RF_DAILY)

    res = pd.DataFrame([
        {"strategy": "ARMA net", **{k: sA[k] for k in
         ("cum_return", "ann_return", "ann_vol", "sharpe", "n_trades", "max_dd")}},
        {"strategy": "ARMA-GARCH net", **{k: sB[k] for k in
         ("cum_return", "ann_return", "ann_vol", "sharpe", "n_trades", "max_dd")}},
        {"strategy": f"MA({ma_best['S']},{ma_best['L']}) net", **{k: sMA[k] for k in
         ("cum_return", "ann_return", "ann_vol", "sharpe", "n_trades", "max_dd")}},
        {"strategy": "Buy-and-Hold", **{k: sBH[k] for k in
         ("cum_return", "ann_return", "ann_vol", "sharpe", "n_trades", "max_dd")}},
    ])
    print("\nOUT-OF-SAMPLE TEST PERFORMANCE")
    print(res.round(4).to_string(index=False))

    hitA = 100 * np.mean(np.sign(fc_a) == np.sign(test_ret))
    hitB = 100 * np.mean(np.sign(fc_b) == np.sign(test_ret))
    print(f"\nDirectional hit rates - A: {hitA:.2f}%  B: {hitB:.2f}%")

    # 6. pairwise tests -----------------------------------------------------
    pairs = {"ARMA vs ARMA-GARCH": (netA, netB),
             "ARMA vs MA": (netA, net_ma),
             "ARMA vs B&H": (netA, test_ret),
             "ARMA-GARCH vs MA": (netB, net_ma),
             "ARMA-GARCH vs B&H": (netB, test_ret),
             "MA vs B&H": (net_ma, test_ret)}
    print("\nPAIRWISE TESTS (HAC mean diff / Ledoit-Wolf Sharpe)")
    for name, (r1, r2) in pairs.items():
        h, hp = hac_mean_diff_test(np.asarray(r1) - np.asarray(r2))
        l, lp = ledoit_wolf_sharpe_test(r1, r2)
        print(f"{name:22s} HAC={h:7.4f} p={hp:.4f}   LW={l:7.4f} p={lp:.4f}")

    # 7. transaction-cost sensitivity ---------------------------------------
    print("\nTRANSACTION-COST SENSITIVITY (annualised net Sharpe)")
    print(f"{'TC(bps)':>8s} {'ARMA':>9s} {'ARMA-G':>9s} {'MA':>9s}")
    for tcj in (0.0, 0.0005, 0.001, 0.0015, 0.002, 0.003, 0.005):
        trA, trB = trade_units(sigA), trade_units(sigB)
        a = summarise_returns(sigA * test_ret - tcj * trA, trA, RF_DAILY)
        b = summarise_returns(sigB * test_ret - tcj * trB, trB, RF_DAILY)
        m = summarise_returns(sig_ma * test_ret - tcj * tr_ma, tr_ma, RF_DAILY)
        print(f"{tcj*1e4:8.1f} {a['sharpe']:9.4f} {b['sharpe']:9.4f} "
              f"{m['sharpe']:9.4f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/mrk_close.csv")
    ap.add_argument("--tc", type=float, default=TC_DEFAULT)
    ap.add_argument("--fast", action="store_true",
                    help="smaller grids for a quick run")
    ap.add_argument("--outdir", default="outputs")
    args = ap.parse_args()
    run(args.data, args.tc, args.fast, args.outdir)
