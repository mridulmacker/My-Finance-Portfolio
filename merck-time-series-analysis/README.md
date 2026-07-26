# Time Series Analysis of Merck & Co. (MRK): Returns, Trading Rules, and Volatility Forecasting

Empirical time-series econometrics on 26 years of daily MRK equity data (2000–2025), covering three questions:

1. **What do the returns look like?** Stationarity, normality, and white-noise testing.
2. **Is the conditional mean tradeable?** Rolling-window ARMA / ARMA-GARCH trading strategies evaluated out-of-sample against technical and passive benchmarks, net of transaction costs.
3. **Is the conditional variance forecastable?** One-day-ahead forecasting of the high–low log price spread (a range-based volatility proxy) with Log-AR(1), Log-HAR (Corsi, 2009), and Log-ARMA models.

**Headline finding:** the classic asymmetry of financial econometrics reproduced in a single security — the *mean* of daily returns is unpredictable (no strategy beats buy-and-hold, even at zero transaction costs), while the *variance* is highly predictable (HAR and ARMA spread forecasts beat an AR(1) benchmark with Diebold–Mariano p < 0.0001 under both MSPE and QLIKE loss).

Full write-up with all figures and detailed results: [Merck_Time_Series_Analysis.pdf](./Merck_Time_Series_Analysis.pdf)

## Background

This project was completed individually as coursework for **BMAN 71122 Time Series Econometrics** at Alliance Manchester Business School. All analysis, code, and write-up are my own work.

## Key results (MRK 2000–2025)

### Task 1 — Return properties
| Property | Finding |
|---|---|
| Stationarity (ADF) | Log prices non-stationary (p = 0.348); log returns stationary (p = 0.001) |
| Normality (Jarque–Bera) | Rejected decisively (JB ≈ 178,582); kurtosis ≈ 28, mild negative skew |
| White noise (Ljung–Box) | No linear autocorrelation for lags 1–14; weak dependence beyond lag 15 |
| Independence | Rejected: squared residuals show strong ARCH effects (Q = 43.26, p < 0.0001) |

### Task 2 — Trading strategies (test sample: Oct 2020 – Dec 2025, net of 10 bps costs)
| Strategy | Cum. return | Sharpe | Trades |
|---|---:|---:|---:|
| ARMA(2,2) | −54.2% | −0.32 | 524 |
| ARMA(2,2)-GARCH(1,1) | −9.8% | −0.23 | 142 |
| MA(20,100) crossover | +10.9% | −0.09 | 18 |
| **Buy-and-hold** | **+37.1%** | **+0.21** | 0 |

Hit rates are indistinguishable from a coin toss (~50%), and zero-cost Sharpe ratios remain negative — the underperformance reflects a genuine absence of exploitable mean predictability, not frictions. Consistent with weak-form market efficiency (Fama, 1970).

### Task 3 — Spread forecasting (rolling one-day-ahead, ~4,600 OOS forecasts)
| Model | MSPE (×10⁻⁴) | QLIKE | MZ R² |
|---|---:|---:|---:|
| Log-HAR | 1.0021 | 0.0824 | 0.416 |
| Log-ARMA(1,2) | 1.0074 | 0.0829 | 0.416 |
| Log-AR(1) benchmark | 1.26 | 0.1044 | 0.309 |

HAR and ARMA are statistically indistinguishable from each other (DM p ≈ 0.28) but both dominate the AR(1) benchmark decisively (DM p < 0.0001 under both losses).

## Repository structure

```
├── download_data.py             # fetch MRK OHLC via yfinance
├── task1_return_properties.py   # ADF, Jarque–Bera, Ljung–Box, ACF/PACF, PSD
├── task2_trading_strategies.py  # rolling ARMA/ARMA-GARCH rules + benchmarks + HAC/LW tests
├── task3_spread_forecasting.py  # Log-AR(1)/HAR/ARMA spread forecasts + MZ + DM tests
├── utils.py                     # performance, signals, HAC & Ledoit–Wolf inference
└── requirements.txt
```
## Setup & usage

```bash
pip install -r requirements.txt

# 1. get the data
python download_data.py --start 2000-01-01 --end 2025-12-31

# 2. run the three analyses
python task1_return_properties.py --data data/mrk_close.csv
python task2_trading_strategies.py --data data/mrk_close.csv        # add --fast for a quick pass
python task3_spread_forecasting.py --data data/mrk_ohlc.csv
```

Task 2's full grid search with per-day rolling forecasts takes a while on 26 years of data; `--fast` shrinks the threshold grids for a quicker run. Task 3 accepts `--refit-arma` to re-estimate the AR(1)/ARMA coefficients periodically, removing the frozen-coefficient limitation of the fixed-coefficient default.

## Methodological highlights

- **Strict anti-snooping protocol** (Task 2): 60/20/20 train/development/test split; all hyperparameter search on the development sample; the test sample touched exactly once (White, 2000).
- **Honest cost accounting**: proportional transaction costs on position *changes*, a 10-day minimum holding period, and a full cost-sensitivity sweep with numerical break-even costs.
- **HAC inference throughout**: Newey–West mean-differential tests and Ledoit–Wolf (2008) Sharpe-ratio tests on overlapping daily return differentials.
- **Proper log-scale forecasting** (Task 3): log-normal bias correction `exp(ŷ + σ²/2)` on back-transformation; evaluation with both MSPE and the proxy-robust, asymmetric QLIKE loss; Mincer–Zarnowitz efficiency regressions with Wald tests.

## References

Key sources: Corsi (2009, *J. Financial Econometrics*); Parkinson (1980, *J. Business*); Taylor (2007, *Asset Price Dynamics, Volatility, and Prediction*); Hansen & Lunde (2005, *J. Applied Econometrics*); White (2000, *Econometrica*); Ledoit & Wolf (2008, *J. Empirical Finance*); Diebold & Mariano (1995, *JBES*); Fama (1970, *J. Finance*).

## Disclaimer

Academic research code. Nothing here is investment advice — indeed, the central result of Task 2 is that these signals *lose* money.
