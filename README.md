# My-Finance-Portfolio

Financial modelling, valuation and portfolio risk projects, built using Python, Excel and Bloomberg Terminal.

## Projects

### 🔢 Derivatives Pricing
Priced European and American equity options using seven numerical methods in Python, including binomial trees, finite-difference PDE methods, Monte Carlo simulation, and Longstaff-Schwartz. Validated every method against the Black-Scholes closed-form solution to within 10⁻⁴, and diagnosed sources of model bias including discrete-monitoring bias and volatility miscalibration.

[View project](./derivatives-pricing/)

### 💼 Portfolio Risk & Active Management (MSCI Europe)
Used Python to calculate VaR, Sharpe ratio and correlation matrix across FTSE stocks, then constructed and managed a €1bn active European equity portfolio against the MSCI Europe Index using Bloomberg Terminal, applying multi-factor risk decomposition, Brinson attribution and Monte Carlo VaR. The portfolio delivered 8.47% against a 9.09% benchmark, at lower total risk and VaR than benchmark; Brinson attribution diagnosed the shortfall as allocation-driven, concentrated in Novo Nordisk, rather than a failure of risk control.

[View project](./Active%20Portfolio%20Management.pdf) · [Risk model](./Portfolio%20Risk%20Model.pdf)

### 📊 Rolls Royce DCF Valuation
Built a discounted cash flow model projecting free cash flow growth from £3.3bn to £5.1bn over five years, calculating WACC via CAPM at 8.5%. Discounted projected cash flows and terminal value back to present value, arriving at an intrinsic value of 1,320p against a market price of 1,110p, a 19% upside supporting a Buy recommendation.

[View project](./Rolls-Royce_DCF_Model.xlsx)

### 📉 Merck (MRK) Time Series Analysis
Empirical time-series econometrics on 26 years of daily MRK equity data (2000–2025): tested return properties (stationarity, normality, white-noise), evaluated rolling-window ARMA/ARMA-GARCH trading strategies out-of-sample against passive and technical benchmarks net of transaction costs, and forecast the high-low price spread using Log-HAR and Log-ARMA models. Found returns unpredictable at short horizons but volatility strongly persistent and forecastable, consistent with weak-form market efficiency.

[View project folder](./merck-time-series-analysis/)

## About me

MSc Finance student at the University of Manchester, CFA Level 1 candidate. Previously a Data & Insights Analyst at Axtria, working with SQL and Python on financial and commercial analysis for global healthcare clients.

[LinkedIn](https://www.linkedin.com/in/mridulmacker/)
