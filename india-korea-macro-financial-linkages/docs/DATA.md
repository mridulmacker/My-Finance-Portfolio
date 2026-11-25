# Data sources and local inputs

## Scope

The dissertation studies monthly observations from 2000–2025. Korea's M3 series starts in November 2003, so its specifications have a shorter sample. The financial-openness comparison covers 2000–2023.

The original `master_dataset.csv` was not supplied with this publication. Local MATLAB dataset and result snapshots are not published. In particular, the transformed MATLAB dataset contains source-derived observations and is not a substitute for documenting and obtaining the source data. Use the dissertation's Appendix A, especially Table A2, for the recorded sources and definitions.

## Source map

| Series | Source recorded in the dissertation |
| --- | --- |
| Nifty 50, KOSPI, MSCI World | Bloomberg Terminal |
| USD/INR and USD/KRW | Bloomberg Terminal |
| Industrial production, CPI and M3 | OECD |
| India call money rate | Reserve Bank of India |
| Korea short-term interest rate | OECD |
| US three-month Treasury-bill rate | US Department of the Treasury |
| Capital-account openness | Chinn–Ito index, 2023 dataset |

Source naming and transformations follow the submitted dissertation. No source data are downloaded by the supplied scripts; the scripts expect a preassembled local input file. Access to the original Bloomberg series may require an appropriate subscription.

## Expected master CSV columns

Create `code/master_dataset.csv` with these column names. The spelling `usd_kwr` is retained from the original code, even though the currency code is KRW.

| Column | Meaning / input form |
| --- | --- |
| `date` | Monthly dates that MATLAB `datetime` can parse |
| `nifty`, `kospi`, `msci_world` | Equity index levels |
| `india_ip`, `korea_ip` | Industrial-production levels |
| `usd_inr`, `usd_kwr` | Exchange-rate levels |
| `india_cpi`, `korea_cpi` | Consumer-price index levels |
| `india_m3`, `korea_m3` | Broad-money levels |
| `india_call_rate`, `korea_ir` | Domestic interest-rate levels |
| `us3m_tbills` | US three-month Treasury-bill rate |
| `india_wpi` | Additional wholesale-price series referenced by the original preparation scripts |

Rows should be sorted chronologically with one observation per month. Use missing numeric entries for unavailable data, rather than filling the early Korea M3 sample with invented values. Log-transformed input series must be positive wherever observed.

The preparation scripts create `ln_` levels, `dln_` log differences and the differenced rate columns `d_indiacallrate`, `d_koreair` and `d_us3mtbills`. Model scripts select complete observations separately for each specification.

## Financial-openness input

`code/kaopen.m` expects `code/kaopen_2023.xls`, with columns `country_name`, `year` and `ka_open`. Its original computer-specific absolute path has been replaced by this relative filename. The plotting logic is otherwise unchanged.

## Generated files

With inputs and dependencies available, scripts can write `master_dataset_transformed.mat`, individual `india_*_results.mat` and `korea_*_results.mat` files, `vecm_*.mat` models, and figure exports into the current working directory. These local data and generated files are ignored by Git. The four figures under `figures/` are saved outputs supplied by the author, rather than outputs regenerated during publication.
