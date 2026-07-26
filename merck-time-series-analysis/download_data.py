"""Download daily MRK OHLC data and write the two CSVs used by the tasks.

Usage:
    python src/download_data.py [--start 2000-01-01] [--end 2025-12-31]

Writes:
    data/mrk_close.csv   (Date, Close)            -> tasks 1 and 2
    data/mrk_ohlc.csv    (Date, High, Low, Close) -> task 3
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yfinance as yf


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ticker", default="MRK")
    ap.add_argument("--start", default="2000-01-01")
    ap.add_argument("--end", default="2025-12-31")
    ap.add_argument("--outdir", default="data")
    args = ap.parse_args()

    df = yf.download(args.ticker, start=args.start, end=args.end,
                     auto_adjust=True, progress=False)
    if df.empty:
        raise SystemExit("Download returned no data - check ticker/dates.")
    if hasattr(df.columns, "droplevel") and df.columns.nlevels > 1:
        df.columns = df.columns.droplevel(1)
    df = df.reset_index()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    df[["Date", "Close"]].to_csv(out / "mrk_close.csv", index=False)
    df[["Date", "High", "Low", "Close"]].to_csv(out / "mrk_ohlc.csv", index=False)
    print(f"Wrote {len(df)} rows to {out}/mrk_close.csv and {out}/mrk_ohlc.csv")


if __name__ == "__main__":
    main()
