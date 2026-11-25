# Publication and provenance

This project accompanies Mridul Macker's 2026 MSc Finance dissertation at Alliance Manchester Business School, University of Manchester.

## Public PDF

The public PDF contains all 53 pages of the supplied dissertation. On the cover, the student ID has been replaced with the author's name. Descriptive title and author metadata have been added. Text on pages 2–53 was checked against the supplied PDF and is unchanged. The original uploaded dissertation remains separate from this public copy.

## Research code

All 41 supplied `.m` files are included. Numerical code has not been rewritten, and original filenames are retained. One portability edit was made: the absolute OneDrive path in `kaopen.m` was replaced by `kaopen_2023.xls`, to be placed in the local working folder. Other scripts are byte-for-byte copies of the supplied files.

Two exact duplicate pairs are retained:

| Preferred file | Historical alias |
| --- | --- |
| `step10_korea_d4.m` | `stw.m` |
| `step32_results_51.m` | `results5.1.m` |

The guide documents intermediate versus corrected analyses, including the later date-aligned crisis-robustness script `recheck.m`. The archive is not presented as a validated one-command reproduction of the final dissertation. The source CSV was not supplied. No missing results or code have been invented.

## Outputs and data

The four PNG figures are unchanged copies of the author's uploaded outputs. Filenames reflect the research working directory and do not match the numbering of every figure in the final dissertation.

The 24 supplied `.mat` data/result snapshots are excluded from the public package, along with screenshots and unrelated portfolio files. Source definitions and local input requirements are documented in `DATA.md`.

## Checks performed

- Verified all supplied MATLAB source files are represented, with the sole path edit recorded above.
- Checked the public PDF page count, cover text and unchanged body text.
- Checked documentation links against the publication file manifest.
- Verified published binary file hashes against the prepared local files.
- MATLAB model execution, toolbox compatibility and numerical agreement with the dissertation have not been revalidated.
