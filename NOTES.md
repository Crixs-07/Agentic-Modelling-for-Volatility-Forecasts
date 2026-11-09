# Working notes: MS-GARCH systemic risk analysis

Running log of decisions, verified package behaviour and open issues for
`R/ms_garch_systemic_risk.R`. Read this before changing the script. Update it
whenever a decision is taken or an issue is closed.

## How to run

From the repository root:

```bash
Rscript R/ms_garch_systemic_risk.R > output/run.log 2>&1
```

- Inputs in `data/`:
  - `prices_raw_2017-01-01_2026-09-30.rds`: frozen Yahoo Finance snapshot
    (list of 21 `xts` objects). Yahoo rejects requests from cloud servers
    (HTTP 429), so cloud runs depend on this file. Yahoo also revises adjusted
    prices after every dividend, so the snapshot is what makes runs reproducible.
    Source: object `prices_raw` in `MS_Data.RData` (private repository
    Crixs-07/Details-on-MS-project), saved unchanged with `saveRDS`. Checked on
    2026-11-09: same common sample (2017-01-04 to 2026-09-29, T = 2489) and same
    excluded dates as the user's run of 2026-11-03.
  - `Fear Index Europe __ VSTOXX.xlsx` (main script, Phase 12) and `h_v2tx.txt`
    (`R/revision_tests.R`): stoxx.com is blocked by the cloud proxy (HTTP 403),
    so in the cloud both files were rebuilt from the object `VSTOXX` in
    `MS_Data.RData` (columns `Date`, `Symbol`, `Indexvalue`, 1999-01-04 to
    2025-11-07): enough for the sample to October 2025, not for the main
    script's sample to September 2026. Phase 12 is skipped with a warning if
    the `.xlsx` is missing.
  - Both are third-party data: **the repository must stay private.** Do not
    commit `data/` while the repository is public.
- Outputs in `output/` (not versioned): `tables/*.csv`, `figures/*.png`, the
  rolling cache `rolling_*.rds` and `sessionInfo.txt`.
- Caches: the price snapshot and the rolling results are reused while
  `USE_CACHE = TRUE`. The rolling cache file name encodes the methodological
  parameters, so changing them forces a new estimate. Changing the code
  without changing those parameters does not: delete the `.rds` in `output/`
  or set `USE_CACHE <- FALSE`.
- Versions used so far: R 4.3.3 (cloud) and the user's local R; MSGARCH 2.51,
  rugarch 1.5.6. The MSGARCH behaviour below was verified on 2.51 only.
  In the cloud (2026-11-09): rugarch 1.5-6 needs Rcpp >= 1.1.1 (installed
  1.1.2), FinTS 0.4-9, moments 0.14.1. The environment's setup script does not
  install MSGARCH, rugarch, FinTS or moments yet.
- Run time: the rolling stage took 4.9 min on the user's Mac without warm
  start; warm start roughly doubles the number of fits. Cloud, 4 cores
  (2026-11-09): main script 15 min (rolling 12.8 min); `R/revision_tests.R`
  13.5 min (rolling 12.8 min) plus 6.6 min for P10.

## Model and methodological choices

Parameters live in Phase 0 of the script.

| Item | Current setting | Status |
|---|---|---|
| Sample | 2017-01-01 to 2026-09-30 (`to` is exclusive: last date 2026-09-29) | original |
| Banks / weights | 21 largest euro-area banks; static market caps of June 2025 (USD bn) | original; look-ahead and survivorship to be acknowledged in the paper |
| `HOLIDAY_METHOD` | `"interp"` (linear interpolation of prices on single-market closures) | original; **open (Q1)**, reviewer recommends `"common"` |
| `MIN_QUOTED_SHARE` | `0.5`: drop dates on which fewer than half of the banks are quoted | decided by the user |
| Non-trading days | dates with zero volume for every quoted ticker are always dropped (2019-05-01, 2019-12-25) | decided |
| `ROLL_SAMPLE` | `"ticker"`: rolling on each bank's own sample | original; **open (Q2)** |
| `RET_SCALE` | `100` (percent returns) | package calibration; *not* chosen for a higher likelihood (see findings) |
| MS model | MS(2)-GJR-GARCH(1,1), Student-t, `nu` shared across regimes, zero mean | original |
| Benchmarks | sGARCH(1,1)-t with constant mean (rugarch, original); GJR-GARCH(1,1)-t zero mean, K = 1 (MSGARCH, like-for-like) | **open (Q3)**: which one the paper reports |
| High regime | regime with the larger closed-form unconditional variance (always regime 2, see MSGARCH notes) | decided |
| Rolling | window 504, step 10, horizon 10; filtered probability at the window's last observation × Pʰ (no look-ahead) | original |
| Full-sample starts | `MS_MULTISTART`: default start + 3 persistent starts built from the K = 1 GJR fit; keep the highest LL | decided |
| Rolling starts | `ROLL_WARM_START`: default start + previous window's estimate; keep the highest LL | decided, not yet run |
| Sector indices | equal-weighted mean; cap-weighted mean renormalised over the banks available on each date | decided |
| Information criteria | totals on the decimal scale, over observations 2..T for every model | decided |
| Events (plots, Table 1) | 2020-03-18 COVID/PEPP, 2022-02-24 Ukraine, 2023-03-10 SVB, 2023-03-15 Credit Suisse, 2025-02-04 US tariffs, 2025-04-02 "Liberation Day" tariffs | decided by the user on 2026-11-09 (all events mentioned in the paper, the test script and `benchmark_prelim`) |

## Bugs fixed in the original script

1. `full_join` appended unmatched dates at the bottom: prices were not in date
   order before interpolation and before the rolling returns. Now sorted.
2. `UncVol(fit)` returns one number, so `which.max()` always picked regime 1,
   which `FitML` makes the *low* regime. The old `Prob_High` was Pr(low) in
   every window, and every Phase 11 result was inverted.
3. Rolling AIC/BIC were always `NA` (`summary()` returns `$AIC`/`$BIC`, not `$IC`).
4. rugarch reports AIC/BIC divided by T and MSGARCH reports totals, yet the two
   were compared directly.
5. `auto.arima(max.p = 0, max.q = 0)` always returned ARMA(0,0); removed with no
   change in results.
6. The cap-weighted index was not renormalised when banks were missing.
7. An arbitrary deletion of the last row of `daily_high_prob` was removed.
8. VSTOXX was read from two different absolute paths; it now comes from `data/`.
9. The end date `2026-09-31` did not exist.

## MSGARCH 2.51 behaviour (verified in the package source)

1. **`UncVol()`** (`R/UncVol.R`) returns a single simulated value for the whole
   MS process, not one per regime. The simulation starts from `data = c(1, 1)`
   and averages steps 1000–5000 (`tmp[nburn:nahead]` with nburn 5000 and
   nahead 1000), so it is biased upwards for highly persistent regimes.
   Use the closed form α0 / (1 − α1 − 0.5·α2 − β) (`src/gjrGARCH.h`,
   `set_vol`), as `regime_uncond_vol()` does.
2. **Regime ordering**: `FitML` sorts regimes by that closed-form unconditional
   variance (`f_sort_par` in `R/Utils.R`). With K = 2 and identical specs,
   regime 1 is always low and regime 2 always high.
3. **The likelihood excludes the first observation** (`calc_lndMat` in
   `src/MSgarch.h`, `eval_model` in `src/SingleRegime.h`): n = T − 1 in every
   scale adjustment and information criterion.
4. **`ctr$par0` combined with `regime.const`** is broken. `FitML` calls
   `f_remove_regimeconstpar()` on the already transformed parameters, and that
   function applies `abs()`, so ν ≈ 5 becomes a starting value of ν ≈ 97 (ν is
   bounded in [2.1, 100]). `fit_ms()` works around this by transforming the
   start itself and passing it straight to the optimiser.
5. **Optimiser**: `optim` BFGS (maxit 100) from a single start. The convergence
   code is discarded, and an `optim` error surfaces only as "$ operator is
   invalid for atomic vectors". `make_optim_bfgs()` records both.
6. **Standard errors** (`f_InferenceFun`, eigen of the Hessian) can fail after a
   successful optimisation, and the whole fit is then lost. `fit_ms()` re-fits
   without SE.
7. The Hamilton filter starts from the ergodic distribution, so the likelihood
   does not depend on regime labels.
8. `TransMat()` is row-stochastic: rows `t|k`, columns `t+1|k`;
   `P_2_1` = Pr(s_{t+1} = 1 | s_t = 2). `State()$SmoothProb` has T + 1 rows.
9. Bounds: ν in [2.1, 100]; GJR lower (1e-7, 1e-6, 1e-4, 0), upper
   (100, 0.9999, 10, 0.9999); the stationarity constraint is
   α1 + 0.5·α2 + β < 0.99999999.
10. **Main cause of failed fits**: when P approaches the identity matrix,
    `loadparam` (`src/MSgarch.h`) cannot invert the matrix it uses for the
    ergodic distribution ("matrix multiplication: problem with matrix
    inverse"). The error escapes `f_nll`, aborts `optim` and the window is
    lost. Proposal P8 in `R/revision_tests.R` treats the point as
    inadmissible (penalty 1e10). Fits that never hit the error are unchanged.
    On synthetic data this removed every failure (4 → 0).
11. **`State(fit, newdata = x)` appends `x` to the estimation data**
    (`State.MSGARCH_ML_FIT`: `data <- c(object$data, newdata)`). Passing the
    estimation window again as `newdata` duplicates it. To filter new days
    with fixed parameters, call `State(object = spec, par = fit$par, data = ...)`
    on the full series, as `R/revision_tests.R` does (verified 2026-11-09).

## Results of the cloud runs of 2026-11-09 (current code)

Logs in `output/run.log`, `output/revision_run.log` and
`output/revision_run2.log` (not versioned). Specification unchanged; only the
events were updated. These numbers supersede the section below.

**Main script** (sample 2017-01-04 to 2026-09-29, T = 2489):

- Full sample (Table 4): `SE_OK` and convergence for all 21 banks. Most banks
  now have persistent regimes (P11, P22 between 0.95 and 0.998). Two degenerate
  chains that switch almost every day: BPE.MI (P11 = 0.06, P22 = 0.11) and,
  new, BBVA.MC (P11 = 0.46, P22 = 0.08), chosen by the multi-start because it
  gains 5.8 in LL over the default start. The 2026-11-03 extreme-day regimes
  of ABN and RBI are gone.
- Multi-start gains (Table 4b): CBK 8.8, BBVA 5.8, INGA 4.9, SAN 3.3, ABN 2.0,
  KBC 1.2; below 1 for the others.
- Information criteria on the full sample: the K = 1 GJR has a lower BIC than
  the MS-GJR for most banks.
- Rolling: 0 failed windows out of 4158 (65 on 2026-11-03), up to 7
  iteration-limit stops per bank; warm start chosen in 65–80% of windows.
  Mean durations 1.6–25 days, long-run Pr(high) 0.19–0.73, mean cross-bank
  Spearman correlation of `Prob_High` from −0.002 (RBI) to 0.29 (BNP).
- Sector index: the first date (2018-12-21) has 20 banks and the last
  (2026-09-22) only 1, because each bank's rolling runs on its own sample
  (`ROLL_SAMPLE = "ticker"`): the index drops to 0.11 on that day.
- Diagnostics unchanged: Ljung–Box on residuals rejects for INGA, BNP, ACA, GLE
  (p < 0.05; INGA and ACA p < 0.01), on squared residuals for INGA and ABN;
  KPSS rejects for BBVA, SAB, SAN, UCG, BPE, ABN; ADF and PP reject a unit root
  for all. Table 1b shows nothing anomalous for BPE (zero returns 0.6%).

**`R/revision_tests.R`** (sample to 2025-10-31, T = 2258; `PROB_METHOD =
"filtered"`, `HIGH_RULE = "condvol"`):

- 0 failed windows (P8), 175 windows per series.
- P3: the two high-state rules disagree in 0–22% of windows (portfolio 21%,
  BIRG 22%).
- Table B.1 (BIC, same sample): the K = 1 GJR is best for 18 banks, sGARCH for
  BKT, SAB, ABN; the MS-GJR is never best.
- Table 1: Pearson correlation with RV22 0.32 (SPI_EW), 0.25 (SPI_cap), 0.12
  (portfolio); with VSTOXX 0.20, 0.16, 0.10; Spearman with VSTOXX close to 0.
- P9: SPI_EW − portfolio Prob(High) has mean 0.04 and sd 0.27; SPI > 0.5 with
  portfolio ≤ 0.5 on 13% of days, the reverse on 19%. Yearly means do not
  follow stress: the portfolio's mean Prob(High) is 0.26 in 2020 and 0.53 in
  2019, because "high" is defined relative to each 504-day window.
- P10: removing warm start changes the bank-level series a lot (Pearson
  0.48–0.81 with the warm-start series; high/low classification differs on
  10–26% of days); warm start is kept only where it has the higher LL.
  SPI_EW warm vs single start: Pearson 0.91.
  The single-start SPI correlates more with RV22 (0.43 vs 0.32). The
  probabilities therefore depend on the optimiser's starting point (local optima).
- Reproducibility: re-estimating 3 default-start windows per series
  reproduces the cache exactly (max difference 0).
- PNG figures of the test script have a transparent background (`ggsave`
  with `theme_minimal`).

## Findings from the run of 2026-11-03 (superseded)

From the user's local run of 2026-11-03 (T = 2489, `MIN_QUOTED_SHARE = 0.5`),
*before* the fixes for points 3–4 above and before warm start.

- **Multiple local optima.** In the scale check (Phase 9b), the true LL
  differences between percent and decimal fits have mixed signs (−0.65 to
  +0.32) while the transition matrices differ a lot (e.g. INGA window 801:
  P22 = 0.94 vs 0.13). No scale is better; single-start BFGS is unreliable.
  An earlier claim of "+4.9 in favour of percent" was an artefact of point 3.
- **Full sample (Table 4).** Regime types differ across banks:
  - isolated extreme-day regimes: ABN, RBI, BBVA, CBK (high regime lasting 1–4
    days, UncVol_R2 10–18% per day);
  - single breaks: SAN (one switch), MB;
  - weak classification: RCM 65–87 for CABK, SAB, MB, ACA, UCG and others.
- **BPE.MI** converges to a degenerate chain that switches almost every day
  (P11 ≈ 0.07, P22 ≈ 0.10) from every start. Table 1b checks whether this is a
  data problem.
- **Multi-start gains** (Table 4b, with the ν bug still present): INGA +4.9,
  CBK +3.8, CABK +0.8.
- **Rolling** (single default start): 65 failed windows out of 4158 (0–7 per
  bank), up to 17 iteration-limit stops (DBK). Mean durations are only 2–20
  days and the long-run Pr(high) is 0.23–0.52. Cross-bank Spearman correlation
  of `Prob_High` is low: the per-bank mean ranges from 0.04 to 0.27, and is
  close to 0 for ABN, RBI and BPE.
- **Diagnostics.** Ljung–Box shows residual mean autocorrelation for INGA and
  ACA (p < 0.01) and for BNP and GLE (p < 0.05), plus remaining ARCH effects
  for INGA and ABN. KPSS rejects stationarity at 5% for BBVA, SAB, SAN, UCG,
  ABN and BPE, while ADF and PP reject a unit root for every series.

## Paper revision (November 2026): answers from the code

`R/revision_tests.R` implements the proposals below without touching the main
script (outputs in `output/revision/`). Tested on synthetic data only so far.

- **Cap-weighted weights**: fixed for the whole sample. The vector `m_caps`
  (USD bn, "June 2025") comes from the original script. The code does not say
  where the figures were taken from: the source must be documented.
- **Rolling probabilities**: the main script estimates once per window and
  then propagates the filtered probability of the window's last day,
  p_T Pʰ for h = 1..10, *without* using the returns of those 10 days. The
  daily-looking series is therefore an h-step-ahead prediction that decays
  towards the ergodic probability and jumps at each re-estimation.
  Proposal P2: filter the 10 following days with fixed parameters,
  Pr(s_t | I_t), which still involves no look-ahead.
- **High regime when φ ≈ 1**: MSGARCH cannot return φ = α1 + 0.5·α2 + β = 1
  (constraint φ < 0.99999999), so a reported "1" is a rounded value.
  The rule α0 / (1 − φ) then labels the near-integrated regime "high"
  almost automatically, whatever its volatility level in the window.
  Proposal P3: rank regimes by their conditional volatility averaged with the
  smoothed probabilities. The test script records how often the two rules disagree.
- **Interpolation**: every internal gap is interpolated, whatever its cause.
  Gaps are at most 2 consecutive days, but not all are holidays (e.g.
  2018-09-20, no price for any of the 5 Spanish banks). The test script
  classifies every gap.
- **IQR outliers**: there is no IQR rule anywhere in the code. The test
  script produces a list (k = 3, flag only) to check by hand before the paper
  can say they were "checked and kept".
- **BIC**: the main script already compares sGARCH(1,1)-t (rugarch), GJR-t
  K = 1 and MS-GJR (Table 5). Numbers from the original script are not
  comparable (per-observation vs total IC, one-observation offset), so
  Table B.1 must be recomputed.
- **756-day window**: never run in this pipeline; `RUN_756` in the test script.
- **MSGARCH `State()` indices**: `FiltProb[t]` = Pr(s_t | I_t), with row 1 =
  ergodic distribution. `PredProb` has T + 1 rows: `PredProb[t]` =
  Pr(s_t | I_{t−1}) for t ≥ 2, and the out-of-sample one-step forecast is
  `PredProb[T + 1]`, not `PredProb[T]`. `SmoothProb` has T + 1 rows.
  `UncVol()` is not per regime (see point 1 of the MSGARCH section).
- **VSTOXX**: the official STOXX file `h_v2tx.txt` (daily, from 1999,
  semicolon-separated `Date;Symbol;Indexvalue`) is downloaded by the test script.

## `benchmark_prelim` (Crixs-07/Details-on-MS-project): ideas, not results

Values, tables and code in that folder are unverified and are not copied.
Each idea is re-implemented in `R/revision_tests.R` and recomputed.

| Idea | Problems found in the prelim code | Status here |
|---|---|---|
| Rolling on the equal-weighted portfolio | `UncVol(fit)[k]` (point 1: "high" = low regime); `State(fit, newdata = window)` duplicates the window (point 11); mean of log returns | P4 |
| Table 1: correlations, event maxima, lead–lag | — | P4 |
| SPI − portfolio gap and disagreement days | — | P9 (added 2026-11-09) |
| Yearly means | — | P9 (added 2026-11-09) |
| Warm vs single start | — | P10 (added 2026-11-09) |
| Random multi-start in the rolling | — | not adopted (open issue) |
| One-step-ahead predicted probability | indices of the R code | not adopted: the user keeps P2 for now |
| Cap on GARCH persistence (0.9995) | arbitrary threshold; Python model, not MSGARCH | not adopted: it changes the specification |
| Two more events (2023-03-15, 2025-04-02) | — | adopted |

## Open issues and next steps

1. Repository visibility: on 2026-11-09 the GitHub page of this repository
   still answered without authentication (public). `data/` stays uncommitted
   until it is private.
2. BPE.MI and now BBVA.MC: degenerate full-sample chains. Data artefact or
   genuine feature? Decide on cleaning or exclusion.
3. P10: the rolling probabilities depend on the starting values. Decide how
   to report this (more starts, single start as robustness, or a different
   specification).
4. Table B.1: the MS-GJR is never preferred by BIC on the paper's sample. The
   paper's motivation for the MS model must be revisited.
5. Ticker-specific rolling samples leave dates with 1 or 20 banks at the ends
   of the sector indices (see Q2): trim those dates or use a common sample.
6. PNG figures of `R/revision_tests.R`: add `bg = "white"` to `ggsave`.
7. Cloud setup script: install MSGARCH 2.51, rugarch 1.5-6 (with Rcpp >= 1.1.1),
   FinTS, moments.
8. Q1 holiday handling (`interp` vs `common`) and Q2 rolling sample: user
   decisions, both change every estimate.
9. Q3 benchmark for the paper (sGARCH with mean vs GJR K = 1).
10. **Main methodological issue**: "high" means different things across banks
   (outlier day, long crisis period, break), and in the rolling stage regimes
   last only days. The count and the sector indices add up states that are not
   comparable. Options to discuss: classification based on conditional
   volatility, persistence restrictions, longer windows, smoothing.
11. Mean specification: residual autocorrelation for INGA, ACA, BNP and GLE.
12. Paper caveats: static end-of-sample weights (look-ahead), selection of
   banks listed today (survivorship), KPSS rejections.
13. Rolling robustness: warm start (main script) and P8 (test script) removed
   every failure; see issue 3 for the dependence on starting values.
