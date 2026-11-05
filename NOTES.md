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
  - `Fear Index Europe __ VSTOXX.xlsx`: columns `Date`, `Indexvalue`. Phase 12
    is skipped with a warning if the file is missing.
  - Both are third-party data: **the repository must stay private.**
- Outputs in `output/` (not versioned): `tables/*.csv`, `figures/*.png`, the
  rolling cache `rolling_*.rds` and `sessionInfo.txt`.
- Caches: the price snapshot and the rolling results are reused while
  `USE_CACHE = TRUE`. The rolling cache file name encodes the methodological
  parameters, so changing them forces a new estimate. Changing the code
  without changing those parameters does not: delete the `.rds` in `output/`
  or set `USE_CACHE <- FALSE`.
- Versions used so far: R 4.3.3 (cloud) and the user's local R; MSGARCH 2.51,
  rugarch 1.5.6. The MSGARCH behaviour below was verified on 2.51 only.
- Run time: the rolling stage took 4.9 min on the user's Mac without warm
  start; warm start roughly doubles the number of fits.

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

## Empirical findings so far

From the user's local run of 2026-11-03 (T = 2489, `MIN_QUOTED_SHARE = 0.5`),
*before* the fixes for points 3–4 above and before warm start. Re-check after
the next run.

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

## Open issues and next steps

1. Re-run Phases 7–10 with the latest fixes. Check `SE_OK`, Table 4b, the
   `FitML` messages, `roll_checks` (`N_failed`, `N_maxit`, `Share_warm`) and
   Table 1b for BPE.
2. BPE.MI: data artefact or genuine feature? Decide on cleaning or exclusion.
3. Q1 holiday handling (`interp` vs `common`) and Q2 rolling sample: user
   decisions, both change every estimate.
4. Q3 benchmark for the paper (sGARCH with mean vs GJR K = 1).
5. **Main methodological issue**: "high" means different things across banks
   (outlier day, long crisis period, break), and in the rolling stage regimes
   last only days. The count and the sector indices add up states that are not
   comparable. Options to discuss: classification based on conditional
   volatility, persistence restrictions, longer windows, smoothing.
6. Mean specification: residual autocorrelation for INGA, ACA, BNP and GLE.
7. Paper caveats: static end-of-sample weights (look-ahead), selection of
   banks listed today (survivorship), KPSS rejections.
8. Rolling robustness: if warm start does not remove most failures, consider
   more starts or a higher `maxit`, and compare with the cached single-start
   results.
