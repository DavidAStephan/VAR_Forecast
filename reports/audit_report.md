# SOE-BVAR Suite — Audit Synthesis Report

## 1. Executive summary

The SOE-BVAR forecasting suite is, on the whole, **correctly implemented and well-engineered**. The core econometrics — block-recursive conjugate NIW, independent Normal-IW Gibbs, Villani steady-state, and stochastic-volatility engines; the Minnesota/SOC/DIO priors and GLP lambda selection; the iterated predictive simulation; the density-combination schemes; and the expanding-window OOS harness — were independently verified and found mathematically sound, with strong empirical confirmation (block exogeneity exact, no-look-ahead enforced, weights normalised, transform round-trips exact).

Directly addressing the three concerns:

- **(i) Are the VAR models (incl. COVID corrections) implemented properly?** Mostly yes. The four engines and the Lenza-Primiceri COVID volatility scaling are correct, internally consistent, and no-look-ahead-safe. There is **one genuine, high-severity correctness defect**: the Student-t predictive draws in the **SV engine and the UCSV benchmark** are not standardised — they draw raw `rt(df=nu)` (variance `nu/(nu-2)`) and scale by `exp(h/2)`, whereas estimation calibrates `exp(h)` to be the conditional variance. This systematically over-disperses every SV/UCSV predictive density by `sqrt(nu/(nu-2))` (~4–12% in tranquil periods, up to ~30% in COVID windows), miscalibrating their density scores and combination weights.
- **(ii) Is the data correct?** Yes, with caveats. Ingestion is fail-loud (no silent substitution), quarter stamping is consistent, and the synthetic DGP is a genuine block-exogeneity ground truth. Two real but **latent/low-severity** issues exist: a fragile `expected = max(n)` coverage heuristic that could collapse a quarterly-native series on a duplicate row, and a config start-date/panel-stamp off-by-one-quarter that is currently masked by the data-bound 1997Q4 start.
- **(iii) Is the forecast averaging working correctly?** Yes. Equal, log-score, and optimal-pool schemes are correct, with weights non-negative and summing to exactly 1 across every cell, and strict no-look-ahead. The only confirmed issues are in the **BMA scheme, which is documented diagnostic-only and not used in production**: a mislabelled "one-step" horizon for medium/far buckets and a `1e-12` weight-floor that distorts BMA's pooled log density.

**Overall confidence: High.** The production forecast path (VAR engines + COVID corrections + averaging) is sound except for the SV/UCSV t-variance bug, which is the one finding that materially affects reported production output and should be fixed promptly. Most remaining findings are nits, latent fragilities, or confined to diagnostic-only code.

---

## 2. Severity-ranked findings

| Sev | Status | Subsystem | File:lines | Issue |
|---|---|---|---|---|
| High | confirmed | SV engine | R/forecast.R:145–147 | Unstandardised `rt(df=nu)` inflates SV predictive variance by `nu/(nu-2)` |
| High | confirmed | UCSV benchmark | R/benchmarks.R:224–225 | Identical unstandardised-t variance bug in UCSV transitory shock |
| High | confirmed | OOS scoring/reporting | R/evaluate.R:222–238 (via R/report.R:81–89) | ex-COVID exclusion is a no-op (exact-date vs quarter-stamp); -Inf log scores contaminate headline tables |
| Medium | confirmed | Data: panel construction | R/data_sources.R:165–175 | `expected = max(n)` coverage heuristic can NA-collapse a quarterly-native series on a duplicate obs |
| Medium | disputed | Engine: Villani ss | R/engines.R:279–289 | Foreign steady-state `psi` contaminated by domestic data through full Sigma (breaks exact block exogeneity for the means) |
| Low | confirmed | Data: panel construction | R/data_sources.R:252–253 | Config mid-quarter start date vs start-of-quarter panel stamps → off-by-one-quarter (latent) |
| Low | confirmed | Combination (BMA, diagnostic) | R/combine.R:83–92 | BMA labelled "one-step" but uses h=5/h=9 for medium/far buckets |
| Low | confirmed | Combination (BMA, diagnostic) | R/combine.R:140–141 | Pooled log density floors zero weights at 1e-12 (distorts BMA only) |
| Low | confirmed | Predictive simulation | R/forecast.R (gibbs/ss/sv storage) | Forecast draws use only first 1000 of 1500 stored posterior draws; benchmarks recycle 1..200 |
| Low | confirmed | OOS harness | R/evaluate.R:280–287, 244–263 | DM Newey-West assumes contiguous origins; non-contiguous common sets misalign lags (latent) |
| Low | confirmed | OOS self-tests | R/evaluate.R:304–331 | No-look-ahead/reproducibility self-tests cover only one member at one origin |
| Low | confirmed | Data acquisition | R/data_sources.R:178–179 | Misleading docstring claims a FRED fallback that does not exist |
| Low | disputed | Benchmark: rw | R/benchmarks.R:10–12 | `sd(diff(y)*w)` is not a proper weighted sd |
| Low | disputed | COVID corrections | R/evaluate.R:51–58, 78–80 | Benchmarks skip the two-pass COVID-scale refinement applied to VARs |
| Nit | confirmed | Engine: conj_br | R/engines.R:120–131 | `conj_br` assumes nd≥1; degenerate nf==M partition would mis-index (never triggers) |
| Nit | confirmed | Transforms | R/transforms.R:61–67 | Dead code: `year_ended()` helper never called (ye computed inline) |
| Nit/disputed | disputed | Priors / conj_br | R/engines.R:130 | Domestic-block omega hardcodes lag_decay=1 (latent inconsistency) |
| Nit/disputed | disputed | Engine: conj_br | R/engines.R:144–151 | Domestic SOC keep-filter conflates zero-ybar SOC dummy with DIO branch |

---

## 3. Detailed findings

### Concern (i): VAR models and COVID corrections

#### [HIGH, confirmed] SV engine: unstandardised t-draw inflates predictive variance — R/forecast.R:145–147
**What is wrong.** In estimation, `fit_sv` whitens residuals with `w <- exp(-h/2)/sqrt(mix)` (R/engines.R:393–394), where `mix` is stochvol's t-mixing weight `tau` with `E[tau]=1`. This calibrates `exp(h)` to be the conditional variance. But `simulate_paths.post_sv` draws `eps <- rt(1, df=nu_i)` and forms `exp(hs[i]/2)*eps`. Since `Var[rt(nu)] = nu/(nu-2) > 1`, each predictive draw has variance `exp(h)*nu/(nu-2)`, not the calibrated `exp(h)`.

**Why it matters.** Predictive SD is inflated by `sqrt(nu/(nu-2))` at every horizon for both SV VAR members (`ndraw=1500`, live with `covid.sv_t_errors=TRUE`). This over-dispersion miscalibrates CRPS, log score, and PIT, and biases the optimal-pool / log-score combination weights against otherwise sound models. Magnitude ~4–12% SD at `nu≈10–25`; up to ~30% in COVID windows where `nu` falls.

**Evidence.** Monte Carlo (2e6 draws): `Var(sqrt(tau)*z)=0.999` while `Var(rt(5))=1.664=nu/(nu-2)`; stochvol's own `svsim` docs confirm conditional volatility is `sqrt(nu/(nu-2))*exp(h/2)`. The Gaussian branch (`rnorm` when `nu` non-finite) is correct; only the finite-`nu` t-branch is defective. No README/DECISIONS entry justifies it.

**Fix.** `eps <- if (is.finite(nu_i) && nu_i > 2) rt(1, df=nu_i)*sqrt((nu_i-2)/nu_i) else rnorm(1)`. `nu>2` is guaranteed by the `nu-2 ~ Exp` prior.

#### [HIGH, confirmed] UCSV benchmark: identical unstandardised-t variance bug — R/benchmarks.R:224–225
**What is wrong.** `fit_ucsv` estimates the transitory conditional variance as `exp(he)*mix_e` with stochvol `tau` (mean 1), so `exp(he)` is the calibrated transitory variance. `simulate_paths.post_ucsv` then draws `eps <- rt(1, df=nu)` and forms `tau + exp(he/2)*eps`, over-inflating the transitory variance by `nu/(nu-2)`.

**Why it matters.** UCSV is both a scoring benchmark every VAR must beat **and** a pool member. Over-dispersing it systematically widens its intervals, biasing CRPS/log-score, making the benchmark artificially easy to beat, and contaminating pool weights.

**Evidence.** The suite's own predictive-SD diagnostic at R/benchmarks.R:189 computes `pred_sd = sqrt(var(tauT) + mean(exp(heT)) + mean(s2uT))` — treating `exp(heT)` as the transitory variance with no `nu/(nu-2)` factor — confirming developer intent and the simulator's inconsistency. Inflation `nu/(nu-2)` is ~1.67 at `nu=5`, ~1.15 at `nu=15` on the transitory component (one of three variance contributors).

**Fix.** Same as above: standardise the t-draw via `*sqrt((nu-2)/nu)` for finite `nu`.

#### [MEDIUM, disputed] Villani ss: foreign steady-state `psi` contaminated by domestic data — R/engines.R:279–289
See Section 4 (Disputed). Summary: the `psi` Gibbs step uses the full M×M `Sigma_inv`, so foreign steady states are nudged by domestic residuals, breaking exact block exogeneity *for the means*. Verifiers split on whether this is a defect or the correct Villani full conditional. Affects only one of seven members and only the steady-state means (lag dynamics remain exactly block-exogenous).

#### [NIT, confirmed] conj_br assumes nd≥1 — R/engines.R:120–131
A degenerate `nf == M` partition would mis-index via `(nf+1):M` (which yields a descending `c(M+1, M)` in R, not an empty range) and produce a malformed domestic block. Never triggers — every configured set includes domestic variables. Fix: add `stopifnot(nd >= 1, nf >= 1)`.

---

### Concern (ii): Data

#### [MEDIUM, confirmed] Fragile `expected = max(n)` coverage heuristic — R/data_sources.R:165–175
**What is wrong.** Coverage is judged against `expected <- max(n)` (count in the fullest quarter). For the quarterly-native series that actually bind this panel (ABS national accounts, RBA trimmed-mean CPI, GDPC1, WPI, FRERTWI), every quarter normally has n=1. If any single quarter ever carries 2 observations (a republished print, a vintage/revision row, or two dates in one quarter), `expected` jumps to 2 and every genuine single-obs quarter gets coverage `1/2 = 0.5 < 0.8 = min_coverage` and is silently NA'd — collapsing essentially the whole series. There is no date de-duplication before binning.

**Why it matters.** A silent, severe failure mode for vintage provider data. Verifiers split on severity: one rated medium (silent, severe), one rated low (downstream NA guards in transforms.R:32–34 and data_sources.R:243–244 would fail loud rather than silent, and the trigger requires an upstream anomaly). Either way it is a genuine fragility.

**Evidence.** Reproduced: a quarterly series with one 2-obs quarter NAs all subsequent single-obs quarters. `median(n)` is robust in both directions (keeps single-obs quarters; still flags a 2-of-3-month trailing partial). The existing test (test-data-real.R:21–27) does not cover the duplicate-row case.

**Fix.** Replace `expected <- max(n)` with `expected <- stats::median(as.integer(n))` (or the mode), and/or de-duplicate per quarter before counting.

#### [LOW, confirmed] Config start-date vs panel-stamp off-by-one-quarter — R/data_sources.R:252–253
`config$data$start = "1990-03-01"` (intent "1990 Q1", mid-quarter) but the panel is stamped at quarter START (`cut(date,"quarter")` → 1990-01-01). The filter `panel$date >= as.Date("1990-03-01")` drops the intended 1990Q1 and starts at 1990Q2. Latent today because the balanced trim binds the start at 1997Q4 (README documents start:1990 as "aspirational"). The end filter is only coincidentally safe. Fix: floor config boundary dates to start-of-quarter before filtering, or standardise config on start-of-quarter dates.

#### [LOW, confirmed] Misleading docstring: nonexistent FRED fallback — R/data_sources.R:178–179
The `download_real_data` header claims it will "fall back per-series to FRED if a key is available." No such fallback exists; `fetch_series` dispatches strictly on the configured provider and a failed fetch errors loudly (lines 216–217). Comment-only; could mislead a maintainer into assuming silent substitution is possible (the opposite of reality). Fix: delete/correct the clause.

---

### Concern (iii): Forecast averaging

#### [LOW, confirmed] BMA mislabelled "one-step" for medium/far buckets — R/combine.R:83–92
`tr1 <- tr[tr$h == min(tr$h), ]` gives `min(tr$h) = 1` only for the near bucket; medium/far use h=5/h=9 predictive likelihoods. The in-code comment and README D10 call this "one-step." **Impact limited: BMA is explicitly diagnostic-only and not used in production forecasts** — this is a label/documentation mismatch. Fix: filter on `tr$h == 1` (fallback to equal when empty), or correct the label.

#### [LOW, confirmed] Pooled log density floors zero weights at 1e-12 — R/combine.R:140–141
`lw <- log(pmax(w, 1e-12))` lets zero-weight members contribute spurious mass to the pooled log density, whereas the CRPS path drops them (rmultinom) and PIT zeros them (`sum(w*pit)`). With `shrink_kappa=0.30`, logscore/pool never produce exact zeros, so **only BMA (unshrunk, diagnostic-only) is affected** — but there the distortion is material (one verifier found worst-case ~140+ log points where a zero-weight well-forecasting member dominates the floored background). Real defect, but confined to diagnostic output. Fix: `pos <- w > 0; logsumexp(log(w[pos]) + ss$logdens[pos])`.

#### [LOW, confirmed] Forecast draws never sample the last third of the posterior — R/forecast.R (gibbs/ss/sv)
`forecast_draws=1000` but gibbs/ss/sv store `ndraw=1500`; the recycling index `k <- ((d-1) %% post$ndraw)+1` only ever uses stored draws 1..1000, discarding 1001..1500. Benchmarks (`bench_ndraw=800 < 1000`) recycle draws 1..200 twice. No inferential bias (post-burn-in draws are exchangeable; shocks redrawn per path), but wasted compute and mildly inflated benchmark draw correlation. Note: **conj_br is exempt** (it stores `ndraw = forecast_draws = 1000`), so the finding's blanket "VAR engines" phrasing and its cited lines 99–100 are inaccurate for that engine. Fix: set `forecast_draws == ndraw`, or systematically subsample posterior indices.

#### [NIT, confirmed] Dead code: `year_ended()` never called — R/transforms.R:61–67
The tested rolling-4-quarter-sum helper is invoked nowhere in production; ye scoring is computed inline in `score_member` (R/evaluate.R:166–178). Two equivalent implementations that could drift. Fix: route the inline path through the helper, or delete the helper.

---

### Cross-cutting (scoring / harness)

#### [HIGH, confirmed] ex-COVID score exclusion is a no-op; -Inf log scores contaminate headline tables — R/evaluate.R:222–238 (triggered from R/report.R:81–89)
**What is wrong.** `summarise_scores()` excludes COVID quarters with `as.Date(scores$date) %in% as.Date(exclude_dates)` — an **exact** date match. Real-data realizations are stamped at quarter START (2020-04-01, 2020-07-01) while `cfg$covid$quarters` are mid-quarter (2020-03-01/06-01/09-01). These never coincide, so the exclusion removes nothing on real data. This is precisely the quarter-stamp gotcha that `covid.R` deliberately neutralises with `.qidx()` (R/covid.R:24–31), whose own comment warns that exact-date matching "would silently disable the treatment."

**Why it matters.** The 2020Q2/Q3 realizations fall outside the predictive kernel-density support, producing -Inf log densities. `mean(..., na.rm=TRUE)` does **not** drop -Inf, so ~4.4% of summary cells (62/1410) report -Inf mean log density in **both** `scores_by_horizon.csv` and the supposedly clean `scores_by_horizon_excovid.csv`, across 10 members. The headline log-score comparison is unusable and the ex-COVID safeguard is silently dead. (Note: the underlying draws, posteriors, and DM losses are unaffected — DM already filters non-finite at R/evaluate.R:246 — so this is a reporting/summary defect, not corruption of the forecasts themselves.)

**Evidence.** Reproduced on the cached store: `summarise_scores(allscores)` and `summarise_scores(allscores, exclude_dates=excl)` produce byte-identical logdens columns; applying the `.qidx`-based exclusion drops the -Inf cell count from 62 to 0.

**Fix.** Match at quarter granularity: `!(.qidx(scores$date) %in% .qidx(exclude_dates))`. Additionally guard log-score aggregation against -Inf (treat non-finite as NA before averaging, or report a trimmed mean / median log score).

#### [LOW, confirmed] DM Newey-West assumes contiguous origins — R/evaluate.R:280–287, 244–263
`dm_test` computes lag-l autocovariances positionally; `dm_vs_reference` passes the differential over `common <- intersect(base$origin, alt$origin)`. If an interior origin were missing (e.g. an SV failure) the common set becomes non-contiguous and the NW correction is applied at the wrong calendar lag, biasing the DM variance. Not triggered in the current run (all members share contiguous origins; the harness runs every member over identical origins with no interior skips), and only affects a diagnostic test's variance, not any forecast. Fix: assert contiguity, or index autocovariances by actual origin spacing.

#### [LOW, confirmed] Self-tests cover only one member at one origin — R/evaluate.R:304–331
`test_no_lookahead()` and `test_reproducibility()` exercise only `all_members(cfg)[[1]]` (small_minn, gibbs) at the first origin. The shared slice point in `harness_forecast` means the no-look-ahead structural guarantee holds for all engines, but per-engine RNG/seeding paths (SV, benchmarks) are never checked, so a reproducibility regression there could pass CI. Coverage gap, not a demonstrated bug. Fix: loop both tests over representative members and a COVID-spanning origin.

---

## 4. Disputed / needs-human-review

#### Foreign steady-state contamination in the ss engine — R/engines.R:279–289 (claimed medium)
- **Confirmed view:** The `psi` full conditional uses the full M×M `Sigma_inv`, so the foreign sub-vector of `psi` depends on domestic residuals. Isolation test (A and Sigma fixed at posterior means, block-exog verified at ~7e-6): shifting domestic GDP by 1 sd moves foreign `psi` by up to ~1 posterior sd; a counterfactual foreign-only update (using `Sigma_ff/U_ff`) is invariant. README D3/D5 frame block exogeneity and the inherited foreign steady state without noting this contamination, and the gated `block_exog_metric` checks only A. → Genuine defect, medium.
- **Refuted view:** This *is* the canonical, correct Villani (2009) steady-state full conditional. Block exogeneity restricts the AR coefficient matrix (Granger non-causality), not the error covariance; with a non-block-diagonal Sigma (correlations up to 0.83 here), GLS legitimately uses cross-equation residual information for the foreign means — that is the true posterior, not a leak. The finder's suggested foreign-only fix would make the Gibbs step an *incorrect* full conditional. At most a documentation clarification.
- **Adjudication needed:** Is exact block exogeneity *of the means* a required contract (README intent) or only of the lag dynamics? If the former, the foreign-only update or a documentation change is warranted; if the latter, the code is correct as-is.
- **RESOLVED (2026-06-22): follow Villani.** Block exogeneity is a contract on the *lag dynamics only*, not the steady-state means. The full-Σ ψ full conditional is the correct Villani (2009) posterior and is kept unchanged; a foreign-only update would be an incorrect full conditional. Documented in README D3 ("Scope — dynamics, not the steady-state means") and a code note in `fit_ss` (R/engines.R). No code change to the sampler.

#### RW weighted sd: `sd(diff(y)*w)` — R/benchmarks.R:10–12 (claimed low)
- **Confirmed view:** Not a proper weighted sd — `sd()` re-estimates the mean of scaled increments and uses n-1; `fit_ucmean` (R/benchmarks.R:235–237) implements the textbook weighted form correctly, proving it was available. Crude but inert in outcome.
- **Refuted view:** With `w = 1/s_t`, `diff(y)*w` is the *whitened* (homoskedastic) increment, for which an ordinary sd with n-1 is the correct GLS-style estimator — mirroring the engines' row-weighting. Monte Carlo shows the code estimator is unbiased and marginally *better* than the proposed "proper weighted sd." The downward-bias rationale is methodologically backwards.
- **Adjudication needed:** Both agree the output is approximately correct. This is a statistical-style question, not a forecast-affecting defect.

#### Benchmarks skip two-pass COVID-scale refinement — R/evaluate.R:51–58, 78–80 (claimed low)
- **Confirmed view:** VAR members re-estimate COVID scales with weighted `ar_sigmas`; benchmarks do a single unweighted-sigma pass. Demonstrated to differ at 1 of 24 COVID-active origins (scale 4→8, s_future differing by ~0.11), contradicting README D17.1's intent that benchmarks share the same outlier-robust calibration.
- **Refuted view:** The coarse scale/decay grids make the marginal-likelihood argmax robust; recomputing both passes on the cached panel gave *zero* difference in scales/rho/s_future at every COVID origin. Skipping a redundant pass for univariate targets is defensible.
- **Adjudication needed:** Whether the (at most one-origin, modest) asymmetry matters for benchmark comparability is a judgement call; cheap to resolve by adding the refinement pass for symmetry.

#### Two minor conj_br/priors disputes (nit-level)
- **Domestic SOC keep-filter zero-ybar edge case** (R/engines.R:144–151): both verifiers agree the all-zero retained row is inert (contributes nothing to the cross-products); split only on whether to flag it. Trigger requires a p-observation mean to hit exactly 0.0 — unreachable on real data. Cosmetic.
- **Domestic-block omega hardcodes lag_decay=1** (R/engines.R:130): code observation is accurate (`idx$lag` vs `idx$lag^lag_decay` elsewhere), but `lag_decay` is never plumbed into conj_br at all (neither foreign nor domestic block), so the claimed *inconsistency between blocks* cannot occur and nothing miscomputes today. Latent fragility / clarity nit only.

---

## 5. Per-subsystem assessment

- **Engine: conj_br (block-recursive conjugate NIW)** — ✅ sound. Block exogeneity exact by construction (foreign VAR carries no domestic regressors); verified forecast-level exogeneity (perturbing domestic history left foreign draws bit-identical), posterior-mean = B1 over 20k draws, SOC/DIO splitting correct. Only nit-level fragility (no nd≥1 guard).
- **Engine: gibbs (independent Normal-IW)** — ✅ sound. vec(B)|Sigma reproduces OLS under diffuse prior; Sigma|B uses correct IW; block exogeneity via asymmetric prior (leakage ~7e-6). The ill-conditioning concern was refuted: Cholesky is backward-stable for the penalty-driven structure (solve residual ~1e-9 even far tighter than config). Minor: priors.R default 1e-5 vs config 1e-4 (only tests use the default).
- **Engine: ss (Villani steady-state)** — ⚠️ issues. Core math correct (reparametrisation, three full conditionals, unit-test collapse, forecasts converge to psi). One disputed medium: foreign `psi` informed by domestic data via full Sigma. The near-singular-U concern was refuted (the ss prior provably keeps Pp positive definite).
- **Engine: sv (stochastic volatility)** — ⚠️ issues. Triangular CCM factorisation, block exogeneity, lag/contemporaneous reconstruction, and AR(1) log-vol propagation all correct. One confirmed **high**-severity defect: unstandardised t predictive draws.
- **Priors (Minnesota, SOC/DIO, GLP)** — ✅ sound. Minnesota moments, conjugate implied SDs, and GLP log-ML all verified to machine precision against independent computations; lambda selection strictly up-to-origin. Only latent/nit items (lag_decay hardcode, minimal-nu0 convention refuted as standard).
- **Predictive simulation** — ✅ sound (one low item). Parameter+shock uncertainty propagation, companion recursion, lag stacking, conj_br block append, and reproducibility all verified; h=1 means reproduce posterior-mean predictions. Only the wasted last-third-of-posterior draws (no bias).
- **COVID corrections (LP scaling + t-errors + dummy)** — ✅ sound. Quarter-matching via `.qidx` robust to the start vs mid-quarter stamp; recursive scale build-up (2→2/8→2/8/8) with no pre-2020 leakage; engine switch correct (SV/UCSV ignore weights/shock_scale). All findings here refuted/disputed as documented or benign.
- **Data acquisition + panel construction** — ⚠️ issues. Fail-loud ingestion, mtime-keyed cache, consistent quarter stamping, genuine synthetic ground truth all verified. Real issues: fragile `max(n)` coverage heuristic (medium), config start off-by-one (low, latent), stale FRED-fallback docstring (low).
- **Transforms + spec + level reconstruction** — ✅ sound. dlog/level/loglevel transforms exact; cum = `100*(log lvl_{t+h} − log lvl_t)` verified to machine precision; ye index arithmetic correct and no-look-ahead. Only dead-code nit (`year_ended()`).
- **Forecast averaging (combination)** — ✅ sound. Weights non-negative and sum to exactly 1 across every cell; equal/logscore/pool correct and non-degenerate; strict no-look-ahead verified. Confirmed issues confined to BMA (diagnostic-only).
- **OOS harness / no-look-ahead / scoring / DM** — ⚠️ issues. No-look-ahead structurally enforced and empirically verified; reproducibility deterministic; q/ye/cum alignment and CRPS/log/PIT correct; DM has correct Bartlett/NW + Harvey correction. One confirmed **high** reporting defect (ex-COVID no-op + -Inf contamination) plus two low items (DM contiguity, thin self-tests).
- **Benchmarks (rw, ar4, ucsv, ucmean)** — ⚠️ issues. RW linear-variance growth, AR(4) NIG posterior + stationarity deflation, and ucmean GLS-mean all verified. One confirmed **high** defect (UCSV t-variance, shared with SV); rw-sd and two-pass items disputed.

---

## 6. What was verified clean

The following load-bearing properties were checked and held, scoping the assurance:

- **Block exogeneity (lag dynamics):** exact for conj_br (forecast draws bit-identical under domestic perturbation) and sv; ~6e-6–7e-6 leakage for gibbs/ss. Caveat: ss steady-state *means* are contaminated (disputed).
- **No-look-ahead:** structurally enforced at a single slice chokepoint (`td[first:t]`); COVID scales, AR sigmas, and GLP lambda all consume only sliced data; corrupting post-origin data through the entry point is caught by `test_no_lookahead`. Combination training strictly `origin+h <= t`.
- **Weight normalisation:** all combo and final-forecast weights ≥0 and sum to exactly 1 across every scheme×variable×bucket×origin; equal=1/n; logscore/pool non-degenerate and pool beats equal at tested cells.
- **Posterior correctness:** gibbs vec(B)|Sigma reproduces OLS under diffuse prior; conjugate posterior mean = B1; GLP log-ML and Minnesota moments match independent references to machine precision.
- **Transform round-trips:** dlog/level/loglevel exact; cum telescopes to the log-level difference; ye=cum=Σq at h=4 to machine precision.
- **Quarter alignment:** `.qidx` maps both start-stamp (real) and mid-quarter (synthetic) 2020 dates to identical indices — no off-by-one in COVID treatment (the reporting layer is the exception, see the high finding).
- **Reproducibility:** `derive_seed` deterministic, collision-free, scheduling-independent; identical draws on re-run for gibbs and sv.
- **DM construction:** correct Bartlett/Newey-West truncation (lag h−1, h+3 for overlapping ye) with Harvey small-sample correction and consistent loss orientation.
- **config_hash:** covers the nine estimation-relevant source files so the OOS draw cache cannot serve stale results.

---

## 7. Prioritised recommendations

1. **Fix the SV/UCSV t-variance bug (high, production-affecting).** In `simulate_paths.post_sv` (R/forecast.R:146) and `simulate_paths.post_ucsv` (R/benchmarks.R:224), standardise the t-innovation: `eps <- if (is.finite(nu) && nu > 2) rt(1, df=nu)*sqrt((nu-2)/nu) else rnorm(1)`. This is the single change that materially corrects reported production density forecasts and combination weights for the SV members and the UCSV benchmark. Re-run OOS afterward (it busts the cache hash by design).
2. **Fix the ex-COVID exclusion + -Inf contamination (high, reporting).** In `summarise_scores` (R/evaluate.R:224) match at quarter granularity via `.qidx` instead of exact dates, and guard log-score aggregation against non-finite values (NA-out -Inf, or report median/trimmed-mean log score). This makes the headline and ex-COVID scorecards usable.
3. **Harden the panel coverage heuristic (medium).** Replace `expected <- max(n)` with `median`/mode and de-duplicate observations per quarter before counting (R/data_sources.R:171), to remove the silent series-collapse risk on duplicated provider rows.
4. **Adjudicate the ss steady-state contamination (disputed medium).** Decide whether block exogeneity of the *means* is a required contract. If yes, update the foreign `psi` sub-vector using `Sigma_ff/U_ff` (and add a test); if no, document in README D3/D5 that ss block exogeneity applies only to lag dynamics.
5. **Low-effort correctness/clarity fixes:** normalise config boundary dates to start-of-quarter (R/data_sources.R:252); align `forecast_draws` with `ndraw` or subsample posterior indices (R/forecast.R); fix the BMA "one-step" label/horizon and zero-weight floor (R/combine.R:86, 140) — even though diagnostic-only; correct the stale FRED-fallback docstring (R/data_sources.R:178).
6. **Robustness hardening (low/nit):** add `stopifnot(nd >= 1, nf >= 1)` to `fit_conj_br`; assert contiguous origins (or index NW by origin spacing) in DM; extend `test_no_lookahead`/`test_reproducibility` to all engines + a COVID-spanning origin; plumb `lag_decay` into conj_br or document it as fixed at 1; remove or wire up the dead `year_ended()` helper.
7. **Consider** adding the weighted-sigma COVID refinement pass for benchmarks (or documenting its omission) for symmetry with VAR members, and aligning the `minnesota_prior` default `block_exog_sd` (1e-5) with config (1e-4) — both low priority, no demonstrated production impact.

---

## 8. Resolution log (2026-06-22)

All confirmed correctness findings and the low/nit cleanups were implemented; the
**disputed** items were deliberately left for adjudication (one verifier judged
some of those "fixes" would make the code *incorrect*). Verified afterward: all
unit tests pass; the standardised t-draw has unit variance (MC: 0.999 vs the old
ν/(ν−2)); the ex-COVID summary now drops the 62 −∞ cells → 0 (the full table
keeps them by design as pooling evidence).

| Finding | Status | Action |
|---|---|---|
| SV t-variance (forecast.R) | ✅ fixed | standardise: `rt(ν)*sqrt((ν−2)/ν)` for finite ν>2 |
| UCSV t-variance (benchmarks.R) | ✅ fixed | same standardisation on the transitory shock |
| ex-COVID no-op + −∞ (evaluate.R) | ✅ fixed | exclude at quarter granularity via `.qidx`; full-table −∞ kept intentionally |
| Panel coverage `max(n)` (data_sources.R) | ✅ fixed | `expected <- median(n)` (robust to duplicate rows; passes all `to_quarterly` tests) |
| Config start off-by-one (data_sources.R) | ✅ fixed | floor config bounds to quarter start before filtering |
| Forecast draws use only 1st 1000/1500 (forecast.R/benchmarks.R) | ✅ fixed | `k <- floor((d−1)*post$ndraw/ndraw)+1` — spread across the full posterior |
| BMA "one-step" mislabel (combine.R) | ✅ fixed | comment corrected (shortest-horizon-in-bucket; diagnostic-only) |
| BMA zero-weight floor (combine.R) | ✅ fixed | pool over positive-weight members only |
| Stale FRED-fallback docstring (data_sources.R) | ✅ fixed | corrected to "no silent substitution" |
| conj_br nd≥1/nf≥1 (engines.R) | ✅ fixed | added `stopifnot` |
| DM contiguity (evaluate.R) | ✅ fixed | warn once if origins non-contiguous |
| Thin self-tests (evaluate.R) | ✅ fixed | no-look-ahead + reproducibility now loop one member per engine; NLA adds a COVID-spanning origin |
| lag_decay hardcode (engines.R) | ✅ documented | comment: decay exponent fixed at 1 in conj_br |
| Dead `year_ended()` (transforms.R) | ✅ documented | kept (tested reference); comment notes the inline scorer must match it |
| ss foreign-ψ contamination | ✅ resolved (follow Villani) | code kept (correct full conditional); documented that block-exog is a lag-dynamics contract only (README D3 + fit_ss note) |
| rw weighted sd | ⚠️ deferred (disputed) | verifier MC showed the current estimator is unbiased/fine |
| benchmark two-pass COVID refinement | ⚠️ deferred (disputed) | verifier found zero difference at every COVID origin |
| conj_br SOC zero-ybar edge; minnesota default 1e-5 | ⚠️ deferred | cosmetic/unreachable; default change could perturb tests, no production impact |

After the code fixes the full OOS pipeline was re-run (the edits bust the cache
hash by design) and the scorecard regenerated.