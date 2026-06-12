# DECISIONS.md — modelling choices and rationale

Running log, in the spirit of "decide, document, proceed". Each entry: the
choice, the rationale, the rejected alternative, and the config key that
controls it. Dated 2026-06-12 unless noted.

---

## D1. Variable roster and blocks

**Choice.** Foreign block: world activity (`f_act`, US industrial production
proxy), commodity prices (`f_comm`, RBA Index of Commodity Prices), foreign
policy rate (`f_rate`, fed funds), plus an alternative trade-weighted activity
measure (`f_act_tw`, G20 GDP) reserved for an alt-foreign variant. Domestic
small core: GDP growth, trimmed-mean CPI inflation, unemployment rate, cash
rate, real TWI. Medium set adds terms of trade, WPI growth, employment growth,
consumption growth, 10-year AGS yield (13 variables total).

**Why.** The small core covers the four policy targets plus the exchange rate
(the SOE adjustment channel). The medium set follows Bańbura–Giannone–Reichlin
(2010): medium systems often forecast better with enough shrinkage. The US
proxy is the standard caveat for Australia (Asia matters more now); the
trade-weighted alternative is wired into the spec (`sets: [alt_foreign]`) as a
documented suite-diversity hook.

**Rejected.** A large (20+) system — recursive re-estimation cost grows
quadratically and the medium set already exercises the shrinkage machinery.
**Config:** `variables`, each member's `set`.

## D2. Transforms and Minnesota deltas

**Choice.** Growth-type series enter as quarterly 100·Δlog with `delta = 0`;
interest/unemployment rates as levels with `delta = 1`; real TWI and ToT as
100·log levels with `delta = 1`. Trimmed-mean CPI (published as a quarterly %
change) is cumulated to an index in the data layer so all series arrive as raw
levels and one transform layer serves both data paths.

**Why.** `delta = 0` for mean-reverting growth (white-noise prior centre) and
`delta = 1` for persistent levels is the textbook Minnesota assignment
(Litterman 1986); modelling rates in levels keeps the cash rate and the policy
question in interpretable units, with the I(1)-discipline handled by priors
(D5) rather than differencing away level information.

**Rejected.** Differencing everything (loses cointegration-ish level
relations); modelling log CPI level (long-memory inflation drifts wreck a
fixed-mean VAR over a 35-year window). **Config:** `variables.*.transform`,
`variables.*.delta`.

## D3. Block exogeneity as a first-class property

**Choice.** Every VAR member imposes the foreign/domestic partition
(Cushman–Zha 1997; Zha 1999): foreign equations exclude domestic information.
Three mechanisms, by engine: (i) `gibbs`/`ss`: near-zero prior mean **and**
variance (`mcmc.block_exog_prior_sd: 1e-4`) on domestic-lag coefficients in
foreign equations — possible only because these engines use an independent
(non-Kronecker) prior; (ii) `conj_br`: exact, by estimating the foreign block
as its own VAR (block-recursive, the Bloor–Matheson RBNZ approach) with the
domestic block conditioned on contemporaneous foreign values; (iii) `sv`:
exact, by dropping domestic columns from foreign equations in the
triangularised system. The max |posterior-mean domestic coefficient in a
foreign equation| is logged per fit, checked against synthetic ground truth,
and gated in `assert_diagnostics()`.

**Why.** Without it the model implies Australia moves world demand and
commodity prices, which contaminates medium-horizon dynamics; RBA RDP 2013-06
shows the restriction materially raises the role of foreign shocks.

**Rejected.** A pure Kronecker-conjugate suite (cannot represent asymmetric
prior variances — this is exactly why the spec requires the Gibbs engine).
**Config:** `variables.*.block`, `mcmc.block_exog_prior_sd`.

## D4. Tightness: GLP-style marginal-likelihood selection

**Choice.** Overall tightness λ is selected by maximising the closed-form
conjugate NIW marginal likelihood over a grid (`glp.lambda_grid`), re-run *at
every forecast origin inside the recursive loop*, and shared by the members
whose `prior.lambda` is `auto`. Two members keep fixed λ (0.05 tight, 0.4
loose) as suite diversity.

**Why.** Giannone–Lenza–Primiceri (2015): treat shrinkage hyperparameters as
parameters, not folklore. The conjugate marginal likelihood is cheap and
closed-form; a full hierarchical MCMC over λ inside every origin × member
would multiply compute for little forecast gain. Re-selecting per origin keeps
the no-look-ahead property exact.

**Rejected.** Hand-set λ everywhere (indefensible); full hierarchical NUTS/MH
over hyperparameters (cost). **Config:** `glp.*`, member `prior.lambda`.

## D5. Long-run priors: SOC + DIO, and a steady-state member

**Choice.** Sum-of-coefficients (Doan–Litterman–Sims) and
dummy-initial-observation (Sims–Zha) artificial observations on the
`delta = 1` level variables in selected members
(`prior.soc/dio`); and a Villani (2009) steady-state member with informative
priors on the unconditional means: trimmed-mean inflation 2.5 % p.a. (target
midpoint, sd 0.6pp), GDP growth 2.8 % p.a., NAIRU 4.5 %, neutral cash rate
3.5 %, foreign anchors likewise; loose sample-mean anchors for the TWI/ToT.

**Why.** Iterated forecasts at h = 8–12 converge to the model's unconditional
mean, which plain Minnesota pins down poorly; the steady-state
reparametrisation lets policy-relevant long-run information (and the foreign
block's steady states, which the domestic forecast inherits) enter directly.

**Rejected.** SOC/DIO on growth-rate variables (their sample means are already
well-estimated; the dummies would just bias short-run dynamics).
**Config:** member `prior.soc/dio/soc_mu/dio_delta`, `variables.*.ss_mean/ss_sd`.

## D6. Stochastic volatility via exact triangular equation-by-equation

**Choice.** The SV engine factorises the reduced form as a recursive system
(foreign block ordered first), estimating each conditional equation
independently: Bayesian regression with SV errors, β drawn by weighted least
squares given the log-variance path, the path updated with
`stochvol::svsample_fast_cpp`. This is the Carriero–Clark–Marcellino (2019)
equation-by-equation idea; with independent per-equation priors the posterior
factorises *exactly*, so per-equation estimation is not an approximation.

**Why.** Keeps the recursive OOS loop feasible (the medium SV member estimates
in seconds) while delivering the density-calibration benefits of SV
(Clark 2011) and protecting against the 2020 outliers.

**Rejected.** Joint VAR-SV with a full covariance Gibbs (cost explodes in the
recursive loop); Lenza–Primiceri outlier dummies (SV with the COVID quarter in
the synthetic DGP already absorbs it; noted as an extension).
**Config:** members with `engine: sv`.

## D7. Iterated multi-step forecasts

**Choice.** All members iterate the one-step model forward, simulating shocks
and parameter draws jointly (h = 1…12); no direct/local-projection members.

**Why.** Marcellino–Stock–Watson (2006): iterated wins at long horizons when
the one-step model is adequate, and only iterated forecasts give internally
coherent joint densities for year-ended transforms and fan charts.

**Rejected.** A direct-forecast robustness member — 12 extra estimations per
member-origin for a known-weaker density forecaster; revisit if point RMSE at
h ≥ 8 looks broken. **Config:** none (structural).

## D8. The suite roster

Seven VARs + four univariate anchors (see `suite` in config):
`small_minn` (Gibbs, GLP-λ, SOC+DIO), `small_ss` (steady-state),
`small_sv` (SV), `small_loose_p5` (λ=0.4, p=5), `medium_minn` (13-var SV,
p=2), `medium_conj` (13-var block-recursive conjugate, SOC+DIO),
`small_tight` (λ=0.05 conjugate), plus `rw`, `ar4`, `ucsv`, `ucmean`.
Diversity axes per the brief: prior family, size, tightness, lag length,
volatility model. The members genuinely fail differently: the tight small
model is robust at long horizons, the medium SV member wins short-horizon
density, the steady-state member anchors h ≥ 8, the benchmarks discipline
everything.

**Rejected.** An explicit US-proxy vs trade-weighted pair of otherwise
identical members in the default roster (the alt foreign series is a DBnomics
pull that fails without network; the variant is wired in config via the
`alt_foreign` set but not in the default suite, so the offline run stays
green). **Config:** `suite`, `benchmarks`.

## D9. UCSV applied per target variable, outlier-robust specification

**Choice.** A UC-SV model (local level + FFBS) is run for *each* target
variable, not just inflation, with three deliberate deviations from the
canonical twin-SV Stock–Watson form, all forced by COVID-sized outliers in the
evaluation window:
(i) **t-distributed transitory errors** (`nu ~ Exp(0.1)` via stochvol) — the
Gaussian version goes sticky (trend ESS ≈ 2) at origins where the outlier is
the sample endpoint;
(ii) **constant trend-shock variance with a conjugate IG update** instead of a
second SV process — the trend-vol dimension is where the endpoint
trend-vs-noise bimodality piles up, and a direct Gibbs draw removes it;
(iii) **tight vol-of-vol prior** (`sigma2 ~ Gamma(2, 50)`, mean 0.04) — the
Stock–Watson tradition *fixes* the vol smoothing parameter near 0.2; a loose
prior lets an endpoint outlier inflate vol-of-vol and the 12-step-ahead
variance forecast explodes through exp(h/2) compounding (observed: 99.5%
unemployment quantiles in the hundreds). Updated 2026-06-12 after the first
full synthetic run failed the convergence gate at origins 117–120.

**Why per target.** It is the canonical inflation benchmark, and as a generic
trend-plus-noise-with-SV model it is a sensible univariate density anchor for
the other targets too. **Config:** `benchmarks`, `mcmc.bench_*`.

**Convergence gating.** The trend/noise *split* of a near-white series is
weakly identified in any UC model (the classic variance pile-up): chains for
the decomposition parameters (trend variance, SV level) and for the trend
endpoint mix slowly even in very long runs, while the predictive *sum* —
the only thing the forecasting suite consumes — is well identified. UCSV
convergence is therefore gated on a Monte-Carlo-precision criterion for the
forecast-relevant quantity: the MCSE of the trend endpoint must be below 15 %
of the one-step predictive sd, with an adaptive triple-thinning retry before
failure. Raw ESS values are still computed and reported (trend-endpoint ESS
appears as `ess_min` in the diagnostics table; it can be low at origins where
an outlier sits at the boundary, which is expected and documented rather than
hidden).

## D10. Combination: linear pools, per variable × horizon bucket, shrunk

**Choice.** Density-level linear pools with weights per target variable and
per horizon bucket (1–4, 5–8, 9–12): equal; recursive log-score weights
(Jore–Mitchell–Vahey 2010) with forgetting 0.95; optimal pool
(Hall–Mitchell / Geweke–Amisano) maximising the discounted historical pool log
score, softmax-parametrised BFGS; BMA via one-step predictive likelihood,
**reported as a diagnostic only**. Non-BMA weights are shrunk toward equal
with κ = 0.3; weights at origin t use only scores with s + h ≤ t; equal
weights until 8 training origins exist.

**Why.** Bucketing by horizon recognises that the best member changes with
horizon (the design premise of the suite); shrinkage + forgetting are the
standard guards against weight overfitting (the combination puzzle — equal
weights are the benchmark to beat, so they are also a scheme). Geweke–Amisano:
the optimal pool does not degenerate to one model; BMA does, which is why it
is reported but not recommended. Pool log scores use log-sum-exp; pool CRPS
uses draws resampled from the mixture; pool PIT is the weighted member PIT
(exact mixture CDF).

**Rejected as default.** Del Negro–Hasegawa–Schorfheide dynamic pools — the
forgetting-factor log-score scheme already provides time variation; a full
dynamic pool is the natural next extension. **Config:** `combination.*`.

## D11. Evaluation design

**Choice.** Expanding window (rolling available via `evaluation.window`),
origins = the last `max_origins` (36) quarters ending at T−1, final-vintage
data with the caveat stated loudly in the report. Scoring: RMSE, log score
(kernel density from draws), CRPS, PIT, by horizon and variable, quarterly and
year-ended. DM tests with Harvey correction and Newey–West (h−1) variance vs
the AR(4) and RW anchors.

**Why.** Expanding windows match how a central bank actually re-estimates;
36 origins balances test power against compute. Real-time vintages for
Australia are not reliably available key-free; pretending otherwise would be
worse than documenting the caveat. DM rather than Giacomini–White: with an
expanding window and these sample sizes the conditional GW adds machinery
without changing conclusions; noted as an extension.
**Config:** `evaluation.*`.

## D12. Synthetic DGP

**Choice.** A stationary block-exogenous VAR(2) in the transformed units with
a structured (not random) coefficient pattern — foreign→domestic activity
loadings, a Taylor-rule-ish cash-rate equation, monetary transmission with the
right signs — Gaussian shocks whose contemporaneous loadings run strictly
foreign→domestic, a COVID-style outlier quarter with partial rebound, and
integration back to raw levels so the transform layer runs identically to the
real path.

**Why.** Random off-diagonal coefficients compound badly under high
persistence (rate series wandered to ±10 %); a structured DGP keeps every
series in realistic ranges, gives the block-exogeneity diagnostic a known
ground truth, and exercises the same code paths as real data.
**Config:** `synthetic.*`, `data.source`.

## D13. Seeds and parallelism

**Choice.** One master seed; every stochastic task (member × origin, final
forecasts, mixture resampling) derives its own seed by hashing the master seed
with a task key, so results are identical regardless of worker scheduling.
Parallelism is over origins inside each member target (`furrr`), because the
recursive loop, not any single fit, is the bottleneck.

**Rejected.** future's built-in L'Ecuyer streams alone (reproducible only for
a fixed worker count). **Config:** `master_seed`, `parallel.workers`.

## D14. Caching

**Choice.** Raw downloads cached in `data/raw/` with a manifest and freshness
window; per-(member, origin) OOS results cached in `cache/oos_<confighash>/`
keyed by a hash of the estimation-relevant config sections, so combination
experiments re-use estimation; `targets` provides DAG-level caching on top.

**Why.** The recursive loop is the expensive stage; weight-scheme tweaks should
cost seconds, not an hour. **Config:** `data.cache_max_age_days`; hash covers
`master_seed, data, synthetic, variables, horizons, mcmc, glp, suite,
benchmarks, evaluation`.

## D15. Conditional forecasts

**Choice.** Hard conditioning by substitution (overwrite the conditioned
variable's path each step as it feeds the recursion), exposed as an optional
argument, documented as an approximation.

**Rejected.** Waggoner–Zha conditional simulation — correct but heavy, and
conditional forecasting is an optional capability in this spec, not a scored
deliverable.

## D16. Real-data details

**Choice.** RBA series via `readrba`, ABS via `readabs`, US/world via DBnomics
(`rdbnomics`, key-free) with a FRED fallback only if `FRED_API_KEY` is set.
Mixed-frequency series are quarterly-averaged. Failures stop with a message
pointing at the synthetic path (always-runnable guarantee lives there, not in
fragile retries).

**Caveat recorded.** Trimmed-mean CPI is published as a quarterly rate and is
cumulated to a synthetic index purely so the uniform dlog transform
reproduces it; the round trip is exact up to float error. The same `pre:
pct_change` mechanism handles any source published as a % change.

**Freshness caveat (updated after the real-data run).** The key-free foreign
activity series (IMF IFS US industrial production via DBnomics) lags by about
a year, and the balanced panel is trimmed to the stalest series — the 2026-06
real run therefore ends 2024 Q4. Setting `FRED_API_KEY` switches `f_act` to
FRED INDPRO and recovers the missing quarters. The long-history RBA commodity
index is `GRCPAISDR` (from 1982); the bulk-spot variant `GRCPAISAD` only
starts in 2009 and silently truncated the panel — found and fixed in the
first real-data run.
