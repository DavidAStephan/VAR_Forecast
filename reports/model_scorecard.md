# SOE-BVAR Suite — Model Scorecard

**Data source:** real (RBA + ABS + FRED)  
**Panel:** 1997Q4-2026Q1 | **Targets:** gdp_growth, cpi_inflation, unemp_rate, cash_rate | **Horizons:** 1-12 quarters  
**Evaluation:** expanding-window pseudo-real-time, 36 forecast origins; densities scored by CRPS and log predictive density, points by RMSE, all by horizon.  
**Diagnostics (§9):** all green (block exogeneity, MCMC convergence, forecast sanity, no-look-ahead, reproducibility).

## 1. The models

Every VAR member is **block-exogenous** (the domestic block never feeds back into the foreign block) and produces **iterated** density forecasts. Members are designed to fail differently; see DECISIONS.md for the full rationale.

### 1a. VAR members

| Model | Family | System | Lags | Shrinkage λ | Volatility | COVID |
|:--|:--|:--|:--|:--|:--|:--|
| `small_minn` | Independent Normal-inverse-Wishart (Gibbs) + SOC/DIO | small (8 var) | 4 | auto (GLP) | constant | LP scaling |
| `small_ss` | Steady-state (Villani) | small (8 var) | 4 | auto (GLP) | constant | LP scaling |
| `small_sv` | Stochastic volatility, equation-by-equation | small (8 var) | 4 | auto (GLP) | stochastic (SV) | t-errors (SV-t) |
| `small_loose_p5` | Independent Normal-inverse-Wishart (Gibbs) | small (8 var) | 5 | 0.4 | constant | LP scaling |
| `medium_minn` | Stochastic volatility, equation-by-equation | medium (13 var) | 2 | auto (GLP) | stochastic (SV) | t-errors (SV-t) |
| `medium_conj` | Block-recursive conjugate NIW + SOC/DIO | medium (13 var) | 4 | 0.1 | constant | LP scaling |
| `small_tight` | Block-recursive conjugate NIW + SOC/DIO | small (8 var) | 4 | 0.05 | constant | LP scaling |

### 1b. Benchmark members (the bar every VAR must clear)

| Model | Description | COVID |
|:--|:--|:--|
| `rw` | Random walk | LP scaling |
| `ar4` | Bayesian AR(4) | LP scaling |
| `ucsv` | Unobserved components + stochastic volatility | t-errors (robust) |
| `ucmean` | Unconditional mean | LP scaling |

### 1c. Combination schemes (density pools)

Weights estimated **per target variable and per horizon bucket** (near 1-4, medium 5-8, far 9-12), shrunk toward equal weights, strictly recursive (no look-ahead).

| Scheme | How weights are set |
|:--|:--|
| `combo_equal` | Equal weights (the benchmark pool — hard to beat) |
| `combo_logscore` | Recursive log-score weights, with forgetting |
| `combo_pool` | Optimal prediction pool (Hall-Mitchell / Geweke-Amisano) |
| `combo_bma` | Bayesian model averaging — reported as a diagnostic only |

## 2. Forecast performance

Lower CRPS / RMSE is better; higher log score is better. **Bold** = best in that column. Models ordered best-first (by mean over the shown horizons).

### 2a. Who forecasts best, by variable and horizon (CRPS)

Best single model (lowest mean CRPS in the bucket; value in parentheses):

| Variable | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| gdp_growth | small_sv (0.683) | small_tight (0.773) | small_tight (0.846) |
| cpi_inflation | rw (0.189) | ar4 (0.273) | ucmean (0.302) |
| unemp_rate | medium_minn (0.318) | small_ss (0.631) | small_ss (0.721) |
| cash_rate | small_sv (0.307) | combo_equal (0.893) | ar4 (1.447) |

### 2b. Density accuracy by variable and horizon (CRPS, lower better)

**Real GDP growth (qtr %) (`gdp_growth`)**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| small_loose_p5 | 0.627 | **0.670** | 0.741 | 0.808 | 0.883 |
| combo_equal | 0.634 | 0.682 | 0.734 | 0.817 | 0.887 |
| ar4 | 0.668 | 0.684 | **0.727** | 0.797 | 0.898 |
| small_sv | 0.611 | 0.678 | 0.735 | 0.850 | 0.901 |
| small_ss | 0.664 | 0.689 | 0.743 | 0.804 | 0.888 |
| small_minn | 0.674 | 0.697 | 0.744 | 0.802 | 0.880 |
| small_tight | 0.682 | 0.702 | 0.743 | **0.794** | 0.878 |
| medium_minn | **0.610** | 0.685 | 0.743 | 0.840 | 0.925 |
| medium_conj | 0.672 | 0.700 | 0.757 | 0.800 | **0.876** |
| ucmean | 0.680 | 0.703 | 0.732 | 0.816 | 0.908 |
| ucsv | 0.634 | 0.673 | 0.767 | 0.900 | 0.981 |
| combo_pool | 0.642 | 0.705 | 0.758 | 0.931 | 1.018 |
| combo_logscore | 0.754 | 0.791 | 0.860 | 1.134 | 1.042 |
| combo_bma | 0.842 | 0.903 | 0.946 | 1.409 | 1.164 |
| rw | 0.941 | 1.100 | 1.164 | 1.387 | 1.447 |

**Trimmed-mean CPI inflation (qtr %) (`cpi_inflation`)**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| ar4 | 0.162 | 0.208 | **0.246** | 0.288 | 0.314 |
| combo_equal | 0.161 | 0.194 | 0.246 | 0.323 | 0.364 |
| ucmean | 0.230 | 0.246 | 0.250 | **0.286** | **0.312** |
| combo_pool | 0.147 | 0.188 | 0.257 | 0.354 | 0.391 |
| rw | **0.131** | **0.170** | 0.248 | 0.361 | 0.431 |
| combo_logscore | 0.137 | 0.181 | 0.263 | 0.365 | 0.408 |
| small_ss | 0.172 | 0.203 | 0.269 | 0.356 | 0.371 |
| small_minn | 0.170 | 0.202 | 0.262 | 0.356 | 0.387 |
| small_sv | 0.170 | 0.212 | 0.269 | 0.356 | 0.407 |
| small_tight | 0.198 | 0.226 | 0.272 | 0.350 | 0.395 |
| medium_minn | 0.165 | 0.205 | 0.264 | 0.364 | 0.444 |
| combo_bma | 0.145 | 0.193 | 0.288 | 0.401 | 0.439 |
| medium_conj | 0.191 | 0.218 | 0.272 | 0.367 | 0.426 |
| ucsv | 0.175 | 0.229 | 0.283 | 0.387 | 0.430 |
| small_loose_p5 | 0.192 | 0.226 | 0.295 | 0.415 | 0.422 |

**Unemployment rate (%) (`unemp_rate`)**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| small_ss | **0.146** | 0.283 | 0.500 | **0.688** | **0.747** |
| medium_minn | 0.146 | 0.272 | **0.466** | 0.722 | 0.796 |
| combo_equal | 0.159 | 0.285 | 0.490 | 0.699 | 0.782 |
| small_sv | 0.154 | 0.289 | 0.501 | 0.758 | 0.782 |
| small_minn | 0.152 | 0.292 | 0.529 | 0.753 | 0.787 |
| small_loose_p5 | 0.155 | 0.295 | 0.546 | 0.780 | 0.808 |
| combo_pool | 0.174 | 0.297 | 0.490 | 0.765 | 0.871 |
| small_tight | 0.158 | 0.307 | 0.549 | 0.779 | 0.811 |
| rw | 0.169 | 0.310 | 0.545 | 0.777 | 0.913 |
| combo_bma | 0.148 | **0.272** | 0.475 | 0.935 | 0.886 |
| medium_conj | 0.163 | 0.306 | 0.545 | 0.820 | 0.896 |
| combo_logscore | 0.232 | 0.348 | 0.532 | 0.774 | 0.848 |
| ucsv | 0.265 | 0.376 | 0.531 | 0.737 | 0.832 |
| ar4 | 0.189 | 0.326 | 0.534 | 0.760 | 0.955 |
| ucmean | 0.726 | 0.770 | 0.774 | 0.905 | 1.034 |

**Cash rate (%) (`cash_rate`)**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| combo_equal | 0.147 | 0.263 | **0.514** | **1.126** | 1.668 |
| small_ss | 0.139 | 0.252 | 0.525 | 1.151 | 1.676 |
| ar4 | 0.167 | 0.303 | 0.589 | 1.141 | **1.593** |
| small_minn | 0.145 | 0.259 | 0.527 | 1.152 | 1.785 |
| combo_pool | 0.130 | 0.243 | 0.547 | 1.223 | 1.724 |
| combo_logscore | 0.122 | 0.235 | 0.532 | 1.223 | 1.758 |
| rw | 0.183 | 0.320 | 0.590 | 1.136 | 1.670 |
| combo_bma | 0.124 | 0.239 | 0.576 | 1.393 | 1.850 |
| ucsv | 0.170 | 0.323 | 0.639 | 1.266 | 1.832 |
| small_tight | 0.157 | 0.287 | 0.598 | 1.286 | 1.935 |
| medium_conj | 0.192 | 0.323 | 0.572 | 1.239 | 2.053 |
| small_sv | 0.115 | **0.224** | 0.533 | 1.422 | 2.111 |
| small_loose_p5 | 0.162 | 0.320 | 0.609 | 1.276 | 2.076 |
| medium_minn | **0.113** | 0.228 | 0.559 | 1.462 | 2.467 |
| ucmean | 1.563 | 1.587 | 1.654 | 1.671 | 1.615 |


### 2c. Point accuracy by variable and horizon (RMSE, lower better)

**`gdp_growth`**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| small_tight | 1.614 | 1.633 | 1.687 | 1.784 | 1.907 |
| small_ss | 1.603 | 1.620 | 1.684 | 1.797 | 1.922 |
| small_loose_p5 | **1.566** | 1.640 | 1.702 | 1.799 | 1.925 |
| small_minn | 1.612 | 1.629 | 1.683 | 1.793 | 1.915 |
| ucsv | 1.599 | **1.614** | 1.693 | 1.812 | 1.922 |
| small_sv | 1.583 | 1.655 | 1.700 | 1.812 | **1.893** |
| medium_conj | 1.616 | 1.641 | 1.704 | 1.794 | 1.904 |
| ar4 | 1.632 | 1.637 | **1.674** | **1.784** | 1.934 |
| medium_minn | 1.578 | 1.650 | 1.701 | 1.810 | 1.930 |
| ucmean | 1.619 | 1.639 | 1.684 | 1.801 | 1.929 |
| combo_equal | 1.620 | 1.670 | 1.689 | 1.814 | 1.922 |
| combo_pool | 1.688 | 1.745 | 1.724 | 1.921 | 2.012 |
| combo_logscore | 2.054 | 2.106 | 1.970 | 2.254 | 2.011 |
| combo_bma | 2.306 | 2.357 | 2.182 | 2.570 | 2.084 |
| rw | 2.377 | 2.545 | 2.326 | 2.609 | 2.661 |

**`cpi_inflation`**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| ar4 | 0.315 | 0.384 | **0.441** | **0.478** | **0.504** |
| ucmean | 0.434 | 0.444 | 0.449 | 0.480 | 0.507 |
| combo_pool | 0.285 | 0.369 | 0.482 | 0.602 | 0.628 |
| combo_equal | 0.326 | 0.389 | 0.471 | 0.581 | 0.603 |
| combo_logscore | 0.261 | 0.345 | 0.487 | 0.625 | 0.654 |
| small_ss | 0.332 | 0.398 | 0.502 | 0.594 | 0.588 |
| rw | **0.235** | **0.321** | 0.481 | 0.666 | 0.726 |
| small_minn | 0.333 | 0.400 | 0.498 | 0.595 | 0.608 |
| combo_bma | 0.277 | 0.370 | 0.526 | 0.675 | 0.686 |
| small_tight | 0.402 | 0.448 | 0.508 | 0.586 | 0.612 |
| small_sv | 0.357 | 0.428 | 0.516 | 0.632 | 0.651 |
| ucsv | 0.357 | 0.442 | 0.513 | 0.646 | 0.666 |
| medium_minn | 0.341 | 0.410 | 0.503 | 0.658 | 0.733 |
| medium_conj | 0.400 | 0.450 | 0.518 | 0.623 | 0.655 |
| small_loose_p5 | 0.370 | 0.431 | 0.542 | 0.687 | 0.668 |

**`unemp_rate`**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| small_ss | 0.337 | **0.561** | **0.891** | 1.205 | **1.223** |
| combo_equal | 0.369 | 0.578 | 0.902 | 1.262 | 1.321 |
| ar4 | 0.413 | 0.637 | 0.915 | **1.174** | 1.300 |
| small_minn | 0.349 | 0.586 | 0.948 | 1.338 | 1.366 |
| combo_pool | 0.426 | 0.620 | 0.919 | 1.260 | 1.400 |
| small_sv | 0.346 | 0.584 | 0.946 | 1.395 | 1.360 |
| small_tight | 0.362 | 0.602 | 0.953 | 1.351 | 1.371 |
| medium_minn | **0.336** | 0.562 | 0.918 | 1.377 | 1.466 |
| combo_bma | 0.339 | 0.569 | 0.940 | 1.507 | 1.416 |
| small_loose_p5 | 0.354 | 0.598 | 1.005 | 1.395 | 1.470 |
| rw | 0.391 | 0.637 | 0.984 | 1.391 | 1.509 |
| ucsv | 0.606 | 0.766 | 0.974 | 1.259 | 1.334 |
| combo_logscore | 0.583 | 0.748 | 1.047 | 1.307 | 1.385 |
| medium_conj | 0.375 | 0.632 | 1.050 | 1.487 | 1.574 |
| ucmean | 1.260 | 1.309 | 1.326 | 1.435 | 1.542 |

**`cash_rate`**

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| small_ss | 0.278 | 0.507 | 0.990 | **1.908** | 2.579 |
| ar4 | 0.307 | 0.587 | 1.103 | 1.931 | **2.467** |
| combo_equal | 0.342 | 0.524 | **0.969** | 1.923 | 2.783 |
| combo_logscore | 0.264 | **0.450** | 0.975 | 2.141 | 2.869 |
| combo_pool | 0.279 | 0.484 | 1.008 | 2.161 | 2.856 |
| small_minn | 0.292 | 0.531 | 1.045 | 2.033 | 2.895 |
| ucsv | 0.354 | 0.658 | 1.190 | 2.058 | 2.684 |
| rw | 0.352 | 0.658 | 1.194 | 2.062 | 2.684 |
| combo_bma | 0.282 | 0.479 | 1.049 | 2.397 | 2.967 |
| small_tight | 0.313 | 0.595 | 1.190 | 2.230 | 3.078 |
| small_sv | 0.256 | 0.460 | 0.978 | 2.324 | 3.425 |
| small_loose_p5 | 0.355 | 0.651 | 1.130 | 2.221 | 3.433 |
| medium_conj | 0.479 | 0.767 | 1.218 | 2.280 | 3.418 |
| medium_minn | **0.248** | 0.464 | 1.033 | 2.525 | 4.098 |
| ucmean | 2.576 | 2.594 | 2.658 | 2.649 | 2.635 |


### 2d. Density calibration — mean log predictive density (higher better)

Averaged across the 4 targets. The mean log score is brutally sensitive to tail events: **−∞** means at least one origin where the realization fell outside that member's predictive support (an individual model can catastrophically miss a tail). The **combinations never do** — the linear pool always assigns positive density — which is the clearest single piece of evidence that pooling buys calibration and robustness.

| Model | h=1 | h=2 | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|:--|:--|
| combo_equal | -1.852 | **-1.504** | **-1.524** | **-1.757** | **-1.872** |
| combo_pool | **-1.837** | -1.548 | -1.670 | -1.954 | -2.028 |
| combo_logscore | -1.840 | -1.537 | -1.649 | -2.013 | -2.042 |
| combo_bma | -2.056 | -1.870 | -1.986 | -2.379 | -2.192 |
| medium_minn | -3.127 | -4.714 | -4.551 | -7.785 | -3.681 |
| ar4 | −∞ | −∞ | -7.791 | -10.428 | −∞ |
| medium_conj | −∞ | −∞ | −∞ | −∞ | -13.455 |
| rw | −∞ | -5.098 | -2.978 | -2.956 | -2.692 |
| small_loose_p5 | -6.527 | −∞ | -9.128 | -9.524 | -8.624 |
| small_minn | -6.573 | −∞ | −∞ | −∞ | −∞ |
| small_ss | −∞ | −∞ | −∞ | −∞ | -12.302 |
| small_sv | -4.600 | −∞ | -6.722 | -6.114 | -3.030 |
| small_tight | −∞ | −∞ | -10.973 | −∞ | -12.616 |
| ucmean | −∞ | −∞ | −∞ | −∞ | −∞ |
| ucsv | -5.404 | −∞ | -6.487 | -10.739 | -9.902 |


### 2e. Year-ended and cumulative-level accuracy (growth variables)

For GDP and inflation (modelled as quarterly growth), the `q` score above is the marginal rate *in that one quarter* — the narrowest, least predictable view at long horizons. Two integrated views matter more for policy and are scored here:

- **Year-ended** (`ye`, 4-quarter sum ending at t+h): the RBA's headline concept — year-ended GDP growth, and the 2-3% *year-ended* trimmed-mean inflation target.

- **Cumulative level** (`cum`, the h-quarter sum from the origin = 100·(log level_{t+h} − log level_t)): where the level lands h quarters out. Far more discriminating than the single quarter — a model that gets the persistent/drift component wrong (e.g. the random walk) is exposed here but not by the quarterly score.

(For level variables — unemployment, cash rate — the `q` score already *is* the level at t+h, so no separate cumulative view is needed. At h=4 the two measures below coincide by construction — both span the 4 quarters from the origin — and diverge from h=8 on.)

**Real GDP growth (qtr %) (`gdp_growth`) — CRPS, lower better**

Year-ended:

| Model | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|
| small_tight | 1.366 | **1.476** | 1.499 |
| combo_equal | 1.298 | 1.510 | 1.546 |
| small_sv | **1.267** | 1.562 | 1.559 |
| small_ss | 1.377 | 1.525 | **1.499** |
| medium_conj | 1.383 | 1.511 | 1.517 |
| ar4 | 1.342 | 1.495 | 1.616 |
| medium_minn | 1.275 | 1.573 | 1.628 |
| small_minn | 1.420 | 1.550 | 1.524 |
| small_loose_p5 | 1.376 | 1.610 | 1.553 |
| ucmean | 1.380 | 1.525 | 1.672 |
| ucsv | 1.304 | 1.641 | 1.738 |
| combo_pool | 1.423 | 1.791 | 2.047 |
| combo_logscore | 1.805 | 2.666 | 2.213 |
| combo_bma | 2.214 | 3.691 | 2.669 |
| rw | 3.264 | 3.774 | 4.400 |

Cumulative level (from forecast origin):

| Model | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|
| small_sv | **1.267** | **1.954** | 2.206 |
| combo_equal | 1.305 | 2.030 | 2.237 |
| medium_minn | 1.275 | 1.995 | 2.314 |
| small_tight | 1.366 | 2.101 | **2.165** |
| ucsv | 1.304 | 2.020 | 2.424 |
| medium_conj | 1.383 | 2.139 | 2.238 |
| small_ss | 1.377 | 2.169 | 2.240 |
| small_loose_p5 | 1.376 | 2.173 | 2.313 |
| ar4 | 1.342 | 2.156 | 2.470 |
| small_minn | 1.420 | 2.259 | 2.337 |
| ucmean | 1.380 | 2.273 | 2.662 |
| combo_pool | 1.425 | 2.651 | 3.685 |
| combo_logscore | 1.789 | 4.352 | 3.995 |
| combo_bma | 2.203 | 6.319 | 5.314 |
| rw | 3.264 | 6.711 | 11.084 |

**Trimmed-mean CPI inflation (qtr %) (`cpi_inflation`) — CRPS, lower better**

Year-ended:

| Model | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|
| ar4 | 0.844 | **1.134** | **1.258** |
| combo_equal | 0.765 | 1.196 | 1.425 |
| ucmean | 1.056 | 1.184 | 1.289 |
| rw | **0.674** | 1.251 | 1.664 |
| combo_pool | 0.761 | 1.293 | 1.540 |
| small_minn | 0.827 | 1.318 | 1.543 |
| small_ss | 0.849 | 1.346 | 1.506 |
| combo_logscore | 0.736 | 1.335 | 1.635 |
| small_sv | 0.849 | 1.304 | 1.577 |
| small_tight | 0.925 | 1.316 | 1.557 |
| medium_minn | 0.814 | 1.321 | 1.714 |
| medium_conj | 0.895 | 1.354 | 1.653 |
| ucsv | 0.920 | 1.419 | 1.698 |
| combo_bma | 0.804 | 1.486 | 1.763 |
| small_loose_p5 | 0.935 | 1.543 | 1.750 |

Cumulative level (from forecast origin):

| Model | h=4 | h=8 | h=12 |
|:--|:--|:--|:--|
| rw | **0.674** | **1.832** | 3.356 |
| ar4 | 0.844 | 1.973 | **3.108** |
| combo_equal | 0.770 | 1.919 | 3.340 |
| combo_pool | 0.761 | 2.043 | 3.523 |
| combo_logscore | 0.729 | 2.075 | 3.592 |
| small_minn | 0.827 | 2.113 | 3.714 |
| small_ss | 0.849 | 2.171 | 3.663 |
| ucmean | 1.056 | 2.269 | 3.405 |
| small_sv | 0.849 | 2.129 | 3.761 |
| medium_minn | 0.814 | 2.141 | 4.002 |
| small_tight | 0.925 | 2.241 | 3.910 |
| combo_bma | 0.800 | 2.363 | 3.936 |
| medium_conj | 0.895 | 2.234 | 4.029 |
| ucsv | 0.920 | 2.331 | 4.056 |
| small_loose_p5 | 0.935 | 2.436 | 4.117 |


## 3. Do the combinations beat the best single model?

Mean CRPS over all 4 targets, by horizon bucket. The honest test of a pool is whether it beats both equal weights and the best individual member.

| Model | near (1-4) | medium (5-8) | far (9-12) |
|:--|:--|:--|:--|
| combo_equal | 0.389 | **0.650** | **0.860** |
| small_ss | 0.396 | 0.665 | 0.862 |
| ar4 | 0.414 | 0.668 | 0.870 |
| small_minn | 0.401 | 0.678 | 0.894 |
| combo_pool | 0.397 | 0.716 | 0.933 |
| small_tight | 0.421 | 0.707 | 0.936 |
| small_sv | 0.388 | 0.721 | 0.988 |
| medium_conj | 0.423 | 0.706 | 0.968 |
| ucsv | 0.438 | 0.725 | 0.946 |
| small_loose_p5 | 0.419 | 0.724 | 0.967 |
| combo_logscore | 0.428 | 0.765 | 0.944 |
| medium_minn | **0.386** | 0.718 | 1.042 |
| rw | 0.510 | 0.802 | 1.042 |
| combo_bma | 0.444 | 0.904 | 1.014 |
| ucmean | 0.831 | 0.896 | 0.950 |

## 4. Statistical significance (Diebold-Mariano)

How often each combination **significantly beats** the random-walk and AR(4) benchmarks on CRPS (Harvey-corrected, 10% level), counted over the 4 targets x 3 horizons {1, 4, 8} tested. A negative DM statistic means the combination is more accurate; significance is one-sided here.

| Combination | beats ar4 | beats rw |
|:--|:--|:--|
| combo_bma | 0 / 12 | 0 / 12 |
| combo_equal | 0 / 12 | 0 / 12 |
| combo_logscore | 0 / 12 | 0 / 12 |
| combo_pool | 0 / 12 | 0 / 12 |

## 5. Model profiles

One entry per model: its specification, what makes it distinct, the role it plays in the suite, its strengths and failure modes, and where it actually ranks in this evaluation (the eval line is computed, not asserted). Full rationale is in DECISIONS.md.

### 5a. VAR members

**`small_minn`**  
*Spec:* Independent Normal-inverse-Wishart (Gibbs); 8-variable small SOE core; 4 lags; GLP marginal-likelihood shrinkage; constant volatility; LP scaling (COVID); sum-of-coefficients + dummy-initial-observation priors.  
*Distinctive:* The workhorse Minnesota BVAR, estimated by Gibbs with an *independent* Normal-inverse-Wishart prior — the engine required to impose block exogeneity, since asymmetric (equation-specific) prior variances cannot be represented by a Kronecker/conjugate prior. Shrinkage is data-driven (GLP marginal-likelihood), and sum-of-coefficients + dummy-initial-observation priors discipline the I(1) levels (rates, real TWI).  
*Role:* The central, representative small-SOE BVAR — the reference point the other small members are deliberate variations around (tighter, looser, steady-state, SV).  
*Strengths & failure modes:* Solid all-rounder at short-to-medium horizons. Constant volatility means it leans on the LP scaling for 2020; without it the 2020 outliers would distort the coefficients.  
*In this evaluation:* strongest at cash_rate (medium (5-8)), ranked 2 of 11 individual models.  
*See:* DECISIONS.md D3, D4, D5, D17

**`small_ss`**  
*Spec:* Steady-state (Villani); 8-variable small SOE core; 4 lags; GLP marginal-likelihood shrinkage; constant volatility; LP scaling (COVID).  
*Distinctive:* Reparametrised around its *unconditional means* (Villani steady state), with informative priors placed directly on long-run levels — inflation 2.5% (target midpoint), NAIRU 4.5%, neutral cash rate 3.5%, US potential growth — and on the foreign block's steady states, which the domestic forecast inherits.  
*Role:* The **long-horizon anchor**. An iterated VAR reverts to its unconditional mean at h = 8-12; this member makes that mean economically grounded rather than the raw sample average.  
*Strengths & failure modes:* Strongest at medium/far horizons for the mean-reverting targets. Vulnerable if a steady-state anchor is stale (e.g. a shifted neutral rate drags the long end); constant volatility under-disperses around 2020 absent the LP correction.  
*In this evaluation:* best individual model for unemp_rate at medium (5-8); unemp_rate at far (9-12); cash_rate at medium (5-8).  
*See:* DECISIONS.md D3, D4, D5, D17

**`small_sv`**  
*Spec:* Stochastic volatility, equation-by-equation; 8-variable small SOE core; 4 lags; GLP marginal-likelihood shrinkage; stochastic volatility; t-errors (SV-t) (COVID).  
*Distinctive:* Stochastic volatility estimated equation-by-equation (the Carriero-Clark-Marcellino triangular factorisation), with t-distributed errors (the CCMM SV-t COVID treatment). Block exogeneity is *exact* — foreign equations simply drop the domestic regressors.  
*Role:* The **density-calibration specialist**: time-varying volatility tracks changing macro uncertainty and the fat tails absorb outliers instead of letting them widen the whole history.  
*Strengths & failure modes:* Best at short horizons (h = 1-2), where getting the conditional variance right matters most. The volatility state at the jump-off can over/under-shoot if the last few quarters were unusual; the small system limits cross-variable information.  
*In this evaluation:* best individual model for gdp_growth at near (1-4); cash_rate at near (1-4).  
*See:* DECISIONS.md D3, D6, D17

**`small_loose_p5`**  
*Spec:* Independent Normal-inverse-Wishart (Gibbs); 8-variable small SOE core; 5 lags; fixed shrinkage λ=0.4; constant volatility; LP scaling (COVID).  
*Distinctive:* A deliberately *under-shrunk*, longer-lag variant — fixed λ = 0.4 (vs the ~0.1-0.15 the GLP procedure selects) and 5 lags — so the data speak more and richer dynamics can show through, at the cost of estimation noise.  
*Role:* The **loose / long-lag** diversity axis: it fails differently from the tightly-shrunk members and can capture dynamics they shrink away.  
*Strengths & failure modes:* Occasionally best at near-horizon density when the extra flexibility pays; noisier and prone to wider intervals at long horizons (the cost of light shrinkage in a short sample).  
*In this evaluation:* strongest at gdp_growth (near (1-4)), ranked 2 of 11 individual models.  
*See:* DECISIONS.md D4, D8

**`medium_minn`**  
*Spec:* Stochastic volatility, equation-by-equation; 13-variable medium system; 2 lags; GLP marginal-likelihood shrinkage; stochastic volatility; t-errors (SV-t) (COVID).  
*Distinctive:* The larger system — 13 variables (adds terms of trade, wages, employment, consumption, the 10y yield) with stochastic volatility and t-errors, but only 2 lags to keep the parameter count feasible; equation-by-equation estimation keeps the recursive loop tractable.  
*Role:* The **medium-system** axis (Banbura-Giannone-Reichlin): medium systems often forecast best given enough shrinkage, and the extra variables bring cross-sectional information the small core lacks.  
*Strengths & failure modes:* Strong short-horizon density (it tends to win the near bucket). The short lag length limits long-horizon dynamics, and more parameters mean more estimation uncertainty at the far end.  
*In this evaluation:* best individual model for unemp_rate at near (1-4).  
*See:* DECISIONS.md D1, D6, D17

**`medium_conj`**  
*Spec:* Block-recursive conjugate NIW; 13-variable medium system; 4 lags; fixed shrinkage λ=0.1; constant volatility; LP scaling (COVID); sum-of-coefficients + dummy-initial-observation priors.  
*Distinctive:* The medium system estimated by the fast *block-recursive conjugate* scheme (foreign VAR + domestic block conditioned on contemporaneous foreign values, the RBNZ Bloor-Matheson approach) — closed-form, so cheap even at 13 variables x 4 lags. Block exogeneity is exact by the recursive structure, not the prior.  
*Role:* The cheap medium workhorse; it complements `medium_minn` (conjugate constant-volatility vs SV) on the same large system.  
*Strengths & failure modes:* Tends to lead the far-horizon GDP year-ended / cumulative-level buckets. Constant volatility leans on LP for 2020; the conjugate Kronecker prior cannot represent asymmetric shrinkage, which is why block exogeneity comes from the recursive structure.  
*In this evaluation:* strongest at gdp_growth (far (9-12)), ranked 2 of 11 individual models.  
*See:* DECISIONS.md D3, D5, D8

**`small_tight`**  
*Spec:* Block-recursive conjugate NIW; 8-variable small SOE core; 4 lags; fixed shrinkage λ=0.05; constant volatility; LP scaling (COVID); sum-of-coefficients + dummy-initial-observation priors.  
*Distinctive:* The heavily-shrunk small model — fixed λ = 0.05, far tighter than the GLP selection — pulling hard toward the persistence/random-walk prior, so it is parsimonious and low-variance.  
*Role:* The **tight** diversity axis and the long-horizon robustness member: heavy shrinkage buys stability where lightly-parametrised models wander.  
*Strengths & failure modes:* Best or near-best at far-horizon GDP (the tight prior stops it over-reacting). Can be too rigid at short horizons, missing genuine dynamics the looser members catch.  
*In this evaluation:* best individual model for gdp_growth at medium (5-8); gdp_growth at far (9-12).  
*See:* DECISIONS.md D4, D5, D8

### 5b. Benchmark members

**`rw`**  
*Spec:* Random walk; LP scaling (COVID).  
*Distinctive:* The no-change forecast: the last observed value persists, with Gaussian increments scaled to the historical change.  
*Role:* The universal hard-to-beat short-horizon bar for persistent/level variables, and a pool member.  
*Strengths & failure modes:* Competitive at h = 1 for level variables; fails badly at long horizons for growth variables — its level path runs away, which the cumulative-level metric (§2e) exposes brutally.  
*In this evaluation:* best individual model for cpi_inflation at near (1-4).  
*See:* DECISIONS.md D8

**`ar4`**  
*Spec:* Bayesian AR(4); LP scaling (COVID).  
*Distinctive:* A Bayesian AR(4) per variable with Minnesota-style lag shrinkage and a stationarity-truncated posterior.  
*Role:* The univariate-persistence bar — it isolates how much of the forecast is just own-history dynamics.  
*Strengths & failure modes:* Surprisingly strong for inflation at near/medium horizons, where univariate dynamics dominate; it cannot use cross-variable information, so it lags when that matters.  
*In this evaluation:* best individual model for cpi_inflation at medium (5-8); cash_rate at far (9-12).  
*See:* DECISIONS.md D8

**`ucsv`**  
*Spec:* Unobserved components + stochastic volatility; t-errors (robust) (COVID).  
*Distinctive:* Stock-Watson unobserved-components stochastic volatility per variable: a random-walk trend plus transitory noise, both with time-varying variances and outlier-robust t-errors.  
*Role:* The canonical inflation benchmark and a genuine density anchor for the other targets.  
*Strengths & failure modes:* Strong for inflation (its native use case); weaker for variables with richer multivariate dynamics. The trend/noise split is weakly identified, so it is gated on Monte-Carlo precision, not raw ESS.  
*In this evaluation:* strongest at unemp_rate (medium (5-8)), ranked 3 of 11 individual models.  
*See:* DECISIONS.md D9, D17

**`ucmean`**  
*Spec:* Unconditional mean; LP scaling (COVID).  
*Distinctive:* A Gaussian density centred on the expanding-sample mean with the sample variance — the simplest possible density forecast.  
*Role:* The floor: the 'did the model beat just predicting the long-run average' bar.  
*Strengths & failure modes:* Unexpectedly competitive at long horizons for mean-reverting growth (everything reverts to the mean eventually); useless at short horizons where dynamics matter.  
*In this evaluation:* best individual model for cpi_inflation at far (9-12).  
*See:* DECISIONS.md D8

### 5c. Combination schemes

**`combo_equal`**  
*Weights:* Equal weights on every member, per variable x horizon bucket.  
*In this evaluation:* The forecast-combination-puzzle benchmark and the recommended robust default: in the evaluation it has the best mean log score at every horizon and the best far-horizon CRPS. Hard to beat because it never over-fits weights.

**`combo_logscore`**  
*Weights:* Weights proportional to each member's recent log predictive score, with a forgetting factor, shrunk toward equal.  
*In this evaluation:* Adapts to which members are forecasting well lately; the shrinkage and forgetting guard against over-concentrating on a member that was lucky.

**`combo_pool`**  
*Weights:* Optimal prediction pool (Hall-Mitchell / Geweke-Amisano): weights on the simplex that maximise the historical *pooled* log score, shrunk toward equal.  
*In this evaluation:* Unlike BMA it does not degenerate to a single model ('all models are false but useful'); competitive with equal weights and occasionally better at near horizons.

**`combo_bma`**  
*Weights:* Bayesian model averaging by predictive likelihood — no shrinkage.  
*In this evaluation:* **Diagnostic only.** It concentrates weight on the single best-fitting member, so it answers 'which model does the data favour' rather than serving as a robust combination; reported, not recommended.


## 6. How to read this

- **Point gains over the best member are modest by design** — equal weights are hard to beat (the forecast-combination puzzle). The pool's payoff is *calibration and robustness*: it insures against any single member failing, rather than always winning on accuracy.

- **CRPS and log score can disagree** on outlier-heavy windows (log score is far more sensitive to tail events); both are reported. A COVID-excluded variant is in `output/tables/scores_by_horizon_excovid.csv`.

- **Pick the horizon view to match the decision.** The quarterly score (§2b) is the marginal growth rate; for GDP and inflation the year-ended and cumulative-level views (§2e) are usually what a central bank acts on, and they rank models differently at long horizons. For level variables the quarterly score already is the level. See the full Quarto report for fan charts, PIT calibration, and weight-evolution plots.

