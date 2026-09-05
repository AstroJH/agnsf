# AGN Structure Function (AGNSF)

**AGNSF** is a C++ library with Python bindings for calculating structure functions of astronomical light curves.

> 🚧 **Work in progress.**  
> The project is currently under active development and is not yet ready for general use.
> 
> ⚠️ **Note**  
> AGNSF is currently tested on macOS with a C++17 compiler and Python 3.

## Planned Features

- [x] Structure function (SF) calculation
- [x] Ensemble structure function
  - [x] Pooled structure function
  - [x] Aggregated structure function
- [x] Uncertainty estimation *(currently for SF/ESF only)*
  - [x] Measurement uncertainty
  - [x] Sampling uncertainty
  - [x] Bootstrap and jackknife methods
- [x] Time-delay analysis
  - [x] DCF and ICCF
  - [x] Peak and centroid lag estimates
  - [x] FR/RSS uncertainty estimation
  - [x] Transfer-function fitting
- [x] Variability amplitude estimation
- [ ] Power spectral density (PSD) analysis *(still a long way off)*
- [x] Python interface through NumPy and pybind11
- [x] Efficient NumPy interface with zero-copy access for compatible arrays;
  automatic conversion for other inputs
- [x] CLI for calculating SF/ESF directly from light-curve files and writing results to FITS or CSV

## Installation and Build
### Build CLI

```sh
cmake -S . -B build
cmake --build build -j
```

### Install with pip

```sh
python -m pip install git+https://github.com/AstroJH/agnsf.git
```

This installs both the Python bindings and the AGNSF CLI.

## Structure Function Estimators

The structure function can be computed with several estimators, selected through
`SFMethod` (`sf`, `pooled_sf`, and the per-curve step of `ensemble_sf`):

- **SecondOrder** (default) — noise-corrected second-order estimator:

  $$\mathrm{SF}^2(\tau)=
  \left\langle \Delta^2 \right\rangle -
  \left\langle \sigma_i^2 + \sigma_j^2 \right\rangle
  $$

- **SecondOrderNoNoise** — second-order estimator without noise subtraction:

  $$\mathrm{SF}^2(\tau)=
  \left\langle \Delta^2 \right\rangle
  $$

- **MeanAbsoluteDeviation** — mean-absolute-deviation based estimator, noise-corrected:

  $$\mathrm{SF}^2(\tau)=
  \frac{\pi}{2}\left\langle |\Delta| \right\rangle^2 -
  \left\langle \sigma_i^2 + \sigma_j^2 \right\rangle
  $$

- **MeanAbsoluteDeviationNoNoise** — without noise subtraction:

  $$\mathrm{SF}^2(\tau)=
  \frac{\pi}{2}\left\langle |\Delta| \right\rangle^2
  $$

where $\Delta = m_j - m_i$ runs over all pairs in the lag bin. In every case
$\mathrm{SF}(\tau)=\sqrt{\mathrm{SF}^2(\tau)}$ when the result is finite and
non-negative, otherwise it is NaN.

## Ensemble Structure Functions

AGNSF supports two approaches for calculating ensemble structure functions:

- **Pooled ESF** — combines all valid light-curve pairs across the ensemble and calculates the SF from the pooled pairs. Each pair contributes equally.

- **Aggregated ESF** — calculates the SF independently for each light curve and combines the individual SFs using one of two methods:
  - *Root-mean-square (default)*:

    $$\mathrm{ESF}(\tau)=
    \sqrt{\left\langle \mathrm{SF}_k^2(\tau) \right\rangle_k}
    $$

  - *Mean*:

    $$\mathrm{ESF}(\tau)=
    \left\langle \mathrm{SF}_k(\tau) \right\rangle_k.
    $$

  Each contributing light curve is weighted equally within each lag bin. Only light curves with a finite contribution are included in each bin (a finite SF² for the root-mean-square method, a finite SF for the mean method).

## Uncertainty Estimation

AGNSF can estimate uncertainties on the structure function for every lag
bin. Uncertainty is reported as an **asymmetric interval** `[lower, upper]`
(`Uncertainty`), where NaN means "not estimated". Symmetric estimates are
represented by `lower == upper`.

### Configuration

`UncertaintyConfig` controls estimation:

- `measurement` (`Off` / `Analytic`): per-bin measurement uncertainty.
- `sampling` (`Off` / `Analytic` / `Jackknife` / `Bootstrap`): source-to-source
  sampling uncertainty **(ESF only)**.
- `n_bootstrap`, `bootstrap_seed`: bootstrap replicates and RNG seed (fixed
  seed gives reproducible results).

```python
import agnsf

cfg = agnsf.UncertaintyConfig()
cfg.measurement = agnsf.UncertaintyMethod.Analytic
cfg.sampling    = agnsf.UncertaintyMethod.Bootstrap
cfg.n_bootstrap = 500
cfg.bootstrap_seed = 42

r = agnsf.ensemble_sf(times, values, errors, bins, uncertainty=cfg)
bin0 = r.bins[0]

# point estimate
bin0.sf

# propagated measurement uncertainty
(bin0.measurement.lower, bin0.measurement.upper)

# bootstrap 16/84 percentile interval
(bin0.sampling.lower, bin0.sampling.upper)
```

### Measurement uncertainty (SF and ESF)

The analytic estimator propagates the **within-bin standard error of the
mean** of the per-pair quantities to SF:

- Second-order estimators:
 
  $$\mathrm{SE}=
  \frac{\mathrm{std}(\Delta^2 - \mathrm{noise})}{\sqrt{N}}
  $$

  or, for the no-noise variants,

  $$\mathrm{SE}=
  \frac{\mathrm{std}(\Delta^2)}{\sqrt{N}}.
  $$

- Mean-absolute-deviation estimators: the standard error of
  $\langle|\Delta|\rangle$ is propagated through
  $f(x) = \frac{\pi}{2}x^2$
  using the delta method.

For ESF, per-curve measurement uncertainties are combined in quadrature
(independence assumption) and propagated through the aggregation.

> Note: pairs sharing points are treated as independent; this is an
> approximation. A Monte Carlo perturbation option is planned as an
> alternative.

### Sampling uncertainty (ESF only)

Sampling uncertainty exploits the independence of different light curves:

- `Analytic` (aggregated ESF, default):

  $$\mathrm{SE} =
  \frac{\mathrm{std}(x_1,\ldots,x_n)}{\sqrt{n}}
  $$

  where $x_i$ are the per-curve quantities used by the aggregation.

- `Jackknife`: leave-one-curve-out recomputation of the ESF.
- `Bootstrap`: curve-level resampling with replacement; the reported interval
  is the 16th/84th percentile (naturally asymmetric).

`Analytic` sampling is not defined for the pooled ESF (no per-curve values are
averaged); use `Jackknife` or `Bootstrap` there.

A single light curve never estimates sampling uncertainty (`sampling` stays
NaN), and measurement/sampling are kept as separate fields by design.

## Rest-frame lags

All Python computation and file functions accept a `redshift` keyword.
When given, observed times are divided by `1 + z`, so lag bins are
interpreted in the source rest frame:

```python
r = agnsf.sf(time, value, error, bins, redshift=0.5)  # dt_rest = dt_obs / (1+z)
```

The redshift is applied inside the native kernels (i.e., C++ lib; per pair, without
copying the time arrays), and is also available in the CLI
(`--redshift <z>`) and in path-list files (per-curve, optional second
column; a missing column means `z = 0`).

## Time-delay Analysis

AGNSF provides a lightweight interface for measuring time delays between two
astronomical light curves. It currently supports the **discrete correlation
function (DCF)** and **interpolated cross-correlation function (ICCF)**, with
either peak or centroid lag estimates.

The main interface is `agnsf.timedelay.lag()`:

```python
import agnsf

result = agnsf.timedelay.lag(
    time1, value1, error1,
    time2, value2, error2,
    lag_range=(-50, 50),
    step=1.0,
    method="dcf",
    estimate="centroid",
    uncertainty="fr_rss",
    n_realizations=500,
    seed=0,
)
```

Here, the first light curve is treated as the continuum (driving) light curve
and the second as the response light curve. A positive lag means that the
response lags the continuum.

`method` can be `"dcf"` or `"iccf"`, and `estimate` can be `"peak"` or
`"centroid"`. The centroid is calculated from the CCF points above a specified
fraction of the maximum correlation coefficient.

FR/RSS uncertainty estimation is available through
`uncertainty="fr_rss"`. The resulting interval is reported as the 16th and
84th percentiles of the lag distribution from the Monte Carlo realizations.
Set `uncertainty=None` to disable uncertainty estimation.

The returned `LagResult` contains the selected lag estimate, both peak and
centroid estimates, the FR/RSS interval when requested, and the full
cross-correlation curve (`tau`, `ccf`, and `count`).

## Command-line Interface (CLI)

The AGNSF command-line program computes the structure function of a
single light curve or the ensemble structure function of multiple
light curves.

```sh
# SF of one light curve
agnsf --input lc.csv --output sf.txt

# Pooled ESF from a path-list file
agnsf --input @list.txt --output esf.txt

# Aggregated ESF with bootstrap sampling uncertainty
agnsf -i @list.txt -o esf.txt \
      --esf aggregated \
      --sampling bootstrap --n-bootstrap 500 --bootstrap-seed 42
```

When the input path is prefixed with `@`, the file is interpreted as a
path-list file containing one light-curve file path per line. Each
line may optionally carry a second column with the source redshift
(per-curve rest-frame correction):

```
data/lc1.csv
data/lc2.csv 0.5
data/lc3.fits
```

For single light-curve input, pass `--redshift <z>` for a global
rest-frame correction.


For all available options, SF/ESF methods, bin specifications, column
settings, uncertainty estimation, and output format, see:

```sh
agnsf --help
```

