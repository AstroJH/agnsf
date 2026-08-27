# AGN Structure Function (AGNSF)

**AGNSF** is a C++ library with Python bindings for calculating structure functions of astronomical light curves.

> 🚧 **Work in progress.**  
> The project is currently under active development and is not yet ready for general use.

## Planned features

- [x] Structure function (SF) calculation

- [x] Ensemble structure function
  - [x] Pooled structure function
  - [x] Aggregated structure function

- [ ] Monte Carlo uncertainty estimation

- [x] Python interface through NumPy and pybind11

- [x] Efficient NumPy interface with zero-copy access for compatible arrays;
  automatic conversion for other inputs

- [x] Command-line interface for calculating structure functions directly from
  light-curve files and writing results to FITS or CSV

- [x] Python interface supporting both in-memory light curves and file paths

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

## Uncertainty estimation

AGNSF can estimate uncertainties on the structure function for every lag
bin. Uncertainty is reported as an **asymmetric interval** `[lower, upper]`
(`SFUncertainty`), where NaN means "not estimated". Symmetric estimates are
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

## Command-line interface

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
path-list file containing one light-curve file path per line. For example:

```
data/lc1.csv
data/lc2.csv
data/lc3.fits
```

For all available options, SF/ESF methods, bin specifications, column
settings, uncertainty estimation, and output format, see:

```sh
agnsf --help
```

