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

- [ ] Command-line interface for calculating structure functions directly from
  light-curve files and writing results to FITS or CSV

- [ ] Python interface supporting both in-memory light curves and file paths

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