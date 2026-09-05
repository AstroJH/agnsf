"""Single-curve variability statistics (user-facing Python layer).

All numbers are computed by the C++ core (`agnsf::variability`); this
module only wraps the result:

    result = agnsf.variability.measure(t, v, e)
    print(result.fvar, result.fvar_uncertainty)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ._agnsf import (
    Uncertainty,
    variability_measure as _variability_measure,
)

__all__ = ["measure", "VariabilityResult"]


@dataclass
class VariabilityResult:
    """Variability statistics of one light curve.

    Definitions follow Vaughan et al. (2003); ``F_var`` (fractional
    variability amplitude), ``nxs`` (normalized excess variance) and
    ``xs`` (excess variance) carry analytic uncertainties from the
    same reference (``None`` when not defined).
    """

    n: int
    valid: bool
    mean: float
    weighted_mean: float
    weighted_mean_error: float
    stddev: float
    stddev_error: float
    peak_to_peak: float
    peak_to_peak_noise_corrected: float
    sigma_m: float
    fvar: float
    fvar_uncertainty: Uncertainty | None
    nxs: float
    nxs_uncertainty: Uncertainty | None
    xs: float
    xs_uncertainty: Uncertainty | None
    chi2: float
    chi2_dof: float
    chi2_q: float
    von_neumann: float

    @property
    def fvar_error(self) -> float:
        """Symmetric uncertainty of ``fvar`` (NaN when not estimated)."""
        if self.fvar_uncertainty is None:
            return float("nan")
        return 0.5 * (self.fvar_uncertainty.upper - self.fvar_uncertainty.lower)

    @property
    def nxs_error(self) -> float:
        if self.nxs_uncertainty is None:
            return float("nan")
        return 0.5 * (self.nxs_uncertainty.upper - self.nxs_uncertainty.lower)

    @property
    def xs_error(self) -> float:
        if self.xs_uncertainty is None:
            return float("nan")
        return 0.5 * (self.xs_uncertainty.upper - self.xs_uncertainty.lower)


def measure(
    time: Sequence[float],
    value: Sequence[float],
    error: Sequence[float],
    *,
    err_sys: float = 0.0,
    weighted: bool = False,
) -> VariabilityResult:
    """Compute single-curve variability statistics.

    Parameters
    ----------
    time, value, error:
        Observed time, flux, and per-point uncertainty arrays.
    err_sys:
        Systematic error floor added in quadrature when subtracting
        the noise contribution.
    weighted:
        Use the inverse-variance weighted mean inside the variance /
        intrinsic-amplitude estimates.
    """
    raw = _variability_measure(
        time,
        value,
        error,
        float(err_sys),
        bool(weighted),
    )

    return VariabilityResult(
        n=raw.n,
        valid=raw.valid,
        mean=raw.mean,
        weighted_mean=raw.weighted_mean,
        weighted_mean_error=raw.weighted_mean_error,
        stddev=raw.stddev,
        stddev_error=raw.stddev_error,
        peak_to_peak=raw.peak_to_peak,
        peak_to_peak_noise_corrected=raw.peak_to_peak_noise_corrected,
        sigma_m=raw.sigma_m,
        fvar=raw.fvar,
        fvar_uncertainty=raw.fvar_uncertainty,
        nxs=raw.nxs,
        nxs_uncertainty=raw.nxs_uncertainty,
        xs=raw.xs,
        xs_uncertainty=raw.xs_uncertainty,
        chi2=raw.chi2,
        chi2_dof=raw.chi2_dof,
        chi2_q=raw.chi2_q,
        von_neumann=raw.von_neumann,
    )
