"""Time-delay (lag) analysis — user-facing Python layer.

Wraps the low-level cross-correlation bindings with a Pythonic API for
AGN reverberation / torus lag measurements.

Example:

    result = agnsf.timedelay.lag(
        t_opt, f_opt, e_opt,      # continuum (optical)
        t_nir, f_nir, e_nir,      # response (near-IR)
        lag_range=(-50, 50), step=1.0,
        method="iccf", estimate="centroid",
        uncertainty="fr_rss", n_realizations=1000, seed=42,
    )
    print(result.lag, result.lower, result.upper)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from ._agnsf import (
    CrossCorrelationMethod,
    LagEstimate,
    Uncertainty,
    cross_correlate as _cross_correlate,
    lag_uncertainty as _lag_uncertainty,
)

__all__ = ["lag", "LagResult", "CrossCorrelationMethod", "LagEstimate"]

_METHODS = {
    "dcf": CrossCorrelationMethod.Dcf,
    "iccf": CrossCorrelationMethod.Iccf,
}

_ESTIMATES = {
    "peak": LagEstimate.Peak,
    "centroid": LagEstimate.Centroid,
}


def _coerce_method(method: str | CrossCorrelationMethod) -> CrossCorrelationMethod:
    if isinstance(method, CrossCorrelationMethod):
        return method
    key = str(method).strip().lower()
    if key in _METHODS:
        return _METHODS[key]
    raise ValueError(f"method must be 'dcf' or 'iccf', got {method!r}")


def _coerce_estimate(estimate: str | LagEstimate) -> LagEstimate:
    if isinstance(estimate, LagEstimate):
        return estimate
    key = str(estimate).strip().lower()
    if key in _ESTIMATES:
        return _ESTIMATES[key]
    raise ValueError(f"estimate must be 'peak' or 'centroid', got {estimate!r}")


@dataclass
class LagResult:
    """Result of a time-delay analysis.

    Attributes
    ----------
    lag:
        The selected lag estimate (peak or centroid).
    lag_peak, lag_centroid:
        Both estimates from the cross-correlation curve.
    lower, upper:
        FR/RSS uncertainty interval on the selected lag (NaN if not
        estimated).
    peak_value:
        Maximum correlation value.
    tau, ccf, count:
        Trial-lag grid, correlation values, and per-bin counts.
    """

    lag: float
    lag_peak: float
    lag_centroid: float
    lower: float
    upper: float
    peak_value: float
    tau: np.ndarray
    ccf: np.ndarray
    count: np.ndarray


def lag(
    time1: Sequence[float],
    value1: Sequence[float],
    error1: Sequence[float],
    time2: Sequence[float],
    value2: Sequence[float],
    error2: Sequence[float],
    *,
    lag_range: tuple[float, float] = (-50.0, 50.0),
    step: float = 1.0,
    method: str | CrossCorrelationMethod = "dcf",
    estimate: str | LagEstimate = "peak",
    uncertainty: str | bool = "fr_rss",
    n_realizations: int = 1000,
    seed: int = 0,
    dcf_bin_width: float | None = None,
    centroid_threshold: float = 0.8,
    min_overlap: int = 3,
) -> LagResult:
    """Estimate the time delay between a continuum and a response light curve.

    Positive lag means the response lags the continuum
    (response(t + lag) correlates with continuum(t)).

    Parameters
    ----------
    time1, value1, error1:
        Continuum (driving) light curve.
    time2, value2, error2:
        Response light curve.
    lag_range:
        (min, max) trial-lag grid in the same time units.
    step:
        Trial-lag spacing.
    method:
        "dcf" (discrete correlation function) or "iccf" (interpolated
        cross-correlation function).
    estimate:
        "peak" or "centroid".
    uncertainty:
        "fr_rss" (or True) enables FR/RSS Monte Carlo; None/False disables.
    n_realizations, seed:
        FR/RSS realizations and RNG seed.
    dcf_bin_width:
        DCF bin width (defaults to ``step``).
    centroid_threshold:
        Fraction of the peak used as the centroid threshold.
    min_overlap:
        Minimum overlapping points for a valid correlation.
    """
    lo, hi = lag_range
    method_enum = _coerce_method(method)
    estimate_enum = _coerce_estimate(estimate)
    bin_width = step if dcf_bin_width is None else float(dcf_bin_width)

    raw = _cross_correlate(
        time1, value1, error1,
        time2, value2, error2,
        float(lo), float(hi), float(step),
        method_enum,
        bin_width,
        float(centroid_threshold),
        int(min_overlap),
    )

    if uncertainty:
        interval = _lag_uncertainty(
            time1, value1, error1,
            time2, value2, error2,
            float(lo), float(hi), float(step),
            estimate_enum,
            method_enum,
            bin_width,
            float(centroid_threshold),
            int(min_overlap),
            int(n_realizations),
            int(seed),
            True,
            True,
        )
        lower = interval.lower
        upper = interval.upper
    else:
        lower = float("nan")
        upper = float("nan")

    lag_value = (
        raw.lag_centroid
        if estimate_enum == LagEstimate.Centroid
        else raw.lag_peak
    )

    return LagResult(
        lag=lag_value,
        lag_peak=raw.lag_peak,
        lag_centroid=raw.lag_centroid,
        lower=lower,
        upper=upper,
        peak_value=raw.peak_value,
        tau=np.asarray(raw.tau),
        ccf=np.asarray(raw.ccf),
        count=np.asarray(raw.count),
    )
