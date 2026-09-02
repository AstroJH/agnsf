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
from typing import Mapping, Sequence

import numpy as np

from ._agnsf import (
    CrossCorrelationMethod,
    FitParameter,
    LagEstimate,
    OptimizationAlgorithm,
    TransferFunctionShape,
    Uncertainty,
    cross_correlate as _cross_correlate,
    default_transfer_function_parameters as _default_transfer_function_parameters,
    fit_transfer_function as _fit_transfer_function,
    lag_uncertainty as _lag_uncertainty,
    transfer_function_curve as _transfer_function_curve,
    transfer_function_model_response as _transfer_function_model_response,
)

__all__ = [
    "lag",
    "LagResult",
    "CrossCorrelationMethod",
    "LagEstimate",
    "fit_transfer_function",
    "TransferFunctionResult",
    "TransferFunctionShape",
    "transfer_function_curve",
    "model_response",
]

_METHODS = {
    "dcf": CrossCorrelationMethod.Dcf,
    "iccf": CrossCorrelationMethod.Iccf,
}

_ESTIMATES = {
    "peak": LagEstimate.Peak,
    "centroid": LagEstimate.Centroid,
}

_SHAPES = {
    "gaussian": TransferFunctionShape.Gaussian,
    "top_hat": TransferFunctionShape.TopHat,
    "tophat": TransferFunctionShape.TopHat,
}

_ALGORITHMS = {
    "bobyqa": OptimizationAlgorithm.BoundedBobyqa,
    "bounded_bobyqa": OptimizationAlgorithm.BoundedBobyqa,
    "nelder_mead": OptimizationAlgorithm.NelderMead,
    "nelder_mead_simplex": OptimizationAlgorithm.NelderMead,
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


def _coerce_shape(shape: str | TransferFunctionShape) -> TransferFunctionShape:
    if isinstance(shape, TransferFunctionShape):
        return shape
    key = str(shape).strip().lower().replace("-", "_")
    if key in _SHAPES:
        return _SHAPES[key]
    raise ValueError(f"shape must be 'gaussian' or 'top_hat', got {shape!r}")


def _coerce_algorithm(
    algorithm: str | OptimizationAlgorithm,
) -> OptimizationAlgorithm:
    if isinstance(algorithm, OptimizationAlgorithm):
        return algorithm
    key = str(algorithm).strip().lower().replace("-", "_")
    if key in _ALGORITHMS:
        return _ALGORITHMS[key]
    raise ValueError(
        f"algorithm must be 'bobyqa' or 'nelder_mead', got {algorithm!r}"
    )


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



@dataclass
class TransferFunctionResult:
    """Result of a transfer-function fit.

    Attributes
    ----------
    converged:
        Whether the optimizer converged within the tolerances.
    offset, amplitude:
        Best linear parameters (background level and integrated
        response amplitude); the linear parameters are profiled
        analytically during the fit.
    lag, width:
        Best non-linear parameters (transfer-function center and
        width).
    chi2:
        Minimum chi^2.
    evaluations:
        Objective evaluations used by the final local optimization.
    n_valid_points:
        Response points whose convolution window lies fully inside
        the continuum coverage at the best fit.
    message:
        Optimizer status message.
    """

    converged: bool
    offset: float
    amplitude: float
    lag: float
    width: float
    chi2: float
    evaluations: int
    n_valid_points: int
    message: str


def fit_transfer_function(
    time1: Sequence[float],
    value1: Sequence[float],
    error1: Sequence[float],
    time2: Sequence[float],
    value2: Sequence[float],
    error2: Sequence[float],
    *,
    shape: str | TransferFunctionShape = "gaussian",
    offset: float | None = None,
    amplitude: float | None = None,
    lag: float | None = None,
    width: float | None = None,
    bounds: Mapping[str, tuple[float, float]] | None = None,
    grid_step: float | None = None,
    lag_restarts: int = 0,
    algorithm: str | OptimizationAlgorithm = "bobyqa",
    max_evaluations: int = 10000,
    xtol_rel: float = 1e-6,
    ftol_rel: float = 1e-8,
) -> TransferFunctionResult:
    """Fit a parametric transfer function Psi(tau) by chi^2 minimization.

    Model:

        R(t) = offset + amplitude * (C * Psi)(t)

    where C is the continuum light curve and Psi is a normalized
    Gaussian or top-hat centered at ``lag`` with ``width``. The linear
    parameters (offset, amplitude) are profiled analytically; only
    (lag, width) are searched by the bounded optimizer (default:
    NLopt BOBYQA). A lag scan (``lag_restarts``) reduces the risk of
    landing in a secondary minimum of the multimodal chi^2 vs lag
    surface.

    Parameters
    ----------
    time1, value1, error1:
        Continuum (driving) light curve.
    time2, value2, error2:
        Response light curve.
    shape:
        ``"gaussian"`` (default) or ``"top_hat"``.
    offset, amplitude, lag, width:
        Optional initial values; ``None`` uses data-driven defaults.
    bounds:
        Optional dict mapping parameter names (``"offset"``,
        ``"amplitude"``, ``"lag"``, ``"width"``) to ``(lower, upper)``.
        A parameter with ``lower == upper`` is pinned.
    grid_step:
        Continuum interpolation step for the convolution integral;
        ``None`` selects an automatic step.
    lag_restarts:
        Lag scan size: ``0`` = automatic scan, ``1`` = single start,
        ``n > 1`` = explicit scan size.
    algorithm:
        ``"bobyqa"`` (default) or ``"nelder_mead"``.
    """
    shape_enum = _coerce_shape(shape)
    algorithm_enum = _coerce_algorithm(algorithm)
    step = 0.0 if grid_step is None else float(grid_step)

    defaults = _default_transfer_function_parameters(
        time1, value1, error1, time2, value2, error2,
        shape_enum, step,
    )

    names = ("offset", "amplitude", "lag", "width")
    initial = {"offset": offset, "amplitude": amplitude, "lag": lag, "width": width}
    bounds = {} if bounds is None else bounds

    parameters: list[FitParameter] = list(defaults)

    for i, name in enumerate(names):
        if initial[name] is not None:
            parameters[i].value = float(initial[name])

        if name in bounds:
            lo, hi = bounds[name]
            parameters[i].lower = float(lo)
            parameters[i].upper = float(hi)

        if parameters[i].lower > parameters[i].upper:
            raise ValueError(
                f"lower bound exceeds upper bound for {name!r}"
            )

    raw = _fit_transfer_function(
        time1, value1, error1, time2, value2, error2,
        parameters,
        shape_enum,
        step,
        int(lag_restarts),
        algorithm_enum,
        float(xtol_rel),
        float(ftol_rel),
        int(max_evaluations),
    )

    return TransferFunctionResult(
        converged=raw.converged,
        offset=raw.offset,
        amplitude=raw.amplitude,
        lag=raw.lag,
        width=raw.width,
        chi2=raw.chi2,
        evaluations=raw.evaluations,
        n_valid_points=raw.n_valid_points,
        message=raw.message,
    )


def transfer_function_curve(
    taus: Sequence[float],
    *,
    shape: str | TransferFunctionShape = "gaussian",
    lag: float,
    width: float,
) -> np.ndarray:
    """Normalized transfer function Psi(tau) on a lag grid."""
    return np.asarray(
        _transfer_function_curve(
            list(taus),
            _coerce_shape(shape),
            float(lag),
            float(width),
        )
    )


def model_response(
    time1: Sequence[float],
    value1: Sequence[float],
    error1: Sequence[float],
    times: Sequence[float],
    *,
    shape: str | TransferFunctionShape = "gaussian",
    offset: float = 0.0,
    amplitude: float = 1.0,
    lag: float,
    width: float,
    grid_step: float | None = None,
) -> np.ndarray:
    """Model response R(t) = offset + amplitude * (C * Psi)(t)."""
    return np.asarray(
        _transfer_function_model_response(
            time1, value1, error1,
            list(times),
            _coerce_shape(shape),
            float(offset),
            float(amplitude),
            float(lag),
            float(width),
            0.0 if grid_step is None else float(grid_step),
        )
    )
