"""Core structure-function computations (user-facing Python layer).

The in-memory functions accept either enum members or plain strings
for the estimator / aggregation / uncertainty settings, so users do
not need to import the enum classes:

    agnsf.sf(t, v, e, bins, method="mean_absolute_deviation",
             uncertainty="bootstrap")

Redshift: when ``redshift`` is given, lags are converted to the
source rest frame (dt_rest = dt_obs / (1 + z)) inside the native
kernel, so no copy of the time array is made.
"""

from __future__ import annotations

from numbers import Real
from typing import Any, Mapping, Sequence, Union

from ._agnsf import (
    EnsembleMethod,
    LagBins,
    SFMethod,
    SFResult,
    UncertaintyConfig,
    UncertaintyMethod,
    ensemble_sf as _ensemble_sf,
    pooled_sf as _pooled_sf,
    sf as _sf,
)

__all__ = ["sf", "pooled_sf", "ensemble_sf"]

#: accepted spellings of SFMethod (underscores, dashes and spaces
#: are normalized).
_SF_METHODS: dict[str, SFMethod] = {
    "second_order": SFMethod.SecondOrder,
    "second_order_no_noise": SFMethod.SecondOrderNoNoise,
    "mean_absolute_deviation": SFMethod.MeanAbsoluteDeviation,
    "mean_absolute_deviation_no_noise": SFMethod.MeanAbsoluteDeviationNoNoise,
}

_AGGREGATIONS: dict[str, EnsembleMethod] = {
    "sqrt_mean_squared": EnsembleMethod.SqrtMeanSquared,
    "mean_sf": EnsembleMethod.MeanSf,
}

_UNCERTAINTY_METHODS: dict[str, UncertaintyMethod] = {
    "off": UncertaintyMethod.Off,
    "analytic": UncertaintyMethod.Analytic,
    "monte_carlo": UncertaintyMethod.MonteCarlo,
    "jackknife": UncertaintyMethod.Jackknife,
    "bootstrap": UncertaintyMethod.Bootstrap,
}


def _normalize(text: str) -> str:
    return text.strip().lower().replace("-", "_").replace(" ", "_")


def _coerce_sf_method(method: SFMethod | str) -> SFMethod:
    if isinstance(method, SFMethod):
        return method

    if isinstance(method, str):
        key = _normalize(method)

        if key in _SF_METHODS:
            return _SF_METHODS[key]

    raise ValueError(
        f"invalid SF method {method!r}; expected one of "
        f"{sorted(_SF_METHODS)}"
    )


def _coerce_aggregation(aggregation: EnsembleMethod | str) -> EnsembleMethod:
    if isinstance(aggregation, EnsembleMethod):
        return aggregation

    if isinstance(aggregation, str):
        key = _normalize(aggregation)

        if key in _AGGREGATIONS:
            return _AGGREGATIONS[key]

    raise ValueError(
        f"invalid aggregation {aggregation!r}; expected one of "
        f"{sorted(_AGGREGATIONS)}"
    )


def _coerce_uncertainty_method(
    value: UncertaintyMethod | str | bool,
) -> UncertaintyMethod:
    if isinstance(value, UncertaintyMethod):
        return value

    if isinstance(value, bool):
        return UncertaintyMethod.Analytic if value else UncertaintyMethod.Off

    if isinstance(value, str):
        key = _normalize(value)

        if key in _UNCERTAINTY_METHODS:
            return _UNCERTAINTY_METHODS[key]

    raise ValueError(
        f"invalid uncertainty method {value!r}; expected one of "
        f"{sorted(_UNCERTAINTY_METHODS)}"
    )


def _coerce_uncertainty(
    value: UncertaintyConfig | str | bool | Mapping[str, Any] | None,
) -> UncertaintyConfig:
    """Coerce a user-supplied uncertainty setting into an UncertaintyConfig.

    Accepted forms:

    - ``None`` / ``False``: off (default)
    - ``True``: measurement = Analytic
    - ``str``: "analytic" / "measurement" -> measurement = Analytic;
      "jackknife" / "sampling" -> sampling = Jackknife;
      "bootstrap" -> sampling = Bootstrap
    - ``dict``: keys ``measurement``, ``sampling`` (string/enum/bool),
      ``n_bootstrap``, ``bootstrap_seed``
    - :class:`UncertaintyConfig`: passed through unchanged
    """
    if value is None or value is False:
        return UncertaintyConfig()

    if value is True:
        return UncertaintyConfig(measurement=UncertaintyMethod.Analytic)

    if isinstance(value, UncertaintyConfig):
        return value

    if isinstance(value, str):
        key = _normalize(value)

        if key in ("analytic", "measurement"):
            return UncertaintyConfig(measurement=UncertaintyMethod.Analytic)

        if key in ("monte_carlo", "mc"):
            return UncertaintyConfig(measurement=UncertaintyMethod.MonteCarlo)

        if key == "within":
            return UncertaintyConfig(within=UncertaintyMethod.Analytic)

        if key in ("jackknife", "sampling"):
            return UncertaintyConfig(sampling=UncertaintyMethod.Jackknife)

        if key == "bootstrap":
            return UncertaintyConfig(sampling=UncertaintyMethod.Bootstrap)

        raise ValueError(
            f"invalid uncertainty shorthand {value!r}; expected "
            "'analytic' | 'measurement' | 'monte_carlo' | 'within' | "
            "'jackknife' | 'sampling' | 'bootstrap'"
        )

    if isinstance(value, Mapping):
        config = UncertaintyConfig()

        if "measurement" in value:
            config.measurement = _coerce_uncertainty_method(
                value["measurement"]
            )

        if "within" in value:
            config.within = _coerce_uncertainty_method(value["within"])

        if "sampling" in value:
            config.sampling = _coerce_uncertainty_method(value["sampling"])

        if "n_bootstrap" in value:
            config.n_bootstrap = int(value["n_bootstrap"])

        if "bootstrap_seed" in value:
            config.bootstrap_seed = int(value["bootstrap_seed"])

        return config

    raise TypeError(
        f"cannot interpret uncertainty setting {value!r}"
    )


def _check_redshift(redshift: float) -> None:
    if redshift <= -1.0:
        raise ValueError(
            f"redshift must be > -1, got {redshift}"
        )


def _is_scalar_redshift(value: Any) -> bool:
    """True for a single redshift number (including numpy scalars)."""
    return isinstance(value, Real) and not isinstance(value, bool)


def _resolve_redshifts(
    redshift: float | Sequence[float] | None,
    n_curves: int,
) -> list[float] | None:
    """Normalize the ``redshift`` argument to per-curve values.

    A scalar means "the same redshift for every light curve"; a
    sequence must contain one value per light curve. Returns ``None``
    when no correction is needed (all redshifts are zero).
    """
    if redshift is None:
        return None

    if _is_scalar_redshift(redshift):
        _check_redshift(float(redshift))

        if redshift == 0.0:
            return None

        return [float(redshift)] * n_curves

    if isinstance(redshift, str):
        raise TypeError(
            "redshift must be a number or a sequence of numbers, "
            f"got {redshift!r}"
        )

    values = [float(value) for value in redshift]

    if len(values) != n_curves:
        raise ValueError(
            f"redshift must provide one value per light curve "
            f"({n_curves} curves), got {len(values)}"
        )

    for value in values:
        _check_redshift(value)

    if all(value == 0.0 for value in values):
        return None

    return values



def sf(
    time: Sequence[float],
    value: Sequence[float],
    error: Sequence[float],
    bins: LagBins,
    *,
    method: SFMethod | str = SFMethod.SecondOrder,
    uncertainty: UncertaintyConfig | str | bool | Mapping[str, Any] | None = None,
    redshift: float = 0.0,
) -> SFResult:
    """Structure function of a single light curve.

    Parameters
    ----------
    time, value, error:
        Observed time, flux, and per-point uncertainty arrays.
    bins:
        Lag bins (in rest-frame units when ``redshift`` is given).
    method:
        SF estimator, as an :class:`SFMethod` member or a string
        (e.g. ``"mean_absolute_deviation"``).
    uncertainty:
        See :func:`_coerce_uncertainty`; ``None`` disables estimation.
    redshift:
        Source redshift ``z``; lags are converted to the rest frame.
    """
    if not _is_scalar_redshift(redshift):
        raise TypeError(
            "redshift for a single light curve must be a number, "
            f"got {redshift!r}"
        )

    _check_redshift(float(redshift))

    return _sf(
        time,
        value,
        error,
        bins,
        _coerce_sf_method(method),
        _coerce_uncertainty(uncertainty),
        float(redshift),
    )


def pooled_sf(
    times: Sequence[Sequence[float]],
    values: Sequence[Sequence[float]],
    errors: Sequence[Sequence[float]],
    bins: LagBins,
    *,
    method: SFMethod | str = SFMethod.SecondOrder,
    uncertainty: UncertaintyConfig | str | bool | Mapping[str, Any] | None = None,
    redshift: float | Sequence[float] | None = 0.0,
) -> SFResult:
    """Pooled ensemble structure function over several light curves.

    ``redshift`` may be a single number (the same redshift for every
    light curve) or a sequence with one value per light curve
    (per-source rest-frame correction).
    """
    if _is_scalar_redshift(redshift):
        _check_redshift(float(redshift))

        return _pooled_sf(
            times,
            values,
            errors,
            bins,
            _coerce_sf_method(method),
            _coerce_uncertainty(uncertainty),
            float(redshift),
        )

    # Per-curve redshifts are passed to the native kernel, which
    # applies each curve's z per pair (no rest-frame time arrays).
    redshifts = _resolve_redshifts(redshift, len(times))

    return _pooled_sf(
        times,
        values,
        errors,
        bins,
        _coerce_sf_method(method),
        _coerce_uncertainty(uncertainty),
        redshifts if redshifts is not None else 0.0,
    )


def ensemble_sf(
    times: Sequence[Sequence[float]],
    values: Sequence[Sequence[float]],
    errors: Sequence[Sequence[float]],
    bins: LagBins,
    *,
    method: EnsembleMethod | str = EnsembleMethod.SqrtMeanSquared,
    sf_method: SFMethod | str = SFMethod.SecondOrder,
    uncertainty: UncertaintyConfig | str | bool | Mapping[str, Any] | None = None,
    redshift: float | Sequence[float] | None = 0.0,
) -> SFResult:
    """Aggregated ensemble structure function over several light curves.

    ``redshift`` may be a single number (the same redshift for every
    light curve) or a sequence with one value per light curve
    (per-source rest-frame correction).
    """
    if _is_scalar_redshift(redshift):
        _check_redshift(float(redshift))

        return _ensemble_sf(
            times,
            values,
            errors,
            bins,
            _coerce_aggregation(method),
            _coerce_sf_method(sf_method),
            _coerce_uncertainty(uncertainty),
            float(redshift),
        )

    # Per-curve redshifts are passed to the native kernel.
    redshifts = _resolve_redshifts(redshift, len(times))

    return _ensemble_sf(
        times,
        values,
        errors,
        bins,
        _coerce_aggregation(method),
        _coerce_sf_method(sf_method),
        _coerce_uncertainty(uncertainty),
        redshifts if redshifts is not None else 0.0,
    )
