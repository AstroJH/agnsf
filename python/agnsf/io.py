"""File-based convenience functions (user-facing Python layer).

The source-based helpers accept a flexible ``source`` argument:

- for :func:`esf_from`, ``source`` may be a path-list file (``str``),
  a list of light-curve paths, or a list of ``(path, redshift)``
  tuples. When a path-list file is used, each line may carry an
  optional second column with the source redshift.
"""

from __future__ import annotations

from typing import Any, Sequence, Union

from ._agnsf import (
    EnsembleMethod,
    LagBins,
    LightCurve,
    SFMethod,
    SFResult,
    UncertaintyConfig,
    read_light_curve as _read_light_curve,
    read_path_list_with_redshift as _read_path_list_with_redshift,
    write_sf_result as _write_sf_result,
    write_table as _write_table,
)
from .analysis import (
    _check_redshift,
    _coerce_aggregation,
    _coerce_sf_method,
    _coerce_uncertainty,
    ensemble_sf,
    pooled_sf,
)

__all__ = [
    "read_light_curve",
    "read_light_curves",
    "sf_file",
    "esf_from",
    "write_table",
    "write_sf_result",
]

Source = Union[
    str,
    Sequence[str],
    Sequence[tuple[str, float]],
]


def read_light_curve(
    path: str,
    *,
    time: str = "time",
    value: str = "value",
    error: str = "error",
) -> LightCurve:
    """Read a light curve from a CSV or FITS file."""
    return _read_light_curve(path, time, value, error)


def _coerce_source_entries(source: Source) -> list[tuple[str, float]]:
    """Normalize a source into a list of ``(path, redshift)`` tuples."""
    if isinstance(source, str):
        entries = _read_path_list_with_redshift(source)

        for _, redshift in entries:
            _check_redshift(redshift)

        return entries

    if isinstance(source, (list, tuple)):
        entries: list[tuple[str, float]] = []

        for item in source:
            if isinstance(item, str):
                entries.append((item, 0.0))
            elif (
                isinstance(item, (tuple, list))
                and len(item) == 2
                and isinstance(item[0], str)
            ):
                redshift = float(item[1])
                _check_redshift(redshift)
                entries.append((item[0], redshift))
            else:
                raise TypeError(
                    "each source item must be a path (str) or a "
                    "(path, redshift) pair, got "
                    f"{item!r}"
                )

        return entries

    raise TypeError(
        "source must be a path-list file (str), a list of paths, or "
        "a list of (path, redshift) pairs"
    )


def read_light_curves(
    source: Source,
    *,
    time: str = "time",
    value: str = "value",
    error: str = "error",
) -> list[LightCurve]:
    """Read several observed-frame light curves from a source.
    """
    return [
        read_light_curve(path, time=time, value=value, error=error)
        for path, _redshift in _coerce_source_entries(source)
    ]


def sf_file(
    path: str,
    bins: LagBins,
    *,
    method: SFMethod | str = SFMethod.SecondOrder,
    uncertainty: UncertaintyConfig | str | bool | dict[str, Any] | None = None,
    redshift: float = 0.0,
    time: str = "time",
    value: str = "value",
    error: str = "error",
) -> SFResult:
    """Structure function of one light-curve file (CSV or FITS)."""
    curve = read_light_curve(path, time=time, value=value, error=error)

    from .analysis import sf

    return sf(
        curve.time,
        curve.value,
        curve.error,
        bins,
        method=method,
        uncertainty=uncertainty,
        redshift=redshift,
    )


def esf_from(
    source: Source,
    bins: LagBins,
    *,
    kind: str = "pooled",
    method: SFMethod | str = SFMethod.SecondOrder,
    aggregation: EnsembleMethod | str = EnsembleMethod.SqrtMeanSquared,
    uncertainty: UncertaintyConfig | str | bool | dict[str, Any] | None = None,
    time: str = "time",
    value: str = "value",
    error: str = "error",
) -> SFResult:
    """Ensemble structure function from a source (files or path list).

    Parameters
    ----------
    source:
        A path-list file (``str``), a list of light-curve paths, or a
        list of ``(path, redshift)`` pairs.
    kind:
        ``"pooled"`` (default) or ``"aggregated"``.
    method:
        SF estimator (for pooled, the pooled estimator; for aggregated,
        the per-curve estimator).
    aggregation:
        Aggregation for ``kind="aggregated"``:
        ``"sqrt_mean_squared"`` (default) or ``"mean_sf"``.
    """
    entries = _coerce_source_entries(source)

    curves = [
        read_light_curve(path, time=time, value=value, error=error)
        for path, _redshift in entries
    ]

    # Per-source redshifts are passed to the native kernels, which
    # convert lags to the rest frame per pair (no time-array copies).
    redshifts = [redshift for _path, redshift in entries]

    times = [curve.time for curve in curves]
    values = [curve.value for curve in curves]
    errors = [curve.error for curve in curves]

    kind_key = str(kind).strip().lower()

    if kind_key == "pooled":
        return pooled_sf(
            times,
            values,
            errors,
            bins,
            method=_coerce_sf_method(method),
            uncertainty=_coerce_uncertainty(uncertainty),
            redshift=redshifts,
        )

    if kind_key == "aggregated":
        return ensemble_sf(
            times,
            values,
            errors,
            bins,
            method=_coerce_aggregation(aggregation),
            sf_method=_coerce_sf_method(method),
            uncertainty=_coerce_uncertainty(uncertainty),
            redshift=redshifts,
        )

    raise ValueError(
        f"kind must be 'pooled' or 'aggregated', got {kind!r}"
    )


def write_table(
    path: str,
    headers: Sequence[str],
    columns: Sequence[Sequence[float]],
) -> None:
    """Write columns of numbers to a simple text file."""
    _write_table(path, headers, columns)


def write_sf_result(
    path: str,
    bins: LagBins,
    result: SFResult,
) -> None:
    """Write an SFResult to a text file (`# lag count sf_squared sf`)."""
    _write_sf_result(path, bins, result)
