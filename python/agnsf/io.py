"""File-based convenience functions (user-facing Python layer).

The source-based helpers accept a flexible ``source`` argument:

- for :func:`esf_from`, ``source`` may be a path-list file (``str``),
  a list of light-curve paths, or a list of ``(path, redshift)``
  tuples. When a path-list file is used, each line may carry an
  optional second column with the source redshift.
"""

from __future__ import annotations

from typing import Any, Sequence, Union, Literal

from ._agnsf import (
    EnsembleMethod,
    LagBins,
    LightCurve,
    SFMethod,
    SFResult,
    UncertaintyConfig,
    ensemble_sf_from_files as _ensemble_sf_from_files,
    ensemble_sf_from_path_list as _ensemble_sf_from_path_list,
    pooled_sf_from_files as _pooled_sf_from_files,
    pooled_sf_from_path_list as _pooled_sf_from_path_list,
    read_light_curve as _read_light_curve,
    read_path_list_with_redshift as _read_path_list_with_redshift,
    sf_from_file as _sf_from_file,
    write_sf_result as _write_sf_result,
    write_table as _write_table,
)
from .analysis import (
    _check_redshift,
    _coerce_aggregation,
    _coerce_sf_method,
    _coerce_uncertainty,
    _is_scalar_redshift,
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
    """Structure function of one light-curve file (CSV or FITS).

    The file is read by the native layer; the ``time`` / ``value`` /
    ``error`` arguments select the column names.
    """
    if not _is_scalar_redshift(redshift):
        raise TypeError(
            "redshift for a single light curve must be a number, "
            f"got {redshift!r}"
        )

    _check_redshift(float(redshift))

    return _sf_from_file(
        path,
        bins,
        method=_coerce_sf_method(method),
        time=time,
        value=value,
        error=error,
        uncertainty=_coerce_uncertainty(uncertainty),
        redshift=float(redshift),
    )


def esf_from(
    source: Source,
    bins: LagBins,
    *,
    kind: Literal["pooled", "aggregated"] = "pooled",
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
    kind_key = str(kind).strip().lower()

    if kind_key not in ("pooled", "aggregated"):
        raise ValueError(
            f"kind must be 'pooled' or 'aggregated', got {kind!r}"
        )

    uncertainty_config = _coerce_uncertainty(uncertainty)

    if isinstance(source, str):
        # Path-list file: the native layer reads the paths and the
        # optional per-curve redshift column.
        if kind_key == "pooled":
            return _pooled_sf_from_path_list(
                source,
                bins,
                method=_coerce_sf_method(method),
                time=time,
                value=value,
                error=error,
                uncertainty=uncertainty_config,
            )

        return _ensemble_sf_from_path_list(
            source,
            bins,
            method=_coerce_aggregation(aggregation),
            sf_method=_coerce_sf_method(method),
            time=time,
            value=value,
            error=error,
            uncertainty=uncertainty_config,
        )

    # A list of paths and/or (path, redshift) pairs. The native layer
    # reads the files; per-curve redshifts are passed as metadata.
    entries = _coerce_source_entries(source)

    paths = [path for path, _redshift in entries]
    redshifts = [redshift for _path, redshift in entries]

    # Keep the z = 0 fast path: a scalar applies to every light curve.
    redshift_arg: float | list[float] = (
        0.0 if all(z == 0.0 for z in redshifts) else redshifts
    )

    if kind_key == "pooled":
        return _pooled_sf_from_files(
            paths,
            bins,
            method=_coerce_sf_method(method),
            time=time,
            value=value,
            error=error,
            uncertainty=uncertainty_config,
            redshift=redshift_arg,
        )

    return _ensemble_sf_from_files(
        paths,
        bins,
        method=_coerce_aggregation(aggregation),
        sf_method=_coerce_sf_method(method),
        time=time,
        value=value,
        error=error,
        uncertainty=uncertainty_config,
        redshift=redshift_arg,
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
