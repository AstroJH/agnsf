"""Result conversion utilities (user-facing Python layer)."""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from numpy.typing import NDArray

from ._agnsf import LagBins, SFResult

__all__ = ["SFArrays", "to_arrays"]


class SFArrays(NamedTuple):
    """Binned SF result as named NumPy arrays (one entry per lag bin).

    ``measurement_*`` and ``sampling_*`` are NaN when the corresponding
    uncertainty was not estimated.
    """

    lag: NDArray[np.float64]
    sf: NDArray[np.float64]
    count: NDArray[np.int64]
    valid: NDArray[np.bool_]
    measurement_lower: NDArray[np.float64]
    measurement_upper: NDArray[np.float64]
    within_lower: NDArray[np.float64]
    within_upper: NDArray[np.float64]
    sampling_lower: NDArray[np.float64]
    sampling_upper: NDArray[np.float64]


def to_arrays(result: SFResult, bins: LagBins) -> SFArrays:
    """Convert an SFResult into a named set of NumPy arrays.

    ``lag`` is the arithmetic mean of each bin's edges, ``sf`` the
    structure function, ``count`` the number of contributions, and
    ``valid`` marks bins with contributions and a finite ``sf``.
    ``measurement_*`` / ``within_*`` / ``sampling_*`` are NaN when the
    corresponding uncertainty was not estimated.
    """
    lag = np.array([
        np.mean(bins[i])
        for i in range(len(bins))
    ])

    sf = np.array([
        result[i].sf
        for i in range(len(result))
    ])

    count = np.array([
        result[i].count
        for i in range(len(result))
    ])

    valid = (count > 0) & np.isfinite(sf)

    measurement_lower = np.array([
        result[i].measurement.lower
        for i in range(len(result))
    ])

    measurement_upper = np.array([
        result[i].measurement.upper
        for i in range(len(result))
    ])

    within_lower = np.array([
        result[i].within.lower
        for i in range(len(result))
    ])

    within_upper = np.array([
        result[i].within.upper
        for i in range(len(result))
    ])

    sampling_lower = np.array([
        result[i].sampling.lower
        for i in range(len(result))
    ])

    sampling_upper = np.array([
        result[i].sampling.upper
        for i in range(len(result))
    ])

    return SFArrays(
        lag=lag,
        sf=sf,
        count=count,
        valid=valid,
        measurement_lower=measurement_lower,
        measurement_upper=measurement_upper,
        within_lower=within_lower,
        within_upper=within_upper,
        sampling_lower=sampling_lower,
        sampling_upper=sampling_upper,
    )
