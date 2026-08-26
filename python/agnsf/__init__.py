from ._agnsf import (
    EnsembleMethod,
    SFMethod,
    LagBins,
    LightCurve,
    SFBinResult,
    SFResult,
    sf,
    pooled_sf,
    ensemble_sf
)

__all__ = [
    "EnsembleMethod",
    "SFMethod",
    "LagBins",
    "LightCurve",
    "SFBinResult",
    "SFResult",
    "sf",
    "pooled_sf",
    "ensemble_sf",
    "to_numpy"
]

import numpy as np


def to_numpy(result: SFBinResult, bins: LagBins):
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

    valid = (
        (count > 0)
        & np.isfinite(sf)
    )

    return lag, sf, count, valid
