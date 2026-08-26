from ._agnsf import (
    EnsembleMethod,
    SFMethod,
    LagBins,
    LightCurve,
    SFBinResult,
    SFResult,
    sf,
    pooled_sf,
    ensemble_sf,
    read_light_curve,
    read_path_list,
    write_table,
    sf_from_file,
    pooled_sf_from_files,
    pooled_sf_from_path_list,
    ensemble_sf_from_files,
    ensemble_sf_from_path_list,
    write_sf_result
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
    "read_light_curve",
    "read_path_list",
    "write_table",
    "sf_from_file",
    "pooled_sf_from_files",
    "pooled_sf_from_path_list",
    "ensemble_sf_from_files",
    "ensemble_sf_from_path_list",
    "write_sf_result",
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
