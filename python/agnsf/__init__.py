from ._agnsf import (
    EnsembleMethod,
    SFMethod,
    UncertaintyMethod,
    UncertaintyConfig,
    SFUncertainty,
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
    "UncertaintyMethod",
    "UncertaintyConfig",
    "SFUncertainty",
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


def to_numpy(
    result: SFBinResult,
    bins: LagBins,
    uncertainty: bool = False,
):
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

    if not uncertainty:
        return lag, sf, count, valid

    measurement_lower = np.array([
        result[i].measurement.lower
        for i in range(len(result))
    ])

    measurement_upper = np.array([
        result[i].measurement.upper
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

    return (
        lag,
        sf,
        count,
        valid,
        measurement_lower,
        measurement_upper,
        sampling_lower,
        sampling_upper,
    )

