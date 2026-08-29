from typing import Sequence


class SFMethod:
    SecondOrder: SFMethod
    SecondOrderNoNoise: SFMethod
    MeanAbsoluteDeviation: SFMethod
    MeanAbsoluteDeviationNoNoise: SFMethod


class EnsembleMethod:
    SqrtMeanSquared: EnsembleMethod
    MeanSf: EnsembleMethod


class UncertaintyMethod:
    Off: UncertaintyMethod
    Analytic: UncertaintyMethod
    MonteCarlo: UncertaintyMethod
    Jackknife: UncertaintyMethod
    Bootstrap: UncertaintyMethod


class Uncertainty:
    lower: float
    upper: float

    @property
    def estimated(self) -> bool: ...


class UncertaintyConfig:
    measurement: UncertaintyMethod
    within: UncertaintyMethod
    sampling: UncertaintyMethod
    n_bootstrap: int
    bootstrap_seed: int

    def __init__(
        self,
        measurement: UncertaintyMethod = UncertaintyMethod.Off,
        within: UncertaintyMethod = UncertaintyMethod.Off,
        sampling: UncertaintyMethod = UncertaintyMethod.Off,
        n_bootstrap: int = 100,
        bootstrap_seed: int = 0,
    ) -> None: ...


class SFBinResult:
    count: int
    sf_squared: float
    sf: float
    measurement: Uncertainty
    within: Uncertainty
    sampling: Uncertainty


class SFResult:
    @property
    def bins(self) -> list[SFBinResult]: ...

    def __len__(self) -> int: ...

    def __getitem__(self, index: int) -> SFBinResult: ...


class LagBins:
    class GridType:
        Custom: GridType
        Linear: GridType
        Logarithmic: GridType

    def __init__(
        self,
        edges: Sequence[float],
    ) -> None: ...

    @staticmethod
    def linear(
        min: float,
        max: float,
        step: float,
    ) -> LagBins: ...

    @staticmethod
    def logarithmic(
        min: float,
        max: float,
        step: float,
    ) -> LagBins: ...

    def __len__(self) -> int: ...

    def __getitem__(self, index: int) -> tuple[float, float]: ...

    @property
    def grid_type(self) -> GridType: ...

    @property
    def min(self) -> float: ...

    @property
    def max(self) -> float: ...

    @property
    def size(self) -> int: ...

    @property
    def edges(self) -> list[float]: ...

    def contains(
        self,
        lag: float,
    ) -> bool: ...

    def index(
        self,
        lag: float,
    ) -> int: ...


class LightCurve:
    def __init__(
        self,
        time: Sequence[float],
        value: Sequence[float],
        error: Sequence[float],
    ) -> None: ...

    @property
    def size(self) -> int: ...

    @property
    def time(self) -> list[float]: ...

    @property
    def value(self) -> list[float]: ...

    @property
    def error(self) -> list[float]: ...


class LightCurveView:
    @property
    def size(self) -> int: ...

    @property
    def time_address(self) -> int: ...

    @property
    def value_address(self) -> int: ...

    @property
    def error_address(self) -> int: ...


def inspect(
    time: object,
    value: object,
    error: object,
) -> dict[str, object]:
    """Report whether the inputs would be copied by the sf() binding."""


def sf(
    time: Sequence[float],
    value: Sequence[float],
    error: Sequence[float],
    bins: LagBins,
    method: SFMethod = SFMethod.SecondOrder,
    uncertainty: UncertaintyConfig = UncertaintyConfig(),
    redshift: float = 0.0,
) -> SFResult: ...


def pooled_sf(
    time: Sequence[Sequence[float]],
    value: Sequence[Sequence[float]],
    error: Sequence[Sequence[float]],
    bins: LagBins,
    method: SFMethod = SFMethod.SecondOrder,
    uncertainty: UncertaintyConfig = UncertaintyConfig(),
    redshift: float | Sequence[float] = 0.0,
) -> SFResult: ...


def ensemble_sf(
    time: Sequence[Sequence[float]],
    value: Sequence[Sequence[float]],
    error: Sequence[Sequence[float]],
    bins: LagBins,
    method: EnsembleMethod = EnsembleMethod.SqrtMeanSquared,
    sf_method: SFMethod = SFMethod.SecondOrder,
    uncertainty: UncertaintyConfig = UncertaintyConfig(),
    redshift: float | Sequence[float] = 0.0,
) -> SFResult: ...


def read_light_curve(
    path: str,
    time: str = "time",
    value: str = "value",
    error: str = "error",
) -> LightCurve: ...


def read_path_list(
    path: str,
) -> list[str]: ...


def read_path_list_with_redshift(
    path: str,
) -> list[tuple[str, float]]: ...


def write_table(
    path: str,
    headers: Sequence[str],
    columns: Sequence[Sequence[float]],
) -> None: ...


def sf_from_file(
    path: str,
    bins: LagBins,
    method: SFMethod = SFMethod.SecondOrder,
    time: str = "time",
    value: str = "value",
    error: str = "error",
    uncertainty: UncertaintyConfig = UncertaintyConfig(),
    redshift: float = 0.0,
) -> SFResult: ...


def pooled_sf_from_files(
    paths: Sequence[str],
    bins: LagBins,
    method: SFMethod = SFMethod.SecondOrder,
    time: str = "time",
    value: str = "value",
    error: str = "error",
    uncertainty: UncertaintyConfig = UncertaintyConfig(),
    redshift: float | Sequence[float] = 0.0,
) -> SFResult: ...


def pooled_sf_from_path_list(
    path_list_file: str,
    bins: LagBins,
    method: SFMethod = SFMethod.SecondOrder,
    time: str = "time",
    value: str = "value",
    error: str = "error",
    uncertainty: UncertaintyConfig = UncertaintyConfig(),
) -> SFResult: ...


def ensemble_sf_from_files(
    paths: Sequence[str],
    bins: LagBins,
    method: EnsembleMethod = EnsembleMethod.SqrtMeanSquared,
    sf_method: SFMethod = SFMethod.SecondOrder,
    time: str = "time",
    value: str = "value",
    error: str = "error",
    uncertainty: UncertaintyConfig = UncertaintyConfig(),
    redshift: float | Sequence[float] = 0.0,
) -> SFResult: ...


def ensemble_sf_from_path_list(
    path_list_file: str,
    bins: LagBins,
    method: EnsembleMethod = EnsembleMethod.SqrtMeanSquared,
    sf_method: SFMethod = SFMethod.SecondOrder,
    time: str = "time",
    value: str = "value",
    error: str = "error",
    uncertainty: UncertaintyConfig = UncertaintyConfig(),
) -> SFResult: ...


def write_sf_result(
    path: str,
    bins: LagBins,
    result: SFResult,
) -> None: ...


class CrossCorrelationMethod:
    Dcf: CrossCorrelationMethod
    Iccf: CrossCorrelationMethod


class LagEstimate:
    Peak: LagEstimate
    Centroid: LagEstimate


class LagResult:
    lag_peak: float
    lag_centroid: float
    peak_value: float
    tau: list[float]
    ccf: list[float]
    count: list[int]


def cross_correlate(
    time1: Sequence[float],
    value1: Sequence[float],
    error1: Sequence[float],
    time2: Sequence[float],
    value2: Sequence[float],
    error2: Sequence[float],
    grid_min: float,
    grid_max: float,
    grid_step: float,
    method: CrossCorrelationMethod = CrossCorrelationMethod.Dcf,
    dcf_bin_width: float = 1.0,
    centroid_threshold: float = 0.8,
    min_overlap: int = 3,
) -> LagResult: ...


def lag_uncertainty(
    time1: Sequence[float],
    value1: Sequence[float],
    error1: Sequence[float],
    time2: Sequence[float],
    value2: Sequence[float],
    error2: Sequence[float],
    grid_min: float,
    grid_max: float,
    grid_step: float,
    estimate: LagEstimate = LagEstimate.Peak,
    method: CrossCorrelationMethod = CrossCorrelationMethod.Dcf,
    dcf_bin_width: float = 1.0,
    centroid_threshold: float = 0.8,
    min_overlap: int = 3,
    n_realizations: int = 1000,
    seed: int = 0,
    flux_randomization: bool = True,
    random_subset: bool = True,
) -> Uncertainty: ...
