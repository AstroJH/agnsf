from typing import Sequence


class SFBinResult:
    count: int
    sf_squared: float
    sf: float


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


def sf(
    time: Sequence[float],
    value: Sequence[float],
    error: Sequence[float],
    bins: LagBins,
) -> SFResult: ...


def pooled_sf(
    time: Sequence[Sequence[float]],
    value: Sequence[Sequence[float]],
    error: Sequence[Sequence[float]],
    bins: LagBins,
) -> SFResult:
    ...
