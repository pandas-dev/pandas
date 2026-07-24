from typing import Literal, overload

from pymeeus.Angle import Angle
from pymeeus.Epoch import Epoch

class JupiterMoons:
    @staticmethod
    def jupiter_system_angles(epoch: Epoch) -> tuple[float, float]: ...
    @staticmethod
    def rectangular_positions_jovian_equatorial(
        epoch: Epoch, tofk5: bool = True, solar: bool = False, do_correction: bool = True
    ) -> tuple[
        tuple[float, float, float], tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]
    ]: ...
    @overload
    @staticmethod
    def apparent_rectangular_coordinates(
        epoch: Epoch,
        X: float,
        Y: float,
        Z: float,
        OMEGA: float,
        psi: float,
        i: float,
        lambda_0: float,
        beta_0: float,
        D: float = 0,
        isFictional: Literal[True] = ...,
    ) -> float: ...
    @overload
    @staticmethod
    def apparent_rectangular_coordinates(
        epoch: Epoch,
        X: float,
        Y: float,
        Z: float,
        OMEGA: float,
        psi: float,
        i: float,
        lambda_0: float,
        beta_0: float,
        D: float = 0,
        isFictional: Literal[False] | None = False,
    ) -> tuple[float, float, float]: ...
    @staticmethod
    def calculate_delta(
        epoch: Epoch,
    ) -> tuple[float, float, Angle, Angle, float] | tuple[float, float, Literal[0], Literal[0], Literal[0]]: ...
    @overload
    @staticmethod
    def correct_rectangular_positions(
        R: float, i_sat: int, DELTA: float, X_coordinate: list[float] | tuple[float, float, float]
    ) -> tuple[float, float, float]: ...
    @overload
    @staticmethod
    def correct_rectangular_positions(
        R: float, i_sat: int, DELTA: float, X_coordinate: float, Y_coordinate: float = 0, Z_coordinate: float = 0
    ) -> tuple[float, float, float]: ...
    @overload
    @staticmethod
    def check_phenomena(epoch: Epoch, check_all: Literal[True] = True, i_sat: int = 0) -> list[list[float]]: ...
    @overload
    @staticmethod
    def check_phenomena(epoch: Epoch, check_all: Literal[False] | None, i_sat: int = 0) -> tuple[float, float]: ...
    @staticmethod
    def is_phenomena(epoch: Epoch) -> list[list[bool]]: ...
    @staticmethod
    def check_coordinates(X: float, Y: float) -> float: ...
    @staticmethod
    def check_occultation(
        X: float = 0, Y: float = 0, Z: float = 0, epoch: Epoch | None = None, i_sat: int | None = None
    ) -> float: ...
    @staticmethod
    def check_eclipse(
        X_0: float = 0, Y_0: float = 0, Z_0: float = 0, epoch: Epoch | None = None, i_sat: int | None = None
    ) -> float: ...

def main() -> None: ...
