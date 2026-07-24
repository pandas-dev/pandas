from typing import Final

from pymeeus.Angle import Angle
from pymeeus.Epoch import Epoch

PLUTO_ARGUMENT: Final[list[tuple[float, float, float]]]
PLUTO_LONGITUDE: Final[list[tuple[float, float]]]
PLUTO_LATITUDE: Final[list[tuple[float, float]]]
PLUTO_RADIUS_VECTOR: Final[list[tuple[float, float]]]

class Pluto:
    @staticmethod
    def geometric_heliocentric_position(epoch: Epoch) -> tuple[Angle, Angle, float]: ...
    @staticmethod
    def geocentric_position(epoch: Epoch) -> tuple[Angle, Angle]: ...

def main() -> None: ...
