import datetime

from .core import UnitedStates

class Illinois(UnitedStates): ...

class ChicagoIllinois(Illinois):
    def get_pulaski_day(self, year: int) -> tuple[datetime.date, str]: ...
