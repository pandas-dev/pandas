import datetime

from .core import UnitedStates

class NorthCarolina(UnitedStates):
    def get_christmas_shifts(self, year: int) -> list[tuple[datetime.date, str]]: ...
