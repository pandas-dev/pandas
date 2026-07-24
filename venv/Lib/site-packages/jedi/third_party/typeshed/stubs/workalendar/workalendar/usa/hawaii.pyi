import datetime

from .core import UnitedStates

class Hawaii(UnitedStates):
    def get_statehood_day(self, year: int) -> tuple[datetime.date, str]: ...
