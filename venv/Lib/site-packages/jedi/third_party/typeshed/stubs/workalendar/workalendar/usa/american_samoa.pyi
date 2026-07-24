import datetime

from .core import UnitedStates

class AmericanSamoa(UnitedStates):
    def get_flag_day(self, year: int) -> tuple[datetime.date, str]: ...
