import datetime
from typing import ClassVar

from .core import UnitedStates

class Indiana(UnitedStates):
    label_washington_birthday_december: ClassVar[str]
    def get_washington_birthday_december(self, year: int) -> tuple[datetime.date, str]: ...
    def get_primary_election_day(self, year: int) -> tuple[datetime.date, str]: ...
