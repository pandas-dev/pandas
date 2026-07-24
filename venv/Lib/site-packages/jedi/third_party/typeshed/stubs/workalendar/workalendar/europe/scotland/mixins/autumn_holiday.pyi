import datetime
from typing import ClassVar

class AutumHoliday:
    include_autumn_holiday: ClassVar[bool]
    autumn_holiday_label: ClassVar[str]

class AutumnHolidayLastMondaySeptember(AutumHoliday):
    def get_autumn_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class AutumnHolidayFirstMondayOctober(AutumHoliday):
    def get_autumn_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class AutumnHolidaySecondMondayOctober(AutumHoliday):
    def get_autumn_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class AutumnHolidayThirdMondayOctober(AutumHoliday):
    def get_autumn_holiday(self, year: int) -> tuple[datetime.date, str]: ...
