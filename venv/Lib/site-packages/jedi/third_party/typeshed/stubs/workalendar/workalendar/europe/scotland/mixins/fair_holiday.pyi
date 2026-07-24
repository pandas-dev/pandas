import datetime
from typing import ClassVar

class FairHoliday:
    include_fair_holiday: ClassVar[bool]
    fair_holiday_label: ClassVar[str]

class FairHolidayLastMondayJune(FairHoliday):
    def get_fair_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class FairHolidayFirstMondayJuly(FairHoliday):
    def get_fair_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class FairHolidaySecondMondayJuly(FairHoliday):
    def get_fair_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class FairHolidayThirdMondayJuly(FairHoliday):
    def get_fair_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class FairHolidayLastMondayJuly(FairHoliday):
    def get_fair_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class FairHolidayFourthFridayJuly(FairHoliday):
    def get_fair_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class FairHolidayFirstMondayAugust(FairHoliday):
    def get_fair_holiday(self, year: int) -> tuple[datetime.date, str]: ...
