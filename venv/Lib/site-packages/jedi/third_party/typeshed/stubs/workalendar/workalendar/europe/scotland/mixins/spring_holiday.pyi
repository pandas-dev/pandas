import datetime
from typing import ClassVar

class SpringHoliday:
    include_spring_holiday: ClassVar[bool]

class SpringHolidayFirstMondayApril(SpringHoliday):
    def get_spring_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class SpringHolidaySecondMondayApril(SpringHoliday):
    def get_spring_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class SpringHolidayTuesdayAfterFirstMondayMay(SpringHoliday):
    def get_spring_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class SpringHolidayLastMondayMay(SpringHoliday):
    def get_spring_holiday(self, year: int) -> tuple[datetime.date, str]: ...

class SpringHolidayFirstMondayJune(SpringHoliday):
    def get_spring_holiday(self, year: int) -> tuple[datetime.date, str]: ...
