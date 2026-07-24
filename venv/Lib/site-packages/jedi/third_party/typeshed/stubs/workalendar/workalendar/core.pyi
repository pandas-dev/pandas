import datetime
from _typeshed import Incomplete, StrPath
from collections.abc import Generator, Iterable
from typing import ClassVar, Final, Literal, TypeVar, overload

_D = TypeVar("_D", datetime.date, datetime.datetime)
_DT = TypeVar("_DT", datetime.date, datetime.datetime, datetime.timedelta)

MON: Final = 0
TUE: Final = 1
WED: Final = 2
THU: Final = 3
FRI: Final = 4
SAT: Final = 5
SUN: Final = 6
ISO_MON: Final = 1
ISO_TUE: Final = 2
ISO_WED: Final = 3
ISO_THU: Final = 4
ISO_FRI: Final = 5
ISO_SAT: Final = 6
ISO_SUN: Final = 7

@overload
def cleaned_date(day: _D, keep_datetime: Literal[True]) -> _D: ...
@overload
def cleaned_date(day: datetime.date | datetime.datetime, keep_datetime: Literal[False] | None = False) -> datetime.date: ...
def daterange(start: _DT, end: _DT) -> Generator[_DT]: ...

class ChristianMixin:
    EASTER_METHOD: ClassVar[Literal[1, 2, 3] | None]
    include_epiphany: ClassVar[bool]
    include_clean_monday: ClassVar[bool]
    include_annunciation: ClassVar[bool]
    include_fat_tuesday: ClassVar[bool]
    fat_tuesday_label: ClassVar[str | None]
    include_ash_wednesday: ClassVar[bool]
    ash_wednesday_label: ClassVar[str]
    include_palm_sunday: ClassVar[bool]
    include_holy_thursday: ClassVar[bool]
    holy_thursday_label: ClassVar[str]
    include_good_friday: ClassVar[bool]
    good_friday_label: ClassVar[str]
    include_easter_monday: ClassVar[bool]
    include_easter_saturday: ClassVar[bool]
    easter_saturday_label: ClassVar[str]
    include_easter_sunday: ClassVar[bool]
    include_all_saints: ClassVar[bool]
    include_immaculate_conception: ClassVar[bool]
    immaculate_conception_label: ClassVar[str]
    include_christmas: ClassVar[bool]
    christmas_day_label: ClassVar[str]
    include_christmas_eve: ClassVar[bool]
    include_ascension: ClassVar[bool]
    include_assumption: ClassVar[bool]
    include_whit_sunday: ClassVar[bool]
    whit_sunday_label: ClassVar[str]
    include_whit_monday: ClassVar[bool]
    whit_monday_label: ClassVar[str]
    include_corpus_christi: ClassVar[bool]
    include_boxing_day: ClassVar[bool]
    boxing_day_label: ClassVar[str]
    include_all_souls: ClassVar[bool]
    def get_fat_tuesday(self, year: int) -> datetime.date: ...
    def get_ash_wednesday(self, year: int) -> datetime.date: ...
    def get_palm_sunday(self, year: int) -> datetime.date: ...
    def get_holy_thursday(self, year: int) -> datetime.date: ...
    def get_good_friday(self, year: int) -> datetime.date: ...
    def get_clean_monday(self, year: int) -> datetime.date: ...
    def get_easter_saturday(self, year: int) -> datetime.date: ...
    def get_easter_sunday(self, year: int) -> datetime.date: ...
    def get_easter_monday(self, year: int) -> datetime.date: ...
    def get_ascension_thursday(self, year: int) -> datetime.date: ...
    def get_whit_monday(self, year: int) -> datetime.date: ...
    def get_whit_sunday(self, year: int) -> datetime.date: ...
    def get_corpus_christi(self, year: int) -> datetime.date: ...
    def shift_christmas_boxing_days(self, year: int) -> list[tuple[datetime.date, str]]: ...
    def get_variable_days(self, year: int) -> list[tuple[datetime.date, str]]: ...

class WesternMixin(ChristianMixin):
    EASTER_METHOD: ClassVar[Literal[1, 2, 3]] = 3
    WEEKEND_DAYS: ClassVar[tuple[int, ...]]

class OrthodoxMixin(ChristianMixin):
    EASTER_METHOD: ClassVar[Literal[1, 2, 3]] = 2
    WEEKEND_DAYS: ClassVar[tuple[int, ...]]
    include_orthodox_christmas: ClassVar[bool]
    orthodox_christmas_day_label: ClassVar[str]
    def get_fixed_holidays(self, year: int) -> list[tuple[datetime.date, str]]: ...

class LunarMixin:
    @staticmethod
    def lunar(year: int, month: int, day: int) -> datetime.date: ...

class ChineseNewYearMixin(LunarMixin):
    include_chinese_new_year_eve: ClassVar[bool]
    chinese_new_year_eve_label: ClassVar[str]
    include_chinese_new_year: ClassVar[bool]
    chinese_new_year_label: ClassVar[str]
    include_chinese_second_day: ClassVar[bool]
    chinese_second_day_label: ClassVar[str]
    include_chinese_third_day: ClassVar[bool]
    chinese_third_day_label: ClassVar[str]
    shift_sunday_holidays: ClassVar[bool]
    shift_start_cny_sunday: ClassVar[bool]
    def get_chinese_new_year(self, year: int) -> list[tuple[datetime.date, str]]: ...
    def get_variable_days(self, year: int) -> list[tuple[datetime.date, str]]: ...
    def get_shifted_holidays(self, dates: Iterable[tuple[_D, str]]) -> Generator[tuple[_D, str]]: ...
    def get_calendar_holidays(self, year: int) -> list[tuple[datetime.date, str]]: ...

class CalverterMixin:
    conversion_method: Incomplete
    ISLAMIC_HOLIDAYS: ClassVar[tuple[tuple[int, int, str], ...]]
    def __init__(self, *args, **kwargs) -> None: ...
    def converted(self, year: int) -> list[tuple[int, int, int]]: ...
    def calverted_years(self, year: int) -> list[int]: ...
    def get_islamic_holidays(self) -> tuple[tuple[int, int, str], ...]: ...
    def get_delta_islamic_holidays(self, year: int) -> datetime.timedelta | None: ...
    def get_variable_days(self, year: int) -> list[tuple[datetime.date, str]]: ...

class IslamicMixin(CalverterMixin):
    WEEKEND_DAYS: ClassVar[tuple[int, ...]]
    conversion_method: Incomplete
    include_prophet_birthday: ClassVar[bool]
    include_day_after_prophet_birthday: ClassVar[bool]
    include_start_ramadan: ClassVar[bool]
    include_eid_al_fitr: ClassVar[bool]
    length_eid_al_fitr: int
    eid_al_fitr_label: ClassVar[str]
    include_eid_al_adha: ClassVar[bool]
    eid_al_adha_label: ClassVar[str]
    length_eid_al_adha: int
    include_day_of_sacrifice: ClassVar[bool]
    day_of_sacrifice_label: ClassVar[str]
    include_islamic_new_year: ClassVar[bool]
    include_laylat_al_qadr: ClassVar[bool]
    include_nuzul_al_quran: ClassVar[bool]
    def get_islamic_holidays(self) -> tuple[tuple[int, int, str]]: ...

class CoreCalendar:
    FIXED_HOLIDAYS: ClassVar[tuple[tuple[int, int, str], ...]]
    WEEKEND_DAYS: ClassVar[tuple[int, ...]]
    def __init__(self) -> None: ...
    def name(cls) -> str: ...
    def get_fixed_holidays(self, year: int) -> list[tuple[datetime.date, str]]: ...
    def get_variable_days(self, year: int) -> list[tuple[datetime.date, str]]: ...
    def get_calendar_holidays(self, year: int) -> list[tuple[datetime.date, str]]: ...
    def holidays(self, year: int | None = None) -> list[tuple[datetime.date, str]]: ...
    def get_holiday_label(self, day: datetime.date | datetime.datetime) -> str | None: ...
    def holidays_set(self, year: int | None = None) -> set[datetime.date]: ...
    def get_weekend_days(self) -> tuple[int, ...]: ...
    def is_working_day(
        self,
        day: datetime.date | datetime.datetime,
        extra_working_days: Iterable[datetime.date | datetime.datetime] | None = None,
        extra_holidays: Iterable[datetime.date | datetime.datetime] | None = None,
    ) -> bool: ...
    def is_holiday(
        self, day: datetime.date | datetime.datetime, extra_holidays: Iterable[datetime.date | datetime.datetime] | None = None
    ) -> bool: ...
    def add_working_days(
        self,
        day: datetime.date | datetime.datetime,
        delta: int,
        extra_working_days: Iterable[datetime.date | datetime.datetime] | None = None,
        extra_holidays: Iterable[datetime.date | datetime.datetime] | None = None,
        keep_datetime: bool = False,
    ) -> datetime.date: ...
    def sub_working_days(
        self,
        day: datetime.date | datetime.datetime,
        delta: int,
        extra_working_days: Iterable[datetime.date | datetime.datetime] | None = None,
        extra_holidays: Iterable[datetime.date | datetime.datetime] | None = None,
        keep_datetime: bool = False,
    ) -> datetime.date: ...
    def find_following_working_day(self, day: datetime.date) -> datetime.date: ...
    @staticmethod
    def get_nth_weekday_in_month(
        year: int, month: int, weekday: int, n: int = 1, start: datetime.date | Literal[False] | None = None
    ) -> datetime.date | None: ...
    @staticmethod
    def get_last_weekday_in_month(year: int, month: int, weekday: int) -> datetime.date: ...
    @staticmethod
    def get_iso_week_date(year: int, week_nb: int, weekday: int = 1) -> datetime.date: ...
    @staticmethod
    def get_first_weekday_after(day: _D, weekday: int) -> _D: ...
    def get_working_days_delta(
        self,
        start: datetime.date | datetime.datetime,
        end: datetime.date | datetime.datetime,
        include_start: bool = False,
        extra_working_days: Iterable[datetime.date | datetime.datetime] | None = None,
        extra_holidays: Iterable[datetime.date | datetime.datetime] | None = None,
    ) -> int: ...
    def export_to_ical(
        self, period: tuple[int, int] | list[int] = [2000, 2030], target_path: StrPath | None = None
    ) -> str | None: ...

class Calendar(CoreCalendar):
    include_new_years_day: ClassVar[bool]
    include_new_years_eve: ClassVar[bool]
    shift_new_years_day: ClassVar[bool]
    include_labour_day: ClassVar[bool]
    labour_day_label: ClassVar[str]
    def __init__(self, **kwargs) -> None: ...
    def get_fixed_holidays(self, year: int) -> list[tuple[datetime.date, str]]: ...
    def get_variable_days(self, year: int) -> list[tuple[datetime.date, str]]: ...

class WesternCalendar(WesternMixin, Calendar): ...
class OrthodoxCalendar(OrthodoxMixin, Calendar): ...

class ChineseNewYearCalendar(ChineseNewYearMixin, Calendar):
    WEEKEND_DAYS: ClassVar[tuple[int, ...]]

class IslamicCalendar(IslamicMixin, Calendar): ...

class IslamoWesternCalendar(IslamicMixin, WesternMixin, Calendar):
    FIXED_HOLIDAYS: ClassVar[tuple[tuple[int, int, str], ...]]
