from collections.abc import Generator, Iterable

from ..core import _D, OrthodoxCalendar

class Bulgaria(OrthodoxCalendar):
    def get_shifted_holidays(self, days: Iterable[tuple[_D, str]]) -> Generator[tuple[_D, str]]: ...
