from typing import Final

EPOCH: Final = 2394646.5
EPOCH_GREGORIAN_YEAR: Final = 1844
TEHRAN: Final = 51.4215, 35.6944
WEEKDAYS: Final = ("Jamál", "Kamál", "Fidál", "Idál", "Istijlál", "Istiqlál", "Jalál")
MONTHS: Final[tuple[str, ...]]
ENGLISH_MONTHS: Final[tuple[str, ...]]
BAHA: Final = 1
JALAL: Final = 2
JAMAL: Final = 3
AZAMAT: Final = 4
NUR: Final = 5
RAHMAT: Final = 6
KALIMAT: Final = 7
KAMAL: Final = 8
ASMA: Final = 9
IZZAT: Final = 10
MASHIYYAT: Final = 11
ILM: Final = 12
QUDRAT: Final = 13
QAWL: Final = 14
MASAIL: Final = 15
SHARAF: Final = 16
SULTAN: Final = 17
MULK: Final = 18
AYYAMIHA: Final = 19
ALA: Final = 20

def gregorian_nawruz(year: int) -> tuple[int, int]: ...
def to_jd(year: int, month: int, day: int) -> float: ...
def from_jd(jd: float) -> tuple[int, int, int]: ...
def from_gregorian(year: int, month: int, day: int) -> tuple[int, int, int]: ...
def to_gregorian(year: int, month: int, day: int) -> tuple[int, int, int]: ...
def month_length(year: int, month: int) -> int: ...
def monthcalendar(year: int, month: int) -> list[list[int | None]]: ...
def format(year: int, month: int, day: int, lang: str | None = None) -> str: ...
