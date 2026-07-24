from enum import Enum
from typing import Final

class _Days(Enum):
    ROSH_HASHANA = "Rosh Hashana"
    YOM_KIPPUR = "Yom Kippur"
    SUCCOS = "Succos"
    SHMINI_ATZERES = "Shmini Atzeres"
    SIMCHAS_TORAH = "Simchas Torah"
    CHANUKA = "Chanuka"
    TU_BSHVAT = "Tu B'shvat"
    PURIM_KATAN = "Purim Katan"
    PURIM = "Purim"
    SHUSHAN_PURIM = "Shushan Purim"
    PESACH = "Pesach"
    PESACH_SHENI = "Pesach Sheni"
    LAG_BAOMER = "Lag Ba'omer"
    SHAVUOS = "Shavuos"
    TU_BAV = "Tu B'av"
    TZOM_GEDALIA = "Tzom Gedalia"
    TENTH_OF_TEVES = "10 of Teves"
    TAANIS_ESTHER = "Taanis Esther"
    SEVENTEENTH_OF_TAMUZ = "17 of Tamuz"
    NINTH_OF_AV = "9 of Av"

MONTH_NAMES: Final[list[str]]
MONTH_NAMES_HEBREW: Final[list[str]]
FAST_DAYS: Final[list[str]]
FAST_DAYS_HEBREW: Final[list[str]]
FESTIVALS: Final[list[str]]
FESTIVALS_HEBREW: Final[list[str]]
WEEKDAYS: Final[dict[int, str]]
