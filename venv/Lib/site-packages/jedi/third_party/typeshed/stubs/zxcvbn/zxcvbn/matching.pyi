from collections.abc import Iterable, Mapping
from decimal import Decimal
from re import Pattern
from typing import Any, Literal, TypedDict, type_check_only
from typing_extensions import NotRequired

from .adjacency_graphs import _Graph

@type_check_only
class _Match(TypedDict):
    pattern: Literal["dictionary", "spatial", "repeat", "sequence", "regex", "date"]
    token: str
    i: int
    j: int
    guesses: NotRequired[int]  # all patterns except 'date'
    guesses_log10: NotRequired[float]  # all patterns except 'date'

    # pattern == 'date'
    separator: NotRequired[str]
    year: NotRequired[int]
    month: NotRequired[int]
    day: NotRequired[int]

    # pattern == 'dictionary'
    matched_word: NotRequired[str]
    dictionary_name: NotRequired[str]
    l33t: NotRequired[bool]
    reversed: NotRequired[bool]
    rank: NotRequired[int]
    base_guesses: NotRequired[int | Decimal]  # Decimal for 'repeat', see below
    uppercase_variations: NotRequired[int]
    l33t_variations: NotRequired[int]

    # pattern == 'spatial'
    turns: NotRequired[int]

    # pattern == 'repeat'
    base_token: NotRequired[str]
    # base_guesses: NotRequired[Decimal]
    base_matches: NotRequired[list[Any]]  # Any = _Match, https://github.com/python/mypy/issues/731
    repeat_count: NotRequired[float]

    # pattern == 'regex'
    regex_name: NotRequired[str]

def build_ranked_dict(ordered_list: Iterable[str]) -> dict[str, int]: ...

RANKED_DICTIONARIES: dict[str, dict[str, int]]

def add_frequency_lists(frequency_lists_: Mapping[str, Iterable[str]]) -> None: ...

GRAPHS: dict[str, dict[str, list[str | None]]]
L33T_TABLE: dict[str, list[str]]
REGEXEN: dict[str, Pattern[str]]
DATE_MAX_YEAR: int
DATE_MIN_YEAR: int
DATE_SPLITS: dict[int, list[list[int]]]

def omnimatch(password: str, _ranked_dictionaries: dict[str, dict[str, int]] = ...) -> list[_Match]: ...
def dictionary_match(password: str, _ranked_dictionaries: dict[str, dict[str, int]] = ...) -> list[_Match]: ...
def reverse_dictionary_match(password: str, _ranked_dictionaries: dict[str, dict[str, int]] = ...) -> list[_Match]: ...
def relevant_l33t_subtable(password: str, table: Mapping[str, Iterable[str]]) -> dict[str, list[str]]: ...
def enumerate_l33t_subs(table: Mapping[str, Iterable[str]]) -> list[dict[str, str]]: ...
def translate(string: str, chr_map: Mapping[str, str]) -> str: ...
def l33t_match(
    password: str, _ranked_dictionaries: dict[str, dict[str, int]] = ..., _l33t_table: dict[str, list[str]] = ...
) -> list[_Match]: ...
def repeat_match(password: str, _ranked_dictionaries: dict[str, dict[str, int]] = ...) -> list[_Match]: ...
def spatial_match(
    password: str, _graphs: dict[str, _Graph] = ..., _ranked_dictionaries: dict[str, dict[str, int]] = ...
) -> list[_Match]: ...

SHIFTED_RX: Pattern[str]

def spatial_match_helper(password: str, graph: _Graph, graph_name: str) -> list[_Match]: ...

MAX_DELTA: int

def sequence_match(password: str, _ranked_dictionaries: dict[str, dict[str, int]] = ...) -> list[_Match]: ...
def regex_match(
    password: str, _regexen: dict[str, Pattern[str]] = ..., _ranked_dictionaries: dict[str, dict[str, int]] = ...
) -> list[_Match]: ...
def date_match(password: str, _ranked_dictionaries: dict[str, dict[str, int]] = ...) -> list[_Match]: ...
@type_check_only
class _DM(TypedDict):
    month: int
    day: int

@type_check_only
class _DMY(TypedDict):
    year: int
    month: int
    day: int

def map_ints_to_dmy(ints: tuple[int, int, int]) -> _DMY | None: ...
def map_ints_to_dm(ints: tuple[int, int]) -> _DM | None: ...
def two_to_four_digit_year(year: int) -> int: ...
