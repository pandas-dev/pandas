from collections.abc import Iterable
from decimal import Decimal
from re import Pattern
from typing import TypedDict, type_check_only

from .adjacency_graphs import _Graph
from .matching import _Match

def calc_average_degree(graph: _Graph) -> float: ...

BRUTEFORCE_CARDINALITY: int
MIN_GUESSES_BEFORE_GROWING_SEQUENCE: int
MIN_SUBMATCH_GUESSES_SINGLE_CHAR: int
MIN_SUBMATCH_GUESSES_MULTI_CHAR: int
MIN_YEAR_SPACE: int
REFERENCE_YEAR: int

@type_check_only
class _GuessesResult(TypedDict):
    password: str
    guesses: int
    guesses_log10: float
    sequence: list[_Match]

def nCk(n: int, k: int) -> float: ...
def most_guessable_match_sequence(
    password: str, matches: Iterable[_Match], _exclude_additive: bool = False
) -> _GuessesResult: ...
def estimate_guesses(match: _Match, password: str) -> Decimal: ...
def bruteforce_guesses(match: _Match) -> int: ...
def dictionary_guesses(match: _Match) -> int: ...
def repeat_guesses(match: _Match) -> Decimal: ...
def sequence_guesses(match: _Match) -> int: ...
def regex_guesses(match: _Match) -> int | None: ...
def date_guesses(match: _Match) -> int: ...

KEYBOARD_AVERAGE_DEGREE: float
KEYPAD_AVERAGE_DEGREE: float
KEYBOARD_STARTING_POSITIONS: int
KEYPAD_STARTING_POSITIONS: int

def spatial_guesses(match: _Match) -> int: ...

START_UPPER: Pattern[str]
END_UPPER: Pattern[str]
ALL_UPPER: Pattern[str]
ALL_LOWER: Pattern[str]

def uppercase_variations(match: _Match) -> int: ...
def l33t_variations(match: _Match) -> int: ...
