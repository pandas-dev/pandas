from decimal import Decimal
from typing import Literal, TypedDict, type_check_only

@type_check_only
class _TimeEstimate(TypedDict):
    crack_times_seconds: _CrackTimeSeconds
    crack_times_display: _CrackTimesDisplay
    score: Literal[0, 1, 2, 3, 4]

@type_check_only
class _CrackTimeSeconds(TypedDict):
    online_throttling_100_per_hour: Decimal
    online_no_throttling_10_per_second: Decimal
    offline_slow_hashing_1e4_per_second: Decimal
    offline_fast_hashing_1e10_per_second: Decimal

@type_check_only
class _CrackTimesDisplay(TypedDict):
    online_throttling_100_per_hour: str
    online_no_throttling_10_per_second: str
    offline_slow_hashing_1e4_per_second: str
    offline_fast_hashing_1e10_per_second: str

def estimate_attack_times(guesses: Decimal | float) -> _TimeEstimate: ...
def guesses_to_score(guesses: Decimal) -> Literal[0, 1, 2, 3, 4]: ...
def display_time(seconds: float) -> str: ...
def float_to_decimal(f: float) -> Decimal: ...
