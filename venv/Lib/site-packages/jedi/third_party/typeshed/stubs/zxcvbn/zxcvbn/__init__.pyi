import datetime
from collections.abc import Iterable
from decimal import Decimal
from typing import TypedDict, type_check_only

from .feedback import _Feedback
from .matching import _Match
from .time_estimates import _TimeEstimate

@type_check_only
class _Result(_TimeEstimate, TypedDict):
    password: str
    guesses: Decimal
    guesses_log10: float
    sequence: list[_Match]
    calc_time: datetime.timedelta
    feedback: _Feedback

def zxcvbn(password: str, user_inputs: Iterable[object] | None = None, max_length: int = 72) -> _Result: ...
