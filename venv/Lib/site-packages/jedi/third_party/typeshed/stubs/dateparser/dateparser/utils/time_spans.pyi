import datetime
from typing import TypedDict, TypeVar, type_check_only
from typing_extensions import NotRequired

from dateparser.conf import Settings

_DateT = TypeVar("_DateT", bound=datetime.date)

@type_check_only
class _SpanInformation(TypedDict):
    type: str
    direction: str
    matched_text: str
    start_pos: int
    end_pos: int
    number: NotRequired[int]

def get_week_start(date: _DateT, start_of_week: str = "monday") -> _DateT: ...
def get_week_end(date: _DateT, start_of_week: str = "monday") -> _DateT: ...
def detect_time_span(text: str) -> _SpanInformation | None: ...
def generate_time_span(
    span_info: _SpanInformation, base_date: _DateT | None = None, settings: Settings | None = None
) -> tuple[_DateT, _DateT]: ...
