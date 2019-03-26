from datetime import datetime, timedelta
from pathlib import Path
from typing import AnyStr, IO, Union

from pandas._libs.tslibs import Period, NaT, Timedelta, Timestamp


DateTimeLike = Union[datetime, timedelta, Period, Timedelta, Timestamp]


NullableDateTimeLike = Union[NaT, DateTimeLike]


FilePathOrBuffer = Union[str, Path, IO[AnyStr]]
