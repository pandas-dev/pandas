from datetime import (
    datetime,
    tzinfo,
)
from typing import (
    Callable,
    Optional,
    Union,
)

import numpy as np

# imported from dateutil.tz
dateutil_gettz: Callable[[str], tzinfo]


def tz_standardize(tz: tzinfo) -> tzinfo: ...

def tz_compare(start: Optional[tzinfo], end: Optional[tzinfo]) -> bool: ...

def infer_tzinfo(
    start: Optional[datetime], end: Optional[datetime],
) -> Optional[tzinfo]: ...

# ndarrays returned are both int64_t
def get_dst_info(tz: tzinfo) -> tuple[np.ndarray, np.ndarray, str]: ...

def maybe_get_tz(tz: Optional[Union[str, int, np.int64, tzinfo]]) -> Optional[tzinfo]: ...

def get_timezone(tz: tzinfo) -> Union[tzinfo, str]: ...

def is_utc(tz: Optional[tzinfo]) -> bool: ...
