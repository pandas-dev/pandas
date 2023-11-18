from datetime import (
    datetime,
    tzinfo,
)

import numpy as np

DT64NS_DTYPE: np.dtype
TD64NS_DTYPE: np.dtype

def localize_pydatetime(dt: datetime, tz: tzinfo | None) -> datetime: ...
