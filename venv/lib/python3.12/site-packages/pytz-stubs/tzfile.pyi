from typing import IO

from pytz.tzinfo import DstTzInfo

def build_tzinfo(zone: str, fp: IO[bytes]) -> DstTzInfo: ...
