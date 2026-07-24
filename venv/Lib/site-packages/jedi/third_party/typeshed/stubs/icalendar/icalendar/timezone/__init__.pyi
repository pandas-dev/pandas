from ..timezone.tzp import TZP as TZP  # to prevent "tzp" from being defined here
from .tzid import tzid_from_dt as tzid_from_dt, tzid_from_tzinfo as tzid_from_tzinfo, tzids_from_tzinfo as tzids_from_tzinfo

__all__ = ["TZP", "tzp", "use_pytz", "use_zoneinfo", "tzid_from_tzinfo", "tzid_from_dt", "tzids_from_tzinfo"]

tzp: TZP

def use_pytz() -> None: ...
def use_zoneinfo() -> None: ...
