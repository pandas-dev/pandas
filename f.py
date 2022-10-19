from typing import NamedTuple


class ISO8601Info(NamedTuple):
    format: str = b""
    date_sep: str = b""
    time_sep: str = b""
    micro_or_tz: str = b""
    year: bool = False
    month: bool = False
    day: bool = False
    hour: bool = False
    minute: bool = False
    second: bool = False


def format_is_iso(f: str):
    """
    Does format match the iso8601 set that can be handled by the C parser?
    Generally of form YYYY-MM-DDTHH:MM:SS - date separator can be different
    but must be consistent.  Leading 0s in dates and times are optional.

    no, needs doing in c. eff...
    """
    excluded_formats = ["%Y%m%d", "%Y%m", "%Y"]

    if f in excluded_formats:
        return ISO8601Info()
    for date_sep in [" ", "/", "\\", "-", ".", ""]:
        for time_sep in [" ", "T"]:
            for micro_or_tz in ["", "%z", "%Z", ".%f", ".%f%z", ".%f%Z"]:
                if f"%Y{date_sep}%m{date_sep}%d{time_sep}%H:%M:%S{micro_or_tz}" == f:
                    return ISO8601Info(
                        format=f.encode("utf-8"),
                        date_sep=date_sep.encode("utf-8"),
                        time_sep=time_sep.encode("utf-8"),
                        micro_or_tz=micro_or_tz.encode("utf-8"),
                        year=True,
                        month=True,
                        day=True,
                        hour=True,
                        minute=True,
                        second=True,
                    )
                elif f"%Y{date_sep}%m{date_sep}%d{time_sep}%H:%M" == f:
                    return ISO8601Info(
                        format=f.encode("utf-8"),
                        date_sep=date_sep.encode("utf-8"),
                        time_sep=time_sep.encode("utf-8"),
                        micro_or_tz=micro_or_tz.encode("utf-8"),
                        year=True,
                        month=True,
                        day=True,
                        hour=True,
                        minute=True,
                    )
                elif f"%Y{date_sep}%m{date_sep}%d{time_sep}%H" == f:
                    return ISO8601Info(
                        format=f.encode("utf-8"),
                        date_sep=date_sep.encode("utf-8"),
                        time_sep=time_sep.encode("utf-8"),
                        micro_or_tz=micro_or_tz.encode("utf-8"),
                        year=True,
                        month=True,
                        day=True,
                        hour=True,
                    )
                elif f"%Y{date_sep}%m{date_sep}%d" == f:
                    return ISO8601Info(
                        format=f.encode("utf-8"),
                        date_sep=date_sep.encode("utf-8"),
                        time_sep=time_sep.encode("utf-8"),
                        micro_or_tz=micro_or_tz.encode("utf-8"),
                        year=True,
                        month=True,
                        day=True,
                    )
                elif f"%Y{date_sep}%m" == f:
                    return ISO8601Info(
                        format=f.encode("utf-8"),
                        date_sep=date_sep.encode("utf-8"),
                        time_sep=time_sep.encode("utf-8"),
                        micro_or_tz=micro_or_tz.encode("utf-8"),
                        year=True,
                        month=True,
                    )
    return ISO8601Info()


if __name__ == "__main__":
    print(format_is_iso("%Y-%m-%d %H:%M:%S%z"))
# print(format_is_iso('%Y%m%d %H'))
