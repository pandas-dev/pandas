import math
from typing import Union

INF = float("inf")
MINUS_INF = float("-inf")
NaN = float("NaN")


def floatToGoString(d):
    d = float(d)
    if d == INF:
        return '+Inf'
    elif d == MINUS_INF:
        return '-Inf'
    elif math.isnan(d):
        return 'NaN'
    else:
        s = repr(d)
        dot = s.find('.')
        # Go switches to exponents sooner than Python.
        # We only need to care about positive values for le/quantile.
        if d > 0 and dot > 6:
            mantissa = f'{s[0]}.{s[1:dot]}{s[dot + 1:]}'.rstrip('0.')
            return f'{mantissa}e+0{dot - 1}'
        return s


def parse_version(version_str: str) -> tuple[Union[int, str], ...]:
    version: list[Union[int, str]] = []
    for part in version_str.split('.'):
        try:
            version.append(int(part))
        except ValueError:
            version.append(part)

    return tuple(version)
