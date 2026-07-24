from builtins import bytes as bytes, str as str
from collections import OrderedDict as OrderedDict

# If simplejson is installed, JSONDecodeError is actually imported from there.
from json import JSONDecodeError as JSONDecodeError
from typing import Literal
from typing_extensions import TypeAlias
from urllib.parse import (
    quote as quote,
    quote_plus as quote_plus,
    unquote as unquote,
    unquote_plus as unquote_plus,
    urldefrag as urldefrag,
    urlencode as urlencode,
    urljoin as urljoin,
    urlparse as urlparse,
    urlsplit as urlsplit,
    urlunparse as urlunparse,
)
from urllib.request import getproxies as getproxies, parse_http_list as parse_http_list, proxy_bypass as proxy_bypass

is_urllib3_1: bool
is_py2: Literal[False]
is_py3: Literal[True]
has_simplejson: bool

builtin_str: TypeAlias = str  # noqa: Y042
basestring: tuple[type, ...]
numeric_types: tuple[type, ...]
integer_types: tuple[type, ...]
