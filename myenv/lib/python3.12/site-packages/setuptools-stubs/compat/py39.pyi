import sys
from typing import Final

if sys.version_info >= (3, 10):
    LOCALE_ENCODING: Final = "locale"
else:
    LOCALE_ENCODING: Final = None
