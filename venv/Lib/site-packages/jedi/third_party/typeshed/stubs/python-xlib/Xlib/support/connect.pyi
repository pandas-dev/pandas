# Ignore OpenVMS in typeshed
from typing import Final

from Xlib.support.unix_connect import get_auth as get_auth, get_display as get_display, get_socket as get_socket

platform: Final[str]
