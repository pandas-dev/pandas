from collections.abc import Collection, Iterator, KeysView
from enum import Enum
from http.cookiejar import Cookie, CookiePolicy, MozillaCookieJar
from http.cookies import SimpleCookie
from typing import Final, TextIO, TypeVar

from . import _LoggerProtocol
from .minicurses import MultilinePrinter
from .utils._utils import YoutubeDLError
from .YoutubeDL import YoutubeDL

CHROMIUM_BASED_BROWSERS: Final[set[str]]
SUPPORTED_BROWSERS: Final[set[str]]

class _LinuxKeyring(Enum):
    BASICTEXT = 5
    GNOMEKEYRING = 4
    KWALLET = 1
    KWALLET5 = 2
    KWALLET6 = 3

SUPPORTED_KEYRINGS: Final[KeysView[str]]

class YDLLogger(_LoggerProtocol):
    def warning(self, message: str, only_once: bool = False) -> None: ...  # type: ignore[override]

    class ProgressBar(MultilinePrinter):
        def print(self, message: str) -> None: ...

    def progress_bar(self) -> ProgressBar: ...

class CookieLoadError(YoutubeDLError): ...

class YoutubeDLCookieJar(MozillaCookieJar):
    def __init__(self, filename: str | None = None, delayload: bool = False, policy: CookiePolicy | None = None) -> None: ...
    def open(self, file: str, *, write: bool = False) -> Iterator[TextIO]: ...
    def get_cookie_header(self, url: str) -> str: ...
    def get_cookies_for_url(self, url: str) -> list[Cookie]: ...
    def load(self, filename: str | None = None, ignore_discard: bool = True, ignore_expires: bool = True) -> None: ...
    def save(self, filename: str | None = None, ignore_discard: bool = True, ignore_expires: bool = True) -> None: ...

def load_cookies(cookie_file: str, browser_specification: str | None, ydl: YoutubeDL) -> YoutubeDLCookieJar: ...
def extract_cookies_from_browser(
    browser_name: str,
    profile: str | None = None,
    logger: _LoggerProtocol = ...,
    *,
    keyring: _LinuxKeyring | None = None,
    container: str | None = None,
) -> YoutubeDLCookieJar: ...

_T = TypeVar("_T", bound=MozillaCookieJar)

def parse_safari_cookies(data: bytes, jar: _T | None = None, logger: _LoggerProtocol = ...) -> _T: ...

class ChromeCookieDecryptor:
    def decrypt(self, encrypted_value: bytes) -> str: ...

class LinuxChromeCookieDecryptor(ChromeCookieDecryptor):
    def __init__(
        self,
        browser_keyring_name: str,
        logger: _LoggerProtocol,
        *,
        keyring: _LinuxKeyring | None = None,
        meta_version: int | None = None,
    ) -> None: ...
    @staticmethod
    def derive_key(password: bytes) -> bytes: ...

class MacChromeCookieDecryptor(ChromeCookieDecryptor):
    def __init__(self, browser_keyring_name: str, logger: YDLLogger, meta_version: int | None = None) -> None: ...
    @staticmethod
    def derive_key(password: bytes) -> bytes: ...

class WindowsChromeCookieDecryptor(ChromeCookieDecryptor):
    def __init__(self, browser_root: str, logger: YDLLogger, meta_version: int | None = None) -> None: ...

def get_cookie_decryptor(
    browser_root: str,
    browser_keyring_name: str,
    logger: _LoggerProtocol,
    *,
    keyring: _LinuxKeyring | None = None,
    meta_version: int | None = None,
) -> ChromeCookieDecryptor: ...

class ParserError(Exception): ...

class DataParser:
    def __init__(self, data: bytes, logger: YDLLogger) -> None: ...
    def read_bytes(self, num_bytes: int) -> bytes: ...
    def expect_bytes(self, expected_value: bytes, message: str) -> None: ...
    def read_uint(self, big_endian: bool = False) -> int: ...
    def read_double(self, big_endian: bool = False) -> float: ...
    def read_cstring(self) -> bytes: ...
    def skip(self, num_bytes: int, description: str = "unknown") -> None: ...
    def skip_to(self, offset: int, description: str = "unknown") -> None: ...
    def skip_to_end(self, description: str = "unknown") -> None: ...

def pbkdf2_sha1(password: bytes, salt: bytes, iterations: int, key_length: int) -> bytes: ...

class LenientSimpleCookie(SimpleCookie):
    def load(self, data: str | Collection[str]) -> None: ...
