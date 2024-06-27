# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
A set of utilities for writing output to the console.
"""

import contextlib
import locale
import logging
import os
import sys
import textwrap
import time

from asv_runner import util

WIN = os.name == "nt"


def isatty(file):
    """
    Determines if a file is a tty.

    #### Parameters
    **file** (`file-like object`)
    : The file-like object to check.

    #### Returns
    **isatty** (`bool`)
    : Returns `True` if the file is a tty, `False` otherwise.

    #### Notes
    Most built-in Python file-like objects have an `isatty` member,
    but some user-defined types may not. In such cases, this function
    assumes those are not ttys.
    """
    return file.isatty() if hasattr(file, "isatty") else False


def _color_text(text, color):
    """
    Returns a string wrapped in ANSI color codes for coloring the text in a terminal.

    #### Parameters
    **text** (`str`)
    : The string to colorize.

    **color** (`str`)
    : An ANSI terminal color name. Must be one of the following:
      'black', 'red', 'green', 'brown', 'blue', 'magenta', 'cyan', 'lightgrey',
      'default', 'darkgrey', 'lightred', 'lightgreen', 'yellow', 'lightblue',
      'lightmagenta', 'lightcyan', 'white', or '' (the empty string).

    #### Returns
    **colored_text** (`str`)
    : The input string, bounded by the appropriate ANSI color codes.

    #### Notes
    This function wraps the input text with ANSI color codes based on the given color.
    It won't actually affect the text until it is printed to the terminal.
    """
    color_mapping = {
        "black": "0;30",
        "red": "0;31",
        "green": "0;32",
        "brown": "0;33",
        "blue": "0;34",
        "magenta": "0;35",
        "cyan": "0;36",
        "lightgrey": "0;37",
        "default": "0;39",
        "darkgrey": "1;30",
        "lightred": "1;31",
        "lightgreen": "1;32",
        "yellow": "1;33",
        "lightblue": "1;34",
        "lightmagenta": "1;35",
        "lightcyan": "1;36",
        "white": "1;37",
    }

    color_code = color_mapping.get(color, "0;39")
    return f"\033[{color_code}m{text}\033[0m"


# A dictionary of Unicode characters that have reasonable representations in ASCII.
# This dictionary contains Unicode characters as keys and their corresponding ASCII
# representations as values. This allows for convenient replacement of these specific
# Unicode characters with ASCII ones to prevent them from being replaced by '?'.
#
# The mapping currently includes:
# - 'μ' maps to 'u'
# - '·' maps to '-'
# - '±' maps to '~'
#
# You can find additional characters that might need an entry using:
# `grep -P  -n '[^\x00-\x7F]' -r *`
# in the `asv` source directory.
_unicode_translations = {ord("μ"): "u", ord("·"): "-", ord("±"): "~"}


def _write_with_fallback(s, fileobj):
    """
    Writes the supplied string to the given file-like object, handling potential
    UnicodeEncodeErrors by falling back to the locale's preferred encoding.

    #### Parameters
    `s` (`str`):
    The Unicode string to be written to the file-like object. Raises a `ValueError`
    if `s` is not a Unicode string.

    `fileobj` (file-like object):
    The file-like object to which the string `s` is to be written. On Python 3,
    this must be a text stream. On Python 2, this must be a `file` byte stream.

    #### Notes
    This function first tries to write the input string `s` to the file object
    `fileobj`. If a `UnicodeError` occurs during this process (indicating that the
    string contains characters not representable in the file's encoding), the function
    falls back to encoding the string in the locale's preferred encoding before writing.

    If the string `s` still cannot be encoded in the locale's preferred encoding, the
    function translates the string to replace problematic Unicode characters with
    ASCII ones using the `_unicode_translations` dictionary, and then encodes and
    writes the resulting string to `fileobj` using the "replace" error handling scheme
    (which replaces any non-encodable characters with a suitable replacement marker).

    After the write operation, the function flushes the file object's output buffer to
    ensure that the written data is actually saved to the file.
    """
    if not isinstance(s, str):
        raise ValueError("Input string is not a Unicode string")

    with contextlib.suppress(UnicodeError):
        fileobj.write(s)
        return
    # Fall back to writing bytes
    enc = locale.getpreferredencoding()
    try:
        b = s.encode(enc)
    except UnicodeError:
        s = s.translate(_unicode_translations)
        b = s.encode(enc, errors="replace")

    fileobj.flush()
    fileobj.buffer.write(b)


def color_print(*args, **kwargs):
    """
    Prints colored and styled text to the terminal using ANSI escape sequences.

    #### Parameters
    *args (`tuple` of `str`):
    The positional arguments should come in pairs (`msg`, `color`), where `msg`
    is the string to display and `color` is the color to display it in. `color`
    is an ANSI terminal color name. Must be one of: black, red, green, brown,
    blue, magenta, cyan, lightgrey, default, darkgrey, lightred, lightgreen,
    yellow, lightblue, lightmagenta, lightcyan, white, or '' (the empty string).

    `file` (writable file-like object, optional):
    Where to write to. Defaults to `sys.stdout`. If `file` is not a tty (as determined
    by calling its `isatty` member, if one exists), no coloring will be included. It's
    passed as a keyword argument.

    `end` (`str`, optional):
    The ending of the message. Defaults to "\n". The `end` will be printed after
    resetting any color or font state. It's passed as a keyword argument.

    #### Notes
    This function allows you to print text in various colors to the console, which can
    be helpful for distinguishing different kinds of output or for drawing attention to
    particular messages.

    It works by applying ANSI escape sequences to the input strings according to the
    specified colors. These escape sequences are interpreted by the terminal emulator
    to apply the specified colors and styles.

    #### Example
    ```{code-block} python
    color_print('This is the color ', 'default', 'GREEN', 'green')
    ```
    """
    file = kwargs.get("file", sys.stdout)
    end = kwargs.get("end", "\n")

    if isatty(file) and not WIN:
        for i in range(0, len(args), 2):
            msg = args[i]
            color = "" if i + 1 == len(args) else args[i + 1]
            if color:
                msg = _color_text(msg, color)
            _write_with_fallback(msg, file)

    else:
        for i in range(0, len(args), 2):
            msg = args[i]
            _write_with_fallback(msg, file)

    _write_with_fallback(end, file)


def get_answer_default(prompt, default, use_defaults=False):
    """
    Prompts the user for input and returns the entered value or a default.

    #### Parameters
    `prompt` (`str`):
    The string that is presented to the user.

    `default` (any):
    The value returned if the user doesn't enter anything and just hits Enter. This
    value is also shown in the prompt to indicate to the user what the default is.

    `use_defaults` (`bool`, optional):
    If True, the function will immediately return the default value without prompting
    the user for input. Defaults to False.

    #### Returns
    The user's input, or the provided default value if the user didn't enter anything.

    #### Notes
    This function enhances the built-in `input` function by allowing a default value
    to be specified, which is returned if the user doesn't enter anything.
    """
    color_print(f"{prompt} [{default}]: ", end="")

    if use_defaults:
        return default

    x = input()
    return default if x.strip() == "" else x


def truncate_left(s, l):
    return f"...{s[-(l - 3):]}" if len(s) > l else s


class Log:
    def __init__(self):
        self._indent = 1
        self._total = 0
        self._count = 0
        self._logger = logging.getLogger()
        self._needs_newline = False
        self._last_dot = time.time()
        self._colorama = False
        if sys.platform in {"win32", "cli"}:
            try:
                import colorama

                colorama.init()
                self._colorama = True
            except Exception as exc:
                print(f"On Windows or cli, colorama is suggested, but got {exc}")

    def _stream_formatter(self, record):
        """
        The formatter for standard output
        """
        if self._needs_newline:
            color_print("")
        parts = record.msg.split("\n", 1)
        first_line = parts[0]
        rest = None if len(parts) == 1 else parts[1]
        indent = self._indent + 1
        continued = getattr(record, "continued", False)

        if self._total:
            progress_msg = f"[{self._count / self._total:6.02%}] "
            if not continued:
                color_print(progress_msg, end="")
            indent += len(progress_msg)

        if not continued:
            color_print("·" * self._indent, end="")
            color_print(" ", end="")
        else:
            color_print(" " * indent, end="")

        if hasattr(record, "color"):
            color = record.color
        elif record.levelno < logging.DEBUG:
            color = "default"
        elif record.levelno < logging.INFO:
            color = "default"
        elif record.levelno < logging.WARN:
            if self._indent == 1:
                color = "green"
            elif self._indent == 2:
                color = "blue"
            else:
                color = "default"
        elif record.levelno < logging.ERROR:
            color = "brown"
        else:
            color = "red"

        color_print(first_line, color, end="")
        if rest is not None:
            color_print("")
            detail = textwrap.dedent(rest)
            spaces = " " * indent
            for line in detail.split("\n"):
                color_print(spaces, end="")
                color_print(line)

        self._needs_newline = True
        sys.stdout.flush()

    @contextlib.contextmanager
    def indent(self):
        """
        A context manager to increase the indentation level.
        """
        self._indent += 1
        yield
        self._indent -= 1

    def dot(self):
        if isatty(sys.stdout):
            if time.time() > self._last_dot + 1.0:
                color_print(".", "darkgrey", end="")
                sys.stdout.flush()
                self._last_dot = time.time()

    def set_nitems(self, n):
        """
        Set the number of remaining items to process.  Each of these
        steps should be incremented through using `step`.

        Can be called multiple times. The progress percentage is ensured
        to be non-decreasing, except if 100% was already reached in which
        case it is restarted from 0%.
        """
        try:
            # Ensure count/total is nondecreasing
            self._total = util.ceildiv(n * self._total, self._total - self._count)
            self._count = self._total - n
        except ZeroDivisionError:
            # Reset counting from start
            self._total = n
            self._count = 0

    def step(self):
        """
        Write that a step has been completed.  A percentage is
        displayed along with it.

        If we are stepping beyond the number of items, stop counting.
        """
        self._count = min(self._total, self._count + 1)

    def enable(self, verbose=False):
        sh = logging.StreamHandler()
        sh.emit = self._stream_formatter
        self._logger.addHandler(sh)
        if verbose:
            self._logger.setLevel(logging.DEBUG)
        else:
            self._logger.setLevel(logging.INFO)

    @contextlib.contextmanager
    def set_level(self, level):
        orig_level = self._logger.level
        if not self.is_debug_enabled():
            self._logger.setLevel(level)
        try:
            yield
        finally:
            self._logger.setLevel(orig_level)

    def is_debug_enabled(self):
        return self._logger.getEffectiveLevel() <= logging.DEBUG

    def _message(
        self, routine, message, reserve_space=False, color=None, continued=False
    ):
        kwargs = {}
        extra = {}
        if color is not None:
            extra["color"] = color
        if continued:
            extra["continued"] = True
        if extra:
            kwargs["extra"] = extra

        if reserve_space:
            max_width = max(16, util.terminal_width - 33)
            message = truncate_left(message, max_width)
            self._prev_message = message

        routine(message, **kwargs)

    def info(self, *args, **kwargs):
        self._message(self._logger.info, *args, **kwargs)

    def warning(self, *args, **kwargs):
        self._message(self._logger.warning, *args, **kwargs)

    def debug(self, *args, **kwargs):
        self._message(self._logger.debug, *args, **kwargs)

    def error(self, *args, **kwargs):
        self._message(self._logger.error, *args, **kwargs)

    def add(self, msg):
        if self._needs_newline:
            _write_with_fallback(msg, sys.stdout)
            sys.stdout.flush()
        else:
            self.info(msg)

    def add_padded(self, msg):
        """
        Final part of two-part info message.
        Should be preceded by a call to info/warn/...(msg, reserve_space=True)
        """
        if self._prev_message is None:
            # No previous part: print as an info message
            self.info(msg)
            return

        padding_length = (
            util.terminal_width - len(self._prev_message) - 14 - 1 - len(msg)
        )
        if WIN:
            padding_length -= 1
        padding = " " * padding_length

        self._prev_message = None
        self.add(f" {padding}{msg}")

    def flush(self):
        """
        Flush any trailing newlines. Needs to be called before printing
        to stdout via other means, after using Log.
        """
        if self._needs_newline:
            color_print("")
            self._needs_newline = False
        sys.stdout.flush()
