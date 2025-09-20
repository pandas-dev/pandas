import functools
import inspect
import pydoc
import sys
import types
import warnings
from types import TracebackType
from typing import Any, Callable, Optional, Tuple

import stack_data
from pygments.token import Token

from IPython import get_ipython
from IPython.core import debugger
from IPython.utils import path as util_path
from IPython.utils import py3compat
from IPython.utils.PyColorize import Theme, TokenStream, theme_table

_sentinel = object()
INDENT_SIZE = 8


@functools.lru_cache
def count_lines_in_py_file(filename: str) -> int:
    """
    Given a filename, returns the number of lines in the file
    if it ends with the extension ".py". Otherwise, returns 0.
    """
    if not filename.endswith(".py"):
        return 0
    else:
        try:
            with open(filename, "r") as file:
                s = sum(1 for line in file)
        except UnicodeError:
            return 0
    return s


def get_line_number_of_frame(frame: types.FrameType) -> int:
    """
    Given a frame object, returns the total number of lines in the file
    containing the frame's code object, or the number of lines in the
    frame's source code if the file is not available.

    Parameters
    ----------
    frame : FrameType
        The frame object whose line number is to be determined.

    Returns
    -------
    int
        The total number of lines in the file containing the frame's
        code object, or the number of lines in the frame's source code
        if the file is not available.
    """
    filename = frame.f_code.co_filename
    if filename is None:
        print("No file....")
        lines, first = inspect.getsourcelines(frame)
        return first + len(lines)
    return count_lines_in_py_file(filename)


def _safe_string(value: Any, what: Any, func: Any = str) -> str:
    # Copied from cpython/Lib/traceback.py
    try:
        return func(value)
    except:
        return f"<{what} {func.__name__}() failed>"


def _format_traceback_lines(
    lines: list[stack_data.Line],
    theme: Theme,
    has_colors: bool,
    lvals_toks: list[TokenStream],
) -> TokenStream:
    """
    Format tracebacks lines with pointing arrow, leading numbers,
    this assumes the stack have been extracted using stackdata.


    Parameters
    ----------
    lines : list[Line]
    """
    numbers_width = INDENT_SIZE - 1
    tokens: TokenStream = []

    for stack_line in lines:
        if stack_line is stack_data.LINE_GAP:
            toks = [(Token.LinenoEm, "   (...)")]
            tokens.extend(toks)
            continue

        lineno = stack_line.lineno
        line = stack_line.render(pygmented=has_colors).rstrip("\n") + "\n"
        if stack_line.is_current:
            # This is the line with the error
            pad = numbers_width - len(str(lineno))
            toks = [
                (Token.LinenoEm, theme.make_arrow(pad)),
                (Token.LinenoEm, str(lineno)),
                (Token, " "),
                (Token, line),
            ]
        else:
            num = "%*s" % (numbers_width, lineno)
            toks = [
                (Token.LinenoEm, str(num)),
                (Token, " "),
                (Token, line),
            ]

        tokens.extend(toks)
        if lvals_toks and stack_line.is_current:
            for lv in lvals_toks:
                tokens.append((Token, " " * INDENT_SIZE))
                tokens.extend(lv)
                tokens.append((Token, "\n"))
            # strip the last newline
            tokens = tokens[:-1]

    return tokens


# some internal-use functions
def text_repr(value: Any) -> str:
    """Hopefully pretty robust repr equivalent."""
    # this is pretty horrible but should always return *something*
    try:
        return pydoc.text.repr(value)  # type: ignore[call-arg]
    except KeyboardInterrupt:
        raise
    except:
        try:
            return repr(value)
        except KeyboardInterrupt:
            raise
        except:
            try:
                # all still in an except block so we catch
                # getattr raising
                name = getattr(value, "__name__", None)
                if name:
                    # ick, recursion
                    return text_repr(name)
                klass = getattr(value, "__class__", None)
                if klass:
                    return "%s instance" % text_repr(klass)
                return "UNRECOVERABLE REPR FAILURE"
            except KeyboardInterrupt:
                raise
            except:
                return "UNRECOVERABLE REPR FAILURE"


def eqrepr(value: Any, repr: Callable[[Any], str] = text_repr) -> str:
    return "=%s" % repr(value)


def nullrepr(value: Any, repr: Callable[[Any], str] = text_repr) -> str:
    return ""


def _tokens_filename(
    em: bool,
    file: str | None,
    *,
    lineno: int | None = None,
) -> TokenStream:
    """
    Format filename lines with custom formatting from caching compiler or `File *.py` by default

    Parameters
    ----------
    em: wether bold or not
    file : str
    """
    Normal = Token.NormalEm if em else Token.Normal
    Filename = Token.FilenameEm if em else Token.Filename
    ipinst = get_ipython()
    if (
        ipinst is not None
        and (data := ipinst.compile.format_code_name(file)) is not None
    ):
        label, name = data
        if lineno is None:
            return [
                (Normal, label),
                (Normal, " "),
                (Filename, name),
            ]
        else:
            return [
                (Normal, label),
                (Normal, " "),
                (Filename, name),
                (Filename, f", line {lineno}"),
            ]
    else:
        name = util_path.compress_user(
            py3compat.cast_unicode(file, util_path.fs_encoding)
        )
        if lineno is None:
            return [
                (Normal, "File "),
                (Filename, name),
            ]
        else:
            return [
                (Normal, "File "),
                (Filename, f"{name}:{lineno}"),
            ]


def _simple_format_traceback_lines(
    lnum: int,
    index: int,
    lines: list[tuple[str, tuple[str, bool]]],
    lvals_toks: list[TokenStream],
    theme: Theme,
) -> TokenStream:
    """
    Format tracebacks lines with pointing arrow, leading numbers

    This should be equivalent to _format_traceback_lines, but does not rely on stackdata
    to format the lines

    This is due to the fact that stackdata may be slow on super long and complex files.

    Parameters
    ==========

    lnum: int
        number of the target line of code.
    index: int
        which line in the list should be highlighted.
    lines: list[string]
    lvals_toks: pairs of token type and str
        Values of local variables, already colored, to inject just after the error line.
    """
    for item in lvals_toks:
        assert isinstance(item, list)
        for subit in item:
            assert isinstance(subit[1], str)

    numbers_width = INDENT_SIZE - 1
    res_toks: TokenStream = []
    for i, (line, (new_line, err)) in enumerate(lines, lnum - index):
        if not err:
            line = new_line

        colored_line = line
        if i == lnum:
            # This is the line with the error
            pad = numbers_width - len(str(i))
            line_toks = [
                (Token.LinenoEm, theme.make_arrow(pad)),
                (Token.LinenoEm, str(lnum)),
                (Token, " "),
                (Token, colored_line),
            ]
        else:
            padding_num = "%*s" % (numbers_width, i)

            line_toks = [
                (Token.LinenoEm, padding_num),
                (Token, " "),
                (Token, colored_line),
            ]
        res_toks.extend(line_toks)

        if lvals_toks and i == lnum:
            for lv in lvals_toks:
                res_toks.extend(lv)
            # res_toks.extend(lvals_toks)
    return res_toks


class FrameInfo:
    """
    Mirror of stack data's FrameInfo, but so that we can bypass highlighting on
    really long frames.
    """

    description: Optional[str]
    filename: Optional[str]
    lineno: int
    # number of context lines to use
    context: Optional[int]
    raw_lines: list[str]
    _sd: stack_data.core.FrameInfo
    frame: Any

    @classmethod
    def _from_stack_data_FrameInfo(
        cls, frame_info: stack_data.core.FrameInfo | stack_data.core.RepeatedFrames
    ) -> "FrameInfo":
        return cls(
            getattr(frame_info, "description", None),
            getattr(frame_info, "filename", None),  # type: ignore[arg-type]
            getattr(frame_info, "lineno", None),  # type: ignore[arg-type]
            getattr(frame_info, "frame", None),
            getattr(frame_info, "code", None),
            sd=frame_info,
            context=None,
        )

    def __init__(
        self,
        description: Optional[str],
        filename: str,
        lineno: int,
        frame: Any,
        code: Optional[types.CodeType],
        *,
        sd: Any = None,
        context: int | None = None,
    ):
        assert isinstance(lineno, (int, type(None))), lineno
        self.description = description
        self.filename = filename
        self.lineno = lineno
        self.frame = frame
        self.code = code
        self._sd = sd
        self.context = context

        # self.lines = []
        if sd is None:
            try:
                # return a list of source lines and a starting line number
                self.raw_lines = inspect.getsourcelines(frame)[0]
            except OSError:
                self.raw_lines = [
                    "'Could not get source, probably due dynamically evaluated source code.'"
                ]

    @property
    def variables_in_executing_piece(self) -> list[Any]:
        if self._sd is not None:
            return self._sd.variables_in_executing_piece  # type:ignore[misc]
        else:
            return []

    @property
    def lines(self) -> list[Any]:
        from executing.executing import NotOneValueFound

        assert self._sd is not None
        try:
            return self._sd.lines  # type: ignore[misc]
        except NotOneValueFound:

            class Dummy:
                lineno = 0
                is_current = False

                def render(self, *, pygmented: bool) -> str:
                    return "<Error retrieving source code with stack_data see ipython/ipython#13598>"

            return [Dummy()]

    @property
    def executing(self) -> Any:
        if self._sd is not None:
            return self._sd.executing
        else:
            return None


class TBTools:
    """Basic tools used by all traceback printer classes."""

    # Number of frames to skip when reporting tracebacks
    tb_offset = 0
    _theme_name: str
    _old_theme_name: str
    call_pdb: bool
    ostream: Any
    debugger_cls: Any
    pdb: Any

    def __init__(
        self,
        color_scheme: Any = _sentinel,
        call_pdb: bool = False,
        ostream: Any = None,
        *,
        debugger_cls: type | None = None,
        theme_name: str = "nocolor",
    ):
        if color_scheme is not _sentinel:
            assert isinstance(color_scheme, str), color_scheme
            warnings.warn(
                "color_scheme is deprecated since IPython 9.0, use theme_name instead, all lowercase",
                DeprecationWarning,
                stacklevel=2,
            )
            theme_name = color_scheme
        if theme_name in ["Linux", "LightBG", "Neutral", "NoColor"]:
            warnings.warn(
                f"Theme names and color schemes are lowercase in IPython 9.0 use {theme_name.lower()} instead",
                DeprecationWarning,
                stacklevel=2,
            )
            theme_name = theme_name.lower()
        # Whether to call the interactive pdb debugger after printing
        # tracebacks or not
        super().__init__()
        self.call_pdb = call_pdb

        # Output stream to write to.  Note that we store the original value in
        # a private attribute and then make the public ostream a property, so
        # that we can delay accessing sys.stdout until runtime.  The way
        # things are written now, the sys.stdout object is dynamically managed
        # so a reference to it should NEVER be stored statically.  This
        # property approach confines this detail to a single location, and all
        # subclasses can simply access self.ostream for writing.
        self._ostream = ostream

        # Create color table
        self.set_theme_name(theme_name)
        self.debugger_cls = debugger_cls or debugger.Pdb

        if call_pdb:
            self.pdb = self.debugger_cls()
        else:
            self.pdb = None

    def _get_ostream(self) -> Any:
        """Output stream that exceptions are written to.

        Valid values are:

        - None: the default, which means that IPython will dynamically resolve
          to sys.stdout.  This ensures compatibility with most tools, including
          Windows (where plain stdout doesn't recognize ANSI escapes).

        - Any object with 'write' and 'flush' attributes.
        """
        return sys.stdout if self._ostream is None else self._ostream

    def _set_ostream(self, val) -> None:  # type:ignore[no-untyped-def]
        assert val is None or (hasattr(val, "write") and hasattr(val, "flush"))
        self._ostream = val

    ostream = property(_get_ostream, _set_ostream)

    @staticmethod
    def _get_chained_exception(exception_value: Any) -> Any:
        cause = getattr(exception_value, "__cause__", None)
        if cause:
            return cause
        if getattr(exception_value, "__suppress_context__", False):
            return None
        return getattr(exception_value, "__context__", None)

    def get_parts_of_chained_exception(
        self, evalue: BaseException | None
    ) -> Optional[Tuple[type, BaseException, TracebackType]]:
        chained_evalue = self._get_chained_exception(evalue)

        if chained_evalue:
            return (
                chained_evalue.__class__,
                chained_evalue,
                chained_evalue.__traceback__,
            )
        return None

    def prepare_chained_exception_message(
        self, cause: BaseException | None
    ) -> list[list[str]]:
        direct_cause = (
            "\nThe above exception was the direct cause of the following exception:\n"
        )
        exception_during_handling = (
            "\nDuring handling of the above exception, another exception occurred:\n"
        )

        if cause:
            message = [[direct_cause]]
        else:
            message = [[exception_during_handling]]
        return message

    @property
    def has_colors(self) -> bool:
        assert self._theme_name == self._theme_name.lower()
        return self._theme_name != "nocolor"

    def set_theme_name(self, name: str) -> None:
        assert name in theme_table
        assert name.lower() == name
        self._theme_name = name
        # Also set colors of debugger
        if hasattr(self, "pdb") and self.pdb is not None:
            self.pdb.set_theme_name(name)

    def set_colors(self, name: str) -> None:
        """Shorthand access to the color table scheme selector method."""

        # todo emit deprecation
        warnings.warn(
            "set_colors is deprecated since IPython 9.0, use set_theme_name instead",
            DeprecationWarning,
            stacklevel=2,
        )
        self.set_theme_name(name)

    def color_toggle(self) -> None:
        """Toggle between the currently active color scheme and nocolor."""
        if self._theme_name == "nocolor":
            self._theme_name = self._old_theme_name
        else:
            self._old_theme_name = self._theme_name
            self._theme_name = "nocolor"

    def stb2text(self, stb: list[str]) -> str:
        """Convert a structured traceback (a list) to a string."""
        return "\n".join(stb)

    def text(
        self,
        etype: type,
        value: BaseException | None,
        tb: TracebackType | None,
        tb_offset: Optional[int] = None,
        context: int = 5,
    ) -> str:
        """Return formatted traceback.

        Subclasses may override this if they add extra arguments.
        """
        tb_list = self.structured_traceback(etype, value, tb, tb_offset, context)
        return self.stb2text(tb_list)

    def structured_traceback(
        self,
        etype: type,
        evalue: BaseException | None,
        etb: Optional[TracebackType] = None,
        tb_offset: Optional[int] = None,
        context: int = 5,
    ) -> list[str]:
        """Return a list of traceback frames.

        Must be implemented by each class.
        """
        raise NotImplementedError()
