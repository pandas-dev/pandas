import inspect
import linecache
import sys
from collections.abc import Sequence
from types import TracebackType
from typing import Any, Optional
from collections.abc import Callable

import stack_data
from pygments.formatters.terminal256 import Terminal256Formatter
from pygments.token import Token

from IPython.utils.PyColorize import Theme, TokenStream, theme_table
from IPython.utils.terminal import get_terminal_size

from .tbtools import (
    FrameInfo,
    TBTools,
    _safe_string,
    _tokens_filename,
    eqrepr,
    get_line_number_of_frame,
    nullrepr,
)

INDENT_SIZE = 8


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
    tokens: TokenStream = [(Token, "\n")]

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
                (Token.Prompt, theme.make_arrow(3)),
                (Token, " "),
                (Token, line),
            ]
        else:
            # num = "%*s" % (numbers_width, lineno)
            toks = [
                # (Token.LinenoEm, str(num)),
                (Token, "..."),
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


class DocTB(TBTools):
    """

    A stripped down version of Verbose TB, simplified to not have too much information when
    running doctests

    """

    tb_highlight = ""
    tb_highlight_style = "default"
    tb_offset: int
    long_header: bool
    include_vars: bool

    _mode: str

    def __init__(
        self,
        # TODO: no default ?
        theme_name: str = "linux",
        call_pdb: bool = False,
        ostream: Any = None,
        tb_offset: int = 0,
        long_header: bool = False,
        include_vars: bool = True,
        check_cache: Callable[[], None] | None = None,
        debugger_cls: type | None = None,
    ):
        """Specify traceback offset, headers and color scheme.

        Define how many frames to drop from the tracebacks. Calling it with
        tb_offset=1 allows use of this handler in interpreters which will have
        their own code at the top of the traceback (VerboseTB will first
        remove that frame before printing the traceback info)."""
        assert isinstance(theme_name, str)
        super().__init__(
            theme_name=theme_name,
            call_pdb=call_pdb,
            ostream=ostream,
            debugger_cls=debugger_cls,
        )
        self.tb_offset = tb_offset
        self.long_header = long_header
        self.include_vars = include_vars
        # By default we use linecache.checkcache, but the user can provide a
        # different check_cache implementation.  This was formerly used by the
        # IPython kernel for interactive code, but is no longer necessary.
        if check_cache is None:
            check_cache = linecache.checkcache
        self.check_cache = check_cache

        self.skip_hidden = True

    def format_record(self, frame_info: FrameInfo) -> str:
        """Format a single stack frame"""
        assert isinstance(frame_info, FrameInfo)

        if isinstance(frame_info._sd, stack_data.RepeatedFrames):
            return theme_table[self._theme_name].format(
                [
                    (Token, "    "),
                    (
                        Token.ExcName,
                        "[... skipping similar frames: %s]" % frame_info.description,
                    ),
                    (Token, "\n"),
                ]
            )

        indent: str = " " * INDENT_SIZE

        assert isinstance(frame_info.lineno, int)
        args, varargs, varkw, locals_ = inspect.getargvalues(frame_info.frame)
        if frame_info.executing is not None:
            func = frame_info.executing.code_qualname()
        else:
            func = "?"
        if func == "<module>":
            call = ""
        else:
            # Decide whether to include variable details or not
            var_repr = eqrepr if self.include_vars else nullrepr
            try:
                scope = inspect.formatargvalues(
                    args, varargs, varkw, locals_, formatvalue=var_repr
                )
                assert isinstance(scope, str)
                call = theme_table[self._theme_name].format(
                    [(Token, "in "), (Token.VName, func), (Token.ValEm, scope)]
                )
            except KeyError:
                # This happens in situations like errors inside generator
                # expressions, where local variables are listed in the
                # line, but can't be extracted from the frame.  I'm not
                # 100% sure this isn't actually a bug in inspect itself,
                # but since there's no info for us to compute with, the
                # best we can do is report the failure and move on.  Here
                # we must *not* call any traceback construction again,
                # because that would mess up use of %debug later on.  So we
                # simply report the failure and move on.  The only
                # limitation will be that this frame won't have locals
                # listed in the call signature.  Quite subtle problem...
                # I can't think of a good way to validate this in a unit
                # test, but running a script consisting of:
                #  dict( (k,v.strip()) for (k,v) in range(10) )
                # will illustrate the error, if this exception catch is
                # disabled.
                call = theme_table[self._theme_name].format(
                    [
                        (Token, "in "),
                        (Token.VName, func),
                        (Token.ValEm, "(***failed resolving arguments***)"),
                    ]
                )

        lvals_toks: list[TokenStream] = []
        if self.include_vars:
            try:
                # we likely want to fix stackdata at some point, but
                # still need a workaround.
                fibp = frame_info.variables_in_executing_piece
                for var in fibp:
                    lvals_toks.append(
                        [
                            (Token, var.name),
                            (Token, " "),
                            (Token.ValEm, "= "),
                            (Token.ValEm, repr(var.value)),
                        ]
                    )
            except Exception:
                lvals_toks.append(
                    [
                        (
                            Token,
                            "Exception trying to inspect frame. No more locals available.",
                        ),
                    ]
                )

        assert frame_info._sd is not None
        result = theme_table[self._theme_name].format(
            _tokens_filename(True, frame_info.filename, lineno=frame_info.lineno)
        )
        result += ", " if call else ""
        result += f"{call}\n"
        result += theme_table[self._theme_name].format(
            _format_traceback_lines(
                frame_info.lines,
                theme_table[self._theme_name],
                self.has_colors,
                lvals_toks,
            )
        )
        return result

    def prepare_header(self, etype: str) -> str:
        width = min(75, get_terminal_size()[0])
        head = theme_table[self._theme_name].format(
            [
                (
                    Token,
                    "Traceback (most recent call last):",
                ),
                (Token, " "),
            ]
        )

        return head

    def format_exception(self, etype: Any, evalue: Any) -> Any:
        # Get (safely) a string form of the exception info
        try:
            etype_str, evalue_str = map(str, (etype, evalue))
        except:
            # User exception is improperly defined.
            etype, evalue = str, sys.exc_info()[:2]
            etype_str, evalue_str = map(str, (etype, evalue))

        # PEP-678 notes
        notes = getattr(evalue, "__notes__", [])
        if not isinstance(notes, Sequence) or isinstance(notes, (str, bytes)):
            notes = [_safe_string(notes, "__notes__", func=repr)]

        # ... and format it
        return [
            theme_table[self._theme_name].format(
                [(Token.ExcName, etype_str), (Token, ": "), (Token, evalue_str)]
            ),
            *(
                theme_table[self._theme_name].format([(Token, _safe_string(n, "note"))])
                for n in notes
            ),
        ]

    def format_exception_as_a_whole(
        self,
        etype: type,
        evalue: Optional[BaseException],
        etb: Optional[TracebackType],
        context: int,
        tb_offset: Optional[int],
    ) -> list[list[str]]:
        """Formats the header, traceback and exception message for a single exception.

        This may be called multiple times by Python 3 exception chaining
        (PEP 3134).
        """
        # some locals
        orig_etype = etype
        try:
            etype = etype.__name__  # type: ignore[assignment]
        except AttributeError:
            pass

        tb_offset = self.tb_offset if tb_offset is None else tb_offset
        assert isinstance(tb_offset, int)
        head = self.prepare_header(str(etype))
        assert context == 1, context
        records = self.get_records(etb, context, tb_offset) if etb else []

        frames = []
        skipped = 0
        nskipped = len(records) - 1
        frames.append(self.format_record(records[0]))
        if nskipped:
            frames.append(
                theme_table[self._theme_name].format(
                    [
                        (Token, "\n"),
                        (Token, "    "),
                        (Token, "[... %s skipped frames]" % nskipped),
                        (Token, "\n"),
                        (Token, "\n"),
                    ]
                )
            )

        formatted_exception = self.format_exception(etype, evalue)
        return [[head] + frames + formatted_exception]

    def get_records(self, etb: TracebackType, context: int, tb_offset: int) -> Any:
        assert context == 1, context
        assert etb is not None
        context = context - 1
        after = context // 2
        before = context - after
        if self.has_colors:
            base_style = theme_table[self._theme_name].as_pygments_style()
            style = stack_data.style_with_executing_node(base_style, self.tb_highlight)
            formatter = Terminal256Formatter(style=style)
        else:
            formatter = None
        options = stack_data.Options(
            before=before,
            after=after,
            pygments_formatter=formatter,
        )

        # Let's estimate the amount of code we will have to parse/highlight.
        cf: Optional[TracebackType] = etb
        max_len = 0
        tbs = []
        while cf is not None:
            try:
                mod = inspect.getmodule(cf.tb_frame)
                if mod is not None:
                    mod_name = mod.__name__
                    root_name, *_ = mod_name.split(".")
                    if root_name == "IPython":
                        cf = cf.tb_next
                        continue
                max_len = get_line_number_of_frame(cf.tb_frame)

            except OSError:
                max_len = 0
            max_len = max(max_len, max_len)
            tbs.append(cf)
            cf = getattr(cf, "tb_next", None)

        res = list(stack_data.FrameInfo.stack_data(etb, options=options))[tb_offset:]
        res2 = [FrameInfo._from_stack_data_FrameInfo(r) for r in res]
        return res2

    def structured_traceback(
        self,
        etype: type,
        evalue: Optional[BaseException],
        etb: Optional[TracebackType] = None,
        tb_offset: Optional[int] = None,
        context: int = 1,
    ) -> list[str]:
        """Return a nice text document describing the traceback."""
        assert context > 0
        assert context == 1, context
        formatted_exceptions: list[list[str]] = self.format_exception_as_a_whole(
            etype, evalue, etb, context, tb_offset
        )

        termsize = min(75, get_terminal_size()[0])
        theme = theme_table[self._theme_name]
        structured_traceback_parts: list[str] = []
        chained_exceptions_tb_offset = 0
        lines_of_context = 3
        exception = self.get_parts_of_chained_exception(evalue)
        if exception:
            assert evalue is not None
            formatted_exceptions += self.prepare_chained_exception_message(
                evalue.__cause__
            )
            etype, evalue, etb = exception
        else:
            evalue = None
        chained_exc_ids = set()
        while evalue:
            formatted_exceptions += self.format_exception_as_a_whole(
                etype, evalue, etb, lines_of_context, chained_exceptions_tb_offset
            )
            exception = self.get_parts_of_chained_exception(evalue)

            if exception and id(exception[1]) not in chained_exc_ids:
                chained_exc_ids.add(
                    id(exception[1])
                )  # trace exception to avoid infinite 'cause' loop
                formatted_exceptions += self.prepare_chained_exception_message(
                    evalue.__cause__
                )
                etype, evalue, etb = exception
            else:
                evalue = None

        # we want to see exceptions in a reversed order:
        # the first exception should be on top
        for fx in reversed(formatted_exceptions):
            structured_traceback_parts += fx

        return structured_traceback_parts

    def debugger(self, force: bool = False) -> None:
        raise RuntimeError("canot rundebugger in Docs mode")

    def handler(self, info: tuple[Any, Any, Any] | None = None) -> None:
        (etype, evalue, etb) = info or sys.exc_info()
        self.tb = etb
        ostream = self.ostream
        ostream.flush()
        ostream.write(self.text(etype, evalue, etb))  # type:ignore[arg-type]
        ostream.write("\n")
        ostream.flush()

    # Changed so an instance can just be called as VerboseTB_inst() and print
    # out the right info on its own.
    def __call__(self, etype: Any = None, evalue: Any = None, etb: Any = None) -> None:
        """This hook can replace sys.excepthook (for Python 2.1 or higher)."""
        if etb is None:
            self.handler()
        else:
            self.handler((etype, evalue, etb))
        try:
            self.debugger()
        except KeyboardInterrupt:
            print("\nKeyboardInterrupt")
