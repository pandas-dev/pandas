from __future__ import annotations

import os.path
import sys
import traceback
from collections import defaultdict
from collections.abc import Callable, Iterable, Iterator
from itertools import chain
from typing import Final, Literal, NoReturn, TextIO, TypeVar
from typing_extensions import Self

from librt.internal import (
    ReadBuffer,
    WriteBuffer,
    read_bool,
    read_int as read_int_bare,
    write_bool,
    write_int as write_int_bare,
)

from mypy import errorcodes as codes
from mypy.cache import (
    ErrorTuple,
    read_int,
    read_int_list,
    read_str,
    read_str_opt,
    write_int,
    write_int_list,
    write_str,
    write_str_opt,
)
from mypy.error_formatter import ErrorFormatter
from mypy.errorcodes import IMPORT, IMPORT_NOT_FOUND, IMPORT_UNTYPED, ErrorCode, mypy_error_codes
from mypy.nodes import Context
from mypy.options import Options
from mypy.scope import Scope
from mypy.types import Type
from mypy.util import DEFAULT_SOURCE_OFFSET
from mypy.version import __version__ as mypy_version

T = TypeVar("T")

# Show error codes for some note-level messages (these usually appear alone
# and not as a comment for a previous error-level message).
SHOW_NOTE_CODES: Final = {codes.ANNOTATION_UNCHECKED.code, codes.DEPRECATED.code}

# Do not add notes with links to error code docs to errors with these codes.
# We can tweak this set as we get more experience about what is helpful and what is not.
HIDE_LINK_CODES: Final = {
    # This is a generic error code, so it has no useful docs
    codes.MISC,
    # These are trivial and have some custom notes (e.g. for list being invariant)
    codes.ASSIGNMENT,
    codes.ARG_TYPE,
    codes.RETURN_VALUE,
    # Undefined name/attribute errors are self-explanatory
    codes.ATTR_DEFINED,
    codes.NAME_DEFINED,
    # Overrides have a custom link to docs
    codes.OVERRIDE,
}

BASE_RTD_URL: Final = "https://mypy.rtfd.io/en/stable/_refs.html#code"

# Keep track of the original error code when the error code of a message is changed.
# This is used to give notes about out-of-date "type: ignore" comments.
original_error_codes: Final = {codes.LITERAL_REQ: codes.MISC, codes.TYPE_ABSTRACT: codes.MISC}


class ErrorInfo:
    """Representation of a single error message."""

    # Description of a sequence of imports that refer to the source file
    # related to this error. Each item is a (path, line number) tuple.
    import_ctx: list[tuple[str, int]]
    # Type and function/method where this error occurred. Unqualified, may be None.
    local_ctx: tuple[str | None, str | None]

    # The line number related to this error within file.
    line = 0  # -1 if unknown
    # The column number related to this error with file.
    column = 0  # -1 if unknown
    # The end line number related to this error within file.
    end_line = 0  # -1 if unknown
    # The end column number related to this error with file.
    end_column = 0  # -1 if unknown
    # Either 'error' or 'note'
    severity = ""
    # The error message.
    message = ""
    # The error code.
    code: ErrorCode | None = None

    # If True, we should halt build after the file that generated this error.
    blocker = False
    # Only report this particular messages once per program.
    only_once = False

    # These two are used by the daemon:
    # The fully-qualified id of the source module for this error.
    module: str | None
    # Fine-grained incremental target where this was reported
    target: str | None

    # Lines where `type: ignores` will have effect on this error, for most errors
    # this is just [line]. But sometimes may be custom, e.g. for override errors
    # in methods with multi-line definition.
    origin_span: Iterable[int]
    # For errors on the same line you can use this to customize their sorting
    # (lower value means show first).
    priority: int
    # If True, don't show this message in output, but still record the error.
    hidden = False

    # For notes, specifies (optionally) the error this note is attached to. This is used to
    # simplify error code matching and de-duplication logic for complex multi-line notes.
    parent_error: ErrorInfo | None = None

    def __init__(
        self,
        *,
        import_ctx: list[tuple[str, int]],
        local_ctx: tuple[str | None, str | None],
        line: int,
        column: int,
        end_line: int,
        end_column: int,
        severity: str,
        message: str,
        code: ErrorCode | None,
        blocker: bool,
        only_once: bool,
        module: str | None,
        target: str | None,
        origin_span: Iterable[int] | None = None,
        priority: int = 0,
        parent_error: ErrorInfo | None = None,
    ) -> None:
        self.import_ctx = import_ctx
        self.module = module
        self.local_ctx = local_ctx
        self.line = line
        self.column = column
        self.end_line = end_line
        self.end_column = end_column
        self.severity = severity
        self.message = message
        self.code = code
        self.blocker = blocker
        self.only_once = only_once
        self.origin_span = origin_span or [line]
        self.target = target
        self.priority = priority
        if parent_error is not None:
            assert severity == "note", "Only notes can specify parent errors"
        self.parent_error = parent_error

    def write(self, buf: WriteBuffer) -> None:
        assert self.parent_error is None, "Parent errors not supported yet"
        write_int_bare(buf, len(self.import_ctx))
        for file, line in self.import_ctx:
            write_str(buf, file)
            write_int(buf, line)
        type, function = self.local_ctx
        write_str_opt(buf, type)
        write_str_opt(buf, function)
        write_int(buf, self.line)
        write_int(buf, self.column)
        write_int(buf, self.end_line)
        write_int(buf, self.end_column)
        write_str(buf, self.severity)
        write_str(buf, self.message)
        write_str_opt(buf, self.code.code if self.code else None)
        write_bool(buf, self.blocker)
        write_bool(buf, self.only_once)
        write_str_opt(buf, self.module)
        write_str_opt(buf, self.target)
        write_int_list(buf, list(self.origin_span))
        write_int(buf, self.priority)

    @classmethod
    def read(cls, buf: ReadBuffer) -> ErrorInfo:
        return ErrorInfo(
            import_ctx=[(read_str(buf), read_int(buf)) for _ in range(read_int_bare(buf))],
            local_ctx=(read_str_opt(buf), read_str_opt(buf)),
            line=read_int(buf),
            column=read_int(buf),
            end_line=read_int(buf),
            end_column=read_int(buf),
            severity=read_str(buf),
            message=read_str(buf),
            code=mypy_error_codes[code] if (code := read_str_opt(buf)) else None,
            blocker=read_bool(buf),
            only_once=read_bool(buf),
            module=read_str_opt(buf),
            target=read_str_opt(buf),
            origin_span=read_int_list(buf),
            priority=read_int(buf),
        )


class ErrorWatcher:
    """Context manager that can be used to keep track of new errors recorded
    around a given operation.

    Errors maintain a stack of such watchers. The handler is called starting
    at the top of the stack, and is propagated down the stack unless filtered
    out by one of the ErrorWatcher instances.
    """

    # public attribute for the special treatment of `reveal_type` by
    # `MessageBuilder.reveal_type`:
    filter_revealed_type: bool

    def __init__(
        self,
        errors: Errors,
        *,
        filter_errors: bool | Callable[[str, ErrorInfo], bool] = False,
        save_filtered_errors: bool = False,
        filter_deprecated: bool = False,
        filter_revealed_type: bool = False,
    ) -> None:
        self.errors = errors
        self._has_new_errors = False
        self._filter = filter_errors
        self._filter_deprecated = filter_deprecated
        self.filter_revealed_type = filter_revealed_type
        self._filtered: list[ErrorInfo] | None = [] if save_filtered_errors else None

    def __enter__(self) -> Self:
        self.errors._watchers.append(self)
        return self

    def __exit__(self, exc_type: object, exc_val: object, exc_tb: object) -> Literal[False]:
        last = self.errors._watchers.pop()
        assert last == self
        return False

    def on_error(self, file: str, info: ErrorInfo) -> bool:
        """Handler called when a new error is recorded.

        The default implementation just sets the has_new_errors flag

        Return True to filter out the error, preventing it from being seen by other
        ErrorWatcher further down the stack and from being recorded by Errors
        """
        if info.code == codes.DEPRECATED:
            # Deprecated is not a type error, so it is handled on opt-in basis here.
            if not self._filter_deprecated:
                return False

        self._has_new_errors = True
        if isinstance(self._filter, bool):
            should_filter = self._filter
        elif callable(self._filter):
            should_filter = self._filter(file, info)
        else:
            raise AssertionError(f"invalid error filter: {type(self._filter)}")
        if should_filter and self._filtered is not None:
            self._filtered.append(info)

        return should_filter

    def has_new_errors(self) -> bool:
        return self._has_new_errors

    def filtered_errors(self) -> list[ErrorInfo]:
        assert self._filtered is not None
        return self._filtered


class NonOverlapErrorInfo:
    line: int
    column: int
    end_line: int | None
    end_column: int | None
    kind: str

    def __init__(
        self, *, line: int, column: int, end_line: int | None, end_column: int | None, kind: str
    ) -> None:
        self.line = line
        self.column = column
        self.end_line = end_line
        self.end_column = end_column
        self.kind = kind

    def __eq__(self, other: object) -> bool:
        if isinstance(other, NonOverlapErrorInfo):
            return (
                self.line == other.line
                and self.column == other.column
                and self.end_line == other.end_line
                and self.end_column == other.end_column
                and self.kind == other.kind
            )
        return False

    def __hash__(self) -> int:
        return hash((self.line, self.column, self.end_line, self.end_column, self.kind))


class IterationDependentErrors:
    """An `IterationDependentErrors` instance serves to collect the `unreachable`,
    `redundant-expr`, and `redundant-casts` errors, as well as the revealed types and
    non-overlapping types, handled by the individual `IterationErrorWatcher` instances
    sequentially applied to the same code section."""

    # One set of `unreachable`, `redundant-expr`, and `redundant-casts` errors per
    # iteration step.  Meaning of the tuple items: ErrorCode, message, line, column,
    # end_line, end_column.
    uselessness_errors: list[set[tuple[ErrorCode, str, int, int, int, int]]]

    # One set of unreachable line numbers per iteration step.  Not only the lines where
    # the error report occurs but really all unreachable lines.
    unreachable_lines: list[set[int]]

    # One list of revealed types for each `reveal_type` statement.  Each created list
    # can grow during the iteration.  Meaning of the tuple items: line, column,
    # end_line, end_column:
    revealed_types: dict[tuple[int, int, int | None, int | None], list[Type]]

    # One dictionary of non-overlapping types per iteration step:
    nonoverlapping_types: list[dict[NonOverlapErrorInfo, tuple[Type, Type]]]

    def __init__(self) -> None:
        self.uselessness_errors = []
        self.unreachable_lines = []
        self.nonoverlapping_types = []
        self.revealed_types = defaultdict(list)

    def yield_uselessness_error_infos(self) -> Iterator[tuple[str, Context, ErrorCode]]:
        """Report only those `unreachable`, `redundant-expr`, and `redundant-casts`
        errors that could not be ruled out in any iteration step."""

        persistent_uselessness_errors = set()
        for candidate in set(chain(*self.uselessness_errors)):
            if all(
                (candidate in errors) or (candidate[2] in lines)
                for errors, lines in zip(self.uselessness_errors, self.unreachable_lines)
            ):
                persistent_uselessness_errors.add(candidate)
        for error_info in persistent_uselessness_errors:
            context = Context(line=error_info[2], column=error_info[3])
            context.end_line = error_info[4]
            context.end_column = error_info[5]
            yield error_info[1], context, error_info[0]

    def yield_nonoverlapping_types(
        self,
    ) -> Iterator[tuple[tuple[list[Type], list[Type]], str, Context]]:
        """Report expressions where non-overlapping types were detected for all iterations
        were the expression was reachable."""

        selected = set()
        for candidate in set(chain.from_iterable(self.nonoverlapping_types)):
            if all(
                (candidate in nonoverlap) or (candidate.line in lines)
                for nonoverlap, lines in zip(self.nonoverlapping_types, self.unreachable_lines)
            ):
                selected.add(candidate)

        persistent_nonoverlaps: dict[NonOverlapErrorInfo, tuple[list[Type], list[Type]]] = (
            defaultdict(lambda: ([], []))
        )
        for nonoverlaps in self.nonoverlapping_types:
            for candidate, (left, right) in nonoverlaps.items():
                if candidate in selected:
                    types = persistent_nonoverlaps[candidate]
                    types[0].append(left)
                    types[1].append(right)

        for error_info, types in persistent_nonoverlaps.items():
            context = Context(line=error_info.line, column=error_info.column)
            context.end_line = error_info.end_line
            context.end_column = error_info.end_column
            yield (types[0], types[1]), error_info.kind, context

    def yield_revealed_type_infos(self) -> Iterator[tuple[list[Type], Context]]:
        """Yield all types revealed in at least one iteration step."""

        for note_info, types in self.revealed_types.items():
            context = Context(line=note_info[0], column=note_info[1])
            context.end_line = note_info[2]
            context.end_column = note_info[3]
            yield types, context


class IterationErrorWatcher(ErrorWatcher):
    """Error watcher that filters and separately collects `unreachable` errors,
    `redundant-expr` and `redundant-casts` errors, and revealed types and
    non-overlapping types when analysing code sections iteratively to help avoid
    making too-hasty reports."""

    iteration_dependent_errors: IterationDependentErrors

    def __init__(
        self,
        errors: Errors,
        iteration_dependent_errors: IterationDependentErrors,
        *,
        filter_errors: bool | Callable[[str, ErrorInfo], bool] = False,
        save_filtered_errors: bool = False,
        filter_deprecated: bool = False,
    ) -> None:
        super().__init__(
            errors,
            filter_errors=filter_errors,
            save_filtered_errors=save_filtered_errors,
            filter_deprecated=filter_deprecated,
        )
        self.iteration_dependent_errors = iteration_dependent_errors
        iteration_dependent_errors.uselessness_errors.append(set())
        iteration_dependent_errors.nonoverlapping_types.append({})
        iteration_dependent_errors.unreachable_lines.append(set())

    def on_error(self, file: str, info: ErrorInfo) -> bool:
        """Filter out the "iteration-dependent" errors and notes and store their
        information to handle them after iteration is completed."""

        iter_errors = self.iteration_dependent_errors

        if info.code in (codes.UNREACHABLE, codes.REDUNDANT_EXPR, codes.REDUNDANT_CAST):
            iter_errors.uselessness_errors[-1].add(
                (info.code, info.message, info.line, info.column, info.end_line, info.end_column)
            )
            if info.code == codes.UNREACHABLE:
                iter_errors.unreachable_lines[-1].update(range(info.line, info.end_line + 1))
            return True

        return super().on_error(file, info)


class Errors:
    """Container for compile errors.

    This class generates and keeps tracks of compile errors and the
    current error context (nested imports).
    """

    # Map from files to generated error messages. Is an OrderedDict so
    # that it can be used to order messages based on the order the
    # files were processed.
    error_info_map: dict[str, list[ErrorInfo]]

    # optimization for legacy codebases with many files with errors
    has_blockers: set[str]

    # Files that we have reported the errors for
    flushed_files: set[str]

    # Current error context: nested import context/stack, as a list of (path, line) pairs.
    import_ctx: list[tuple[str, int]]

    # Path name prefix that is removed from all paths, if set.
    ignore_prefix: str | None = None

    # Path to current file.
    file: str = ""

    # Ignore some errors on these lines of each file
    # (path -> line -> error-codes)
    ignored_lines: dict[str, dict[int, list[str]]]

    # Lines that were skipped during semantic analysis e.g. due to ALWAYS_FALSE, MYPY_FALSE,
    # or platform/version checks. Those lines would not be type-checked.
    skipped_lines: dict[str, set[int]]

    # Lines on which an error was actually ignored.
    used_ignored_lines: dict[str, dict[int, list[str]]]

    # Files where all errors should be ignored.
    ignored_files: set[str]

    # Collection of reported only_once messages.
    only_once_messages: set[str]

    # State for keeping track of the current fine-grained incremental mode target.
    # (See mypy.server.update for more about targets.)
    # Current module id.
    target_module: str | None = None
    scope: Scope | None = None

    # Have we seen an import-related error so far? If yes, we filter out other messages
    # in some cases to avoid reporting huge numbers of errors.
    seen_import_error = False

    # Set this flag to record all raw report() calls. Recorded error (per file) can
    # be replayed using by calling set_file() and add_error_info().
    global_watcher = False
    recorded: dict[str, list[ErrorInfo]]

    _watchers: list[ErrorWatcher]

    def __init__(
        self,
        options: Options,
        *,
        read_source: Callable[[str], list[str] | None] | None = None,
        hide_error_codes: bool | None = None,
        error_formatter: ErrorFormatter | None = None,
    ) -> None:
        self.options = options
        self.hide_error_codes = (
            hide_error_codes if hide_error_codes is not None else options.hide_error_codes
        )
        # We use fscache to read source code when showing snippets.
        self.read_source = read_source
        self.error_formatter = error_formatter
        self.initialize()

    def initialize(self) -> None:
        self.error_info_map = {}
        self.flushed_files = set()
        self.import_ctx = []
        self.function_or_member = [None]
        self.ignored_lines = {}
        self.skipped_lines = {}
        self.used_ignored_lines = defaultdict(lambda: defaultdict(list))
        self.ignored_files = set()
        self.only_once_messages = set()
        self.has_blockers = set()
        self.scope = None
        self.target_module = None
        self.seen_import_error = False
        self._watchers = []
        self.global_watcher = False
        self.recorded = defaultdict(list)

    def reset(self) -> None:
        self.initialize()

    def set_ignore_prefix(self, prefix: str) -> None:
        """Set path prefix that will be removed from all paths."""
        prefix = os.path.normpath(prefix)
        # Add separator to the end, if not given.
        if os.path.basename(prefix) != "":
            prefix += os.sep
        self.ignore_prefix = prefix

    def simplify_path(self, file: str) -> str:
        if self.options.show_absolute_path:
            return os.path.abspath(file)
        else:
            file = os.path.normpath(file)
            return remove_path_prefix(file, self.ignore_prefix)

    def set_file(
        self, file: str, module: str | None, options: Options, scope: Scope | None = None
    ) -> None:
        """Set the path and module id of the current file."""
        # The path will be simplified later, in render_messages. That way
        #  * 'file' is always a key that uniquely identifies a source file
        #    that mypy read (simplified paths might not be unique); and
        #  * we only have to simplify in one place, while still supporting
        #    reporting errors for files other than the one currently being
        #    processed.
        self.file = file
        self.target_module = module
        self.scope = scope
        self.options = options

    def set_file_ignored_lines(
        self, file: str, ignored_lines: dict[int, list[str]], ignore_all: bool = False
    ) -> None:
        self.ignored_lines[file] = ignored_lines
        if ignore_all:
            self.ignored_files.add(file)

    def set_skipped_lines(self, file: str, skipped_lines: set[int]) -> None:
        self.skipped_lines[file] = skipped_lines

    def current_target(self) -> str | None:
        """Retrieves the current target from the associated scope.

        If there is no associated scope, use the target module."""
        if self.scope is not None:
            return self.scope.current_target()
        return self.target_module

    def current_module(self) -> str | None:
        return self.target_module

    def import_context(self) -> list[tuple[str, int]]:
        """Return a copy of the import context."""
        return self.import_ctx.copy()

    def set_import_context(self, ctx: list[tuple[str, int]]) -> None:
        """Replace the entire import context with a new value."""
        self.import_ctx = ctx.copy()

    def report(
        self,
        line: int,
        column: int | None,
        message: str,
        code: ErrorCode | None = None,
        *,
        blocker: bool = False,
        severity: str = "error",
        only_once: bool = False,
        origin_span: Iterable[int] | None = None,
        offset: int = 0,
        end_line: int | None = None,
        end_column: int | None = None,
        parent_error: ErrorInfo | None = None,
    ) -> ErrorInfo:
        """Report message at the given line using the current error context.

        Args:
            line: line number of error
            column: column number of error
            message: message to report
            code: error code (defaults to 'misc'; not shown for notes)
            blocker: if True, don't continue analysis after this error
            severity: 'error' or 'note'
            only_once: if True, only report this exact message once per build
            origin_span: lines where `type: ignore`s have effect for this error
                (default is [line])
            offset: number of spaces to prefix this message
            end_line: if known, end line of error location
            end_column: if known, end column of error location
            parent_error: an error this note is attached to (for notes only).
        """
        if self.scope:
            type = self.scope.current_type_name()
            if self.scope.ignored > 0:
                type = None  # Omit type context if nested function
            function = self.scope.current_function_name()
        else:
            type = None
            function = None

        # It looks like there is a bug in how we parse f-strings,
        # we cannot simply assert this yet.
        if end_line is None or end_line < line:
            end_line = line

        if column is None:
            column = -1
        if end_column is None:
            if column == -1:
                end_column = -1
            else:
                end_column = column + 1
        if line == end_line and end_column <= column:
            # Be defensive, similar to the logic for lines above.
            end_column = column + 1

        if offset:
            message = " " * offset + message

        code = code or (parent_error.code if parent_error else None)
        if parent_error is not None:
            assert code == parent_error.code, "Must have same error code as parent"
            assert severity == "note", "Only notes can have parent errors"
        code = code or (codes.MISC if not blocker else None)

        info = ErrorInfo(
            import_ctx=self.import_context(),
            local_ctx=(type, function),
            line=line,
            column=column,
            end_line=end_line,
            end_column=end_column,
            severity=severity,
            message=message,
            code=code,
            blocker=blocker,
            only_once=only_once,
            origin_span=origin_span,
            module=self.current_module(),
            target=self.current_target(),
            parent_error=parent_error,
        )
        if self.global_watcher:
            self.recorded[self.file].append(info)
        self.add_error_info(info)
        return info

    def _add_error_info(self, file: str, info: ErrorInfo) -> None:
        assert file not in self.flushed_files
        # process the stack of ErrorWatchers before modifying any internal state
        # in case we need to filter out the error entirely
        if self._filter_error(file, info):
            return
        if file not in self.error_info_map:
            self.error_info_map[file] = []
        self.error_info_map[file].append(info)
        if info.blocker:
            self.has_blockers.add(file)
        if info.code in (IMPORT, IMPORT_UNTYPED, IMPORT_NOT_FOUND):
            self.seen_import_error = True

    def note_for_info(
        self,
        file: str,
        info: ErrorInfo,
        message: str,
        code: ErrorCode | None,
        *,
        only_once: bool = False,
        priority: int = 0,
    ) -> None:
        """Generate an additional note for an existing ErrorInfo.

        This skip the logic in add_error_info() and goes to _add_error_info().
        """
        info = ErrorInfo(
            import_ctx=info.import_ctx,
            local_ctx=info.local_ctx,
            line=info.line,
            column=info.column,
            end_line=info.end_line,
            end_column=info.end_column,
            severity="note",
            message=message,
            code=code,
            blocker=False,
            only_once=only_once,
            module=info.module,
            target=info.target,
            origin_span=info.origin_span,
            priority=priority,
        )
        self._add_error_info(file, info)

    def report_simple_error(
        self, file: str, line: int, message: str, code: ErrorCode | None
    ) -> None:
        """Generate a simple error in a module.

        This skip the logic in add_error_info() and goes to _add_error_info().
        """
        info = ErrorInfo(
            import_ctx=self.import_context(),
            local_ctx=(None, None),
            line=line,
            column=-1,
            end_line=line,
            end_column=-1,
            severity="error",
            message=message,
            code=code,
            blocker=False,
            only_once=False,
            module=self.current_module(),
            # TODO: can we support more precise targets?
            target=self.target_module,
        )
        self._add_error_info(file, info)

    def get_watchers(self) -> Iterator[ErrorWatcher]:
        """Yield the `ErrorWatcher` stack from top to bottom."""
        i = len(self._watchers)
        while i > 0:
            i -= 1
            yield self._watchers[i]

    def _filter_error(self, file: str, info: ErrorInfo) -> bool:
        """
        process ErrorWatcher stack from top to bottom,
        stopping early if error needs to be filtered out
        """
        return any(w.on_error(file, info) for w in self.get_watchers())

    def add_error_info(self, info: ErrorInfo, *, file: str | None = None) -> None:
        lines = info.origin_span
        file = file or self.file
        # process the stack of ErrorWatchers before modifying any internal state
        # in case we need to filter out the error entirely
        # NB: we need to do this both here and in _add_error_info, otherwise we
        # might incorrectly update the sets of ignored or only_once messages
        if self._filter_error(file, info):
            return
        if not info.blocker:  # Blockers cannot be ignored
            if file in self.ignored_lines:
                # Check each line in this context for "type: ignore" comments.
                # line == end_line for most nodes, so we only loop once.
                for scope_line in lines:
                    if self.is_ignored_error(scope_line, info, self.ignored_lines[file]):
                        err_code = info.code or codes.MISC
                        if not self.is_error_code_enabled(err_code):
                            # Error code is disabled - don't mark the current
                            # "type: ignore" comment as used.
                            return
                        # Annotation requests us to ignore all errors on this line.
                        self.used_ignored_lines[file][scope_line].append(err_code.code)
                        return
            if file in self.ignored_files:
                return
        if info.only_once:
            if info.message in self.only_once_messages:
                return
            self.only_once_messages.add(info.message)
        if (
            self.seen_import_error
            and info.code not in (IMPORT, IMPORT_UNTYPED, IMPORT_NOT_FOUND)
            and self.has_many_errors()
        ):
            # Missing stubs can easily cause thousands of errors about
            # Any types, especially when upgrading to mypy 0.900,
            # which no longer bundles third-party library stubs. Avoid
            # showing too many errors to make it easier to see
            # import-related errors.
            info.hidden = True
            self.report_hidden_errors(file, info)
        self._add_error_info(file, info)
        ignored_codes = self.ignored_lines.get(file, {}).get(info.line, [])
        if ignored_codes and info.code:
            # Something is ignored on the line, but not this error, so maybe the error
            # code is incorrect.
            msg = f"""Error code "{info.code.code}" not covered by "type: ignore[{', '.join(ignored_codes)}]" comment"""
            if info.code in original_error_codes:
                # If there seems to be a "type: ignore" with a stale error
                # code, report a more specific note.
                old_code = original_error_codes[info.code].code
                if old_code in ignored_codes:
                    msg = (
                        f'Error code changed to {info.code.code}; "type: ignore" comment '
                        + "may be out of date"
                    )
            self.note_for_info(file, info, msg, None, only_once=False)
        if (
            self.options.show_error_code_links
            and not self.options.hide_error_codes
            and info.code is not None
            and info.code not in HIDE_LINK_CODES
            and info.code.code in mypy_error_codes
        ):
            message = f"See {BASE_RTD_URL}-{info.code.code} for more info"
            if message in self.only_once_messages:
                return
            self.only_once_messages.add(message)
            self.note_for_info(file, info, message, info.code, only_once=True, priority=20)

    def has_many_errors(self) -> bool:
        if self.options.many_errors_threshold < 0:
            return False
        if len(self.error_info_map) >= self.options.many_errors_threshold:
            return True
        if (
            sum(len(errors) for errors in self.error_info_map.values())
            >= self.options.many_errors_threshold
        ):
            return True
        return False

    def report_hidden_errors(self, file: str, info: ErrorInfo) -> None:
        message = (
            "(Skipping most remaining errors due to unresolved imports or missing stubs; "
            + "fix these first)"
        )
        if message in self.only_once_messages:
            return
        self.only_once_messages.add(message)
        self.note_for_info(file, info, message, None, only_once=True)

    def is_ignored_error(self, line: int, info: ErrorInfo, ignores: dict[int, list[str]]) -> bool:
        if info.blocker:
            # Blocking errors can never be ignored
            return False
        if info.code and not self.is_error_code_enabled(info.code):
            return True
        if line not in ignores:
            return False
        if not ignores[line]:
            # Empty list means that we ignore all errors
            return True
        if info.code and self.is_error_code_enabled(info.code):
            return (
                info.code.code in ignores[line]
                or info.code.sub_code_of is not None
                and info.code.sub_code_of.code in ignores[line]
            )
        return False

    def is_error_code_enabled(self, error_code: ErrorCode) -> bool:
        current_mod_disabled = self.options.disabled_error_codes
        current_mod_enabled = self.options.enabled_error_codes

        if error_code in current_mod_disabled:
            return False
        elif error_code in current_mod_enabled:
            return True
        elif error_code.sub_code_of is not None and error_code.sub_code_of in current_mod_disabled:
            return False
        else:
            return error_code.default_enabled

    def clear_errors_in_targets(self, path: str, targets: set[str]) -> None:
        """Remove errors in specific fine-grained targets within a file."""
        if path in self.error_info_map:
            new_errors = []
            has_blocker = False
            for info in self.error_info_map[path]:
                if info.target not in targets:
                    new_errors.append(info)
                    has_blocker |= info.blocker
                elif info.only_once:
                    self.only_once_messages.remove(info.message)
            self.error_info_map[path] = new_errors
            if not has_blocker and path in self.has_blockers:
                self.has_blockers.remove(path)

    def generate_unused_ignore_errors(self, file: str, is_typeshed: bool = False) -> None:
        if is_typeshed or file in self.ignored_files:
            return
        ignored_lines = self.ignored_lines[file]
        used_ignored_lines = self.used_ignored_lines[file]
        for line, ignored_codes in ignored_lines.items():
            if line in self.skipped_lines[file]:
                continue
            if codes.UNUSED_IGNORE.code in ignored_codes:
                continue
            used_ignored_codes = set(used_ignored_lines[line])
            unused_ignored_codes = [c for c in ignored_codes if c not in used_ignored_codes]
            # `ignore` is used
            if not ignored_codes and used_ignored_codes:
                continue
            # All codes appearing in `ignore[...]` are used
            if ignored_codes and not unused_ignored_codes:
                continue
            # Display detail only when `ignore[...]` specifies more than one error code
            unused_codes_message = ""
            if len(ignored_codes) > 1 and unused_ignored_codes:
                unused_codes_message = f"[{', '.join(unused_ignored_codes)}]"
            message = f'Unused "type: ignore{unused_codes_message}" comment'
            for unused in unused_ignored_codes:
                narrower = set(used_ignored_codes) & codes.sub_code_map[unused]
                if narrower:
                    message += f", use narrower [{', '.join(narrower)}] instead of [{unused}] code"
            # Don't use report() since add_error_info will ignore the error!
            self.report_simple_error(file, line, message, code=codes.UNUSED_IGNORE)

    def generate_ignore_without_code_errors(
        self, file: str, is_warning_unused_ignores: bool, is_typeshed: bool = False
    ) -> None:
        if is_typeshed or file in self.ignored_files:
            return

        used_ignored_lines = self.used_ignored_lines[file]
        for line, ignored_codes in self.ignored_lines[file].items():
            if line in self.skipped_lines[file]:
                continue
            if ignored_codes:
                continue

            # If the `type: ignore` is itself unused and that would be warned about,
            # let that error stand alone
            if is_warning_unused_ignores and not used_ignored_lines[line]:
                continue

            codes_hint = ""
            ignored_codes = sorted(set(used_ignored_lines[line]))
            if ignored_codes:
                codes_hint = f' (consider "type: ignore[{", ".join(ignored_codes)}]" instead)'

            message = f'"type: ignore" comment without error code{codes_hint}'
            # Don't use report() since add_error_info will ignore the error!
            self.report_simple_error(file, line, message, code=codes.IGNORE_WITHOUT_CODE)

    def num_messages(self) -> int:
        """Return the number of generated messages."""
        return sum(len(x) for x in self.error_info_map.values())

    def is_errors(self) -> bool:
        """Are there any generated messages?"""
        return bool(self.error_info_map)

    def is_blockers(self) -> bool:
        """Are the any errors that are blockers?"""
        return bool(self.has_blockers)

    def blocker_module(self) -> str | None:
        """Return the module with a blocking error, or None if not possible."""
        for path in self.has_blockers:
            for err in self.error_info_map[path]:
                if err.blocker:
                    return err.module
        return None

    def is_errors_for_file(self, file: str) -> bool:
        """Are there any errors for the given file?"""
        return file in self.error_info_map and file not in self.ignored_files

    def prefer_simple_messages(self) -> bool:
        """Should we generate simple/fast error messages?

        Return True if errors are not shown to user, i.e. errors are ignored
        or they are collected for internal use only.

        If True, we should prefer to generate a simple message quickly.
        All normal errors should still be reported.
        """
        if self.file in self.ignored_files:
            # Errors ignored, so no point generating fancy messages
            return True
        if self.options.ignore_errors:
            return True
        if self._watchers:
            _watcher = self._watchers[-1]
            if _watcher._filter is True and _watcher._filtered is None:
                # Errors are filtered
                return True
        return False

    def raise_error(self, use_stdout: bool = True) -> NoReturn:
        """Raise a CompileError with the generated messages.

        Render the messages suitable for displaying.
        """
        # self.new_messages() will format all messages that haven't already
        # been returned from a file_messages() call.
        raise CompileError(
            self.new_messages(), use_stdout=use_stdout, module_with_blocker=self.blocker_module()
        )

    def format_messages_default(
        self, error_tuples: list[ErrorTuple], source_lines: list[str] | None
    ) -> list[str]:
        """Return a string list that represents the error messages.

        Use a form suitable for displaying to the user. If self.pretty
        is True also append a relevant trimmed source code line (only for
        severity 'error').
        """
        a: list[str] = []
        for file, line, column, end_line, end_column, severity, message, code in error_tuples:
            s = ""
            if file is not None:
                if self.options.show_column_numbers and line >= 0 and column >= 0:
                    srcloc = f"{file}:{line}:{1 + column}"
                    if self.options.show_error_end and end_line >= 0 and end_column >= 0:
                        srcloc += f":{end_line}:{end_column}"
                elif line >= 0:
                    srcloc = f"{file}:{line}"
                else:
                    srcloc = file
                s = f"{srcloc}: {severity}: {message}"
            else:
                s = message
            if (
                not self.hide_error_codes
                and code
                and (severity != "note" or code in SHOW_NOTE_CODES)
            ):
                # If note has an error code, it is related to a previous error. Avoid
                # displaying duplicate error codes.
                s = f"{s}  [{code}]"
            a.append(s)
            if self.options.pretty:
                # Add source code fragment and a location marker.
                if severity == "error" and source_lines and line > 0:
                    source_line = source_lines[line - 1]
                    source_line_expanded = source_line.expandtabs()
                    min_column = len(source_line) - len(source_line.lstrip())
                    if column < min_column:
                        # Something went wrong, take first non-empty column.
                        column = min_column

                    # Shifts column after tab expansion
                    column = len(source_line[:column].expandtabs())
                    end_column = len(source_line[:end_column].expandtabs())

                    # Note, currently coloring uses the offset to detect source snippets,
                    # so these offsets should not be arbitrary.
                    a.append(" " * DEFAULT_SOURCE_OFFSET + source_line_expanded)
                    marker = "^"
                    if end_line == line and end_column > column:
                        marker = f'^{"~" * (end_column - column - 1)}'
                    elif end_line != line:
                        # just highlight the first line instead
                        marker = f'^{"~" * (len(source_line_expanded) - column - 1)}'
                    a.append(" " * (DEFAULT_SOURCE_OFFSET + column) + marker)
        return a

    def file_messages(self, path: str) -> list[ErrorTuple]:
        """Return an error tuple list of new error messages from a given file."""
        if path not in self.error_info_map:
            return []

        error_info = self.error_info_map[path]
        error_info = [info for info in error_info if not info.hidden]
        error_info = self.remove_duplicates(self.sort_messages(error_info))
        return self.render_messages(path, error_info)

    def format_messages(
        self, path: str, error_tuples: list[ErrorTuple], formatter: ErrorFormatter | None = None
    ) -> list[str]:
        """Return a string list of new error messages from a given file.

        Use a form suitable for displaying to the user.
        """
        self.flushed_files.add(path)
        if formatter is not None:
            errors = create_errors(error_tuples)
            return [formatter.report_error(err) for err in errors]

        source_lines = None
        if self.options.pretty and self.read_source:
            # Find shadow file mapping and read source lines if a shadow file exists for the given path.
            # If shadow file mapping is not found, read source lines
            mapped_path = self.find_shadow_file_mapping(path)
            if mapped_path:
                source_lines = self.read_source(mapped_path)
            else:
                source_lines = self.read_source(path)
        return self.format_messages_default(error_tuples, source_lines)

    def find_shadow_file_mapping(self, path: str) -> str | None:
        """Return the shadow file path for a given source file path or None."""
        if self.options.shadow_file is None:
            return None

        for i in self.options.shadow_file:
            if i[0] == path:
                return i[1]
        return None

    def new_messages(self) -> list[str]:
        """Return a string list of new error messages.

        Use a form suitable for displaying to the user.
        Errors from different files are ordered based on the order in which
        they first generated an error.
        """
        msgs = []
        for path in self.error_info_map.keys():
            if path not in self.flushed_files:
                error_tuples = self.file_messages(path)
                msgs.extend(
                    self.format_messages(path, error_tuples, formatter=self.error_formatter)
                )
        return msgs

    def targets(self) -> set[str]:
        """Return a set of all targets that contain errors."""
        # TODO: Make sure that either target is always defined or that not being defined
        #       is okay for fine-grained incremental checking.
        return {
            info.target for errs in self.error_info_map.values() for info in errs if info.target
        }

    def render_messages(self, file: str, errors: list[ErrorInfo]) -> list[ErrorTuple]:
        """Translate the messages into a sequence of tuples.

        Each tuple is of form (path, line, col, severity, message, code).
        The rendered sequence includes information about error contexts.
        The path item may be None. If the line item is negative, the
        line number is not defined for the tuple.
        """
        file = self.simplify_path(file)
        result: list[ErrorTuple] = []
        prev_import_context: list[tuple[str, int]] = []
        prev_function: str | None = None
        prev_type: str | None = None

        for e in errors:
            # Report module import context, if different from previous message.
            if not self.options.show_error_context:
                pass
            elif e.import_ctx != prev_import_context:
                last = len(e.import_ctx) - 1
                i = last
                while i >= 0:
                    path, line = e.import_ctx[i]
                    fmt = "{}:{}: note: In module imported here"
                    if i < last:
                        fmt = "{}:{}: note: ... from here"
                    if i > 0:
                        fmt += ","
                    else:
                        fmt += ":"
                    # Remove prefix to ignore from path (if present) to
                    # simplify path.
                    path = remove_path_prefix(path, self.ignore_prefix)
                    result.append((None, -1, -1, -1, -1, "note", fmt.format(path, line), None))
                    i -= 1

            # Report context within a source file.
            type, function = e.local_ctx
            if not self.options.show_error_context:
                pass
            elif function != prev_function or type != prev_type:
                if function is None:
                    if type is None:
                        result.append((file, -1, -1, -1, -1, "note", "At top level:", None))
                    else:
                        result.append((file, -1, -1, -1, -1, "note", f'In class "{type}":', None))
                else:

                    if type is None:
                        msg = f'In function "{function}":'
                    else:
                        msg = 'In member "{}" of class "{}":'.format(function, type)
                    result.append((file, -1, -1, -1, -1, "note", msg, None))

            elif type != prev_type:
                if type is None:
                    result.append((file, -1, -1, -1, -1, "note", "At top level:", None))
                else:
                    result.append((file, -1, -1, -1, -1, "note", f'In class "{type}":', None))

            code = e.code.code if e.code is not None else None
            result.append(
                (file, e.line, e.column, e.end_line, e.end_column, e.severity, e.message, code)
            )

            prev_import_context = e.import_ctx
            prev_function = function
            prev_type = type

        return result

    def sort_messages(self, errors: list[ErrorInfo]) -> list[ErrorInfo]:
        """Sort an array of error messages locally by line number.

        I.e., sort a run of consecutive messages with the same
        context by line number, but otherwise retain the general
        ordering of the messages.
        """
        result: list[ErrorInfo] = []
        i = 0
        while i < len(errors):
            i0 = i
            # Find neighbouring errors with the same context and file.
            while i + 1 < len(errors) and errors[i + 1].import_ctx == errors[i].import_ctx:
                i += 1
            i += 1

            # Sort the errors specific to a file according to line number and column.
            a = sorted(errors[i0:i], key=lambda x: (x.line, x.column))
            a = self.sort_within_context(a)
            result.extend(a)
        return result

    def sort_within_context(self, errors: list[ErrorInfo]) -> list[ErrorInfo]:
        """For the same location decide which messages to show first/last.

        Currently, we only compare within the same error code, to decide the
        order of various additional notes.
        """
        result = []
        i = 0
        while i < len(errors):
            i0 = i
            # Find neighbouring errors with the same position and error code.
            while (
                i + 1 < len(errors)
                and errors[i + 1].line == errors[i].line
                and errors[i + 1].column == errors[i].column
                and errors[i + 1].end_line == errors[i].end_line
                and errors[i + 1].end_column == errors[i].end_column
                and errors[i + 1].code == errors[i].code
            ):
                i += 1
            i += 1

            # Sort the messages specific to a given error by priority.
            a = sorted(errors[i0:i], key=lambda x: x.priority)
            result.extend(a)
        return result

    def remove_duplicates(self, errors: list[ErrorInfo]) -> list[ErrorInfo]:
        filtered_errors = []
        seen_by_line: defaultdict[int, set[tuple[str, str]]] = defaultdict(set)
        removed = set()
        for err in errors:
            if err.parent_error is not None:
                # Notes with specified parent are removed together with error below.
                filtered_errors.append(err)
            elif (err.severity, err.message) not in seen_by_line[err.line]:
                filtered_errors.append(err)
                seen_by_line[err.line].add((err.severity, err.message))
            else:
                removed.add(err)
        return [
            err
            for err in filtered_errors
            if err.parent_error is None or err.parent_error not in removed
        ]


class CompileError(Exception):
    """Exception raised when there is a compile error.

    It can be a parse, semantic analysis, type check or other
    compilation-related error.

    CompileErrors raised from an errors object carry all of the
    messages that have not been reported out by error streaming.
    This is patched up by build.build to contain either all error
    messages (if errors were streamed) or none (if they were not).

    """

    messages: list[str]
    use_stdout = False
    # Can be set in case there was a module with a blocking error
    module_with_blocker: str | None = None

    def __init__(
        self, messages: list[str], use_stdout: bool = False, module_with_blocker: str | None = None
    ) -> None:
        super().__init__("\n".join(messages))
        self.messages = messages
        self.use_stdout = use_stdout
        self.module_with_blocker = module_with_blocker


def remove_path_prefix(path: str, prefix: str | None) -> str:
    """If path starts with prefix, return copy of path with the prefix removed.
    Otherwise, return path. If path is None, return None.
    """
    if prefix is not None and path.startswith(prefix):
        return path[len(prefix) :]
    else:
        return path


def report_internal_error(
    err: Exception,
    file: str | None,
    line: int,
    errors: Errors | None,
    options: Options,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
) -> NoReturn:
    """Report internal error and exit.

    This optionally starts pdb or shows a traceback.
    """
    stdout = stdout or sys.stdout
    stderr = stderr or sys.stderr
    # Dump out errors so far, they often provide a clue.
    # But catch unexpected errors rendering them.
    if errors:
        try:
            for msg in errors.new_messages():
                print(msg)
        except Exception as e:
            print("Failed to dump errors:", repr(e), file=stderr)

    # Compute file:line prefix for official-looking error messages.
    if file:
        if line:
            prefix = f"{file}:{line}: "
        else:
            prefix = f"{file}: "
    else:
        prefix = ""

    # Print "INTERNAL ERROR" message.
    print(
        f"{prefix}error: INTERNAL ERROR --",
        "Please try using mypy master on GitHub:\n"
        "https://mypy.readthedocs.io/en/stable/common_issues.html"
        "#using-a-development-mypy-build",
        file=stderr,
    )
    if options.show_traceback:
        print("Please report a bug at https://github.com/python/mypy/issues", file=stderr)
    else:
        print(
            "If this issue continues with mypy master, "
            "please report a bug at https://github.com/python/mypy/issues",
            file=stderr,
        )
    print(f"version: {mypy_version}", file=stderr)

    # If requested, drop into pdb. This overrides show_tb.
    if options.pdb:
        print("Dropping into pdb", file=stderr)
        import pdb

        pdb.post_mortem(sys.exc_info()[2])

    # If requested, print traceback, else print note explaining how to get one.
    if options.raise_exceptions:
        raise err
    if not options.show_traceback:
        if not options.pdb:
            print(
                "{}note: please use --show-traceback to print a traceback "
                "when reporting a bug".format(prefix),
                file=stderr,
            )
    else:
        tberr = traceback.TracebackException.from_exception(err)
        tberr.stack[:0] = traceback.extract_stack()[:-2]
        print("".join(tberr.format()), file=stdout)
        print(f"{prefix}note: use --pdb to drop into pdb", file=stderr)

    # Exit.  The caller has nothing more to say.
    # We use exit code 2 to signal that this is no ordinary error.
    raise SystemExit(2)


class MypyError:
    def __init__(
        self,
        file_path: str,
        line: int,
        column: int,
        end_line: int,
        end_column: int,
        message: str,
        errorcode: str | None,
        severity: Literal["error", "note"],
    ) -> None:
        self.file_path = file_path
        self.line = line
        self.column = column
        self.end_line = end_line
        self.end_column = end_column
        self.message = message
        self.errorcode = errorcode
        self.severity = severity
        self.hints: list[str] = []


# (file_path, line, column)
_ErrorLocation = tuple[str, int, int]


def create_errors(error_tuples: list[ErrorTuple]) -> list[MypyError]:
    errors: list[MypyError] = []
    latest_error_at_location: dict[_ErrorLocation, MypyError] = {}

    for error_tuple in error_tuples:
        file_path, line, column, end_line, end_column, severity, message, errorcode = error_tuple
        if file_path is None:
            continue

        assert severity in ("error", "note")
        if severity == "note":
            error_location = (file_path, line, column)
            error = latest_error_at_location.get(error_location)
            if error is None:
                # This is purely a note, with no error correlated to it
                error = MypyError(
                    file_path,
                    line,
                    column,
                    end_line,
                    end_column,
                    message,
                    errorcode,
                    severity="note",
                )
                errors.append(error)
                continue

            error.hints.append(message)

        else:
            error = MypyError(
                file_path, line, column, end_line, end_column, message, errorcode, severity="error"
            )
            errors.append(error)
            error_location = (file_path, line, column)
            latest_error_at_location[error_location] = error

    return errors
