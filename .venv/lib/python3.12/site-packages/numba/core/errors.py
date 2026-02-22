"""
Numba-specific errors and warnings.
"""


import abc
import contextlib
import os
import warnings
import numba.core.config
import numpy as np
from collections import defaultdict
from functools import wraps
from abc import abstractmethod

# Filled at the end
__all__ = []


def _is_numba_core_config_loaded():
    """
    To detect if numba.core.config has been initialized due to circular imports.
    """
    try:
        numba.core.config
    except AttributeError:
        return False
    else:
        return True


class NumbaWarning(Warning):
    """
    Base category for all Numba compiler warnings.
    """

    def __init__(self, msg, loc=None, highlighting=True, ):
        self.msg = msg
        self.loc = loc

        # If a warning is emitted inside validation of env-vars in
        # numba.core.config. Highlighting will not be available.
        if highlighting and _is_numba_core_config_loaded():
            highlight = termcolor().errmsg
        else:
            def highlight(x):
                return x
        if loc:
            super(NumbaWarning, self).__init__(
                highlight("%s\n%s\n" % (msg, loc.strformat())))
        else:
            super(NumbaWarning, self).__init__(highlight("%s" % (msg,)))


class NumbaPerformanceWarning(NumbaWarning):
    """
    Warning category for when an operation might not be
    as fast as expected.
    """


class NumbaDeprecationWarning(NumbaWarning, DeprecationWarning):
    """
    Warning category for use of a deprecated feature.
    """


class NumbaPendingDeprecationWarning(NumbaWarning, PendingDeprecationWarning):
    """
    Warning category for use of a feature that is pending deprecation.
    """


class NumbaParallelSafetyWarning(NumbaWarning):
    """
    Warning category for when an operation in a prange
    might not have parallel semantics.
    """


class NumbaTypeSafetyWarning(NumbaWarning):
    """
    Warning category for unsafe casting operations.
    """


class NumbaExperimentalFeatureWarning(NumbaWarning):
    """
    Warning category for using an experimental feature.
    """


class NumbaInvalidConfigWarning(NumbaWarning):
    """
    Warning category for using an invalid configuration.
    """


class NumbaPedanticWarning(NumbaWarning):
    """
    Warning category for reporting pedantic messages.
    """
    def __init__(self, msg, **kwargs):
        super().__init__(f"{msg}\n{pedantic_warning_info}")


class NumbaIRAssumptionWarning(NumbaPedanticWarning):
    """
    Warning category for reporting an IR assumption violation.
    """


class NumbaDebugInfoWarning(NumbaWarning):
    """
    Warning category for an issue with the emission of debug information.
    """


class NumbaSystemWarning(NumbaWarning):
    """
    Warning category for an issue with the system configuration.
    """

# These are needed in the color formatting of errors setup


class _ColorScheme(metaclass=abc.ABCMeta):

    @abstractmethod
    def code(self, msg):
        pass

    @abstractmethod
    def errmsg(self, msg):
        pass

    @abstractmethod
    def filename(self, msg):
        pass

    @abstractmethod
    def indicate(self, msg):
        pass

    @abstractmethod
    def highlight(self, msg):
        pass

    @abstractmethod
    def reset(self, msg):
        pass


class _DummyColorScheme(_ColorScheme):

    def __init__(self, theme=None):
        pass

    def code(self, msg):
        pass

    def errmsg(self, msg):
        pass

    def filename(self, msg):
        pass

    def indicate(self, msg):
        pass

    def highlight(self, msg):
        pass

    def reset(self, msg):
        pass


# holds reference to the instance of the terminal color scheme in use
_termcolor_inst = None

try:
    import colorama

    # If the colorama version is < 0.3.9 it can break stdout/stderr in some
    # situations, as a result if this condition is met colorama is disabled and
    # the user is warned. Note that early versions did not have a __version__.
    colorama_version = getattr(colorama, '__version__', '0.0.0')

    if tuple([int(x) for x in colorama_version.split('.')]) < (0, 3, 9):
        msg = ("Insufficiently recent colorama version found. "
               "Numba requires colorama >= 0.3.9")
        # warn the user
        warnings.warn(msg)
        # trip the exception to disable color errors
        raise ImportError

    # If Numba is running in testsuite mode then do not use error message
    # coloring so CI system output is consistently readable without having
    # to read between shell escape characters.
    if os.environ.get('NUMBA_DISABLE_ERROR_MESSAGE_HIGHLIGHTING', None):
        raise ImportError  # just to trigger the exception handler below

except ImportError:

    class NOPColorScheme(_DummyColorScheme):
        def __init__(self, theme=None):
            if theme is not None:
                raise ValueError("specifying a theme has no effect")
            _DummyColorScheme.__init__(self, theme=theme)

        def code(self, msg):
            return msg

        def errmsg(self, msg):
            return msg

        def filename(self, msg):
            return msg

        def indicate(self, msg):
            return msg

        def highlight(self, msg):
            return msg

        def reset(self, msg):
            return msg

    def termcolor():
        global _termcolor_inst
        if _termcolor_inst is None:
            _termcolor_inst = NOPColorScheme()
        return _termcolor_inst

else:

    from colorama import init, reinit, deinit, Fore, Style

    class ColorShell(object):
        _has_initialized = False

        def __init__(self):
            init()
            self._has_initialized = True

        def __enter__(self):
            if self._has_initialized:
                reinit()

        def __exit__(self, *exc_detail):
            Style.RESET_ALL
            deinit()

    class reset_terminal(object):
        def __init__(self):
            self._buf = bytearray(b'')

        def __enter__(self):
            return self._buf

        def __exit__(self, *exc_detail):
            self._buf += bytearray(Style.RESET_ALL.encode('utf-8'))

    # define some default themes, if more are added, update the envvars docs!
    themes = {}

    # No color added, just bold weighting
    themes['no_color'] = {'code': None,
                          'errmsg': None,
                          'filename': None,
                          'indicate': None,
                          'highlight': None,
                          'reset': None, }

    # suitable for terminals with a dark background
    themes['dark_bg'] = {'code': Fore.BLUE,
                         'errmsg': Fore.YELLOW,
                         'filename': Fore.WHITE,
                         'indicate': Fore.GREEN,
                         'highlight': Fore.RED,
                         'reset': Style.RESET_ALL, }

    # suitable for terminals with a light background
    themes['light_bg'] = {'code': Fore.BLUE,
                          'errmsg': Fore.BLACK,
                          'filename': Fore.MAGENTA,
                          'indicate': Fore.BLACK,
                          'highlight': Fore.RED,
                          'reset': Style.RESET_ALL, }

    # suitable for terminals with a blue background
    themes['blue_bg'] = {'code': Fore.WHITE,
                         'errmsg': Fore.YELLOW,
                         'filename': Fore.MAGENTA,
                         'indicate': Fore.CYAN,
                         'highlight': Fore.RED,
                         'reset': Style.RESET_ALL, }

    # suitable for use in jupyter notebooks
    themes['jupyter_nb'] = {'code': Fore.BLACK,
                            'errmsg': Fore.BLACK,
                            'filename': Fore.GREEN,
                            'indicate': Fore.CYAN,
                            'highlight': Fore.RED,
                            'reset': Style.RESET_ALL, }

    default_theme = themes['no_color']

    class HighlightColorScheme(_DummyColorScheme):
        def __init__(self, theme=default_theme):
            self._code = theme['code']
            self._errmsg = theme['errmsg']
            self._filename = theme['filename']
            self._indicate = theme['indicate']
            self._highlight = theme['highlight']
            self._reset = theme['reset']
            _DummyColorScheme.__init__(self, theme=theme)

        def _markup(self, msg, color=None, style=Style.BRIGHT):
            features = ''
            if color:
                features += color
            if style:
                features += style
            with ColorShell():
                with reset_terminal() as mu:
                    mu += features.encode('utf-8')
                    mu += (msg).encode('utf-8')
                return mu.decode('utf-8')

        def code(self, msg):
            return self._markup(msg, self._code)

        def errmsg(self, msg):
            return self._markup(msg, self._errmsg)

        def filename(self, msg):
            return self._markup(msg, self._filename)

        def indicate(self, msg):
            return self._markup(msg, self._indicate)

        def highlight(self, msg):
            return self._markup(msg, self._highlight)

        def reset(self, msg):
            return self._markup(msg, self._reset)

    def termcolor():
        global _termcolor_inst
        if _termcolor_inst is None:
            scheme = themes[numba.core.config.COLOR_SCHEME]
            _termcolor_inst = HighlightColorScheme(scheme)
        return _termcolor_inst


pedantic_warning_info = """
This warning came from an internal pedantic check. Please report the warning
message and traceback, along with a minimal reproducer at:
https://github.com/numba/numba/issues/new?template=bug_report.md
"""

feedback_details = """
Please report the error message and traceback, along with a minimal reproducer
at: https://github.com/numba/numba/issues/new?template=bug_report.md

If more help is needed please feel free to speak to the Numba core developers
directly at: https://gitter.im/numba/numba

Thanks in advance for your help in improving Numba!
"""

unsupported_error_info = """
Unsupported functionality was found in the code Numba was trying to compile.

If this functionality is important to you please file a feature request at:
https://github.com/numba/numba/issues/new?template=feature_request.md
"""

interpreter_error_info = """
Unsupported Python functionality was found in the code Numba was trying to
compile. This error could be due to invalid code, does the code work
without Numba? (To temporarily disable Numba JIT, set the `NUMBA_DISABLE_JIT`
environment variable to non-zero, and then rerun the code).

If the code is valid and the unsupported functionality is important to you
please file a feature request at:
https://github.com/numba/numba/issues/new?template=feature_request.md

To see Python/NumPy features supported by the latest release of Numba visit:
https://numba.readthedocs.io/en/stable/reference/pysupported.html
and
https://numba.readthedocs.io/en/stable/reference/numpysupported.html
"""

constant_inference_info = """
Numba could not make a constant out of something that it decided should be
a constant. This could well be a current limitation in Numba's internals,
however please first check that your code is valid for compilation,
particularly with respect to string interpolation (not supported!) and
the requirement of compile time constants as arguments to exceptions:
https://numba.readthedocs.io/en/stable/reference/pysupported.html?highlight=exceptions#constructs

If the code is valid and the unsupported functionality is important to you
please file a feature request at:
https://github.com/numba/numba/issues/new?template=feature_request.md

If you think your code should work with Numba. %s
""" % feedback_details

typing_error_info = """
This is not usually a problem with Numba itself but instead often caused by
the use of unsupported features or an issue in resolving types.

To see Python/NumPy features supported by the latest release of Numba visit:
https://numba.readthedocs.io/en/stable/reference/pysupported.html
and
https://numba.readthedocs.io/en/stable/reference/numpysupported.html

For more information about typing errors and how to debug them visit:
https://numba.readthedocs.io/en/stable/user/troubleshoot.html#my-code-doesn-t-compile

If you think your code should work with Numba, please report the error message
and traceback, along with a minimal reproducer at:
https://github.com/numba/numba/issues/new?template=bug_report.md
"""

reportable_issue_info = """
-------------------------------------------------------------------------------
This should not have happened, a problem has occurred in Numba's internals.
You are currently using Numba version %s.
%s
""" % (numba.__version__, feedback_details)

error_extras = dict()
error_extras['unsupported_error'] = unsupported_error_info
error_extras['typing'] = typing_error_info
error_extras['reportable'] = reportable_issue_info
error_extras['interpreter'] = interpreter_error_info
error_extras['constant_inference'] = constant_inference_info


def deprecated(arg):
    """Define a deprecation decorator.
    An optional string should refer to the new API to be used instead.

    Example:
      @deprecated
      def old_func(): ...

      @deprecated('new_func')
      def old_func(): ..."""

    subst = arg if isinstance(arg, str) else None

    def decorator(func):
        def wrapper(*args, **kwargs):
            msg = "Call to deprecated function \"{}\"."
            if subst:
                msg += "\n Use \"{}\" instead."
            warnings.warn(msg.format(func.__name__, subst),
                          category=DeprecationWarning, stacklevel=2)
            return func(*args, **kwargs)

        return wraps(func)(wrapper)

    if not subst:
        return decorator(arg)
    else:
        return decorator


class WarningsFixer(object):
    """
    An object "fixing" warnings of a given category caught during
    certain phases.  The warnings can have their filename and lineno fixed,
    and they are deduplicated as well.

    When used as a context manager, any warnings caught by `.catch_warnings()`
    will be flushed at the exit of the context manager.
    """

    def __init__(self, category):
        self._category = category
        # {(filename, lineno, category) -> messages}
        self._warnings = defaultdict(set)

    @contextlib.contextmanager
    def catch_warnings(self, filename=None, lineno=None):
        """
        Store warnings and optionally fix their filename and lineno.
        """
        with warnings.catch_warnings(record=True) as wlist:
            warnings.simplefilter('always', self._category)
            yield

        for w in wlist:
            msg = str(w.message)
            if issubclass(w.category, self._category):
                # Store warnings of this category for deduplication
                filename = filename or w.filename
                lineno = lineno or w.lineno
                self._warnings[filename, lineno, w.category].add(msg)
            else:
                # Simply emit other warnings again
                warnings.warn_explicit(msg, w.category,
                                       w.filename, w.lineno)

    def flush(self):
        """
        Emit all stored warnings.
        """
        def key(arg):
            # It is possible through codegen to create entirely identical
            # warnings, this leads to comparing types when sorting which breaks
            # on Python 3. Key as str() and if the worse happens then `id`
            # creates some uniqueness
            return str(arg) + str(id(arg))

        for (filename, lineno, category), messages in sorted(
                self._warnings.items(), key=key):
            for msg in sorted(messages):
                warnings.warn_explicit(msg, category, filename, lineno)
        self._warnings.clear()

    def __enter__(self):
        return

    def __exit__(self, exc_type, exc_value, traceback):
        self.flush()


class NumbaError(Exception):
    def __init__(self, msg, loc=None, highlighting=True):
        self.msg = msg
        self.loc = loc
        if highlighting:
            highlight = termcolor().errmsg
        else:
            def highlight(x):
                return x

        if loc:
            new_msg = "%s\n%s\n" % (msg, loc.strformat())
        else:
            new_msg = "%s" % (msg,)
        super(NumbaError, self).__init__(highlight(new_msg))

    @property
    def contexts(self):
        try:
            return self._contexts
        except AttributeError:
            self._contexts = lst = []
            return lst

    def add_context(self, msg):
        """
        Add contextual info.  The exception message is expanded with the new
        contextual information.
        """
        if msg in self.contexts:
            # avoid duplicating contexts
            return self
        self.contexts.append(msg)
        f = termcolor().errmsg('{0}\n') + termcolor().filename('During: {1}')
        newmsg = f.format(self, msg)
        self.args = (newmsg,)
        return self

    def patch_message(self, new_message):
        """
        Change the error message to the given new message.
        """
        self.args = (new_message,) + self.args[1:]


class UnsupportedError(NumbaError):
    """
    Numba does not have an implementation for this functionality.
    """


class UnsupportedBytecodeError(Exception):
    """Unsupported bytecode is non-recoverable
    """
    def __init__(self, msg, loc=None):
        super().__init__(f"{msg}. Raised from {loc}")


class UnsupportedRewriteError(UnsupportedError):
    """UnsupportedError from rewrite passes
    """
    pass


class IRError(NumbaError):
    """
    An error occurred during Numba IR generation.
    """
    pass


class RedefinedError(IRError):
    """
    An error occurred during interpretation of IR due to variable redefinition.
    """
    pass


class NotDefinedError(IRError):
    """
    An undefined variable is encountered during interpretation of IR.
    """

    def __init__(self, name, loc=None):
        self.name = name
        msg = ("The compiler failed to analyze the bytecode. "
               "Variable '%s' is not defined." % name)
        super(NotDefinedError, self).__init__(msg, loc=loc)


class VerificationError(IRError):
    """
    An error occurred during IR verification. Once Numba's internal
    representation (IR) is constructed it is then verified to ensure that
    terminators are both present and in the correct places within the IR. If
    it is the case that this condition is not met, a VerificationError is
    raised.
    """
    pass


class DeprecationError(NumbaError):
    """
    Functionality is deprecated.
    """
    pass


class LoweringError(NumbaError):
    """
    An error occurred during lowering.
    """

    def __init__(self, msg, loc=None):
        super(LoweringError, self).__init__(msg, loc=loc)


class UnsupportedParforsError(NumbaError):
    """
    An error occurred because parfors is not supported on the platform.
    """
    pass


class ForbiddenConstruct(LoweringError):
    """
    A forbidden Python construct was encountered (e.g. use of locals()).
    """
    pass


class TypingError(NumbaError):
    """
    A type inference failure.
    """
    pass


class UntypedAttributeError(TypingError):
    def __init__(self, value, attr, loc=None):
        module = getattr(value, 'pymod', None)
        if module is not None and module == np:
            # unsupported numpy feature.
            msg = ("Use of unsupported NumPy function 'numpy.%s' "
                   "or unsupported use of the function.") % attr
        else:
            msg = "Unknown attribute '{attr}' of type {type}"
            msg = msg.format(type=value, attr=attr)
        super(UntypedAttributeError, self).__init__(msg, loc=loc)


class ByteCodeSupportError(NumbaError):
    """
    Failure to extract the bytecode of the user's function.
    """

    def __init__(self, msg, loc=None):
        super(ByteCodeSupportError, self).__init__(msg, loc=loc)


class CompilerError(NumbaError):
    """
    Some high-level error in the compiler.
    """
    pass


class ConstantInferenceError(NumbaError):
    """
    Failure during constant inference.
    """

    def __init__(self, value, loc=None):
        super(ConstantInferenceError, self).__init__(value, loc=loc)


class InternalError(NumbaError):
    """
    For wrapping internal error occurred within the compiler
    """

    def __init__(self, exception):
        super(InternalError, self).__init__(str(exception))
        self.old_exception = exception


class InternalTargetMismatchError(InternalError):
    """For signalling a target mismatch error occurred internally within the
    compiler.
    """
    def __init__(self, kind, target_hw, hw_clazz):
        msg = (f"{kind.title()} being resolved on a target from which it does "
               f"not inherit. Local target is {target_hw}, declared "
               f"target class is {hw_clazz}.")
        super().__init__(msg)


class NonexistentTargetError(InternalError):
    """For signalling that a target that does not exist was requested.
    """
    pass


class RequireLiteralValue(TypingError):
    """
    For signalling that a function's typing requires a constant value for
    some of its arguments.
    """
    pass


class ForceLiteralArg(NumbaError):
    """A Pseudo-exception to signal the dispatcher to type an argument literally

    Attributes
    ----------
    requested_args : frozenset[int]
        requested positions of the arguments.
    """
    def __init__(self, arg_indices, fold_arguments=None, loc=None):
        """
        Parameters
        ----------
        arg_indices : Sequence[int]
            requested positions of the arguments.
        fold_arguments: callable
            A function ``(tuple, dict) -> tuple`` that binds and flattens
            the ``args`` and ``kwargs``.
        loc : numba.ir.Loc or None
        """
        super(ForceLiteralArg, self).__init__(
            "Pseudo-exception to force literal arguments in the dispatcher",
            loc=loc,
        )
        self.requested_args = frozenset(arg_indices)
        self.fold_arguments = fold_arguments

    def bind_fold_arguments(self, fold_arguments):
        """Bind the fold_arguments function
        """
        # to avoid circular import
        from numba.core.utils import chain_exception

        e = ForceLiteralArg(self.requested_args, fold_arguments,
                            loc=self.loc)
        return chain_exception(e, self)

    def combine(self, other):
        """Returns a new instance by or'ing the requested_args.
        """
        if not isinstance(other, ForceLiteralArg):
            m = '*other* must be a {} but got a {} instead'
            raise TypeError(m.format(ForceLiteralArg, type(other)))
        return ForceLiteralArg(self.requested_args | other.requested_args)

    def __or__(self, other):
        """Same as self.combine(other)
        """
        return self.combine(other)


class LiteralTypingError(TypingError):
    """
    Failure in typing a Literal type
    """
    pass


# These Exception classes are just Numba copies of their Python equivalents for
# use internally in cases where we want e.g. type inference to keep on trying.
# Exceptions extending from NumbaError are considered "special" by Numba's
# internals and are treated differently to standard Python exceptions which are
# permitted to just propagate up the stack.

class NumbaValueError(TypingError):
    pass


class NumbaTypeError(TypingError):
    pass


class NumbaAttributeError(TypingError):
    pass


class NumbaAssertionError(TypingError):
    pass


class NumbaNotImplementedError(TypingError):
    pass


class NumbaKeyError(TypingError):
    pass


class NumbaIndexError(TypingError):
    pass


class NumbaRuntimeError(NumbaError):
    pass


def _format_msg(fmt, args, kwargs):
    # If no formatting arguments are supplied, return the string unchanged.
    # This avoids KeyError when fmt contains curly braces, which can be
    # interpreted as format fields.
    if not args and not kwargs:
        return fmt
    return fmt.format(*args, **kwargs)


_numba_path = os.path.dirname(__file__)
loc_info = {}


@contextlib.contextmanager
def new_error_context(fmt_, *args, **kwargs):
    """
    A contextmanager that prepend contextual information to any exception
    raised within.

    The first argument is a message that describes the context.  It can be a
    format string.  If there are additional arguments, it will be used as
    ``fmt_.format(*args, **kwargs)`` to produce the final message string.
    """
    loc = kwargs.get('loc', None)
    if loc is not None and not loc.filename.startswith(_numba_path):
        loc_info.update(kwargs)

    try:
        yield
    except NumbaError as e:
        e.add_context(_format_msg(fmt_, args, kwargs))
        raise


__all__ += [name for (name, value) in globals().items()
            if not name.startswith('_') and isinstance(value, type)
            and issubclass(value, (Exception, Warning))]
