"""
This module contains utils for manipulating target configurations such as
compiler flags.
"""
import re
import zlib
import base64

from types import MappingProxyType
from numba.core import utils


class Option:
    """An option to be used in ``TargetConfig``.
    """
    __slots__ = "_type", "_default", "_doc"

    def __init__(self, type, *, default, doc):
        """
        Parameters
        ----------
        type :
            Type of the option value. It can be a callable.
            The setter always calls ``self._type(value)``.
        default :
            The default value for the option.
        doc : str
            Docstring for the option.
        """
        self._type = type
        self._default = default
        self._doc = doc

    @property
    def type(self):
        return self._type

    @property
    def default(self):
        return self._default

    @property
    def doc(self):
        return self._doc


class _FlagsStack(utils.ThreadLocalStack, stack_name="flags"):
    pass


class ConfigStack:
    """A stack for tracking target configurations in the compiler.

    It stores the stack in a thread-local class attribute. All instances in the
    same thread will see the same stack.
    """
    @classmethod
    def top_or_none(cls):
        """Get the TOS or return None if no config is set.
        """
        self = cls()
        if self:
            flags = self.top()
        else:
            # Note: should this be the default flag for the target instead?
            flags = None
        return flags

    def __init__(self):
        self._stk = _FlagsStack()

    def top(self):
        return self._stk.top()

    def __len__(self):
        return len(self._stk)

    def enter(self, flags):
        """Returns a contextmanager that performs ``push(flags)`` on enter and
        ``pop()`` on exit.
        """
        return self._stk.enter(flags)


class _MetaTargetConfig(type):
    """Metaclass for ``TargetConfig``.

    When a subclass of ``TargetConfig`` is created, all ``Option`` defined
    as class members will be parsed and corresponding getters, setters, and
    delters will be inserted.
    """
    def __init__(cls, name, bases, dct):
        """Invoked when subclass is created.

        Insert properties for each ``Option`` that are class members.
        All the options will be grouped inside the ``.options`` class
        attribute.
        """
        # Gather options from base classes and class dict
        opts = {}
        # Reversed scan into the base classes to follow MRO ordering such that
        # the closest base class is overriding
        for base_cls in reversed(bases):
            opts.update(base_cls.options)
        opts.update(cls.find_options(dct))
        # Store the options into class attribute as a ready-only mapping.
        cls.options = MappingProxyType(opts)

        # Make properties for each of the options
        def make_prop(name, option):
            def getter(self):
                return self._values.get(name, option.default)

            def setter(self, val):
                self._values[name] = option.type(val)

            def delter(self):
                del self._values[name]

            return property(getter, setter, delter, option.doc)

        for name, option in cls.options.items():
            setattr(cls, name, make_prop(name, option))

    def find_options(cls, dct):
        """Returns a new dict with all the items that are a mapping to an
        ``Option``.
        """
        return {k: v for k, v in dct.items() if isinstance(v, Option)}


class _NotSetType:
    def __repr__(self):
        return "<NotSet>"


_NotSet = _NotSetType()


class TargetConfig(metaclass=_MetaTargetConfig):
    """Base class for ``TargetConfig``.

    Subclass should fill class members with ``Option``. For example:

    >>> class MyTargetConfig(TargetConfig):
    >>>     a_bool_option = Option(type=bool, default=False, doc="a bool")
    >>>     an_int_option = Option(type=int, default=0, doc="an int")

    The metaclass will insert properties for each ``Option``. For example:

    >>> tc = MyTargetConfig()
    >>> tc.a_bool_option = True  # invokes the setter
    >>> print(tc.an_int_option)  # print the default
    """

    # Used for compression in mangling.
    # Set to -15 to disable the header and checksum for smallest output.
    _ZLIB_CONFIG = {"wbits": -15}

    def __init__(self, copy_from=None):
        """
        Parameters
        ----------
        copy_from : TargetConfig or None
            if None, creates an empty ``TargetConfig``.
            Otherwise, creates a copy.
        """
        self._values = {}
        if copy_from is not None:
            assert isinstance(copy_from, TargetConfig)
            self._values.update(copy_from._values)

    def __repr__(self):
        # NOTE: default options will be placed at the end and grouped inside
        #       a square bracket; i.e. [optname=optval, ...]
        args = []
        defs = []
        for k in self.options:
            msg = f"{k}={getattr(self, k)}"
            if not self.is_set(k):
                defs.append(msg)
            else:
                args.append(msg)
        clsname = self.__class__.__name__
        return f"{clsname}({', '.join(args)}, [{', '.join(defs)}])"

    def __hash__(self):
        return hash(tuple(sorted(self.values())))

    def __eq__(self, other):
        if isinstance(other, TargetConfig):
            return self.values() == other.values()
        else:
            return NotImplemented

    def values(self):
        """Returns a dict of all the values
        """
        return {k: getattr(self, k) for k in self.options}

    def is_set(self, name):
        """Is the option set?
        """
        self._guard_option(name)
        return name in self._values

    def discard(self, name):
        """Remove the option by name if it is defined.

        After this, the value for the option will be set to its default value.
        """
        self._guard_option(name)
        self._values.pop(name, None)

    def inherit_if_not_set(self, name, default=_NotSet):
        """Inherit flag from ``ConfigStack``.

        Parameters
        ----------
        name : str
            Option name.
        default : optional
            When given, it overrides the default value.
            It is only used when the flag is not defined locally and there is
            no entry in the ``ConfigStack``.
        """
        self._guard_option(name)
        if not self.is_set(name):
            cstk = ConfigStack()
            if cstk:
                # inherit
                top = cstk.top()
                setattr(self, name, getattr(top, name))
            elif default is not _NotSet:
                setattr(self, name, default)

    def copy(self):
        """Clone this instance.
        """
        return type(self)(self)

    def summary(self) -> str:
        """Returns a ``str`` that summarizes this instance.

        In contrast to ``__repr__``, only options that are explicitly set will
        be shown.
        """
        args = [f"{k}={v}" for k, v in self._summary_args()]
        clsname = self.__class__.__name__
        return f"{clsname}({', '.join(args)})"

    def _guard_option(self, name):
        if name not in self.options:
            msg = f"{name!r} is not a valid option for {type(self)}"
            raise ValueError(msg)

    def _summary_args(self):
        """returns a sorted sequence of 2-tuple containing the
        ``(flag_name, flag_value)`` for flag that are set with a non-default
        value.
        """
        args = []
        for k in sorted(self.options):
            opt = self.options[k]
            if self.is_set(k):
                flagval = getattr(self, k)
                if opt.default != flagval:
                    v = (k, flagval)
                    args.append(v)
        return args

    @classmethod
    def _make_compression_dictionary(cls) -> bytes:
        """Returns a ``bytes`` object suitable for use as a dictionary for
        compression.
        """
        buf = []
        # include package name
        buf.append("numba")
        # include class name
        buf.append(cls.__class__.__name__)
        # include common values
        buf.extend(["True", "False"])
        # include all options name and their default value
        for k, opt in cls.options.items():
            buf.append(k)
            buf.append(str(opt.default))
        return ''.join(buf).encode()

    def get_mangle_string(self) -> str:
        """Return a string suitable for symbol mangling.
        """
        zdict = self._make_compression_dictionary()

        comp = zlib.compressobj(zdict=zdict, level=zlib.Z_BEST_COMPRESSION,
                                **self._ZLIB_CONFIG)
        # The mangled string is a compressed and base64 encoded version of the
        # summary
        buf = [comp.compress(self.summary().encode())]
        buf.append(comp.flush())
        return base64.b64encode(b''.join(buf)).decode()

    @classmethod
    def demangle(cls, mangled: str) -> str:
        """Returns the demangled result from ``.get_mangle_string()``
        """
        # unescape _XX sequence
        def repl(x):
            return chr(int('0x' + x.group(0)[1:], 16))
        unescaped = re.sub(r"_[a-zA-Z0-9][a-zA-Z0-9]", repl, mangled)
        # decode base64
        raw = base64.b64decode(unescaped)
        # decompress
        zdict = cls._make_compression_dictionary()
        dc = zlib.decompressobj(zdict=zdict, **cls._ZLIB_CONFIG)
        buf = []
        while raw:
            buf.append(dc.decompress(raw))
            raw = dc.unconsumed_tail
        buf.append(dc.flush())
        return b''.join(buf).decode()
