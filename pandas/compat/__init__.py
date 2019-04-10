"""
compat
======

Cross-compatible functions for Python 2 and 3.

Key items to import for 2/3 compatible code:
* lists: lrange(), lmap(), lzip(), lfilter()
* iterable method compatibility: iteritems, iterkeys, itervalues
  * Uses the original method if available, otherwise uses items, keys, values.
* add_metaclass(metaclass) - class decorator that recreates class with with the
  given metaclass instead (and avoids intermediary class creation)

Other items:
* platform checker
"""
# pylint disable=W0611
# flake8: noqa

import re
from distutils.version import LooseVersion
import sys
import platform
import struct

PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] >= 3
PY36 = sys.version_info >= (3, 6)
PY37 = sys.version_info >= (3, 7)
PYPY = platform.python_implementation() == 'PyPy'


# list-producing versions of the major Python iterating functions
def lrange(*args, **kwargs):
    return list(range(*args, **kwargs))


def lzip(*args, **kwargs):
    return list(zip(*args, **kwargs))


def lmap(*args, **kwargs):
    return list(map(*args, **kwargs))


def lfilter(*args, **kwargs):
    return list(filter(*args, **kwargs))


if PY2:
    def iteritems(obj, **kw):
        return obj.iteritems(**kw)

    def iterkeys(obj, **kw):
        return obj.iterkeys(**kw)

    def itervalues(obj, **kw):
        return obj.itervalues(**kw)

else:
    def iteritems(obj, **kw):
        return iter(obj.items(**kw))

    def iterkeys(obj, **kw):
        return iter(obj.keys(**kw))

    def itervalues(obj, **kw):
        return iter(obj.values(**kw))

# ----------------------------------------------------------------------------
# functions largely based / taken from the six module

# Much of the code in this module comes from Benjamin Peterson's six library.
# The license for this library can be found in LICENSES/SIX and the code can be
# found at https://bitbucket.org/gutworth/six


if PY3:
    def to_str(s):
        """
        Convert bytes and non-string into Python 3 str
        """
        if isinstance(s, bytes):
            s = s.decode('utf-8')
        elif not isinstance(s, str):
            s = str(s)
        return s

    def set_function_name(f, name, cls):
        """ Bind the name/qualname attributes of the function """
        f.__name__ = name
        f.__qualname__ = '{klass}.{name}'.format(
            klass=cls.__name__,
            name=name)
        f.__module__ = cls.__module__
        return f
else:
    def to_str(s):
        """
        Convert unicode and non-string into Python 2 str
        """
        if not isinstance(s, basestring):
            s = str(s)
        return s

    def set_function_name(f, name, cls):
        """ Bind the name attributes of the function """
        f.__name__ = name
        return f


def add_metaclass(metaclass):
    """Class decorator for creating a class with a metaclass."""
    def wrapper(cls):
        orig_vars = cls.__dict__.copy()
        orig_vars.pop('__dict__', None)
        orig_vars.pop('__weakref__', None)
        for slots_var in orig_vars.get('__slots__', ()):
            orig_vars.pop(slots_var)
        return metaclass(cls.__name__, cls.__bases__, orig_vars)
    return wrapper


if PY3:
    def raise_with_traceback(exc, traceback=Ellipsis):
        if traceback == Ellipsis:
            _, _, traceback = sys.exc_info()
        raise exc.with_traceback(traceback)
else:
    # this version of raise is a syntax error in Python 3
    exec("""
def raise_with_traceback(exc, traceback=Ellipsis):
    if traceback == Ellipsis:
        _, _, traceback = sys.exc_info()
    raise exc, None, traceback
""")

raise_with_traceback.__doc__ = """Raise exception with existing traceback.
If traceback is not passed, uses sys.exc_info() to get traceback."""


# dateutil minimum version
import dateutil

if LooseVersion(dateutil.__version__) < LooseVersion('2.5'):
    raise ImportError('dateutil 2.5.0 is the minimum required version')
from dateutil import parser as _date_parser
parse_date = _date_parser.parse


# In Python 3.7, the private re._pattern_type is removed.
# Python 3.5+ have typing.re.Pattern
if PY36:
    import typing
    re_type = typing.re.Pattern
else:
    re_type = type(re.compile(''))


# https://github.com/pandas-dev/pandas/pull/9123
def is_platform_little_endian():
    """ am I little endian """
    return sys.byteorder == 'little'


def is_platform_windows():
    return sys.platform == 'win32' or sys.platform == 'cygwin'


def is_platform_linux():
    return sys.platform == 'linux2'


def is_platform_mac():
    return sys.platform == 'darwin'


def is_platform_32bit():
    return struct.calcsize("P") * 8 < 64
