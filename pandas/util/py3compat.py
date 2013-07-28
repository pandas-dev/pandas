import sys

PY3 = (sys.version_info[0] >= 3)
# import iterator versions of these functions
from six.moves import zip, filter, reduce, map

try:
    import __builtin__ as builtins
    # not writeable when instantiated with string, doesn't handle unicode well
    from cStringIO import StringIO as StringIO
    # always writeable
    from StringIO import StringIO
    BytesIO = StringIO
    import cPickle
except ImportError:
    import builtins
    from io import StringIO, BytesIO
    cStringIO = StringIO
    import pickle as cPickle

if PY3:
    def isidentifier(s):
        return s.isidentifier()

    def str_to_bytes(s, encoding='ascii'):
        return s.encode(encoding)

    def bytes_to_str(b, encoding='utf-8'):
        return b.decode(encoding)

    # list-producing versions of the major Python iterating functions
    def lrange(*args, **kwargs):
        return list(range(*args, **kwargs))

    def lzip(*args, **kwargs):
        return list(zip(*args, **kwargs))

    def lmap(*args, **kwargs):
        return list(map(*args, **kwargs))

    def lfilter(*args, **kwargs):
        return list(filter(*args, **kwargs))

    # need to put range in the namespace
    range = range
    long = int
    unichr = chr
else:
    # Python 2
    import re
    _name_re = re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*$")

    def isidentifier(s, dotted=False):
        return bool(_name_re.match(s))

    def str_to_bytes(s, encoding='ascii'):
        return s

    def bytes_to_str(b, encoding='ascii'):
        return b

    # Python 2-builtin ranges produce lists
    lrange = builtins.range
    lzip = builtins.zip
    lmap = builtins.map
    lfilter = builtins.filter

    # have to explicitly put builtins into the namespace
    range = xrange
    long = long
    unichr = unichr

