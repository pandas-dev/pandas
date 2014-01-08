"""
compat
======

Cross-compatible functions for Python 2 and 3.

Key items to import for 2/3 compatible code:
* iterators: range(), map(), zip(), filter(), reduce()
* lists: lrange(), lmap(), lzip(), lfilter()
* unicode: u() [u"" is a syntax error in Python 3.0-3.2]
* longs: long (int in Python 3)
* callable
* iterable method compatibility: iteritems, iterkeys, itervalues
  * Uses the original method if available, otherwise uses items, keys, values.
* types:
    * text_type: unicode in Python 2, str in Python 3
    * binary_type: str in Python 2, bythes in Python 3
    * string_types: basestring in Python 2, str in Python 3
* bind_method: binds functions to classes
* add_metaclass(metaclass) - class decorator that recreates class with with the
  given metaclass instead (and avoids intermediary class creation)

Python 2.6 compatibility:
* OrderedDict
* Counter

Other items:
* OrderedDefaultDict
"""
# pylint disable=W0611
import functools
import itertools
from distutils.version import LooseVersion
from itertools import product
import sys
import types

PY3 = (sys.version_info[0] >= 3)
PY3_2 = sys.version_info[:2] == (3, 2)

try:
    import __builtin__ as builtins
    # not writeable when instantiated with string, doesn't handle unicode well
    from cStringIO import StringIO as cStringIO
    # always writeable
    from StringIO import StringIO
    BytesIO = StringIO
    import cPickle
    import httplib
except ImportError:
    import builtins
    from io import StringIO, BytesIO
    cStringIO = StringIO
    import pickle as cPickle
    import http.client as httplib


if PY3:
    def isidentifier(s):
        return s.isidentifier()

    def str_to_bytes(s, encoding=None):
        return s.encode(encoding or 'ascii')

    def bytes_to_str(b, encoding=None):
        return b.decode(encoding or 'utf-8')

    # have to explicitly put builtins into the namespace
    range = range
    map = map
    zip = zip
    filter = filter
    reduce = functools.reduce
    long = int
    unichr = chr

    # list-producing versions of the major Python iterating functions
    def lrange(*args, **kwargs):
        return list(range(*args, **kwargs))

    def lzip(*args, **kwargs):
        return list(zip(*args, **kwargs))

    def lmap(*args, **kwargs):
        return list(map(*args, **kwargs))

    def lfilter(*args, **kwargs):
        return list(filter(*args, **kwargs))
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

    # import iterator versions of these functions
    range = xrange
    zip = itertools.izip
    filter = itertools.ifilter
    map = itertools.imap
    reduce = reduce
    long = long
    unichr = unichr

    # Python 2-builtin ranges produce lists
    lrange = builtins.range
    lzip = builtins.zip
    lmap = builtins.map
    lfilter = builtins.filter


def iteritems(obj, **kwargs):
    """replacement for six's iteritems for Python2/3 compat
       uses 'iteritems' if available and otherwise uses 'items'.

       Passes kwargs to method.
    """
    func = getattr(obj, "iteritems", None)
    if not func:
        func = obj.items
    return func(**kwargs)


def iterkeys(obj, **kwargs):
    func = getattr(obj, "iterkeys", None)
    if not func:
        func = obj.keys
    return func(**kwargs)


def itervalues(obj, **kwargs):
    func = getattr(obj, "itervalues", None)
    if not func:
        func = obj.values
    return func(**kwargs)


def bind_method(cls, name, func):
    """Bind a method to class, python 2 and python 3 compatible.

    Parameters
    ----------

    cls : type
        class to receive bound method
    name : basestring
        name of method on class instance
    func : function
        function to be bound as method


    Returns
    -------
    None
    """
    # only python 2 has bound/unbound method issue
    if not PY3:
        setattr(cls, name, types.MethodType(func, None, cls))
    else:
        setattr(cls, name, func)
# ----------------------------------------------------------------------------
# functions largely based / taken from the six module

# Much of the code in this module comes from Benjamin Peterson's six library.
# The license for this library can be found in LICENSES/SIX and the code can be
# found at https://bitbucket.org/gutworth/six

if PY3:
    string_types = str,
    integer_types = int,
    class_types = type,
    text_type = str
    binary_type = bytes

    def u(s):
        return s

    def u_safe(s):
        return s
else:
    string_types = basestring,
    integer_types = (int, long)
    class_types = (type, types.ClassType)
    text_type = unicode
    binary_type = str

    def u(s):
        return unicode(s, "unicode_escape")

    def u_safe(s):
        try:
            return unicode(s, "unicode_escape")
        except:
            return s


string_and_binary_types = string_types + (binary_type,)


try:
    # callable reintroduced in later versions of Python
    callable = callable
except NameError:
    def callable(obj):
        return any("__call__" in klass.__dict__ for klass in type(obj).__mro__)


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


# ----------------------------------------------------------------------------
# Python 2.6 compatibility shims
#

# OrderedDict Shim from  Raymond Hettinger, python core dev
# http://code.activestate.com/recipes/576693-ordered-dictionary-for-py24/
# here to support versions before 2.6
if not PY3:
    # don't need this except in 2.6
    try:
        from thread import get_ident as _get_ident
    except ImportError:
        from dummy_thread import get_ident as _get_ident

try:
    from _abcoll import KeysView, ValuesView, ItemsView
except ImportError:
    pass


class _OrderedDict(dict):

    """Dictionary that remembers insertion order"""
    # An inherited dict maps keys to values.
    # The inherited dict provides __getitem__, __len__, __contains__, and get.
    # The remaining methods are order-aware.
    # Big-O running times for all methods are the same as for regular
    # dictionaries.

    # The internal self.__map dictionary maps keys to links in a doubly linked
    # list.  The circular doubly linked list starts and ends with a sentinel
    # element.  The sentinel element never gets deleted (this simplifies the
    # algorithm).  Each link is stored as a list of length three:  [PREV, NEXT,
    # KEY].

    def __init__(self, *args, **kwds):
        """Initialize an ordered dictionary. Signature is the same as for
        regular dictionaries, but keyword arguments are not recommended
        because their insertion order is arbitrary.
        """
        if len(args) > 1:
            raise TypeError('expected at most 1 arguments, got %d' % len(args))
        try:
            self.__root
        except AttributeError:
            self.__root = root = []                     # sentinel node
            root[:] = [root, root, None]
            self.__map = {}
        self.__update(*args, **kwds)

    def __setitem__(self, key, value, dict_setitem=dict.__setitem__):
        """od.__setitem__(i, y) <==> od[i]=y"""
        # Setting a new item creates a new link which goes at the end of the
        # linked list, and the inherited dictionary is updated with the new
        # key/value pair.
        if key not in self:
            root = self.__root
            last = root[0]
            last[1] = root[0] = self.__map[key] = [last, root, key]
        dict_setitem(self, key, value)

    def __delitem__(self, key, dict_delitem=dict.__delitem__):
        """od.__delitem__(y) <==> del od[y]"""
        # Deleting an existing item uses self.__map to find the link which is
        # then removed by updating the links in the predecessor and successor
        # nodes.
        dict_delitem(self, key)
        link_prev, link_next, key = self.__map.pop(key)
        link_prev[1] = link_next
        link_next[0] = link_prev

    def __iter__(self):
        """od.__iter__() <==> iter(od)"""
        root = self.__root
        curr = root[1]
        while curr is not root:
            yield curr[2]
            curr = curr[1]

    def __reversed__(self):
        """od.__reversed__() <==> reversed(od)"""
        root = self.__root
        curr = root[0]
        while curr is not root:
            yield curr[2]
            curr = curr[0]

    def clear(self):
        """od.clear() -> None.  Remove all items from od."""
        try:
            for node in itervalues(self.__map):
                del node[:]
            root = self.__root
            root[:] = [root, root, None]
            self.__map.clear()
        except AttributeError:
            pass
        dict.clear(self)

    def popitem(self, last=True):
        """od.popitem() -> (k, v), return and remove a (key, value) pair.

        Pairs are returned in LIFO order if last is true or FIFO order if
        false.
        """
        if not self:
            raise KeyError('dictionary is empty')
        root = self.__root
        if last:
            link = root[0]
            link_prev = link[0]
            link_prev[1] = root
            root[0] = link_prev
        else:
            link = root[1]
            link_next = link[1]
            root[1] = link_next
            link_next[0] = root
        key = link[2]
        del self.__map[key]
        value = dict.pop(self, key)
        return key, value

    # -- the following methods do not depend on the internal structure --

    def keys(self):
        """od.keys() -> list of keys in od"""
        return list(self)

    def values(self):
        """od.values() -> list of values in od"""
        return [self[key] for key in self]

    def items(self):
        """od.items() -> list of (key, value) pairs in od"""
        return [(key, self[key]) for key in self]

    def iterkeys(self):
        """od.iterkeys() -> an iterator over the keys in od"""
        return iter(self)

    def itervalues(self):
        """od.itervalues -> an iterator over the values in od"""
        for k in self:
            yield self[k]

    def iteritems(self):
        """od.iteritems -> an iterator over the (key, value) items in od"""
        for k in self:
            yield (k, self[k])

    def update(*args, **kwds):
        """od.update(E, **F) -> None.  Update od from dict/iterable E and F.

        If E is a dict instance, does:        for k in E: od[k] = E[k]
        If E has a .keys() method, does:      for k in E.keys(): od[k] = E[k]
        Or if E is an iterable of items, does:for k, v in E: od[k] = v
        In either case, this is followed by:  for k, v in F.items(): od[k] = v
        """
        if len(args) > 2:
            raise TypeError('update() takes at most 2 positional '
                            'arguments (%d given)' % (len(args),))
        elif not args:
            raise TypeError('update() takes at least 1 argument (0 given)')
        self = args[0]
        # Make progressively weaker assumptions about "other"
        other = ()
        if len(args) == 2:
            other = args[1]
        if isinstance(other, dict):
            for key in other:
                self[key] = other[key]
        elif hasattr(other, 'keys'):
            for key in other.keys():
                self[key] = other[key]
        else:
            for key, value in other:
                self[key] = value
        for key, value in kwds.items():
            self[key] = value
    # let subclasses override update without breaking __init__
    __update = update

    __marker = object()

    def pop(self, key, default=__marker):
        """od.pop(k[,d]) -> v, remove specified key and return the
        corresponding value.  If key is not found, d is returned if given,
        otherwise KeyError is raised.
        """
        if key in self:
            result = self[key]
            del self[key]
            return result
        if default is self.__marker:
            raise KeyError(key)
        return default

    def setdefault(self, key, default=None):
        """od.setdefault(k[,d]) -> od.get(k,d), also set od[k]=d if k not in od
        """
        if key in self:
            return self[key]
        self[key] = default
        return default

    def __repr__(self, _repr_running={}):
        """od.__repr__() <==> repr(od)"""
        call_key = id(self), _get_ident()
        if call_key in _repr_running:
            return '...'
        _repr_running[call_key] = 1
        try:
            if not self:
                return '%s()' % (self.__class__.__name__,)
            return '%s(%r)' % (self.__class__.__name__, list(self.items()))
        finally:
            del _repr_running[call_key]

    def __reduce__(self):
        """Return state information for pickling"""
        items = [[k, self[k]] for k in self]
        inst_dict = vars(self).copy()
        for k in vars(OrderedDict()):
            inst_dict.pop(k, None)
        if inst_dict:
            return (self.__class__, (items,), inst_dict)
        return self.__class__, (items,)

    def copy(self):
        """od.copy() -> a shallow copy of od"""
        return self.__class__(self)

    @classmethod
    def fromkeys(cls, iterable, value=None):
        """OD.fromkeys(S[, v]) -> New ordered dictionary with keys from S and
        values equal to v (which defaults to None).
        """
        d = cls()
        for key in iterable:
            d[key] = value
        return d

    def __eq__(self, other):
        """od.__eq__(y) <==> od==y.  Comparison to another OD is
        order-sensitive while comparison to a regular mapping is
        order-insensitive.
        """
        if isinstance(other, OrderedDict):
            return (len(self) == len(other) and
                    list(self.items()) == list(other.items()))
        return dict.__eq__(self, other)

    def __ne__(self, other):
        return not self == other

    # -- the following methods are only used in Python 2.7 --

    def viewkeys(self):
        """od.viewkeys() -> a set-like object providing a view on od's keys"""
        return KeysView(self)

    def viewvalues(self):
        """od.viewvalues() -> an object providing a view on od's values"""
        return ValuesView(self)

    def viewitems(self):
        """od.viewitems() -> a set-like object providing a view on od's items
        """
        return ItemsView(self)


# {{{ http://code.activestate.com/recipes/576611/ (r11)

try:
    from operator import itemgetter
    from heapq import nlargest
except ImportError:
    pass


class _Counter(dict):

    """Dict subclass for counting hashable objects.  Sometimes called a bag
    or multiset.  Elements are stored as dictionary keys and their counts
    are stored as dictionary values.

    >>> Counter('zyzygy')
    Counter({'y': 3, 'z': 2, 'g': 1})

    """

    def __init__(self, iterable=None, **kwds):
        """Create a new, empty Counter object.  And if given, count elements
        from an input iterable.  Or, initialize the count from another mapping
        of elements to their counts.

        >>> c = Counter()                    # a new, empty counter
        >>> c = Counter('gallahad')          # a new counter from an iterable
        >>> c = Counter({'a': 4, 'b': 2})    # a new counter from a mapping
        >>> c = Counter(a=4, b=2)            # a new counter from keyword args

        """
        self.update(iterable, **kwds)

    def __missing__(self, key):
        return 0

    def most_common(self, n=None):
        """List the n most common elements and their counts from the most
        common to the least.  If n is None, then list all element counts.

        >>> Counter('abracadabra').most_common(3)
        [('a', 5), ('r', 2), ('b', 2)]

        """
        if n is None:
            return sorted(iteritems(self), key=itemgetter(1), reverse=True)
        return nlargest(n, iteritems(self), key=itemgetter(1))

    def elements(self):
        """Iterator over elements repeating each as many times as its count.

        >>> c = Counter('ABCABC')
        >>> sorted(c.elements())
        ['A', 'A', 'B', 'B', 'C', 'C']

        If an element's count has been set to zero or is a negative number,
        elements() will ignore it.

        """
        for elem, count in iteritems(self):
            for _ in range(count):
                yield elem

    # Override dict methods where the meaning changes for Counter objects.

    @classmethod
    def fromkeys(cls, iterable, v=None):
        raise NotImplementedError(
            'Counter.fromkeys() is undefined.  Use Counter(iterable) instead.')

    def update(self, iterable=None, **kwds):
        """Like dict.update() but add counts instead of replacing them.

        Source can be an iterable, a dictionary, or another Counter instance.

        >>> c = Counter('which')
        >>> c.update('witch')           # add elements from another iterable
        >>> d = Counter('watch')
        >>> c.update(d)                 # add elements from another counter
        >>> c['h']                      # four 'h' in which, witch, and watch
        4

        """
        if iterable is not None:
            if hasattr(iterable, 'iteritems'):
                if self:
                    self_get = self.get
                    for elem, count in iteritems(iterable):
                        self[elem] = self_get(elem, 0) + count
                else:
                    dict.update(
                        self, iterable)  # fast path when counter is empty
            else:
                self_get = self.get
                for elem in iterable:
                    self[elem] = self_get(elem, 0) + 1
        if kwds:
            self.update(kwds)

    def copy(self):
        """Like dict.copy() but returns a Counter instance instead of a dict.
        """
        return Counter(self)

    def __delitem__(self, elem):
        """Like dict.__delitem__() but does not raise KeyError for missing
        values.
        """
        if elem in self:
            dict.__delitem__(self, elem)

    def __repr__(self):
        if not self:
            return '%s()' % self.__class__.__name__
        items = ', '.join(map('%r: %r'.__mod__, self.most_common()))
        return '%s({%s})' % (self.__class__.__name__, items)

    # Multiset-style mathematical operations discussed in:
    #       Knuth TAOCP Volume II section 4.6.3 exercise 19
    #       and at http://en.wikipedia.org/wiki/Multiset
    #
    # Outputs guaranteed to only include positive counts.
    #
    # To strip negative and zero counts, add-in an empty counter:
    #       c += Counter()

    def __add__(self, other):
        """Add counts from two counters.

        >>> Counter('abbb') + Counter('bcc')
        Counter({'b': 4, 'c': 2, 'a': 1})

        """
        if not isinstance(other, Counter):
            return NotImplemented
        result = Counter()
        for elem in set(self) | set(other):
            newcount = self[elem] + other[elem]
            if newcount > 0:
                result[elem] = newcount
        return result

    def __sub__(self, other):
        """Subtract count, but keep only results with positive counts.

        >>> Counter('abbbc') - Counter('bccd')
        Counter({'b': 2, 'a': 1})

        """
        if not isinstance(other, Counter):
            return NotImplemented
        result = Counter()
        for elem in set(self) | set(other):
            newcount = self[elem] - other[elem]
            if newcount > 0:
                result[elem] = newcount
        return result

    def __or__(self, other):
        """Union is the maximum of value in either of the input counters.

        >>> Counter('abbb') | Counter('bcc')
        Counter({'b': 3, 'c': 2, 'a': 1})

        """
        if not isinstance(other, Counter):
            return NotImplemented
        _max = max
        result = Counter()
        for elem in set(self) | set(other):
            newcount = _max(self[elem], other[elem])
            if newcount > 0:
                result[elem] = newcount
        return result

    def __and__(self, other):
        """Intersection is the minimum of corresponding counts.

        >>> Counter('abbb') & Counter('bcc')
        Counter({'b': 1})

        """
        if not isinstance(other, Counter):
            return NotImplemented
        _min = min
        result = Counter()
        if len(self) < len(other):
            self, other = other, self
        for elem in filter(self.__contains__, other):
            newcount = _min(self[elem], other[elem])
            if newcount > 0:
                result[elem] = newcount
        return result

if sys.version_info[:2] < (2, 7):
    OrderedDict = _OrderedDict
    Counter = _Counter
else:
    from collections import OrderedDict, Counter

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


# http://stackoverflow.com/questions/4126348
# Thanks to @martineau at SO

from dateutil import parser as _date_parser
import dateutil
if LooseVersion(dateutil.__version__) < '2.0':
    @functools.wraps(_date_parser.parse)
    def parse_date(timestr, *args, **kwargs):
        timestr = bytes(timestr)
        return _date_parser.parse(timestr, *args, **kwargs)
else:
    parse_date = _date_parser.parse


class OrderedDefaultdict(OrderedDict):

    def __init__(self, *args, **kwargs):
        newdefault = None
        newargs = ()
        if args:
            newdefault = args[0]
            if not (newdefault is None or callable(newdefault)):
                raise TypeError('first argument must be callable or None')
            newargs = args[1:]
        self.default_factory = newdefault
        super(self.__class__, self).__init__(*newargs, **kwargs)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):  # optional, for pickle support
        args = self.default_factory if self.default_factory else tuple()
        return type(self), args, None, None, list(self.items())
