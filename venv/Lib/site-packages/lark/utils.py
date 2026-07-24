import unicodedata
import os
from itertools import product
from collections import deque
from typing import Callable, Iterator, List, Optional, Tuple, Type, TypeVar, Union, Dict, Any, Sequence, Iterable, AbstractSet

###{standalone
import sys, re
import logging
from dataclasses import dataclass
from typing import Generic, AnyStr

logger: logging.Logger = logging.getLogger("lark")
logger.addHandler(logging.StreamHandler())
# Set to highest level, since we have some warnings amongst the code
# By default, we should not output any log messages
logger.setLevel(logging.CRITICAL)


NO_VALUE = object()

T = TypeVar("T")


def classify(seq: Iterable, key: Optional[Callable] = None, value: Optional[Callable] = None) -> Dict:
    d: Dict[Any, Any] = {}
    for item in seq:
        k = key(item) if (key is not None) else item
        v = value(item) if (value is not None) else item
        try:
            d[k].append(v)
        except KeyError:
            d[k] = [v]
    return d


def _deserialize(data: Any, namespace: Dict[str, Any], memo: Dict) -> Any:
    if isinstance(data, dict):
        if '__type__' in data:  # Object
            class_ = namespace[data['__type__']]
            return class_.deserialize(data, memo)
        elif '@' in data:
            return memo[data['@']]
        return {key:_deserialize(value, namespace, memo) for key, value in data.items()}
    elif isinstance(data, list):
        return [_deserialize(value, namespace, memo) for value in data]
    return data


_T = TypeVar("_T", bound="Serialize")

class Serialize:
    """Safe-ish serialization interface that doesn't rely on Pickle

    Attributes:
        __serialize_fields__ (List[str]): Fields (aka attributes) to serialize.
        __serialize_namespace__ (list): List of classes that deserialization is allowed to instantiate.
                                        Should include all field types that aren't builtin types.
    """

    def memo_serialize(self, types_to_memoize: List) -> Any:
        memo = SerializeMemoizer(types_to_memoize)
        return self.serialize(memo), memo.serialize()

    def serialize(self, memo = None) -> Dict[str, Any]:
        if memo and memo.in_types(self):
            return {'@': memo.memoized.get(self)}

        fields = getattr(self, '__serialize_fields__')
        res = {f: _serialize(getattr(self, f), memo) for f in fields}
        res['__type__'] = type(self).__name__
        if hasattr(self, '_serialize'):
            self._serialize(res, memo)
        return res

    @classmethod
    def deserialize(cls: Type[_T], data: Dict[str, Any], memo: Dict[int, Any]) -> _T:
        namespace = getattr(cls, '__serialize_namespace__', [])
        namespace = {c.__name__:c for c in namespace}

        fields = getattr(cls, '__serialize_fields__')

        if '@' in data:
            return memo[data['@']]

        inst = cls.__new__(cls)
        for f in fields:
            try:
                setattr(inst, f, _deserialize(data[f], namespace, memo))
            except KeyError as e:
                raise KeyError("Cannot find key for class", cls, e)

        if hasattr(inst, '_deserialize'):
            inst._deserialize()

        return inst


class SerializeMemoizer(Serialize):
    "A version of serialize that memoizes objects to reduce space"

    __serialize_fields__ = 'memoized',

    def __init__(self, types_to_memoize: List) -> None:
        self.types_to_memoize = tuple(types_to_memoize)
        self.memoized = Enumerator()

    def in_types(self, value: Serialize) -> bool:
        return isinstance(value, self.types_to_memoize)

    def serialize(self) -> Dict[int, Any]:  # type: ignore[override]
        return _serialize(self.memoized.reversed(), None)

    @classmethod
    def deserialize(cls, data: Dict[int, Any], namespace: Dict[str, Any], memo: Dict[Any, Any]) -> Dict[int, Any]:  # type: ignore[override]
        return _deserialize(data, namespace, memo)


try:
    import regex
    _has_regex = True
except ImportError:
    _has_regex = False

if sys.version_info >= (3, 11):
    import re._parser as sre_parse
    import re._constants as sre_constants
else:
    import sre_parse
    import sre_constants

categ_pattern = re.compile(r'\\p{[A-Za-z_]+}')

def get_regexp_width(expr: str) -> Union[Tuple[int, int], List[int]]:
    if _has_regex:
        # Since `sre_parse` cannot deal with Unicode categories of the form `\p{Mn}`, we replace these with
        # a simple letter, which makes no difference as we are only trying to get the possible lengths of the regex
        # match here below.
        regexp_final = re.sub(categ_pattern, 'A', expr)
    else:
        if re.search(categ_pattern, expr):
            raise ImportError('`regex` module must be installed in order to use Unicode categories.', expr)
        regexp_final = expr
    try:
        # Fixed in next version (past 0.960) of typeshed
        return [int(x) for x in sre_parse.parse(regexp_final).getwidth()]
    except sre_constants.error:
        if not _has_regex:
            raise ValueError(expr)
        else:
            # sre_parse does not support the new features in regex. To not completely fail in that case,
            # we manually test for the most important info (whether the empty string is matched)
            c = regex.compile(regexp_final)
            # Python 3.11.7 introducded sre_parse.MAXWIDTH that is used instead of MAXREPEAT
            # See lark-parser/lark#1376 and python/cpython#109859
            MAXWIDTH = getattr(sre_parse, "MAXWIDTH", sre_constants.MAXREPEAT)
            if c.match('') is None:
                # MAXREPEAT is a none pickable subclass of int, therefore needs to be converted to enable caching
                return 1, int(MAXWIDTH)
            else:
                return 0, int(MAXWIDTH)


@dataclass(frozen=True)
class TextSlice(Generic[AnyStr]):
    """A view of a string or bytes object, between the start and end indices.

    Never creates a copy.

    Lark accepts instances of TextSlice as input (instead of a string),
    when the lexer is 'basic' or 'contextual'.

    Args:
        text (str or bytes): The text to slice.
        start (int): The start index. Negative indices are supported.
        end (int): The end index. Negative indices are supported.

    Raises:
        TypeError: If `text` is not a `str` or `bytes`.
        AssertionError: If `start` or `end` are out of bounds.

    Examples:
        >>> TextSlice("Hello, World!", 7, -1)
        TextSlice(text='Hello, World!', start=7, end=12)

        >>> TextSlice("Hello, World!", 7, None).count("o")
        1

    """
    text: AnyStr
    start: int
    end: int

    def __post_init__(self):
        if not isinstance(self.text, (str, bytes)):
            raise TypeError("text must be str or bytes")

        if self.start < 0:
            object.__setattr__(self, 'start', self.start + len(self.text))
            assert self.start >=0

        if self.end is None:
            object.__setattr__(self, 'end', len(self.text))
        elif self.end < 0:
            object.__setattr__(self, 'end', self.end + len(self.text))
            assert self.end <= len(self.text)

    @classmethod
    def cast_from(cls, text: 'TextOrSlice') -> 'TextSlice[AnyStr]':
        if isinstance(text, TextSlice):
            return text

        return cls(text, 0, len(text))

    def is_complete_text(self):
        return self.start == 0 and self.end == len(self.text)

    def __len__(self):
        return self.end - self.start

    def count(self, substr: AnyStr):
        return self.text.count(substr, self.start, self.end)

    def rindex(self, substr: AnyStr):
        return self.text.rindex(substr, self.start, self.end)


TextOrSlice = Union[AnyStr, 'TextSlice[AnyStr]']
LarkInput = Union[AnyStr, TextSlice[AnyStr], Any]

###}


_ID_START =    'Lu', 'Ll', 'Lt', 'Lm', 'Lo', 'Mn', 'Mc', 'Pc'
_ID_CONTINUE = _ID_START + ('Nd', 'Nl',)

def _test_unicode_category(s: str, categories: Sequence[str]) -> bool:
    if len(s) != 1:
        return all(_test_unicode_category(char, categories) for char in s)
    return s == '_' or unicodedata.category(s) in categories

def is_id_continue(s: str) -> bool:
    """
    Checks if all characters in `s` are alphanumeric characters (Unicode standard, so diacritics, indian vowels, non-latin
    numbers, etc. all pass). Synonymous with a Python `ID_CONTINUE` identifier. See PEP 3131 for details.
    """
    return _test_unicode_category(s, _ID_CONTINUE)

def is_id_start(s: str) -> bool:
    """
    Checks if all characters in `s` are alphabetic characters (Unicode standard, so diacritics, indian vowels, non-latin
    numbers, etc. all pass). Synonymous with a Python `ID_START` identifier. See PEP 3131 for details.
    """
    return _test_unicode_category(s, _ID_START)


def dedup_list(l: Iterable[T]) -> List[T]:
    """Given a list (l) will removing duplicates from the list,
       preserving the original order of the list. Assumes that
       the list entries are hashable."""
    return list(dict.fromkeys(l))


class Enumerator(Serialize):
    def __init__(self) -> None:
        self.enums: Dict[Any, int] = {}

    def get(self, item) -> int:
        if item not in self.enums:
            self.enums[item] = len(self.enums)
        return self.enums[item]

    def __len__(self):
        return len(self.enums)

    def reversed(self) -> Dict[int, Any]:
        r = {v: k for k, v in self.enums.items()}
        assert len(r) == len(self.enums)
        return r



def combine_alternatives(lists):
    """
    Accepts a list of alternatives, and enumerates all their possible concatenations.

    Examples:
        >>> combine_alternatives([range(2), [4,5]])
        [[0, 4], [0, 5], [1, 4], [1, 5]]

        >>> combine_alternatives(["abc", "xy", '$'])
        [['a', 'x', '$'], ['a', 'y', '$'], ['b', 'x', '$'], ['b', 'y', '$'], ['c', 'x', '$'], ['c', 'y', '$']]

        >>> combine_alternatives([])
        [[]]
    """
    if not lists:
        return [[]]
    assert all(l for l in lists), lists
    return list(product(*lists))

try:
    import atomicwrites
    _has_atomicwrites = True
except ImportError:
    _has_atomicwrites = False

class FS:
    exists = staticmethod(os.path.exists)

    @staticmethod
    def open(name, mode="r", **kwargs):
        if _has_atomicwrites and "w" in mode:
            return atomicwrites.atomic_write(name, mode=mode, overwrite=True, **kwargs)
        else:
            return open(name, mode, **kwargs)


class fzset(frozenset):
    def __repr__(self):
        return '{%s}' % ', '.join(map(repr, self))


def classify_bool(seq: Iterable, pred: Callable) -> Any:
    false_elems = []
    true_elems = [elem for elem in seq if pred(elem) or false_elems.append(elem)]  # type: ignore[func-returns-value]
    return true_elems, false_elems


def bfs(initial: Iterable, expand: Callable) -> Iterator:
    open_q = deque(list(initial))
    visited = set(open_q)
    while open_q:
        node = open_q.popleft()
        yield node
        for next_node in expand(node):
            if next_node not in visited:
                visited.add(next_node)
                open_q.append(next_node)

def bfs_all_unique(initial, expand):
    "bfs, but doesn't keep track of visited (aka seen), because there can be no repetitions"
    open_q = deque(list(initial))
    while open_q:
        node = open_q.popleft()
        yield node
        open_q += expand(node)


def _serialize(value: Any, memo: Optional[SerializeMemoizer]) -> Any:
    if isinstance(value, Serialize):
        return value.serialize(memo)
    elif isinstance(value, list):
        return [_serialize(elem, memo) for elem in value]
    elif isinstance(value, frozenset):
        return list(value)  # TODO reversible?
    elif isinstance(value, dict):
        return {key:_serialize(elem, memo) for key, elem in value.items()}
    # assert value is None or isinstance(value, (int, float, str, tuple)), value
    return value




def small_factors(n: int, max_factor: int) -> List[Tuple[int, int]]:
    """
    Splits n up into smaller factors and summands <= max_factor.
    Returns a list of [(a, b), ...]
    so that the following code returns n:

    n = 1
    for a, b in values:
        n = n * a + b

    Currently, we also keep a + b <= max_factor, but that might change
    """
    assert n >= 0
    assert max_factor > 2
    if n <= max_factor:
        return [(n, 0)]

    for a in range(max_factor, 1, -1):
        r, b = divmod(n, a)
        if a + b <= max_factor:
            return small_factors(r, max_factor) + [(a, b)]
    assert False, "Failed to factorize %s" % n


class OrderedSet(AbstractSet[T]):
    """A minimal OrderedSet implementation, using a dictionary.

    (relies on the dictionary being ordered)
    """
    def __init__(self, items: Iterable[T] =()):
        self.d = dict.fromkeys(items)

    def __contains__(self, item: Any) -> bool:
        return item in self.d

    def add(self, item: T):
        self.d[item] = None

    def __iter__(self) -> Iterator[T]:
        return iter(self.d)

    def remove(self, item: T):
        del self.d[item]

    def __bool__(self):
        return bool(self.d)

    def __len__(self) -> int:
        return len(self.d)

    def __repr__(self):
        return f"{type(self).__name__}({', '.join(map(repr,self))})"
