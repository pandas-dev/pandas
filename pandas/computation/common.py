import collections
from pandas.core.common import is_string


def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not is_string(el):
            for s in flatten(el):
                yield s
        else:
            yield el
