from collections import ChainMap
from typing import List, MutableMapping, TypeVar

_T = TypeVar("_T")


class DeepChainMap(ChainMap):
    maps: List[MutableMapping]  # type: ignore

    def __setitem__(self, key, value):
        for mapping in self.maps:
            if key in mapping:
                mapping[key] = value
                return
        self.maps[0][key] = value

    def __delitem__(self, key):
        for mapping in self.maps:
            if key in mapping:
                del mapping[key]
                return
        raise KeyError(key)

    # FIXME: return type of new_child incorrect in typeshed
    def new_child(self: _T, m) -> _T:  # type: ignore
        return super().new_child(m)  # type: ignore
