from typing import ChainMap, List, MutableMapping, TypeVar

from pandas._typing import T

_KT = TypeVar("_KT")
_VT = TypeVar("_VT")


class DeepChainMap(ChainMap[_KT, _VT]):
    """
    Variant of ChainMap that allows direct updates to inner scopes.

    Only works when all passed mapping are mutable.
    """

    # error: Incompatible types in assignment (expression has type
    #  "List[MutableMapping[_KT, _VT]]", base class "ChainMap" defined the type
    #  as "List[Mapping[_KT, _VT]]")  [assignment]
    maps: List[MutableMapping[_KT, _VT]]  # type: ignore

    def __setitem__(self, key: _KT, value: _VT) -> None:
        for mapping in self.maps:
            if key in mapping:
                mapping[key] = value
                return
        self.maps[0][key] = value

    def __delitem__(self, key: _KT) -> None:
        """
        Raises
        ------
        KeyError
            If `key` doesn't exist.
        """
        for mapping in self.maps:
            if key in mapping:
                del mapping[key]
                return
        raise KeyError(key)

    # FIXME: return type of new_child incorrect in typeshed
    def new_child(self: T, m) -> T:  # type: ignore
        return super().new_child(m)  # type: ignore
