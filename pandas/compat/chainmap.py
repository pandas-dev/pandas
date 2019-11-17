from collections import ChainMap


class DeepChainMap(ChainMap):
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

    def new_child(self, *args, **kwargs) -> "DeepChainMap":
        # ChainMap.new_child returns self.__class__(...) but mypy
        #  doesn't know that, so we annotate it explicitly here.
        return super().new_child(*args, **kwargs)  # type: ignore
