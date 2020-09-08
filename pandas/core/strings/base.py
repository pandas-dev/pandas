import unicodedata


class BaseStringArrayMethods:
    """Base class for it."""

    def __init__(self, array):
        self._array = array

    def upper(self):
        return self._map(lambda x: x.upper())

    def isalnum(self):
        return self._map(str.isalnum, dtype="bool")

    def capitalize(self):
        return self._map(str.capitalize)

    def casefold(self):
        return self._map(str.casefold)

    def title(self):
        return self._map(str.title)

    def swapcase(self):
        return self._map(str.swapcase)

    def lower(self):
        return self._map(str.lower)

    def normalize(self, form):
        f = lambda x: unicodedata.normalize(form, x)
        return self._map(f)
