"""Proxy dictionary for objects stored in a container."""
import weakref


class ProxyDict(dict):
    """A dictionary which uses a container object to store its values."""

    def __init__(self, container):
        self.containerref = weakref.ref(container)
        """A weak reference to the container object.

        .. versionchanged:: 3.0
           The *containerRef* attribute has been renamed into
           *containerref*.

        """

    def __getitem__(self, key):
        if key not in self:
            raise KeyError(key)

        # Values are not actually stored to avoid extra references.
        return self._get_value_from_container(self._get_container(), key)

    def __setitem__(self, key, value):
        # Values are not actually stored to avoid extra references.
        super().__setitem__(key, None)

    def __repr__(self):
        return object.__repr__(self)

    def __str__(self):
        # C implementation does not use `self.__getitem__()`. :(
        return '{' + ", ".join("{k!r}: {v!r}" for k, v in self.items()) + '}'

    def values(self):
        # C implementation does not use `self.__getitem__()`. :(
        return [self[key] for key in self.keys()]

    def itervalues(self):
        # C implementation does not use `self.__getitem__()`. :(
        for key in self.keys():
            yield self[key]

    def items(self):
        # C implementation does not use `self.__getitem__()`. :(
        return [(key, self[key]) for key in self.keys()]

    def iteritems(self):
        # C implementation does not use `self.__getitem__()`. :(
        for key in self.keys():
            yield (key, self[key])

    def _get_container(self):
        container = self.containerref()
        if container is None:
            raise ValueError("the container object does no longer exist")
        return container
