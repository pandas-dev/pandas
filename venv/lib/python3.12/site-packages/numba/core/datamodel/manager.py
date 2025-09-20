import weakref
from collections import ChainMap

from numba.core import types


class DataModelManager(object):
    """Manages mapping of FE types to their corresponding data model
    """

    def __init__(self, handlers=None):
        """
        Parameters
        -----------
        handlers: Mapping[Type, DataModel] or None
            Optionally provide the initial handlers mapping.
        """
        # { numba type class -> model factory }
        self._handlers = handlers or {}
        # { numba type instance -> model instance }
        self._cache = weakref.WeakKeyDictionary()

    def register(self, fetypecls, handler):
        """Register the datamodel factory corresponding to a frontend-type class
        """
        assert issubclass(fetypecls, types.Type)
        self._handlers[fetypecls] = handler

    def lookup(self, fetype):
        """Returns the corresponding datamodel given the frontend-type instance
        """
        try:
            return self._cache[fetype]
        except KeyError:
            pass
        handler = self._handlers[type(fetype)]
        model = self._cache[fetype] = handler(self, fetype)
        return model

    def __getitem__(self, fetype):
        """Shorthand for lookup()
        """
        return self.lookup(fetype)

    def copy(self):
        """
        Make a copy of the manager.
        Use this to inherit from the default data model and specialize it
        for custom target.
        """
        return DataModelManager(self._handlers.copy())

    def chain(self, other_manager):
        """Create a new DataModelManager by chaining the handlers mapping of
        `other_manager` with a fresh handlers mapping.

        Any existing and new handlers inserted to `other_manager` will be
        visible to the new manager. Any handlers inserted to the new manager
        can override existing handlers in `other_manager` without actually
        mutating `other_manager`.

        Parameters
        ----------
        other_manager: DataModelManager
        """
        chained = ChainMap(self._handlers, other_manager._handlers)
        return DataModelManager(chained)

