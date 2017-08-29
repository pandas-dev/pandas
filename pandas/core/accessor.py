# -*- coding: utf-8 -*-
"""

accessor.py contains base classes for implementing accessor properties
that can be mixed into or pinned onto other pandas classes.

"""


class DirNamesMixin(object):
    _accessors = frozenset([])

    def _dir_deletions(self):
        """ delete unwanted __dir__ for this object """
        return self._accessors

    def _dir_additions(self):
        """ add addtional __dir__ for this object """
        rv = set()
        for accessor in self._accessors:
            try:
                getattr(self, accessor)
                rv.add(accessor)
            except AttributeError:
                pass
        return rv

    def __dir__(self):
        """
        Provide method name lookup and completion
        Only provide 'public' methods
        """
        rv = set(dir(type(self)))
        rv = (rv - self._dir_deletions()) | self._dir_additions()
        return sorted(rv)
