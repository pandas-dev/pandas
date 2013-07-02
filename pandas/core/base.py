"""
Base class(es) for all pandas objects.
"""
from pandas.util import py3compat

class StringMixin(object):
    """implements string methods so long as object defines a `__unicode__` method.
    Handles Python2/3 compatibility transparently."""
    # side note - this could be made into a metaclass if more than one object nees
    def __str__(self):
        """
        Return a string representation for a particular object.

        Invoked by str(obj) in both py2/py3.
        Yields Bytestring in Py2, Unicode String in py3.
        """

        if py3compat.PY3:
            return self.__unicode__()
        return self.__bytes__()

    def __bytes__(self):
        """
        Return a string representation for a particular object.

        Invoked by bytes(obj) in py3 only.
        Yields a bytestring in both py2/py3.
        """
        from pandas.core.config import get_option

        encoding = get_option("display.encoding")
        return self.__unicode__().encode(encoding, 'replace')

    def __repr__(self):
        """
        Return a string representation for a particular object.

        Yields Bytestring in Py2, Unicode String in py3.
        """
        return str(self)

class PandasObject(StringMixin):
    """baseclass for various pandas objects"""

    @property
    def _constructor(self):
        """class constructor (for this class it's just `__class__`"""
        return self.__class__

    def __unicode__(self):
        """
        Return a string representation for a particular object.

        Invoked by unicode(obj) in py2 only. Yields a Unicode String in both
        py2/py3.
        """
        # Should be overwritten by base classes
        return object.__repr__(self)
