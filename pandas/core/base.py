# pylint: disable=W0231,E1101

from pandas.core import common as com
from pandas.util import py3compat

class PandasObject(object):
    """ The base class for pandas objects """

    #----------------------------------------------------------------------
    # Reconstruction

    def save(self, path):
        com.save(self, path)

    @classmethod
    def load(cls, path):
        return com.load(path)

    #----------------------------------------------------------------------
    # Formatting

    def __unicode__(self):
        raise NotImplementedError

    def __str__(self):
        """
        Return a string representation for a particular Object

        Invoked by str(df) in both py2/py3.
        Yields Bytestring in Py2, Unicode String in py3.
        """

        if py3compat.PY3:
            return self.__unicode__()
        return self.__bytes__()

    def __bytes__(self):
        """
        Return a string representation for a particular Object

        Invoked by bytes(df) in py3 only.
        Yields a bytestring in both py2/py3.
        """
        encoding = com.get_option("display.encoding")
        return self.__unicode__().encode(encoding, 'replace')

    def __repr__(self):
        """
        Return a string representation for a particular Object

        Yields Bytestring in Py2, Unicode String in py3.
        """
        return str(self)
