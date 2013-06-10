""" io on the clipboard """
from StringIO import StringIO

def read_clipboard(**kwargs):  # pragma: no cover
    """
    Read text from clipboard and pass to read_table. See read_table for the
    full argument list

    Returns
    -------
    parsed : DataFrame
    """
    from pandas.util.clipboard import clipboard_get
    from pandas.io.parsers import table
    text = clipboard_get()
    return read_table(StringIO(text), **kwargs)


def to_clipboard(obj):  # pragma: no cover
    """
    Attempt to write text representation of object to the system clipboard

    Notes
    -----
    Requirements for your platform
      - Linux: xsel command line tool
      - Windows: Python win32 extensions
      - OS X:
    """
    from pandas.util.clipboard import clipboard_set
    clipboard_set(str(obj))


