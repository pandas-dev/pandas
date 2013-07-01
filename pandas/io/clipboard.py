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
    if kwargs.get('sep') is None and kwargs.get('delim_whitespace') is None:
        kwargs['sep'] = '\s+'
    from pandas.util.clipboard import clipboard_get
    from pandas.io.parsers import read_table
    text = clipboard_get()
    return read_table(StringIO(text), **kwargs)


def to_clipboard(obj):  # pragma: no cover
    """
    Attempt to write text representation of object to the system clipboard

    Notes
    -----
    Requirements for your platform
      - Linux: xclip, or xsel (with gtk or PyQt4 modules)
      - Windows:
      - OS X:
    """
    from pandas.util.clipboard import clipboard_set
    clipboard_set(str(obj))


