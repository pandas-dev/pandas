""" io on the clipboard """
from pandas import compat, get_option
from pandas.compat import StringIO

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

    # try to decode (if needed on PY3)
    if compat.PY3:
        try:
            text = compat.bytes_to_str(text,encoding=kwargs.get('encoding') or get_option('display.encoding'))
        except:
            pass
    return read_table(StringIO(text), **kwargs)


def to_clipboard(obj, sep=None, **kwargs):  # pragma: no cover
    """
    Attempt to write text representation of object to the system clipboard
    The clipboard can be then pasted into Excel for example.

    Notes
    -----
    Requirements for your platform
      - Linux: xclip, or xsel (with gtk or PyQt4 modules)
      - Windows:
      - OS X:
    """
    from pandas.util.clipboard import clipboard_set
    try:
        if sep is None:
            sep = '\t'
        buf = StringIO()
        obj.to_csv(buf,sep=sep, **kwargs)
        clipboard_set(buf.getvalue())
    except:
        clipboard_set(str(obj))

