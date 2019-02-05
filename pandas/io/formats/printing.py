"""
printing tools
"""

import sys

from pandas import compat
from pandas.core.config import get_option, pprint_thing


def adjoin(space, *lists, **kwargs):
    """
    Glues together two sets of strings using the amount of space requested.
    The idea is to prettify.

    ----------
    space : int
        number of spaces for padding
    lists : str
        list of str which being joined
    strlen : callable
        function used to calculate the length of each str. Needed for unicode
        handling.
    justfunc : callable
        function used to justify str. Needed for unicode handling.
    """
    strlen = kwargs.pop('strlen', len)
    justfunc = kwargs.pop('justfunc', justify)

    out_lines = []
    newLists = []
    lengths = [max(map(strlen, x)) + space for x in lists[:-1]]
    # not the last one
    lengths.append(max(map(len, lists[-1])))
    maxLen = max(map(len, lists))
    for i, lst in enumerate(lists):
        nl = justfunc(lst, lengths[i], mode='left')
        nl.extend([' ' * lengths[i]] * (maxLen - len(lst)))
        newLists.append(nl)
    toJoin = zip(*newLists)
    for lines in toJoin:
        out_lines.append(_join_unicode(lines))
    return _join_unicode(out_lines, sep='\n')


def justify(texts, max_len, mode='right'):
    """
    Perform ljust, center, rjust against string or list-like
    """
    if mode == 'left':
        return [x.ljust(max_len) for x in texts]
    elif mode == 'center':
        return [x.center(max_len) for x in texts]
    else:
        return [x.rjust(max_len) for x in texts]


def _join_unicode(lines, sep=''):
    try:
        return sep.join(lines)
    except UnicodeDecodeError:
        sep = compat.text_type(sep)
        return sep.join([x.decode('utf-8') if isinstance(x, str) else x
                         for x in lines])


def pprint_thing_encoded(object, encoding='utf-8', errors='replace', **kwds):
    value = pprint_thing(object)  # get unicode representation of object
    return value.encode(encoding, errors, **kwds)


def _enable_data_resource_formatter(enable):
    if 'IPython' not in sys.modules:
        # definitely not in IPython
        return
    from IPython import get_ipython
    ip = get_ipython()
    if ip is None:
        # still not in IPython
        return

    formatters = ip.display_formatter.formatters
    mimetype = "application/vnd.dataresource+json"

    if enable:
        if mimetype not in formatters:
            # define tableschema formatter
            from IPython.core.formatters import BaseFormatter

            class TableSchemaFormatter(BaseFormatter):
                print_method = '_repr_data_resource_'
                _return_type = (dict,)
            # register it:
            formatters[mimetype] = TableSchemaFormatter()
        # enable it if it's been disabled:
        formatters[mimetype].enabled = True
    else:
        # unregister tableschema mime-type
        if mimetype in formatters:
            formatters[mimetype].enabled = False


default_pprint = lambda x, max_seq_items=None: \
    pprint_thing(x, escape_chars=('\t', '\r', '\n'), quote_strings=True,
                 max_seq_items=max_seq_items)


def format_object_summary(obj, formatter, is_justify=True, name=None,
                          indent_for_name=True):
    """
    Return the formatted obj as a unicode string

    Parameters
    ----------
    obj : object
        must be iterable and support __getitem__
    formatter : callable
        string formatter for an element
    is_justify : boolean
        should justify the display
    name : name, optional
        defaults to the class name of the obj
    indent_for_name : bool, default True
        Whether subsequent lines should be be indented to
        align with the name.

    Returns
    -------
    summary string

    """
    from pandas.io.formats.console import get_console_size
    from pandas.io.formats.format import _get_adjustment

    display_width, _ = get_console_size()
    if display_width is None:
        display_width = get_option('display.width') or 80
    if name is None:
        name = obj.__class__.__name__

    if indent_for_name:
        name_len = len(name)
        space1 = "\n%s" % (' ' * (name_len + 1))
        space2 = "\n%s" % (' ' * (name_len + 2))
    else:
        space1 = "\n"
        space2 = "\n "  # space for the opening '['

    n = len(obj)
    sep = ','
    max_seq_items = get_option('display.max_seq_items') or n

    # are we a truncated display
    is_truncated = n > max_seq_items

    # adj can optionally handle unicode eastern asian width
    adj = _get_adjustment()

    def _extend_line(s, line, value, display_width, next_line_prefix):

        if (adj.len(line.rstrip()) + adj.len(value.rstrip()) >=
                display_width):
            s += line.rstrip()
            line = next_line_prefix
        line += value
        return s, line

    def best_len(values):
        if values:
            return max(adj.len(x) for x in values)
        else:
            return 0

    close = u', '

    if n == 0:
        summary = u'[]{}'.format(close)
    elif n == 1:
        first = formatter(obj[0])
        summary = u'[{}]{}'.format(first, close)
    elif n == 2:
        first = formatter(obj[0])
        last = formatter(obj[-1])
        summary = u'[{}, {}]{}'.format(first, last, close)
    else:

        if n > max_seq_items:
            n = min(max_seq_items // 2, 10)
            head = [formatter(x) for x in obj[:n]]
            tail = [formatter(x) for x in obj[-n:]]
        else:
            head = []
            tail = [formatter(x) for x in obj]

        # adjust all values to max length if needed
        if is_justify:

            # however, if we are not truncated and we are only a single
            # line, then don't justify
            if (is_truncated or
                    not (len(', '.join(head)) < display_width and
                         len(', '.join(tail)) < display_width)):
                max_len = max(best_len(head), best_len(tail))
                head = [x.rjust(max_len) for x in head]
                tail = [x.rjust(max_len) for x in tail]

        summary = ""
        line = space2

        for i in range(len(head)):
            word = head[i] + sep + ' '
            summary, line = _extend_line(summary, line, word,
                                         display_width, space2)

        if is_truncated:
            # remove trailing space of last line
            summary += line.rstrip() + space2 + '...'
            line = space2

        for i in range(len(tail) - 1):
            word = tail[i] + sep + ' '
            summary, line = _extend_line(summary, line, word,
                                         display_width, space2)

        # last value: no sep added + 1 space of width used for trailing ','
        summary, line = _extend_line(summary, line, tail[-1],
                                     display_width - 2, space2)
        summary += line

        # right now close is either '' or ', '
        # Now we want to include the ']', but not the maybe space.
        close = ']' + close.rstrip(' ')
        summary += close

        if len(summary) > (display_width):
            summary += space1
        else:  # one row
            summary += ' '

        # remove initial space
        summary = '[' + summary[len(space2):]

    return summary


def format_object_attrs(obj):
    """
    Return a list of tuples of the (attr, formatted_value)
    for common attrs, including dtype, name, length

    Parameters
    ----------
    obj : object
        must be iterable

    Returns
    -------
    list

    """
    attrs = []
    if hasattr(obj, 'dtype'):
        attrs.append(('dtype', "'{}'".format(obj.dtype)))
    if getattr(obj, 'name', None) is not None:
        attrs.append(('name', default_pprint(obj.name)))
    max_seq_items = get_option('display.max_seq_items') or len(obj)
    if len(obj) > max_seq_items:
        attrs.append(('length', len(obj)))
    return attrs
