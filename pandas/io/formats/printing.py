"""
printing tools
"""

import sys

from pandas._config import get_option

from pandas.core.dtypes.inference import is_sequence


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
    strlen = kwargs.pop("strlen", len)
    justfunc = kwargs.pop("justfunc", justify)

    out_lines = []
    newLists = []
    lengths = [max(map(strlen, x)) + space for x in lists[:-1]]
    # not the last one
    lengths.append(max(map(len, lists[-1])))
    maxLen = max(map(len, lists))
    for i, lst in enumerate(lists):
        nl = justfunc(lst, lengths[i], mode="left")
        nl.extend([" " * lengths[i]] * (maxLen - len(lst)))
        newLists.append(nl)
    toJoin = zip(*newLists)
    for lines in toJoin:
        out_lines.append(_join_unicode(lines))
    return _join_unicode(out_lines, sep="\n")


def justify(texts, max_len, mode="right"):
    """
    Perform ljust, center, rjust against string or list-like
    """
    if mode == "left":
        return [x.ljust(max_len) for x in texts]
    elif mode == "center":
        return [x.center(max_len) for x in texts]
    else:
        return [x.rjust(max_len) for x in texts]


def _join_unicode(lines, sep=""):
    try:
        return sep.join(lines)
    except UnicodeDecodeError:
        sep = str(sep)
        return sep.join([x.decode("utf-8") if isinstance(x, str) else x for x in lines])


# Unicode consolidation
# ---------------------
#
# pprinting utility functions for generating Unicode text or
# bytes(3.x)/str(2.x) representations of objects.
# Try to use these as much as possible rather then rolling your own.
#
# When to use
# -----------
#
# 1) If you're writing code internal to pandas (no I/O directly involved),
#    use pprint_thing().
#
#    It will always return unicode text which can handled by other
#    parts of the package without breakage.
#
# 2) if you need to write something out to file, use
#    pprint_thing_encoded(encoding).
#
#    If no encoding is specified, it defaults to utf-8. Since encoding pure
#    ascii with utf-8 is a no-op you can safely use the default utf-8 if you're
#    working with straight ascii.


def _pprint_seq(seq, _nest_lvl=0, max_seq_items=None, **kwds):
    """
    internal. pprinter for iterables. you should probably use pprint_thing()
    rather then calling this directly.

    bounds length of printed sequence, depending on options
    """
    if isinstance(seq, set):
        fmt = "{{{body}}}"
    else:
        fmt = "[{body}]" if hasattr(seq, "__setitem__") else "({body})"

    if max_seq_items is False:
        nitems = len(seq)
    else:
        nitems = max_seq_items or get_option("max_seq_items") or len(seq)

    s = iter(seq)
    # handle sets, no slicing
    r = [
        pprint_thing(next(s), _nest_lvl + 1, max_seq_items=max_seq_items, **kwds)
        for i in range(min(nitems, len(seq)))
    ]
    body = ", ".join(r)

    if nitems < len(seq):
        body += ", ..."
    elif isinstance(seq, tuple) and len(seq) == 1:
        body += ","

    return fmt.format(body=body)


def _pprint_dict(seq, _nest_lvl=0, max_seq_items=None, **kwds):
    """
    internal. pprinter for iterables. you should probably use pprint_thing()
    rather then calling this directly.
    """
    fmt = "{{{things}}}"
    pairs = []

    pfmt = "{key}: {val}"

    if max_seq_items is False:
        nitems = len(seq)
    else:
        nitems = max_seq_items or get_option("max_seq_items") or len(seq)

    for k, v in list(seq.items())[:nitems]:
        pairs.append(
            pfmt.format(
                key=pprint_thing(k, _nest_lvl + 1, max_seq_items=max_seq_items, **kwds),
                val=pprint_thing(v, _nest_lvl + 1, max_seq_items=max_seq_items, **kwds),
            )
        )

    if nitems < len(seq):
        return fmt.format(things=", ".join(pairs) + ", ...")
    else:
        return fmt.format(things=", ".join(pairs))


def pprint_thing(
    thing,
    _nest_lvl=0,
    escape_chars=None,
    default_escapes=False,
    quote_strings=False,
    max_seq_items=None,
):
    """
    This function is the sanctioned way of converting objects
    to a unicode representation.

    properly handles nested sequences containing unicode strings
    (unicode(object) does not)

    Parameters
    ----------
    thing : anything to be formatted
    _nest_lvl : internal use only. pprint_thing() is mutually-recursive
        with pprint_sequence, this argument is used to keep track of the
        current nesting level, and limit it.
    escape_chars : list or dict, optional
        Characters to escape. If a dict is passed the values are the
        replacements
    default_escapes : bool, default False
        Whether the input escape characters replaces or adds to the defaults
    max_seq_items : False, int, default None
        Pass thru to other pretty printers to limit sequence printing

    Returns
    -------
    result - unicode str

    """

    def as_escaped_unicode(thing, escape_chars=escape_chars):
        # Unicode is fine, else we try to decode using utf-8 and 'replace'
        # if that's not it either, we have no way of knowing and the user
        # should deal with it himself.

        try:
            result = str(thing)  # we should try this first
        except UnicodeDecodeError:
            # either utf-8 or we replace errors
            result = str(thing).decode("utf-8", "replace")

        translate = {"\t": r"\t", "\n": r"\n", "\r": r"\r"}
        if isinstance(escape_chars, dict):
            if default_escapes:
                translate.update(escape_chars)
            else:
                translate = escape_chars
            escape_chars = list(escape_chars.keys())
        else:
            escape_chars = escape_chars or tuple()
        for c in escape_chars:
            result = result.replace(c, translate[c])

        return str(result)

    if hasattr(thing, "__next__"):
        return str(thing)
    elif isinstance(thing, dict) and _nest_lvl < get_option(
        "display.pprint_nest_depth"
    ):
        result = _pprint_dict(
            thing, _nest_lvl, quote_strings=True, max_seq_items=max_seq_items
        )
    elif is_sequence(thing) and _nest_lvl < get_option("display.pprint_nest_depth"):
        result = _pprint_seq(
            thing,
            _nest_lvl,
            escape_chars=escape_chars,
            quote_strings=quote_strings,
            max_seq_items=max_seq_items,
        )
    elif isinstance(thing, str) and quote_strings:
        result = "'{thing}'".format(thing=as_escaped_unicode(thing))
    else:
        result = as_escaped_unicode(thing)

    return str(result)  # always unicode


def pprint_thing_encoded(object, encoding="utf-8", errors="replace", **kwds):
    value = pprint_thing(object)  # get unicode representation of object
    return value.encode(encoding, errors, **kwds)


def _enable_data_resource_formatter(enable):
    if "IPython" not in sys.modules:
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
                print_method = "_repr_data_resource_"
                _return_type = (dict,)

            # register it:
            formatters[mimetype] = TableSchemaFormatter()
        # enable it if it's been disabled:
        formatters[mimetype].enabled = True
    else:
        # unregister tableschema mime-type
        if mimetype in formatters:
            formatters[mimetype].enabled = False


default_pprint = lambda x, max_seq_items=None: pprint_thing(
    x, escape_chars=("\t", "\r", "\n"), quote_strings=True, max_seq_items=max_seq_items
)


def format_object_summary(
    obj,
    formatter,
    is_justify=True,
    name=None,
    indent_for_name=True,
    line_break_each_value=False,
):
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
    line_break_each_value : bool, default False
        If True, inserts a line break for each value of ``obj``.
        If False, only break lines when the a line of values gets wider
        than the display width.

        .. versionadded:: 0.25.0

    Returns
    -------
    summary string

    """
    from pandas.io.formats.console import get_console_size
    from pandas.io.formats.format import _get_adjustment

    display_width, _ = get_console_size()
    if display_width is None:
        display_width = get_option("display.width") or 80
    if name is None:
        name = obj.__class__.__name__

    if indent_for_name:
        name_len = len(name)
        space1 = "\n%s" % (" " * (name_len + 1))
        space2 = "\n%s" % (" " * (name_len + 2))
    else:
        space1 = "\n"
        space2 = "\n "  # space for the opening '['

    n = len(obj)
    if line_break_each_value:
        # If we want to vertically align on each value of obj, we need to
        # separate values by a line break and indent the values
        sep = ",\n " + " " * len(name)
    else:
        sep = ","
    max_seq_items = get_option("display.max_seq_items") or n

    # are we a truncated display
    is_truncated = n > max_seq_items

    # adj can optionally handle unicode eastern asian width
    adj = _get_adjustment()

    def _extend_line(s, line, value, display_width, next_line_prefix):

        if adj.len(line.rstrip()) + adj.len(value.rstrip()) >= display_width:
            s += line.rstrip()
            line = next_line_prefix
        line += value
        return s, line

    def best_len(values):
        if values:
            return max(adj.len(x) for x in values)
        else:
            return 0

    close = ", "

    if n == 0:
        summary = "[]{}".format(close)
    elif n == 1 and not line_break_each_value:
        first = formatter(obj[0])
        summary = "[{}]{}".format(first, close)
    elif n == 2 and not line_break_each_value:
        first = formatter(obj[0])
        last = formatter(obj[-1])
        summary = "[{}, {}]{}".format(first, last, close)
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
            if line_break_each_value:
                # Justify each string in the values of head and tail, so the
                # strings will right align when head and tail are stacked
                # vertically.
                head, tail = _justify(head, tail)
            elif is_truncated or not (
                len(", ".join(head)) < display_width
                and len(", ".join(tail)) < display_width
            ):
                # Each string in head and tail should align with each other
                max_length = max(best_len(head), best_len(tail))
                head = [x.rjust(max_length) for x in head]
                tail = [x.rjust(max_length) for x in tail]
            # If we are not truncated and we are only a single
            # line, then don't justify

        if line_break_each_value:
            # Now head and tail are of type List[Tuple[str]]. Below we
            # convert them into List[str], so there will be one string per
            # value. Also truncate items horizontally if wider than
            # max_space
            max_space = display_width - len(space2)
            value = tail[0]
            for max_items in reversed(range(1, len(value) + 1)):
                pprinted_seq = _pprint_seq(value, max_seq_items=max_items)
                if len(pprinted_seq) < max_space:
                    break
            head = [_pprint_seq(x, max_seq_items=max_items) for x in head]
            tail = [_pprint_seq(x, max_seq_items=max_items) for x in tail]

        summary = ""
        line = space2

        for max_items in range(len(head)):
            word = head[max_items] + sep + " "
            summary, line = _extend_line(summary, line, word, display_width, space2)

        if is_truncated:
            # remove trailing space of last line
            summary += line.rstrip() + space2 + "..."
            line = space2

        for max_items in range(len(tail) - 1):
            word = tail[max_items] + sep + " "
            summary, line = _extend_line(summary, line, word, display_width, space2)

        # last value: no sep added + 1 space of width used for trailing ','
        summary, line = _extend_line(summary, line, tail[-1], display_width - 2, space2)
        summary += line

        # right now close is either '' or ', '
        # Now we want to include the ']', but not the maybe space.
        close = "]" + close.rstrip(" ")
        summary += close

        if len(summary) > (display_width) or line_break_each_value:
            summary += space1
        else:  # one row
            summary += " "

        # remove initial space
        summary = "[" + summary[len(space2) :]

    return summary


def _justify(head, tail):
    """
    Justify items in head and tail, so they are right-aligned when stacked.

    Parameters
    ----------
    head : list-like of list-likes of strings
    tail : list-like of list-likes of strings

    Returns
    -------
    tuple of list of tuples of strings
        Same as head and tail, but items are right aligned when stacked
        vertically.

    Examples
    --------
    >>> _justify([['a', 'b']], [['abc', 'abcd']])
    ([('  a', '   b')], [('abc', 'abcd')])
    """
    combined = head + tail

    # For each position for the sequences in ``combined``,
    # find the length of the largest string.
    max_length = [0] * len(combined[0])
    for inner_seq in combined:
        length = [len(item) for item in inner_seq]
        max_length = [max(x, y) for x, y in zip(max_length, length)]

    # justify each item in each list-like in head and tail using max_length
    head = [
        tuple(x.rjust(max_len) for x, max_len in zip(seq, max_length)) for seq in head
    ]
    tail = [
        tuple(x.rjust(max_len) for x, max_len in zip(seq, max_length)) for seq in tail
    ]
    return head, tail


def format_object_attrs(obj, include_dtype=True):
    """
    Return a list of tuples of the (attr, formatted_value)
    for common attrs, including dtype, name, length

    Parameters
    ----------
    obj : object
        must be iterable
    include_dtype : bool
        If False, dtype won't be in the returned list

    Returns
    -------
    list

    """
    attrs = []
    if hasattr(obj, "dtype") and include_dtype:
        attrs.append(("dtype", "'{}'".format(obj.dtype)))
    if getattr(obj, "name", None) is not None:
        attrs.append(("name", default_pprint(obj.name)))
    elif getattr(obj, "names", None) is not None and any(obj.names):
        attrs.append(("names", default_pprint(obj.names)))
    max_seq_items = get_option("display.max_seq_items") or len(obj)
    if len(obj) > max_seq_items:
        attrs.append(("length", len(obj)))
    return attrs
