"""
Utilities for working with strings and text.

Inheritance diagram:

.. inheritance-diagram:: IPython.utils.text
   :parts: 3
"""

import os
import re
import string
import sys
import textwrap
import warnings
from string import Formatter
from pathlib import Path

from typing import (
    List,
    Dict,
    Tuple,
    Optional,
    cast,
    Any,
    Union,
    TypeVar,
)
from collections.abc import Sequence, Mapping, Callable, Iterator

if sys.version_info < (3, 12):
    from typing import Self
else:
    from typing import Self


class LSString(str):
    """String derivative with a special access attributes.

    These are normal strings, but with the special attributes:

        .l (or .list) : value as list (split on newlines).
        .n (or .nlstr): original value (the string itself).
        .s (or .spstr): value as whitespace-separated string.
        .p (or .paths): list of path objects (requires path.py package)

    Any values which require transformations are computed only once and
    cached.

    Such strings are very useful to efficiently interact with the shell, which
    typically only understands whitespace-separated options for commands."""

    __list: List[str]
    __spstr: str
    __paths: List[Path]

    def get_list(self) -> List[str]:
        try:
            return self.__list
        except AttributeError:
            self.__list = self.split('\n')
            return self.__list

    l = list = property(get_list)

    def get_spstr(self) -> str:
        try:
            return self.__spstr
        except AttributeError:
            self.__spstr = self.replace('\n',' ')
            return self.__spstr

    s = spstr = property(get_spstr)

    def get_nlstr(self) -> Self:
        return self

    n = nlstr = property(get_nlstr)

    def get_paths(self) -> List[Path]:
        try:
            return self.__paths
        except AttributeError:
            self.__paths = [Path(p) for p in self.split('\n') if os.path.exists(p)]
            return self.__paths

    p = paths = property(get_paths)

# FIXME: We need to reimplement type specific displayhook and then add this
# back as a custom printer. This should also be moved outside utils into the
# core.

# def print_lsstring(arg):
#     """ Prettier (non-repr-like) and more informative printer for LSString """
#     print("LSString (.p, .n, .l, .s available). Value:")
#     print(arg)
#
#
# print_lsstring = result_display.register(LSString)(print_lsstring)


class SList(list[Any]):
    """List derivative with a special access attributes.

    These are normal lists, but with the special attributes:

    * .l (or .list) : value as list (the list itself).
    * .n (or .nlstr): value as a string, joined on newlines.
    * .s (or .spstr): value as a string, joined on spaces.
    * .p (or .paths): list of path objects (requires path.py package)

    Any values which require transformations are computed only once and
    cached."""

    __spstr: str
    __nlstr: str
    __paths: List[Path]

    def get_list(self) -> Self:
        return self

    l = list = property(get_list)

    def get_spstr(self) -> str:
        try:
            return self.__spstr
        except AttributeError:
            self.__spstr = ' '.join(self)
            return self.__spstr

    s = spstr = property(get_spstr)

    def get_nlstr(self) -> str:
        try:
            return self.__nlstr
        except AttributeError:
            self.__nlstr = '\n'.join(self)
            return self.__nlstr

    n = nlstr = property(get_nlstr)

    def get_paths(self) -> List[Path]:
        try:
            return self.__paths
        except AttributeError:
            self.__paths = [Path(p) for p in self if os.path.exists(p)]
            return self.__paths

    p = paths = property(get_paths)

    def grep(
        self,
        pattern: Union[str, Callable[[Any], re.Match[str] | None]],
        prune: bool = False,
        field: Optional[int] = None,
    ) -> Self:
        """Return all strings matching 'pattern' (a regex or callable)

        This is case-insensitive. If prune is true, return all items
        NOT matching the pattern.

        If field is specified, the match must occur in the specified
        whitespace-separated field.

        Examples::

            a.grep( lambda x: x.startswith('C') )
            a.grep('Cha.*log', prune=1)
            a.grep('chm', field=-1)
        """

        def match_target(s: str) -> str:
            if field is None:
                return s
            parts = s.split()
            try:
                tgt = parts[field]
                return tgt
            except IndexError:
                return ""

        if isinstance(pattern, str):
            pred = lambda x : re.search(pattern, x, re.IGNORECASE)
        else:
            pred = pattern
        if not prune:
            return type(self)([el for el in self if pred(match_target(el))])  # type: ignore [no-untyped-call]
        else:
            return type(self)([el for el in self if not pred(match_target(el))])  # type: ignore [no-untyped-call]

    def fields(self, *fields: List[str]) -> List[List[str]]:
        """Collect whitespace-separated fields from string list

        Allows quick awk-like usage of string lists.

        Example data (in var a, created by 'a = !ls -l')::

            -rwxrwxrwx  1 ville None      18 Dec 14  2006 ChangeLog
            drwxrwxrwx+ 6 ville None       0 Oct 24 18:05 IPython

        * ``a.fields(0)`` is ``['-rwxrwxrwx', 'drwxrwxrwx+']``
        * ``a.fields(1,0)`` is ``['1 -rwxrwxrwx', '6 drwxrwxrwx+']``
          (note the joining by space).
        * ``a.fields(-1)`` is ``['ChangeLog', 'IPython']``

        IndexErrors are ignored.

        Without args, fields() just split()'s the strings.
        """
        if len(fields) == 0:
            return [el.split() for el in self]

        res = SList()
        for el in [f.split() for f in self]:
            lineparts = []

            for fd in fields:
                try:
                    lineparts.append(el[fd])
                except IndexError:
                    pass
            if lineparts:
                res.append(" ".join(lineparts))

        return res

    def sort(  # type:ignore[override]
        self,
        field: Optional[List[str]] = None,
        nums: bool = False,
    ) -> Self:
        """sort by specified fields (see fields())

        Example::

            a.sort(1, nums = True)

        Sorts a by second field, in numerical order (so that 21 > 3)

        """

        #decorate, sort, undecorate
        if field is not None:
            dsu = [[SList([line]).fields(field),  line] for line in self]
        else:
            dsu = [[line,  line] for line in self]
        if nums:
            for i in range(len(dsu)):
                numstr = "".join([ch for ch in dsu[i][0] if ch.isdigit()])
                try:
                    n = int(numstr)
                except ValueError:
                    n = 0
                dsu[i][0] = n


        dsu.sort()
        return type(self)([t[1] for t in dsu])


def indent(instr: str, nspaces: int = 4, ntabs: int = 0, flatten: bool = False) -> str:
    """Indent a string a given number of spaces or tabstops.

    indent(str, nspaces=4, ntabs=0) -> indent str by ntabs+nspaces.

    Parameters
    ----------
    instr : basestring
        The string to be indented.
    nspaces : int (default: 4)
        The number of spaces to be indented.
    ntabs : int (default: 0)
        The number of tabs to be indented.
    flatten : bool (default: False)
        Whether to scrub existing indentation.  If True, all lines will be
        aligned to the same indentation.  If False, existing indentation will
        be strictly increased.

    Returns
    -------
    str : string indented by ntabs and nspaces.

    """
    ind = "\t" * ntabs + " " * nspaces
    if flatten:
        pat = re.compile(r'^\s*', re.MULTILINE)
    else:
        pat = re.compile(r'^', re.MULTILINE)
    outstr = re.sub(pat, ind, instr)
    if outstr.endswith(os.linesep+ind):
        return outstr[:-len(ind)]
    else:
        return outstr


def list_strings(arg: Union[str, List[str]]) -> List[str]:
    """Always return a list of strings, given a string or list of strings
    as input.

    Examples
    --------
    ::

        In [7]: list_strings('A single string')
        Out[7]: ['A single string']

        In [8]: list_strings(['A single string in a list'])
        Out[8]: ['A single string in a list']

        In [9]: list_strings(['A','list','of','strings'])
        Out[9]: ['A', 'list', 'of', 'strings']
    """

    if isinstance(arg, str):
        return [arg]
    else:
        return arg


def marquee(txt: str = "", width: int = 78, mark: str = "*") -> str:
    """Return the input string centered in a 'marquee'.

    Examples
    --------
    ::

        In [16]: marquee('A test',40)
        Out[16]: '**************** A test ****************'

        In [17]: marquee('A test',40,'-')
        Out[17]: '---------------- A test ----------------'

        In [18]: marquee('A test',40,' ')
        Out[18]: '                 A test                 '

    """
    if not txt:
        return (mark*width)[:width]
    nmark = (width-len(txt)-2)//len(mark)//2
    if nmark < 0: nmark =0
    marks = mark*nmark
    return '%s %s %s' % (marks,txt,marks)


def format_screen(strng: str) -> str:
    """Format a string for screen printing.

    This removes some latex-type format codes."""
    # Paragraph continue
    par_re = re.compile(r'\\$',re.MULTILINE)
    strng = par_re.sub('',strng)
    return strng


def dedent(text: str) -> str:
    """Equivalent of textwrap.dedent that ignores unindented first line.

    This means it will still dedent strings like:
    '''foo
    is a bar
    '''

    For use in wrap_paragraphs.
    """

    if text.startswith('\n'):
        # text starts with blank line, don't ignore the first line
        return textwrap.dedent(text)

    # split first line
    splits = text.split('\n',1)
    if len(splits) == 1:
        # only one line
        return textwrap.dedent(text)

    first, rest = splits
    # dedent everything but the first line
    rest = textwrap.dedent(rest)
    return '\n'.join([first, rest])


def strip_email_quotes(text: str) -> str:
    """Strip leading email quotation characters ('>').

    Removes any combination of leading '>' interspersed with whitespace that
    appears *identically* in all lines of the input text.

    Parameters
    ----------
    text : str

    Examples
    --------

    Simple uses::

        In [2]: strip_email_quotes('> > text')
        Out[2]: 'text'

        In [3]: strip_email_quotes('> > text\\n> > more')
        Out[3]: 'text\\nmore'

    Note how only the common prefix that appears in all lines is stripped::

        In [4]: strip_email_quotes('> > text\\n> > more\\n> more...')
        Out[4]: '> text\\n> more\\nmore...'

    So if any line has no quote marks ('>'), then none are stripped from any
    of them ::

        In [5]: strip_email_quotes('> > text\\n> > more\\nlast different')
        Out[5]: '> > text\\n> > more\\nlast different'
    """
    lines = text.splitlines()
    strip_len = 0

    for characters in zip(*lines):
        # Check if all characters in this position are the same
        if len(set(characters)) > 1:
            break
        prefix_char = characters[0]

        if prefix_char in string.whitespace or prefix_char == ">":
            strip_len += 1
        else:
            break

    text = "\n".join([ln[strip_len:] for ln in lines])
    return text


class EvalFormatter(Formatter):
    """A String Formatter that allows evaluation of simple expressions.

    Note that this version interprets a `:`  as specifying a format string (as per
    standard string formatting), so if slicing is required, you must explicitly
    create a slice.

    Note that on Python 3.14+ this version interprets `[]` as indexing operator
    so you need to use generators instead of list comprehensions, for example:
    `list(i for i in range(10))`.

    This is to be used in templating cases, such as the parallel batch
    script templates, where simple arithmetic on arguments is useful.

    Examples
    --------
    ::

        In [1]: f = EvalFormatter()
        In [2]: f.format('{n//4}', n=8)
        Out[2]: '2'

        In [3]: f.format("{greeting[slice(2,4)]}", greeting="Hello")
        Out[3]: 'll'
    """

    def get_field(self, name: str, args: Any, kwargs: Any) -> Tuple[Any, str]:
        v = eval(name, kwargs, kwargs)
        return v, name

#XXX: As of Python 3.4, the format string parsing no longer splits on a colon
# inside [], so EvalFormatter can handle slicing. Once we only support 3.4 and
# above, it should be possible to remove FullEvalFormatter.

class FullEvalFormatter(Formatter):
    """A String Formatter that allows evaluation of simple expressions.
    
    Any time a format key is not found in the kwargs,
    it will be tried as an expression in the kwargs namespace.
    
    Note that this version allows slicing using [1:2], so you cannot specify
    a format string. Use :class:`EvalFormatter` to permit format strings.
    
    Examples
    --------
    ::

        In [1]: f = FullEvalFormatter()
        In [2]: f.format('{n//4}', n=8)
        Out[2]: '2'

        In [3]: f.format('{list(range(5))[2:4]}')
        Out[3]: '[2, 3]'

        In [4]: f.format('{3*2}')
        Out[4]: '6'
    """
    # copied from Formatter._vformat with minor changes to allow eval
    # and replace the format_spec code with slicing
    def vformat(
        self, format_string: str, args: Sequence[Any], kwargs: Mapping[str, Any]
    ) -> str:
        result = []
        conversion: Optional[str]
        for literal_text, field_name, format_spec, conversion in self.parse(
            format_string
        ):
            # output the literal text
            if literal_text:
                result.append(literal_text)

            # if there's a field, output it
            if field_name is not None:
                # this is some markup, find the object and do
                # the formatting

                if format_spec:
                    # override format spec, to allow slicing:
                    field_name = ':'.join([field_name, format_spec])

                # eval the contents of the field for the object
                # to be formatted
                obj = eval(field_name, dict(kwargs))

                # do any conversion on the resulting object
                # type issue in typeshed, fined in https://github.com/python/typeshed/pull/11377
                obj = self.convert_field(obj, conversion)

                # format the object and append to the result
                result.append(self.format_field(obj, ''))

        return ''.join(result)


class DollarFormatter(FullEvalFormatter):
    """Formatter allowing Itpl style $foo replacement, for names and attribute
    access only. Standard {foo} replacement also works, and allows full
    evaluation of its arguments.

    Examples
    --------
    ::

        In [1]: f = DollarFormatter()
        In [2]: f.format('{n//4}', n=8)
        Out[2]: '2'

        In [3]: f.format('23 * 76 is $result', result=23*76)
        Out[3]: '23 * 76 is 1748'

        In [4]: f.format('$a or {b}', a=1, b=2)
        Out[4]: '1 or 2'
    """

    _dollar_pattern_ignore_single_quote = re.compile(
        r"(.*?)\$(\$?[\w\.]+)(?=([^']*'[^']*')*[^']*$)"
    )

    def parse(self, fmt_string: str) -> Iterator[Tuple[Any, Any, Any, Any]]:
        for literal_txt, field_name, format_spec, conversion in Formatter.parse(
            self, fmt_string
        ):
            # Find $foo patterns in the literal text.
            continue_from = 0
            txt = ""
            for m in self._dollar_pattern_ignore_single_quote.finditer(literal_txt):
                new_txt, new_field = m.group(1,2)
                # $$foo --> $foo
                if new_field.startswith("$"):
                    txt += new_txt + new_field
                else:
                    yield (txt + new_txt, new_field, "", None)
                    txt = ""
                continue_from = m.end()
            
            # Re-yield the {foo} style pattern
            yield (txt + literal_txt[continue_from:], field_name, format_spec, conversion)

    def __repr__(self) -> str:
        return "<DollarFormatter>"

#-----------------------------------------------------------------------------
# Utils to columnize a list of string
#-----------------------------------------------------------------------------


def _col_chunks(
    l: List[int], max_rows: int, row_first: bool = False
) -> Iterator[List[int]]:
    """Yield successive max_rows-sized column chunks from l."""
    if row_first:
        ncols = (len(l) // max_rows) + (len(l) % max_rows > 0)
        for i in range(ncols):
            yield [l[j] for j in range(i, len(l), ncols)]
    else:
        for i in range(0, len(l), max_rows):
            yield l[i:(i + max_rows)]


def _find_optimal(
    rlist: List[int], row_first: bool, separator_size: int, displaywidth: int
) -> Dict[str, Any]:
    """Calculate optimal info to columnize a list of string"""
    for max_rows in range(1, len(rlist) + 1):
        col_widths = list(map(max, _col_chunks(rlist, max_rows, row_first)))
        sumlength = sum(col_widths)
        ncols = len(col_widths)
        if sumlength + separator_size * (ncols - 1) <= displaywidth:
            break
    return {'num_columns': ncols,
            'optimal_separator_width': (displaywidth - sumlength) // (ncols - 1) if (ncols - 1) else 0,
            'max_rows': max_rows,
            'column_widths': col_widths
            }


T = TypeVar("T")


def _get_or_default(mylist: List[T], i: int, default: T) -> T:
    """return list item number, or default if don't exist"""
    if i >= len(mylist):
        return default
    else :
        return mylist[i]


def get_text_list(
    list_: List[str], last_sep: str = " and ", sep: str = ", ", wrap_item_with: str = ""
) -> str:
    """
    Return a string with a natural enumeration of items

    >>> get_text_list(['a', 'b', 'c', 'd'])
    'a, b, c and d'
    >>> get_text_list(['a', 'b', 'c'], ' or ')
    'a, b or c'
    >>> get_text_list(['a', 'b', 'c'], ', ')
    'a, b, c'
    >>> get_text_list(['a', 'b'], ' or ')
    'a or b'
    >>> get_text_list(['a'])
    'a'
    >>> get_text_list([])
    ''
    >>> get_text_list(['a', 'b'], wrap_item_with="`")
    '`a` and `b`'
    >>> get_text_list(['a', 'b', 'c', 'd'], " = ", sep=" + ")
    'a + b + c = d'
    """
    if len(list_) == 0:
        return ''
    if wrap_item_with:
        list_ = ['%s%s%s' % (wrap_item_with, item, wrap_item_with) for
                 item in list_]
    if len(list_) == 1:
        return list_[0]
    return '%s%s%s' % (
        sep.join(i for i in list_[:-1]),
        last_sep, list_[-1])
