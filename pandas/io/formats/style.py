"""
Module for applying data-dependent formatting to DataFrames and Series.
"""
from functools import partial
from itertools import product
from contextlib import contextmanager
from uuid import uuid1
import copy
from collections import defaultdict, MutableMapping

try:
    from jinja2 import (
        PackageLoader, Environment, ChoiceLoader, FileSystemLoader
    )
except ImportError:
    msg = "pandas.Styler requires jinja2. "\
          "Please install with `conda install Jinja2`\n"\
          "or `pip install Jinja2`"
    raise ImportError(msg)

from pandas.core.dtypes.common import is_float, is_string_like

import numpy as np
import pandas as pd
from pandas.api.types import is_list_like
from pandas.compat import range
from pandas.core.config import get_option
from pandas.core.generic import _shared_docs
from pandas.core.common import _any_not_none, sentinel_factory
from pandas.core.indexing import _maybe_numeric_slice, _non_reducing_slice
from pandas.util._decorators import Appender
try:
    import matplotlib.pyplot as plt
    from matplotlib import colors
    has_mpl = True
except ImportError:
    has_mpl = False
    no_mpl_message = "{0} requires matplotlib."


@contextmanager
def _mpl(func):
    if has_mpl:
        yield plt, colors
    else:
        raise ImportError(no_mpl_message.format(func.__name__))


class Styler(object):
    """
    A decorator class for a `DataFrame` object to provide an
    HTML renderer configurable by data-dependent parameters.

    Parameters
    ----------
    data : Series or DataFrame
        the object to which the Styler is attached
    precision : int
        precision to round floats to, defaults to pd.options.display.precision
    table_styles : list-like, default None
        list of ``{selector: (attr, value)}`` dicts; see docs for the
        `set_table_styles` method
    uuid : str, default None
        a unique identifier to avoid CSS collisons; generated automatically
    caption : str, default None
        caption to attach to the table

    Attributes
    ----------
    env : Jinja2 Environment
    template : Jinja2 Template
    loader : Jinja2 Loader

    Notes
    -----

    `Styler` renders HTML on a `DataFrame` using its own `Styler.render`
    method; in particular, the Styler class will not affect, and is
    independent from, its `data` object's `DataFrame.to_html` rendering
    method.

    `Styler.render` returns an HTML string, based on attributes which are set
    using methods of `Styler`, and are stored as instance attributes in a
    `Styler` instance.

    Some of these attributes are _data-dependent_ functions: they affect the
    rendered HTML corresponding to an entry, row, or column of `data` based on
    its contents.

    Others are _data-independent_: they affect the rendering of the entire
    table, or only affect entries of the table based on their position.

    Some, like the table values set by `set_table_values`, are something
    in-between: HTML table entries corresponding to empty `DataFrame` entries
    are endowed with CSS selectors which indicate this fact, and which can be
    used to set global styles.

    `Styler.render` first constructs a table of strings corresponding
    to each entry of `data`, possibly with the help of the data-dependent
    function passed into `Styler.format`.  It then constructs dict
    with general formatting determined by the data-independent properties set,
    for instance, by `Styler.set_caption` and `Styler.set_uuid`; and CSS
    attributes applied to cells, rows, and columns by data-dependent functions
    passed into `Styler.apply`, `Styler.applymap`, and `Styler.where`.
    Currently, you cannot apply data-dependent styling to column or row
    headers.  It then passes this dict to Jinja2's render function with a
    fixed template file.

    `Styler also has a `_repr_html_` function which calls
    `Styler.render`, so when a Styler object is the value of the last line
    in a Jupyter cell, the Styler is automatically rendered to the ``Out``
    field of the cell.

    Style data can be transferred between `Styler` instances with the
    `Styler.use` and `Styler.export` methods.

    See Also
    --------
    pandas.DataFrame.style
    """
    loader = PackageLoader("pandas", "io/formats/templates")
    env = Environment(
        loader=loader,
        trim_blocks=True,
    )
    template = env.get_template("html.tpl")

    def __init__(self, data, precision=None, table_styles=None, uuid=None,
                 caption=None, table_attributes=None):
        self.ctx = defaultdict(list)
        self._todo = []

        if not isinstance(data, (pd.Series, pd.DataFrame)):
            raise TypeError("``data`` must be a Series or DataFrame")
        if data.ndim == 1:
            data = data.to_frame()
        if not data.index.is_unique or not data.columns.is_unique:
            raise ValueError("style is not supported for non-unique indices.")

        self.data = data
        self.index = data.index
        self.columns = data.columns

        self.uuid = uuid
        self.table_styles = table_styles
        self.caption = caption
        if precision is None:
            precision = get_option('display.precision')
        self.precision = precision
        self.table_attributes = table_attributes
        self.hidden_index = False
        self.hidden_columns = []

        # display_funcs maps (row, col) -> formatting function

        def default_display_func(x):
            if is_float(x):
                return '{:>.{precision}g}'.format(x, precision=self.precision)
            else:
                return x

        self._display_funcs = defaultdict(lambda: default_display_func)

    def _repr_html_(self):
        """Hooks into Jupyter notebook rich display system."""
        return self.render()

    @Appender(_shared_docs['to_excel'] % dict(
        axes='index, columns', klass='Styler',
        axes_single_arg="{0 or 'index', 1 or 'columns'}",
        optional_by="""
            by : str or list of str
                Name or list of names which refer to the axis items.""",
        versionadded_to_excel='\n    .. versionadded:: 0.20'))
    def to_excel(self, excel_writer, sheet_name='Sheet1', na_rep='',
                 float_format=None, columns=None, header=True, index=True,
                 index_label=None, startrow=0, startcol=0, engine=None,
                 merge_cells=True, encoding=None, inf_rep='inf', verbose=True,
                 freeze_panes=None):
                """
                Returns an Excel representation of the DataFrame.
                Currently, `to_excel` creates a
                `pandas.io.formats.excel.ExcelFormatter` object and passes
                all parameters to (depending on the parameter) its
                constructor, or its `write` method; all keyword arguments
                are described in `ExcelFormatter`'s documentation.

                See Also
                --------
                pandas.io.formats.excel.ExcelFormatter
                """
                from pandas.io.formats.excel import ExcelFormatter
                formatter = ExcelFormatter(self, na_rep=na_rep, cols=columns,
                                           header=header,
                                           float_format=float_format,
                                           index=index,
                                           index_label=index_label,
                                           merge_cells=merge_cells,
                                           inf_rep=inf_rep)
                formatter.write(excel_writer, sheet_name=sheet_name,
                                startrow=startrow,
                                startcol=startcol, freeze_panes=freeze_panes,
                                engine=engine)

    def _translate(self):
        """
        Builds a dictionary representation out of the data processed with
        `_compute()` to be passed to `self.template.render()`.
        """
        table_styles = self.table_styles or []
        caption = self.caption
        ctx = self.ctx
        precision = self.precision
        hidden_index = self.hidden_index
        hidden_columns = self.hidden_columns
        uuid = self.uuid or str(uuid1()).replace("-", "_")
        ROW_HEADING_CLASS = "row_heading"
        COL_HEADING_CLASS = "col_heading"
        INDEX_NAME_CLASS = "index_name"

        DATA_CLASS = "data"
        BLANK_CLASS = "blank"
        BLANK_VALUE = ""

        def format_attr(pair):
            return "{key}={value}".format(**pair)

        # for sparsifying a MultiIndex
        idx_lengths = _get_level_lengths(self.index)
        col_lengths = _get_level_lengths(self.columns, hidden_columns)

        cell_context = dict()

        n_rlvls = self.data.index.nlevels
        n_clvls = self.data.columns.nlevels
        rlabels = self.data.index.tolist()
        clabels = self.data.columns.tolist()

        if n_rlvls == 1:
            rlabels = [[x] for x in rlabels]
        if n_clvls == 1:
            clabels = [[x] for x in clabels]
        clabels = list(zip(*clabels))

        cellstyle = []
        head = []

        for r in range(n_clvls):
            # Blank for Index columns...
            row_es = [{"type": "th",
                       "value": BLANK_VALUE,
                       "display_value": BLANK_VALUE,
                       "is_visible": not hidden_index,
                       "class": " ".join([BLANK_CLASS])}] * (n_rlvls - 1)

            # ... except maybe the last for columns.names
            name = self.data.columns.names[r]
            cs = [BLANK_CLASS if name is None else INDEX_NAME_CLASS,
                  "level{lvl}".format(lvl=r)]
            name = BLANK_VALUE if name is None else name
            row_es.append({"type": "th",
                           "value": name,
                           "display_value": name,
                           "class": " ".join(cs),
                           "is_visible": not hidden_index})

            if clabels:
                for c, value in enumerate(clabels[r]):
                    cs = [COL_HEADING_CLASS, "level{lvl}".format(lvl=r),
                          "col{col}".format(col=c)]
                    cs.extend(cell_context.get(
                        "col_headings", {}).get(r, {}).get(c, []))
                    es = {
                        "type": "th",
                        "value": value,
                        "display_value": value,
                        "class": " ".join(cs),
                        "is_visible": _is_visible(c, r, col_lengths),
                    }
                    colspan = col_lengths.get((r, c), 0)
                    if colspan > 1:
                        es["attributes"] = [
                            format_attr({"key": "colspan", "value": colspan})
                        ]
                    row_es.append(es)
                head.append(row_es)

        if (self.data.index.names and _any_not_none(*self.data.index.names) and
                not hidden_index):
            index_header_row = []

            for c, name in enumerate(self.data.index.names):
                cs = [INDEX_NAME_CLASS,
                      "level{lvl}".format(lvl=c)]
                name = '' if name is None else name
                index_header_row.append({"type": "th", "value": name,
                                         "class": " ".join(cs)})

            index_header_row.extend(
                [{"type": "th",
                  "value": BLANK_VALUE,
                  "class": " ".join([BLANK_CLASS])
                  }] * (len(clabels[0]) - len(hidden_columns)))

            head.append(index_header_row)

        body = []
        for r, idx in enumerate(self.data.index):
            row_es = []
            for c, value in enumerate(rlabels[r]):
                rid = [ROW_HEADING_CLASS, "level{lvl}".format(lvl=c),
                       "row{row}".format(row=r)]
                es = {
                    "type": "th",
                    "is_visible": (_is_visible(r, c, idx_lengths) and
                                   not hidden_index),
                    "value": value,
                    "display_value": value,
                    "id": "_".join(rid[1:]),
                    "class": " ".join(rid)
                }
                rowspan = idx_lengths.get((c, r), 0)
                if rowspan > 1:
                    es["attributes"] = [
                        format_attr({"key": "rowspan", "value": rowspan})
                    ]
                row_es.append(es)

            for c, col in enumerate(self.data.columns):
                cs = [DATA_CLASS, "row{row}".format(row=r),
                      "col{col}".format(col=c)]
                cs.extend(cell_context.get("data", {}).get(r, {}).get(c, []))
                formatter = self._display_funcs[(r, c)]
                value = self.data.iloc[r, c]
                row_es.append({
                    "type": "td",
                    "value": value,
                    "class": " ".join(cs),
                    "id": "_".join(cs[1:]),
                    "display_value": formatter(value),
                    "is_visible": (c not in hidden_columns)
                })
                props = []
                for x in ctx[r, c]:
                    # have to handle empty styles like ['']
                    if x.count(":"):
                        props.append(x.split(":"))
                    else:
                        props.append(['', ''])
                cellstyle.append({'props': props,
                                  'selector': "row{row}_col{col}"
                                  .format(row=r, col=c)})
            body.append(row_es)

        return dict(head=head, cellstyle=cellstyle, body=body, uuid=uuid,
                    precision=precision, table_styles=table_styles,
                    caption=caption, table_attributes=self.table_attributes)

    def format(self, formatter, subset=None):
        """
        Specify the text rendered for each table entry based on a formatter.

        .. versionadded:: 0.18.0

        Parameters
        ----------
        formatter: str, callable, or dict
            A `callable` must return a `str`, and a `dict`'s values
            must be either `str`'s or `callable`'s which return a `str`.
        subset: IndexSlice
            An argument to `DataFrame.loc` that restricts which elements
            `formatter` is applied to.

        Returns
        -------
        self : Styler

        Notes
        -----
        A `DataFrame` entry with contents `x` will be formatted to the
        string representation:
            * `formatter.format(x)` if `formatter` is a `str`.
            * `formatter(x)` if `formatter` is a `callable`.
            * `formatter[<col_name>].format(x)` if `formatter` is a
                `dict`, `x` is in the column defined by `<col_name>`,
                and `formatter[<col_name>]` is a `str`.
            * `formatter[<col_name>].(x)` if `formatter` is a `dict`,
                `x` is in the column defined by `<col_name>`, and
                `formatter[<col_name>]` is a `callable` which returns
                a `str`.

        If no formatter is specified for an entry, the default display value
        for numeric values is the "general" (``g``) format with
        `pd.options.display.precision` precision.

        Unlike the methods which apply CSS to table entries, each `Styler`
        instance can have at most one `formatter` method specified at any
        point; calling `format` a second time with a different `formatter`
        parameter replaces the old `formatter` method with the newly
        specified parameter.

        Note that because `Styler`'s current implementation does not sanitize
        HTML tags, you can insert arbitrary HTML code/javascript into your
        format strings, and they will render/execute.  See examples below.

        Examples
        --------

        >>> df = pd.DataFrame(np.random.randn(4, 2), columns=['a', 'b'])
        >>> df.style.format("{:.2%}")
        >>> df['c'] = ['a', 'b', 'c', 'd']
        >>> df.style.format({'C': str.upper})
        >>> df.style.format("<b>{}</b>") # Make everything bold
        >>> df.style.format(\"\"\"
            <script>
            document.removeChild(document.documentElement))
            </script>
            \"\"\") # You probably don't want to do this
        """
        if subset is None:
            row_locs = range(len(self.data))
            col_locs = range(len(self.data.columns))
        else:
            subset = _non_reducing_slice(subset)
            if len(subset) == 1:
                subset = subset, self.data.columns

            sub_df = self.data.loc[subset]
            row_locs = self.data.index.get_indexer_for(sub_df.index)
            col_locs = self.data.columns.get_indexer_for(sub_df.columns)

        if isinstance(formatter, MutableMapping):
            for col, col_formatter in formatter.items():
                # formatter must be callable, so '{}' are converted to lambdas
                col_formatter = _maybe_wrap_formatter(col_formatter)
                col_num = self.data.columns.get_indexer_for([col])[0]

                for row_num in row_locs:
                    self._display_funcs[(row_num, col_num)] = col_formatter
        else:
            # single scalar to format all cells with
            locs = product(*(row_locs, col_locs))
            for i, j in locs:
                formatter = _maybe_wrap_formatter(formatter)
                self._display_funcs[(i, j)] = formatter
        return self

    def render(self, **kwargs):
        """Render the DataFrame `data` as HTML, as described in the Notes
        for the class.

        Parameters
        ----------
        `**kwargs` :
            Any keyword arguments are passed through to the default renderer
            provided by `jinja2.Environment`. This is useful when you
            need to provide additional variables for a custom
            template.

            .. versionadded:: 0.20

        Returns
        -------
        rendered : str
            the rendered HTML

        Notes
        -----
        `Styler._repr_html_` calls render to do HTML rendering.  This means
        that `Styler.render` is automatically invoked to display the output
        in a Jupyter Notebook cell if the value of the last item is a
        `Styler`.

        In order to render HTML into the Out field of a Notebook cell when
        a `Styler` is not the last value, invoke
        `IPython.display.HTML(self.render())`.

        `Styler` renders HTML by creating a dictionary with the entries:
        * head
        * cellstyle
        * body
        * uuid
        * precision
        * table_styles
        * caption
        * table_attributes
        based on the attributes already assigned to the `Styler` instance,
        and passes this dictionary to the `render` method determined by
        `jinja2.Environment`.

        This dictionary is the penultimate result of the process embodied by
        `Styler` (the final product being, of course, the rendered HTML).
        Dictionary values can be over-ridden by keyword arguments passed to
        `Styler.render`; this accesses the underlying rendering engine
        directly, however, so is unstable: the behavior keywords could change
        if the backend changes. Therefore, unless you have both a really good
        reason, and intimate knowledge of the rendering backend, you probably
        should not pass keyword arguments to `render`.
        """
        self._compute()
        # TODO: namespace all the pandas keys
        d = self._translate()
        # filter out empty styles, every cell will have a class
        # but the list of props may just be [['', '']].
        # so we have the neested anys below
        trimmed = [x for x in d['cellstyle']
                   if any(any(y) for y in x['props'])]
        d['cellstyle'] = trimmed
        d.update(kwargs)
        return self.template.render(**d)

    def _update_ctx(self, attrs):
        """
        Update the state of the Styler. Collects a dict of the form
        {index_label: ['<property>: <value>']}

        Parameters
        ----------
        attrs: Series or DataFrame
            Should contain strings of the form '<property>: <value>; <prop2>:
            <val2>; ...'.  Neither whitespace nor the final semicolon matters.
        """
        for row_label, v in attrs.iterrows():
            for col_label, col in v.iteritems():
                i = self.index.get_indexer([row_label])[0]
                j = self.columns.get_indexer([col_label])[0]
                for pair in col.rstrip(";").split(";"):
                    self.ctx[(i, j)].append(pair)

    def _copy(self, deepcopy=False):
        styler = Styler(self.data, precision=self.precision,
                        caption=self.caption, uuid=self.uuid,
                        table_styles=self.table_styles)
        if deepcopy:
            styler.ctx = copy.deepcopy(self.ctx)
            styler._todo = copy.deepcopy(self._todo)
        else:
            styler.ctx = self.ctx
            styler._todo = self._todo
        return styler

    def __copy__(self):
        """
        Shallow copy by default.
        """
        return self._copy(deepcopy=False)

    def __deepcopy__(self, memo):
        return self._copy(deepcopy=True)

    def clear(self):
        """
        "Reset" the styler, removing any previously applied styles.
        Returns None.
        """
        self.ctx.clear()
        self._todo = []

    def _compute(self):
        """
        Execute the style methods enumerated in `self._todo`,
        which usually got there because they were appended by
        `Styler.apply`, `Styler.applymap`, and `Styler.where`.

        Returns none
        """
        r = self
        for func, args, kwargs in self._todo:
            r = func(self)(*args, **kwargs)
        return r

    def _apply(self, func, axis=0, subset=None, **kwargs):
        subset = slice(None) if subset is None else subset
        subset = _non_reducing_slice(subset)
        data = self.data.loc[subset]
        if axis is not None:
            result = data.apply(func, axis=axis, **kwargs)
        else:
            result = func(data, **kwargs)
            if not isinstance(result, pd.DataFrame):
                raise TypeError(
                    "Function {func!r} must return a DataFrame when "
                    "passed to `Styler.apply` with axis=None"
                    .format(func=func))
            if not (result.index.equals(data.index) and
                    result.columns.equals(data.columns)):
                msg = ('Result of {func!r} must have identical index and '
                       'columns as the input'.format(func=func))
                raise ValueError(msg)

        result_shape = result.shape
        expected_shape = self.data.loc[subset].shape
        if result_shape != expected_shape:
            msg = ("Function {func!r} returned the wrong shape.\n"
                   "Result has shape: {res}\n"
                   "Expected shape:   {expect}".format(func=func,
                                                       res=result.shape,
                                                       expect=expected_shape))
            raise ValueError(msg)
        self._update_ctx(result)
        return self

    def apply(self, func, axis=0, subset=None, **kwargs):
        """
        Assure that when `Styler.render` is called, CSS attributes
        determined by the data-dependent function `func` are added
        to cells, rows or entries.

        Parameters
        ----------
        func : function
            `func` should take a Series or DataFrame (depending
            on `axis`), and return an object with the same shape.
            Must return a DataFrame with identical index and
            column labels to those of `data` when `axis=None`
        axis : int, str or None
            determine whether to add attributes to each row
            (`axis=0` or `'index'`), to each column (`axis=1` or
            `'columns'`) or to the entire DataFrame at once with
            `axis=None`
        subset : IndexSlice
            a valid indexer to limit `data` to *before* applying `func`.
            Consider using a pandas.IndexSlice
        kwargs : dict
            pass along to `func`

        Returns
        -------
        self : Styler

        Notes
        -----
        The output shape of `func` should match the input, i.e., if
        `x` is the input row, column, or table (depending on `axis`),
        then `func(x.shape) == x.shape` should be true.

        In this way, `Styler.apply` handles its parameters similarly to
        `DataFrame.apply`, except that `axis=None` applies `func` to
        the entire DataFrame at once, rather than column-wise or row-wise.

        Examples
        --------
        >>> def highlight_max(x):
        ...     return ['background-color: yellow' if v == x.max() else ''
                        for v in x]
        ...
        >>> df = pd.DataFrame(np.random.randn(5, 2))
        >>> df.style.apply(highlight_max)
        """
        self._todo.append((lambda instance: getattr(instance, '_apply'),
                           (func, axis, subset), kwargs))
        return self

    def _applymap(self, func, subset=None, **kwargs):
        func = partial(func, **kwargs)
        if subset is None:
            subset = pd.IndexSlice[:]
        subset = _non_reducing_slice(subset)
        result = self.data.loc[subset].applymap(func)
        self._update_ctx(result)
        return self

    def applymap(self, func, subset=None, **kwargs):
        """
        Assure that when `Styler.render` is called, CSS attributes
        determined by the data-dependent function parameter `func`
        are added to individual table entries.

        Parameters
        ----------
        func : function
            `func` should take a scalar and return a scalar
        subset : IndexSlice
            determines to which cells `func` will apply
        kwargs : dict
            parameters to pass to `func`

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.where

        """
        self._todo.append((lambda instance: getattr(instance, '_applymap'),
                           (func, subset), kwargs))
        return self

    def where(self, cond, value, other=None, subset=None, **kwargs):
        """
        Assure that when `Styler.render` is called, the CSS attribute
        `value` is added to an entry when `cond` of that entry is True,
        and that the CSS attribute `other` is added when `cond` is False.
        are added to individual table entries

        .. versionadded:: 0.21.0

        Parameters
        ----------
        cond : callable
            `cond` should take a scalar and return a boolean
        value : str
            applied when `cond` returns true
        other : str
            applied when `cond` returns false
        subset : IndexSlice
            a valid indexer to limit `data` to *before* applying `cond`.
            Consider using a pandas.IndexSlice
        kwargs : dict
            pass along to `cond`

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.applymap

        """

        if other is None:
            other = ''

        return self.applymap(lambda val: value if cond(val) else other,
                             subset=subset, **kwargs)

    def set_precision(self, precision):
        """
        Set the precision used to render.

        Parameters
        ----------
        precision: int

        Returns
        -------
        self : Styler
        """
        self.precision = precision
        return self

    def set_table_attributes(self, attributes):
        """
        Set the table attributes. These are the items
        that show up in the opening `<table>` tag in addition
        to to automatic (by default) id. Data-independent attribute.

        Parameters
        ----------
        attributes : string

        Returns
        -------
        self : Styler

        Examples
        --------
        >>> df = pd.DataFrame(np.random.randn(10, 4))
        >>> df.style.set_table_attributes('class="pure-table"')
        # ... <table class="pure-table"> ...
        """
        self.table_attributes = attributes
        return self

    def export(self):
        """
        Export the styles to applied to the current Styler.
        Can be applied to a second style with `Styler.use`.

        Returns
        -------
        styles: list

        See Also
        --------
        Styler.use
        """
        return self._todo

    def use(self, styles):
        """
        Set the styles on the current Styler, possibly using styles
        from `Styler.export`.

        Parameters
        ----------
        styles: list
            list of style functions

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.export
        """
        self._todo.extend(styles)
        return self

    def set_uuid(self, uuid):
        """
        Set the uuid for a Styler.

        Parameters
        ----------
        uuid: str

        Returns
        -------
        self : Styler
        """
        self.uuid = uuid
        return self

    def set_caption(self, caption):
        """
        Set the caption on a Styler

        Parameters
        ----------
        caption: str

        Returns
        -------
        self : Styler
        """
        self.caption = caption
        return self

    def set_table_styles(self, table_styles):
        """
        Set table styles. These are placed in a `<style>` tag before the
        generated HTML table.

        Parameters
        ----------
        table_styles: list
            Each list item should be a dictionary with
            `selector` and `props` keys. `selector` should be a CSS
            selector that the style will be applied to (automatically
            prefixed by the table's UUID) and `props` should be a list of
            tuples with `(attribute, value)`.

        Returns
        -------
        self : Styler

        Notes
        -----
        Entries are rendered with CSS selectors according to gross properties
        of their position and content:
        * Index and Column names include ``index_name`` and ``level<k>``
          where `k` is its level in a MultiIndex
        * Index label cells include
            * ``row_heading``
            * ``row<n>`` where `n` is the numeric position of the row
            * ``level<k>`` where `k` is the level in a MultiIndex
        * Column label cells include
            * ``col_heading``
            * ``col<n>`` where `n` is the numeric position of the column
            * ``level<k>`` where `k` is the level in a MultiIndex

        * Blank cells include ``blank``
        * Data cells include ``data``

        Examples
        --------
        >>> df = pd.DataFrame(np.random.randn(10, 4))
        >>> df.style.set_table_styles(
        ...     [{'selector': 'tr:hover',
        ...       'props': [('background-color', 'yellow')]}]
        ... )
        """
        self.table_styles = table_styles
        return self

    def hide_index(self):
        """
        Hide any indices from rendering.

        .. versionadded:: 0.22.0

        Returns
        -------
        self : Styler
        """
        self.hidden_index = True
        return self

    def hide_columns(self, subset):
        """
        Hide columns from rendering.

        .. versionadded:: 0.22.0

        Parameters
        ----------
        subset: IndexSlice
            An argument to `DataFrame.loc` that identifies which columns
            are hidden.

        Returns
        -------
        self : Styler
        """
        subset = _non_reducing_slice(subset)
        hidden_df = self.data.loc[subset]
        self.hidden_columns = self.columns.get_indexer_for(hidden_df.columns)
        return self

    # -----------------------------------------------------------------------
    # A collection of "builtin" styles
    # -----------------------------------------------------------------------

    @staticmethod
    def _highlight_null(v, null_color):
        return ('background-color: {color}'.format(color=null_color)
                if pd.isna(v) else '')

    def highlight_null(self, null_color='red'):
        """
        Set the background color of table entries with missing values.

        Parameters
        ----------
        null_color: str

        Returns
        -------
        self : Styler
        """
        self.applymap(self._highlight_null, null_color=null_color)
        return self

    def background_gradient(self, cmap='PuBu', low=0, high=0, axis=0,
                            subset=None):
        """
        Color the background of table entries in a gradient according to
        the data in each column (optionally row).
        Requires matplotlib.

        Parameters
        ----------
        cmap: str or colormap
            matplotlib colormap
        low, high: float
            compress the range by these values.
        axis: int or str
            1 or 'columns' for columnwise, 0 or 'index' for rowwise
        subset: IndexSlice
            a valid slice for `data` to limit the style application to

        Returns
        -------
        self : Styler

        Notes
        -----
        Tune `low` and `high` to keep the text legible by
        not using the entire range of the color map. These extend
        the range of the data by `low * (x.max() - x.min())`
        and `high * (x.max() - x.min())` before normalizing.
        """
        subset = _maybe_numeric_slice(self.data, subset)
        subset = _non_reducing_slice(subset)
        self.apply(self._background_gradient, cmap=cmap, subset=subset,
                   axis=axis, low=low, high=high)
        return self

    @staticmethod
    def _background_gradient(s, cmap='PuBu', low=0, high=0):
        """Color background in a range according to the data."""
        with _mpl(Styler.background_gradient) as (plt, colors):
            rng = s.max() - s.min()
            # extend lower / upper bounds, compresses color range
            norm = colors.Normalize(s.min() - (rng * low),
                                    s.max() + (rng * high))
            # matplotlib modifies inplace?
            # https://github.com/matplotlib/matplotlib/issues/5427
            normed = norm(s.values)
            c = [colors.rgb2hex(x) for x in plt.cm.get_cmap(cmap)(normed)]
            return ['background-color: {color}'.format(color=color)
                    for color in c]

    def set_properties(self, subset=None, **kwargs):
        """
        Convenience method for setting one or more non-data dependent
        properties for each cell.

        Parameters
        ----------
        subset: IndexSlice
            a valid slice for `data` to limit the style application to
        kwargs: dict
            property: value pairs to be set for each cell

        Returns
        -------
        self : Styler

        Examples
        --------
        >>> df = pd.DataFrame(np.random.randn(10, 4))
        >>> df.style.set_properties(color="white", align="right")
        >>> df.style.set_properties(**{'background-color': 'yellow'})
        """
        values = ';'.join('{p}: {v}'.format(p=p, v=v)
                          for p, v in kwargs.items())
        f = lambda x: values
        return self.applymap(f, subset=subset)

    @staticmethod
    def _bar_left(s, color, width, base):
        """
        The minimum value is aligned at the left of the cell
        Parameters
        ----------
        color: 2-tuple/list, of [`color_negative`, `color_positive`]
        width: float
            A number between 0 or 100. The largest value will cover `width`
            percent of the cell's width
        base: str
            The base css format of the cell, e.g.:
            `base = 'width: 10em; height: 80%;'`
        Returns
        -------
        self : Styler
        """
        normed = width * (s - s.min()) / (s.max() - s.min())
        zero_normed = width * (0 - s.min()) / (s.max() - s.min())
        attrs = (base + 'background: linear-gradient(90deg,{c} {w:.1f}%, '
                        'transparent 0%)')

        return [base if x == 0 else attrs.format(c=color[0], w=x)
                if x < zero_normed
                else attrs.format(c=color[1], w=x) if x >= zero_normed
                else base for x in normed]

    @staticmethod
    def _bar_center_zero(s, color, width, base):
        """
        Creates a bar chart where the zero is centered in the cell
        Parameters
        ----------
        color: 2-tuple/list, of [`color_negative`, `color_positive`]
        width: float
            A number between 0 or 100. The largest value will cover `width`
            percent of the cell's width
        base: str
            The base css format of the cell, e.g.:
            `base = 'width: 10em; height: 80%;'`
        Returns
        -------
        self : Styler
        """

        # Either the min or the max should reach the edge
        # (50%, centered on zero)
        m = max(abs(s.min()), abs(s.max()))

        normed = s * 50 * width / (100.0 * m)

        attrs_neg = (base + 'background: linear-gradient(90deg, transparent 0%'
                     ', transparent {w:.1f}%, {c} {w:.1f}%, '
                     '{c} 50%, transparent 50%)')

        attrs_pos = (base + 'background: linear-gradient(90deg, transparent 0%'
                     ', transparent 50%, {c} 50%, {c} {w:.1f}%, '
                     'transparent {w:.1f}%)')

        return [attrs_pos.format(c=color[1], w=(50 + x)) if x >= 0
                else attrs_neg.format(c=color[0], w=(50 + x))
                for x in normed]

    @staticmethod
    def _bar_center_mid(s, color, width, base):
        """
        Creates a bar chart where the midpoint is centered in the cell
        Parameters
        ----------
        color: 2-tuple/list, of [`color_negative`, `color_positive`]
        width: float
            A number between 0 or 100. The largest value will cover `width`
            percent of the cell's width
        base: str
            The base css format of the cell, e.g.:
            `base = 'width: 10em; height: 80%;'`
        Returns
        -------
        self : Styler
        """

        if s.min() >= 0:
            # In this case, we place the zero at the left, and the max() should
            # be at width
            zero = 0.0
            slope = width / s.max()
        elif s.max() <= 0:
            # In this case, we place the zero at the right, and the min()
            # should be at 100-width
            zero = 100.0
            slope = width / -s.min()
        else:
            slope = width / (s.max() - s.min())
            zero = (100.0 + width) / 2.0 - slope * s.max()

        normed = zero + slope * s

        attrs_neg = (base + 'background: linear-gradient(90deg, transparent 0%'
                     ', transparent {w:.1f}%, {c} {w:.1f}%, '
                     '{c} {zero:.1f}%, transparent {zero:.1f}%)')

        attrs_pos = (base + 'background: linear-gradient(90deg, transparent 0%'
                     ', transparent {zero:.1f}%, {c} {zero:.1f}%, '
                     '{c} {w:.1f}%, transparent {w:.1f}%)')

        return [attrs_pos.format(c=color[1], zero=zero, w=x) if x > zero
                else attrs_neg.format(c=color[0], zero=zero, w=x)
                for x in normed]

    def bar(self, subset=None, axis=0, color='#d65f5f', width=100,
            align='left'):
        """
        Color the background of each cell with a bar of color `color`,
        and of length proportional to the value of the cell relative to those
        in its column.  Excludes non-numeric data by default.

        Parameters
        ----------
        subset: IndexSlice, default None
            a valid slice for `data` to limit the style application to
        axis: int
        color: str or 2-tuple/list
            If a str is passed, the color is the same for both
            negative and positive numbers. If 2-tuple/list is used, the
            first element is the color_negative and the second is the
            color_positive (eg: ['#d65f5f', '#5fba7d'])
        width: float
            A number between 0 or 100. The largest value will cover `width`
            percent of the cell's width
        align : {'left', 'zero',' mid'}, default 'left'
            - 'left' : the min value starts at the left of the cell
            - 'zero' : a value of zero is located at the center of the cell
            - 'mid' : the center of the cell is at (max-min)/2, or
              if values are all negative (positive) the zero is aligned
              at the right (left) of the cell

              .. versionadded:: 0.20.0

        Returns
        -------
        self : Styler
        """
        subset = _maybe_numeric_slice(self.data, subset)
        subset = _non_reducing_slice(subset)

        base = 'width: 10em; height: 80%;'

        if not(is_list_like(color)):
            color = [color, color]
        elif len(color) == 1:
            color = [color[0], color[0]]
        elif len(color) > 2:
            msg = ("Must pass `color` as string or a list-like"
                   " of length 2: [`color_negative`, `color_positive`]\n"
                   "(eg: color=['#d65f5f', '#5fba7d'])")
            raise ValueError(msg)

        if align == 'left':
            self.apply(self._bar_left, subset=subset, axis=axis, color=color,
                       width=width, base=base)
        elif align == 'zero':
            self.apply(self._bar_center_zero, subset=subset, axis=axis,
                       color=color, width=width, base=base)
        elif align == 'mid':
            self.apply(self._bar_center_mid, subset=subset, axis=axis,
                       color=color, width=width, base=base)
        else:
            msg = ("`align` must be one of {'left', 'zero',' mid'}")
            raise ValueError(msg)

        return self

    def highlight_max(self, subset=None, color='yellow', axis=0):
        """
        Set the color to shade the background of an entry
        containing the maximum value in a column, row, or table.

        Parameters
        ----------
        subset: IndexSlice, default None
            a valid slice for `data` to limit the style application to
        color: str, default 'yellow'
        axis: int, str, or None; default 0
            0 or 'index' for columnwise (default), 1 or 'columns' for rowwise,
            or `None` for tablewise

        Returns
        -------
        self : Styler
        """
        return self._highlight_handler(subset=subset, color=color, axis=axis,
                                       max_=True)

    def highlight_min(self, subset=None, color='yellow', axis=0):
        """
        Set the color to shade the background of an entry
        containing the maximum value in a column, row, or table.

        Parameters
        ----------
        subset: IndexSlice, default None
            a valid slice for `data` to limit the style application to
        color: str, default 'yellow'
        axis: int, str, or None; default 0
            0 or 'index' for columnwise (default), 1 or 'columns' for rowwise,
            or `None` for tablewise

        Returns
        -------
        self : Styler
        """
        return self._highlight_handler(subset=subset, color=color, axis=axis,
                                       max_=False)

    def _highlight_handler(self, subset=None, color='yellow', axis=None,
                           max_=True):
        subset = _non_reducing_slice(_maybe_numeric_slice(self.data, subset))
        self.apply(self._highlight_extrema, color=color, axis=axis,
                   subset=subset, max_=max_)
        return self

    @staticmethod
    def _highlight_extrema(data, color='yellow', max_=True):
        """Highlight the min or max in a Series or DataFrame"""
        attr = 'background-color: {0}'.format(color)
        if data.ndim == 1:  # Series from .apply
            if max_:
                extrema = data == data.max()
            else:
                extrema = data == data.min()
            return [attr if v else '' for v in extrema]
        else:  # DataFrame from .tee
            if max_:
                extrema = data == data.max().max()
            else:
                extrema = data == data.min().min()
            return pd.DataFrame(np.where(extrema, attr, ''),
                                index=data.index, columns=data.columns)

    @classmethod
    def from_custom_template(cls, searchpath, name):
        """
        Factory function for creating a subclass of `Styler`
        with a custom template and Jinja environment.

        Parameters
        ----------
        searchpath : str or list
            Path or paths of directories containing the templates
        name : str
            Name of your custom template to use for rendering

        Returns
        -------
        MyStyler : subclass of Styler
            has the correct `env` and `template` class attributes set.
        """
        loader = ChoiceLoader([
            FileSystemLoader(searchpath),
            cls.loader,
        ])

        class MyStyler(cls):
            env = Environment(loader=loader)
            template = env.get_template(name)

        return MyStyler


def _is_visible(idx_row, idx_col, lengths):
    """
    Index -> {(idx_row, idx_col): bool})
    """
    return (idx_col, idx_row) in lengths


def _get_level_lengths(index, hidden_elements=None):
    """
    Given an index, find the level length for each element.
    Optional argument is a list of index positions which
    should not be visible.

    Returns:
     (level, inital_position): span
    """
    sentinel = sentinel_factory()
    levels = index.format(sparsify=sentinel, adjoin=False, names=False)

    if hidden_elements is None:
        hidden_elements = []

    lengths = {}
    if index.nlevels == 1:
        for i, value in enumerate(levels):
            if(i not in hidden_elements):
                lengths[(0, i)] = 1
        return lengths

    for i, lvl in enumerate(levels):
        for j, row in enumerate(lvl):
            if not get_option('display.multi_sparse'):
                lengths[(i, j)] = 1
            elif (row != sentinel) and (j not in hidden_elements):
                last_label = j
                lengths[(i, last_label)] = 1
            elif (row != sentinel):
                # even if its hidden, keep track of it in case
                # length >1 and later elemens are visible
                last_label = j
                lengths[(i, last_label)] = 0
            elif(j not in hidden_elements):
                lengths[(i, last_label)] += 1

    non_zero_lengths = {}
    for element, length in lengths.items():
        if(length >= 1):
            non_zero_lengths[element] = length

    return non_zero_lengths


def _maybe_wrap_formatter(formatter):
    if is_string_like(formatter):
        return lambda x: formatter.format(x)
    elif callable(formatter):
        return formatter
    else:
        msg = ("Expected a template string or callable, got {formatter} "
               "instead".format(formatter=formatter))
        raise TypeError(msg)
