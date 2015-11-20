"""
Module for applying conditional formatting to
DataFrames and Series.
"""
from functools import partial
from contextlib import contextmanager
from uuid import uuid1
import copy
from collections import defaultdict

try:
    from jinja2 import Template
except ImportError:
    msg = "pandas.Styler requires jinja2. "\
          "Please install with `conda install Jinja2`\n"\
          "or `pip install Jinja2`"
    raise ImportError(msg)

import numpy as np
import pandas as pd
from pandas.compat import lzip
from pandas.core.indexing import _maybe_numeric_slice, _non_reducing_slice
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
    Helps style a DataFrame or Series according to the
    data with HTML and CSS.

    .. versionadded:: 0.17.1

    .. warning::
        This is a new feature and is under active development.
        We'll be adding features and possibly making breaking changes in future
        releases.

    Parameters
    ----------
    data: Series or DataFrame
    precision: int
        precision to round floats to, defaults to pd.options.display.precision
    table_styles: list-like, default None
        list of {selector: (attr, value)} dicts; see Notes
    uuid: str, default None
        a unique identifier to avoid CSS collisons; generated automatically
    caption: str, default None
        caption to attach to the table

    Attributes
    ----------
    tempate: Jinja Template

    Notes
    -----
    Most styling will be done by passing style functions into
    ``Styler.apply`` or ``Styler.applymap``. Style functions should
    return values with strings containing CSS ``'attr: value'`` that will
    be applied to the indicated cells.

    If using in the Jupyter notebook, Styler has defined a ``_repr_html_``
    to automatically render itself. Otherwise call Styler.render to get
    the genterated HTML.

    See Also
    --------
    pandas.DataFrame.style
    """
    template = Template("""
        <style  type="text/css" >
        {% for s in table_styles %}
            #T_{{uuid}} {{s.selector}} {
            {% for p,val in s.props %}
                {{p}}: {{val}};
            {% endfor %}
            }
        {% endfor %}
        {% for s in cellstyle %}
            #T_{{uuid}}{{s.selector}} {
            {% for p,val in s.props %}
                {{p}}: {{val}};
            {% endfor %}
            }
        {% endfor %}
        </style>

        <table id="T_{{uuid}}" {{ table_attributes }}>
        {% if caption %}
            <caption>{{caption}}</caption>
        {% endif %}

        <thead>
            {% for r in head %}
            <tr>
                {% for c in r %}
                <{{c.type}} class="{{c.class}}">{{c.value}}
                {% endfor %}
            </tr>
            {% endfor %}
        </thead>
        <tbody>
            {% for r in body %}
            <tr>
                {% for c in r %}
                <{{c.type}} id="T_{{uuid}}{{c.id}}" class="{{c.class}}">
                    {% if c.value is number %}
                        {{c.value|round(precision)}}
                    {% else %}
                        {{c.value}}
                    {% endif %}
                {% endfor %}
            </tr>
            {% endfor %}
        </tbody>
        </table>
        """)

    def __init__(self, data, precision=None, table_styles=None, uuid=None,
                 caption=None, table_attributes=None):
        self.ctx = defaultdict(list)
        self._todo = []

        if not isinstance(data, (pd.Series, pd.DataFrame)):
            raise TypeError
        if data.ndim == 1:
            data = data.to_frame()
        if not data.index.is_unique or not data.columns.is_unique:
            raise ValueError("style is not supported for non-unique indicies.")

        self.data = data
        self.index = data.index
        self.columns = data.columns

        self.uuid = uuid
        self.table_styles = table_styles
        self.caption = caption
        if precision is None:
            precision = pd.options.display.precision
        self.precision = precision
        self.table_attributes = table_attributes

    def _repr_html_(self):
        '''
        Hooks into Jupyter notebook rich display system.
        '''
        return self.render()

    def _translate(self):
        """
        Convert the DataFrame in `self.data` and the attrs from `_build_styles`
        into a dictionary of {head, body, uuid, cellstyle}
        """
        table_styles = self.table_styles or []
        caption = self.caption
        ctx = self.ctx
        precision = self.precision
        uuid = self.uuid or str(uuid1()).replace("-", "_")
        ROW_HEADING_CLASS = "row_heading"
        COL_HEADING_CLASS = "col_heading"
        DATA_CLASS = "data"
        BLANK_CLASS = "blank"
        BLANK_VALUE = ""

        cell_context = dict()

        n_rlvls = self.data.index.nlevels
        n_clvls = self.data.columns.nlevels
        rlabels = self.data.index.tolist()
        clabels = self.data.columns.tolist()

        idx_values = self.data.index.format(sparsify=False, adjoin=False,
                                            names=False)
        idx_values = lzip(*idx_values)

        if n_rlvls == 1:
            rlabels = [[x] for x in rlabels]
        if n_clvls == 1:
            clabels = [[x] for x in clabels]
        clabels = list(zip(*clabels))

        cellstyle = []
        head = []

        for r in range(n_clvls):
            row_es = [{"type": "th", "value": BLANK_VALUE,
                       "class": " ".join([BLANK_CLASS])}] * n_rlvls
            for c in range(len(clabels[0])):
                cs = [COL_HEADING_CLASS, "level%s" % r, "col%s" % c]
                cs.extend(cell_context.get(
                    "col_headings", {}).get(r, {}).get(c, []))
                row_es.append({"type": "th", "value": clabels[r][c],
                               "class": " ".join(cs)})
            head.append(row_es)

        body = []
        for r, idx in enumerate(self.data.index):
            cs = [ROW_HEADING_CLASS, "level%s" % c, "row%s" % r]
            cs.extend(cell_context.get(
                "row_headings", {}).get(r, {}).get(c, []))
            row_es = [{"type": "th",
                       "value": rlabels[r][c],
                       "class": " ".join(cs)}
                      for c in range(len(rlabels[r]))]

            for c, col in enumerate(self.data.columns):
                cs = [DATA_CLASS, "row%s" % r, "col%s" % c]
                cs.extend(cell_context.get("data", {}).get(r, {}).get(c, []))
                row_es.append({"type": "td", "value": self.data.iloc[r][c],
                               "class": " ".join(cs), "id": "_".join(cs[1:])})
                props = []
                for x in ctx[r, c]:
                    # have to handle empty styles like ['']
                    if x.count(":"):
                        props.append(x.split(":"))
                    else:
                        props.append(['', ''])
                cellstyle.append(
                    {'props': props,
                     'selector': "row%s_col%s" % (r, c)}
                )
            body.append(row_es)

        return dict(head=head, cellstyle=cellstyle, body=body, uuid=uuid,
                    precision=precision, table_styles=table_styles,
                    caption=caption, table_attributes=self.table_attributes)

    def render(self):
        """
        Render the built up styles to HTML

        .. versionadded:: 0.17.1

        Returns
        -------
        rendered: str
            the rendered HTML

        Notes
        -----
        ``Styler`` objects have defined the ``_repr_html_`` method
        which automatically calls ``self.render()`` when it's the
        last item in a Notebook cell. When calling ``Styler.render()``
        directly, wrap the resul in ``IPython.display.HTML`` to view
        the rendered HTML in the notebook.
        """
        self._compute()
        d = self._translate()
        # filter out empty styles, every cell will have a class
        # but the list of props may just be [['', '']].
        # so we have the neested anys below
        trimmed = [x for x in d['cellstyle'] if
                    any(any(y) for y in x['props'])]
        d['cellstyle'] = trimmed
        return self.template.render(**d)

    def _update_ctx(self, attrs):
        """
        update the state of the Styler. Collects a mapping
        of {index_label: ['<property>: <value>']}

        attrs: Series or DataFrame
        should contain strings of '<property>: <value>;<prop2>: <val2>'
        Whitespace shouldn't matter and the final trailing ';' shouldn't
        matter.
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
        Deep copy by default.
        """
        return self._copy(deepcopy=False)

    def __deepcopy__(self, memo):
        return self._copy(deepcopy=True)

    def clear(self):
        '''
        "Reset" the styler, removing any previously applied styles.
        Returns None.
        '''
        self.ctx.clear()
        self._todo = []

    def _compute(self):
        '''
        Execute the style functions built up in `self._todo`.

        Relies on the conventions that all style functions go through
        .apply or .applymap. The append styles to apply as tuples of

        (application method, *args, **kwargs)
        '''
        r = self
        for func, args, kwargs in self._todo:
            r = func(self)(*args, **kwargs)
        return r

    def _apply(self, func, axis=0, subset=None, **kwargs):
        subset = slice(None) if subset is None else subset
        subset = _non_reducing_slice(subset)
        if axis is not None:
            result = self.data.loc[subset].apply(func, axis=axis, **kwargs)
        else:
            # like tee
            result = func(self.data.loc[subset], **kwargs)
        self._update_ctx(result)
        return self

    def apply(self, func, axis=0, subset=None, **kwargs):
        """
        Apply a function column-wise, row-wise, or table-wase,
        updating the HTML representation with the result.

        .. versionadded:: 0.17.1

        Parameters
        ----------
        func: function
        axis: int, str or None
            apply to each column (``axis=0`` or ``'index'``)
            or to each row (``axis=1`` or ``'columns'``) or
            to the entire DataFrame at once with ``axis=None``.
        subset: IndexSlice
            a valid indexer to limit ``data`` to *before* applying the
            function. Consider using a pandas.IndexSlice
        kwargs: dict
            pass along to ``func``

        Returns
        -------
        self

        Notes
        -----
        This is similar to ``DataFrame.apply``, except that ``axis=None``
        applies the function to the entire DataFrame at once,
        rather than column-wise or row-wise.
        """
        self._todo.append((lambda instance: getattr(instance, '_apply'),
                          (func, axis, subset),
                           kwargs))
        return self

    def _applymap(self, func, subset=None, **kwargs):
        func = partial(func, **kwargs)  # applymap doesn't take kwargs?
        if subset is None:
            subset = pd.IndexSlice[:]
        subset = _non_reducing_slice(subset)
        result = self.data.loc[subset].applymap(func)
        self._update_ctx(result)
        return self

    def applymap(self, func, subset=None, **kwargs):
        """
        Apply a function elementwise, updating the HTML
        representation with the result.

        .. versionadded:: 0.17.1

        Parameters
        ----------
        func : function
        subset : IndexSlice
            a valid indexer to limit ``data`` to *before* applying the
            function. Consider using a pandas.IndexSlice
        kwargs : dict
            pass along to ``func``

        Returns
        -------
        self

        """
        self._todo.append((lambda instance: getattr(instance, '_applymap'),
                          (func, subset),
                          kwargs))
        return self

    def set_precision(self, precision):
        """
        Set the precision used to render.

        .. versionadded:: 0.17.1

        Parameters
        ----------
        precision: int

        Returns
        -------
        self
        """
        self.precision = precision
        return self

    def set_table_attributes(self, attributes):
        """
        Set the table attributes. These are the items
        that show up in the opening ``<table>`` tag in addition
        to to automatic (by default) id.

        .. versionadded:: 0.17.1

        Parameters
        ----------
        precision: int

        Returns
        -------
        self
        """
        self.table_attributes = attributes
        return self

    def export(self):
        """
        Export the styles to applied to the current Styler.
        Can be applied to a second style with ``Styler.use``.

        .. versionadded:: 0.17.1

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
        from ``Styler.export``.

        .. versionadded:: 0.17.1

        Parameters
        ----------
        styles: list
            list of style functions

        Returns
        -------
        self

        See Also
        --------
        Styler.export
        """
        self._todo.extend(styles)
        return self

    def set_uuid(self, uuid):
        """
        Set the uuid for a Styler.

        .. versionadded:: 0.17.1

        Parameters
        ----------
        uuid: str

        Returns
        -------
        self
        """
        self.uuid = uuid
        return self

    def set_caption(self, caption):
        """
        Se the caption on a Styler

        .. versionadded:: 0.17.1

        Parameters
        ----------
        caption: str

        Returns
        -------
        self
        """
        self.caption = caption
        return self

    def set_table_styles(self, table_styles):
        """
        Set the table styles on a Styler

        .. versionadded:: 0.17.1

        Parameters
        ----------
        table_styles: list

        Returns
        -------
        self
        """
        self.table_styles = table_styles
        return self

    # -----------------------------------------------------------------------
    # A collection of "builtin" styles
    # -----------------------------------------------------------------------

    @staticmethod
    def _highlight_null(v, null_color):
        return 'background-color: %s' % null_color if pd.isnull(v) else ''

    def highlight_null(self, null_color='red'):
        """
        Shade the background ``null_color`` for missing values.

        .. versionadded:: 0.17.1

        Parameters
        ----------
        null_color: str

        Returns
        -------
        self
        """
        self.applymap(self._highlight_null, null_color=null_color)
        return self

    def background_gradient(self, cmap='PuBu', low=0, high=0,
                            axis=0, subset=None):
        """
        Color the background in a gradient according to
        the data in each column (optionally row).
        Requires matplotlib.

        .. versionadded:: 0.17.1

        Parameters
        ----------
        cmap: str or colormap
            matplotlib colormap
        low, high: float
            compress the range by these values.
        axis: int or str
            1 or 'columns' for colunwise, 0 or 'index' for rowwise
        subset: IndexSlice
            a valid slice for ``data`` to limit the style application to

        Returns
        -------
        self

        Notes
        -----
        Tune ``low`` and ``high`` to keep the text legible by
        not using the entire range of the color map. These extend
        the range of the data by ``low * (x.max() - x.min())``
        and ``high * (x.max() - x.min())`` before normalizing.
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
            return ['background-color: %s' % color for color in c]

    def set_properties(self, subset=None, **kwargs):
        """
        Convience method for setting one or more non-data dependent
        properties or each cell.

        .. versionadded:: 0.17.1

        Parameters
        ----------
        subset: IndexSlice
            a valid slice for ``data`` to limit the style application to
        kwargs: dict
            property: value pairs to be set for each cell

        Returns
        -------
        self : Styler

        Examples
        --------
        >>> df = pd.DataFrame(np.random.randn(10, 4))
        >>> df.style.set_properties(color="white", align="right")
        """
        values = ';'.join('{p}: {v}'.format(p=p, v=v) for p, v in
                          kwargs.items())
        f = lambda x: values
        return self.applymap(f, subset=subset)

    @staticmethod
    def _bar(s, color, width):
        normed = width * (s - s.min()) / (s.max() - s.min())
        attrs = 'width: 10em; height: 80%;'\
            'background: linear-gradient(90deg,'\
            '{c} {w}%, transparent 0%)'
        return [attrs.format(c=color, w=x) for x in normed]

    def bar(self, subset=None, axis=0, color='#d65f5f', width=100):
        """
        Color the background ``color`` proptional to the values in each column.
        Excludes non-numeric data by default.

        .. versionadded:: 0.17.1

        Parameters
        ----------
        subset: IndexSlice, default None
            a valid slice for ``data`` to limit the style application to
        axis: int
        color: str
        width: float
            A number between 0 or 100. The largest value will cover ``width``
            percent of the cell's width

        Returns
        -------
        self
        """
        subset = _maybe_numeric_slice(self.data, subset)
        subset = _non_reducing_slice(subset)
        self.apply(self._bar, subset=subset, axis=axis, color=color,
                   width=width)
        return self

    def highlight_max(self, subset=None, color='yellow', axis=0):
        """
        Highlight the maximum by shading the background

        .. versionadded:: 0.17.1

        Parameters
        ----------
        subset: IndexSlice, default None
            a valid slice for ``data`` to limit the style application to
        color: str, default 'yellow'
        axis: int, str, or None; default None
            0 or 'index' for columnwise, 1 or 'columns' for rowwise
            or ``None`` for tablewise (the default)

        Returns
        -------
        self
        """
        return self._highlight_handler(subset=subset, color=color, axis=axis,
                                       max_=True)

    def highlight_min(self, subset=None, color='yellow', axis=0):
        """
        Highlight the minimum by shading the background

        .. versionadded:: 0.17.1

        Parameters
        ----------
        subset: IndexSlice, default None
            a valid slice for ``data`` to limit the style application to
        color: str, default 'yellow'
        axis: int, str, or None; default None
            0 or 'index' for columnwise, 1 or 'columns' for rowwise
            or ``None`` for tablewise (the default)

        Returns
        -------
        self
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
        '''
        highlight the min or max in a Series or DataFrame
        '''
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



