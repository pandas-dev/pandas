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
    import warnings
    warnings.warn("Importing jinja")
    from jinja2 import Template
    has_jinja = True
except ImportError:
    has_jinja = False

import numpy as np
import pandas as pd
import pandas.core.common as com
from pandas.compat import lzip
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
    data.

    .. versionadded:: 0.17.1

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
    ctx : defaultdict
        maps (row, column) -> [styles] where each style is a
        TODO: <(attribute, value)> paire | <string with 'attribute: value'>

    Notes
    -----
    Most styling will be done by passing style functions into
    Styler.tee, Styler.apply, or Styler.applymap. Style functions should
    return values with strings containing 'attr: value' that will be applied
    to the indicated cells.

    If using in the Jupyter notebook, Styler has defined a _repr_html_
    to automatically render itself. Otherwise call Styler.render to get
    the genterated HTML.
    """
    if has_jinja:
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

            <table id="T_{{uuid}}">
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
                 caption=None):
        if not has_jinja:
            msg = "pandas.Styler requires jinja2. "\
                  "Please install with `conda install Jinja2`\n"\
                  "or `pip install Jinja2`"
            raise ImportError(msg)
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

    def _repr_html_(self):
        '''
        Hooks into Jupyter notebook rich display system.
        '''
        return self.render()

    def _translate(self, table_styles=None, caption=None, uuid=None,
                   precision=None):
        """
        Convert the DataFrame in `self.data` and the attrs from `_build_styles`
        into a dictionary of {head, body, uuid, cellstyle}
        """
        table_styles = table_styles or self.table_styles
        if table_styles is None:
            table_styles = []
        caption = caption or self.caption
        precision = precision or self.precision
        ctx = self.ctx
        uuid = uuid or str(uuid1()).replace("-", "_")
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

        # uuid required to isolate table styling from others
        # in same notebook in ipnb
        return dict(head=head, cellstyle=cellstyle, body=body, uuid=uuid,
                    precision=self.precision, table_styles=table_styles,
                    caption=caption)

    def render(self, table_styles=None, caption=None, uuid=None,
               precision=None, trim=True):
        """
        Render the built up styles to HTML

        Parameters
        ----------
        table_styles: list or None
            defaults to ``self.table_styles``
        caption: str or None
            defaults to ``self.caption``
        uuid: str or None
            defaults to random uuid
        precision: int or None
            defaults to self.precision
        trim: bool
            Small optimization to not create classes for cells that
            have no styles applied to them.

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
        table_styles = table_styles or self.table_styles
        caption = caption or self.caption
        uuid = uuid or self.uuid
        d = self._translate(table_styles=table_styles, caption=caption,
                            uuid=uuid)
        if trim:
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
        else:
            styler.ctx = self.ctx
        return styler

    def __copy__(self):
        """
        Deep copy by default.
        """
        return self._copy(deepcopy=False)

    def __deepcopy__(self, memo):
        return self._copy(deepcopy=True)

    def clear(self):
        self.ctx.clear()

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
        Apply a function, updating the HTML representation with the result.

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
        This is similar to DataFrame.apply, except that axis=None applies
        the function to the entire DataFrame at once, rather tha column
        or rowwise.

        Examples
        --------

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

        Notes
        -----
        Examples
        --------
        """
        self._todo.append((lambda instance: getattr(instance, '_applymap'),
                          (func, subset),
                          kwargs))
        return self

    def set_precision(self, precision):
        self.precision = precision
        return self

    def export(self):
        return self._todo

    def set(self, styles):
        self._todo.extend(styles)
        return self

    def set_uuid(self, uuid):
        self.uuid = uuid
        return self

    def set_caption(self, caption):
        self.caption = caption
        return self

    def set_table_styles(self, table_styles):
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

    def text_shadow(self, color='black'):
        self.set_properties(**{"text-shadow": "-1px 0 {color}, 0 1px {color}, "
                               "1px 0 {color}, 0 -1px {color};".format(
                                   color=color)})
        return self

    def set_properties(self, subset=None, **kwargs):
        """
        Convience method for setting a non-data dependent properties
        on each element.

        .. versionadded:: 0.17.1

        Parameters
        ----------
        subset: IndexSlice
            a valid slice for ``data``
        kwargs: dict
            property: value pairs to be set for each cell

        Returns
        -------
        self : Styler

        Examples
        --------
        df.stle.set_properties(color="white", align="right")
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

    def bar(self, subset=None, axis=0, color='#FFC0CB', width=100):
        """
        Color the background `color` proptional to the values in each column.
        Excludes non-numeric data by default.

        .. versionadded:: 0.17.1

        Parameters
        ----------
        subset: IndexSlice, default None
        axis: int
        color: str
        width: float
            A number between 0 or 100. The largest value will cover this
            percent of the cell's width
        """
        subset = _maybe_numeric_slice(self.data, subset)
        subset = _non_reducing_slice(subset)
        self.apply(self._bar, subset=subset, axis=axis, color=color,
                   width=width)
        return self

    def highlight_max(self, subset=None, color='yellow', axis=None):
        """
        Highlight the maximum by shading the background

        .. versionadded:: 0.17.1

        Parameters
        ----------
        subset: IndexSlice, default None
        color: str, default 'yellow'
        axis: int, str, or None
            0 or 'index' for columnwise, 1 or 'columns' for rowwise
            or ``None`` for tablewise
        """
        return self._highlight_handler(subset=subset, color=color, axis=axis,
                                       max_=True)

    def highlight_min(self, subset=None, color='yellow', axis=None):
        """
        Highlight the minimum by shading the background

        .. versionadded:: 0.17.1

        Parameters
        ----------
        subset: IndexSlice, default None
        color: str, default 'yellow'
        axis: int, str, or None
            0 or 'index' for columnwise, 1 or 'columns' for rowwise
            or ``None`` for tablewise
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


def _non_reducing_slice(slice_):
    """
    Ensurse that a slice doesn't reduce to a Series or Scalar.

    Any user-paseed `subset` should have this called on it
    to make sure we're always working with DataFrames.
    """
    # default to column slice, like DataFrame
    # ['A', 'B'] -> IndexSlices[:, ['A', 'B']]
    kinds = tuple(list(pd.compat.string_types) +
                  [pd.Series, np.ndarray, pd.Index, list])
    if isinstance(slice_, kinds):
        slice_ = pd.IndexSlice[:, slice_]

    def pred(part):
        # true when slice does *not* reduce
        return isinstance(part, slice) or com.is_list_like(part)

    if not com.is_list_like(slice_):
        if not isinstance(slice_, slice):
            # a 1-d slice, like df.loc[1]
            slice_ = [[slice_]]
        else:
            # slice(a, b, c)
            slice_ = [slice_]  # to tuplize later
    else:
        slice_ = [part if pred(part) else [part] for part in slice_]
    return tuple(slice_)


def _maybe_numeric_slice(df, slice_, include_bool=False):
    """
    want nice defaults for background_gradient that don't break
    with non-numeric data. But if slice_ is passed go with that.
    """
    if slice_ is None:
        dtypes = [np.number]
        if include_bool:
            dtypes.append(bool)
        slice_ = pd.IndexSlice[:, df.select_dtypes(include=dtypes).columns]
    return slice_
