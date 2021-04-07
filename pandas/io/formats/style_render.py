from __future__ import annotations

from collections import defaultdict
from functools import partial
from typing import (
    Any,
    Callable,
    DefaultDict,
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)
from uuid import uuid4

import numpy as np

from pandas._config import get_option

from pandas._libs import lib
from pandas._typing import FrameOrSeriesUnion
from pandas.compat._optional import import_optional_dependency

from pandas.core.dtypes.generic import ABCSeries

from pandas import (
    DataFrame,
    Index,
    IndexSlice,
    MultiIndex,
    Series,
    isna,
)
from pandas.api.types import is_list_like
import pandas.core.common as com

jinja2 = import_optional_dependency("jinja2", extra="DataFrame.style requires jinja2.")
from markupsafe import escape as escape_func  # markupsafe is jinja2 dependency

BaseFormatter = Union[str, Callable]
ExtFormatter = Union[BaseFormatter, Dict[Any, Optional[BaseFormatter]]]
CSSPair = Tuple[str, Union[str, int, float]]
CSSList = List[CSSPair]
CSSProperties = Union[str, CSSList]
CSSStyles = List[Dict[str, CSSProperties]]  # = List[CSSDict]
# class CSSDict(TypedDict):  # available when TypedDict is valid in pandas
#     selector: str
#     props: CSSProperties


class StylerRenderer:
    """
    Base class to process rendering a Styler with a specified jinja2 template.
    """

    loader = jinja2.PackageLoader("pandas", "io/formats/templates")
    env = jinja2.Environment(loader=loader, trim_blocks=True)
    template = env.get_template("html.tpl")

    def __init__(
        self,
        data: FrameOrSeriesUnion,
        uuid: str | None = None,
        uuid_len: int = 5,
        table_styles: CSSStyles | None = None,
        table_attributes: str | None = None,
        caption: str | None = None,
        cell_ids: bool = True,
    ):

        # validate ordered args
        if isinstance(data, Series):
            data = data.to_frame()
        if not isinstance(data, DataFrame):
            raise TypeError("``data`` must be a Series or DataFrame")
        if not data.index.is_unique or not data.columns.is_unique:
            raise ValueError("style is not supported for non-unique indices.")
        self.data: DataFrame = data
        self.index: Index = data.index
        self.columns: Index = data.columns
        if not isinstance(uuid_len, int) or not uuid_len >= 0:
            raise TypeError("``uuid_len`` must be an integer in range [0, 32].")
        self.uuid_len = min(32, uuid_len)
        self.uuid = (uuid or uuid4().hex[: self.uuid_len]) + "_"
        self.table_styles = table_styles
        self.table_attributes = table_attributes
        self.caption = caption
        self.cell_ids = cell_ids

        # add rendering variables
        self.hidden_index: bool = False
        self.hidden_columns: Sequence[int] = []
        self.ctx: DefaultDict[tuple[int, int], CSSList] = defaultdict(list)
        self.cell_context: DefaultDict[tuple[int, int], str] = defaultdict(str)
        self._todo: list[tuple[Callable, tuple, dict]] = []
        self.tooltips: Tooltips | None = None
        def_precision = get_option("display.precision")
        self._display_funcs: DefaultDict[  # maps (row, col) -> formatting function
            tuple[int, int], Callable[[Any], str]
        ] = defaultdict(lambda: partial(_default_formatter, precision=def_precision))

    def render(self, **kwargs) -> str:
        """
        Render the ``Styler`` including all applied styles to HTML.

        Parameters
        ----------
        **kwargs
            Any additional keyword arguments are passed
            through to ``self.template.render``.
            This is useful when you need to provide
            additional variables for a custom template.

        Returns
        -------
        rendered : str
            The rendered HTML.

        Notes
        -----
        Styler objects have defined the ``_repr_html_`` method
        which automatically calls ``self.render()`` when it's the
        last item in a Notebook cell. When calling ``Styler.render()``
        directly, wrap the result in ``IPython.display.HTML`` to view
        the rendered HTML in the notebook.

        Pandas uses the following keys in render. Arguments passed
        in ``**kwargs`` take precedence, so think carefully if you want
        to override them:

        * head
        * cellstyle
        * body
        * uuid
        * table_styles
        * caption
        * table_attributes
        """
        self._compute()
        # TODO: namespace all the pandas keys
        d = self._translate()
        d.update(kwargs)
        return self.template.render(**d)

    def _compute(self):
        """
        Execute the style functions built up in `self._todo`.

        Relies on the conventions that all style functions go through
        .apply or .applymap. The append styles to apply as tuples of

        (application method, *args, **kwargs)
        """
        self.ctx.clear()
        r = self
        for func, args, kwargs in self._todo:
            r = func(self)(*args, **kwargs)
        return r

    def _translate(self):
        """
        Convert the DataFrame in `self.data` and the attrs from `_build_styles`
        into a dictionary of {head, body, uuid, cellstyle}.
        """
        ROW_HEADING_CLASS = "row_heading"
        COL_HEADING_CLASS = "col_heading"
        INDEX_NAME_CLASS = "index_name"

        DATA_CLASS = "data"
        BLANK_CLASS = "blank"
        BLANK_VALUE = "&nbsp;"

        # mapping variables
        ctx = self.ctx  # td css styles from apply() and applymap()
        cell_context = self.cell_context  # td css classes from set_td_classes()
        cellstyle_map: DefaultDict[tuple[CSSPair, ...], list[str]] = defaultdict(list)

        # copied attributes
        hidden_index = self.hidden_index
        hidden_columns = self.hidden_columns

        # construct render dict
        d = {
            "uuid": self.uuid,
            "table_styles": _format_table_styles(self.table_styles or []),
            "caption": self.caption,
        }

        # for sparsifying a MultiIndex
        idx_lengths = _get_level_lengths(self.index)
        col_lengths = _get_level_lengths(self.columns, hidden_columns)

        n_rlvls = self.data.index.nlevels
        n_clvls = self.data.columns.nlevels
        rlabels = self.data.index.tolist()
        clabels = self.data.columns.tolist()

        if n_rlvls == 1:
            rlabels = [[x] for x in rlabels]
        if n_clvls == 1:
            clabels = [[x] for x in clabels]
        clabels = list(zip(*clabels))

        head = []
        for r in range(n_clvls):
            # Blank for Index columns...
            row_es = [
                {
                    "type": "th",
                    "value": BLANK_VALUE,
                    "display_value": BLANK_VALUE,
                    "is_visible": not hidden_index,
                    "class": " ".join([BLANK_CLASS]),
                }
            ] * (n_rlvls - 1)

            # ... except maybe the last for columns.names
            name = self.data.columns.names[r]
            cs = [
                BLANK_CLASS if name is None else INDEX_NAME_CLASS,
                f"level{r}",
            ]
            name = BLANK_VALUE if name is None else name
            row_es.append(
                {
                    "type": "th",
                    "value": name,
                    "display_value": name,
                    "class": " ".join(cs),
                    "is_visible": not hidden_index,
                }
            )

            if clabels:
                for c, value in enumerate(clabels[r]):
                    es = {
                        "type": "th",
                        "value": value,
                        "display_value": value,
                        "class": f"{COL_HEADING_CLASS} level{r} col{c}",
                        "is_visible": _is_visible(c, r, col_lengths),
                    }
                    colspan = col_lengths.get((r, c), 0)
                    if colspan > 1:
                        es["attributes"] = f'colspan="{colspan}"'
                    row_es.append(es)
                head.append(row_es)

        if (
            self.data.index.names
            and com.any_not_none(*self.data.index.names)
            and not hidden_index
        ):
            index_header_row = []

            for c, name in enumerate(self.data.index.names):
                cs = [INDEX_NAME_CLASS, f"level{c}"]
                name = "" if name is None else name
                index_header_row.append(
                    {"type": "th", "value": name, "class": " ".join(cs)}
                )

            index_header_row.extend(
                [
                    {
                        "type": "th",
                        "value": BLANK_VALUE,
                        "class": " ".join([BLANK_CLASS, f"col{c}"]),
                    }
                    for c in range(len(clabels[0]))
                    if c not in hidden_columns
                ]
            )

            head.append(index_header_row)
        d.update({"head": head})

        body = []
        for r, row_tup in enumerate(self.data.itertuples()):
            row_es = []
            for c, value in enumerate(rlabels[r]):
                rid = [
                    ROW_HEADING_CLASS,
                    f"level{c}",
                    f"row{r}",
                ]
                es = {
                    "type": "th",
                    "is_visible": (_is_visible(r, c, idx_lengths) and not hidden_index),
                    "value": value,
                    "display_value": value,
                    "id": "_".join(rid[1:]),
                    "class": " ".join(rid),
                }
                rowspan = idx_lengths.get((c, r), 0)
                if rowspan > 1:
                    es["attributes"] = f'rowspan="{rowspan}"'
                row_es.append(es)

            for c, value in enumerate(row_tup[1:]):
                formatter = self._display_funcs[(r, c)]
                row_dict = {
                    "type": "td",
                    "value": value,
                    "display_value": formatter(value),
                    "is_visible": (c not in hidden_columns),
                    "attributes": "",
                }

                # only add an id if the cell has a style
                props: CSSList = []
                if self.cell_ids or (r, c) in ctx:
                    row_dict["id"] = f"row{r}_col{c}"
                    props.extend(ctx[r, c])

                # add custom classes from cell context
                cls = ""
                if (r, c) in cell_context:
                    cls = " " + cell_context[r, c]
                row_dict["class"] = f"{DATA_CLASS} row{r} col{c}{cls}"

                row_es.append(row_dict)
                if props:  # (), [] won't be in cellstyle_map, cellstyle respectively
                    cellstyle_map[tuple(props)].append(f"row{r}_col{c}")
            body.append(row_es)
        d.update({"body": body})

        cellstyle: list[dict[str, CSSList | list[str]]] = [
            {"props": list(props), "selectors": selectors}
            for props, selectors in cellstyle_map.items()
        ]
        d.update({"cellstyle": cellstyle})

        table_attr = self.table_attributes
        use_mathjax = get_option("display.html.use_mathjax")
        if not use_mathjax:
            table_attr = table_attr or ""
            if 'class="' in table_attr:
                table_attr = table_attr.replace('class="', 'class="tex2jax_ignore ')
            else:
                table_attr += ' class="tex2jax_ignore"'
        d.update({"table_attributes": table_attr})

        if self.tooltips:
            d = self.tooltips._translate(self.data, self.uuid, d)

        return d

    def format(
        self,
        formatter: ExtFormatter | None = None,
        subset: slice | Sequence[Any] | None = None,
        na_rep: str | None = None,
        precision: int | None = None,
        escape: bool = False,
    ) -> StylerRenderer:
        """
        Format the text display value of cells.

        Parameters
        ----------
        formatter : str, callable, dict or None
            Object to define how values are displayed. See notes.
        subset : IndexSlice
            An argument to ``DataFrame.loc`` that restricts which elements
            ``formatter`` is applied to.
        na_rep : str, optional
            Representation for missing values.
            If ``na_rep`` is None, no special formatting is applied.

            .. versionadded:: 1.0.0

        precision : int, optional
            Floating point precision to use for display purposes, if not determined by
            the specified ``formatter``.

            .. versionadded:: 1.3.0

        escape : bool, default False
            Replace the characters ``&``, ``<``, ``>``, ``'``, and ``"`` in cell display
            string with HTML-safe sequences. Escaping is done before ``formatter``.

            .. versionadded:: 1.3.0

        Returns
        -------
        self : Styler

        Notes
        -----
        This method assigns a formatting function, ``formatter``, to each cell in the
        DataFrame. If ``formatter`` is ``None``, then the default formatter is used.
        If a callable then that function should take a data value as input and return
        a displayable representation, such as a string. If ``formatter`` is
        given as a string this is assumed to be a valid Python format specification
        and is wrapped to a callable as ``string.format(x)``. If a ``dict`` is given,
        keys should correspond to column names, and values should be string or
        callable, as above.

        The default formatter currently expresses floats and complex numbers with the
        pandas display precision unless using the ``precision`` argument here. The
        default formatter does not adjust the representation of missing values unless
        the ``na_rep`` argument is used.

        The ``subset`` argument defines which region to apply the formatting function
        to. If the ``formatter`` argument is given in dict form but does not include
        all columns within the subset then these columns will have the default formatter
        applied. Any columns in the formatter dict excluded from the subset will
        raise a ``KeyError``.

        When using a ``formatter`` string the dtypes must be compatible, otherwise a
        `ValueError` will be raised.

        Examples
        --------
        Using ``na_rep`` and ``precision`` with the default ``formatter``

        >>> df = pd.DataFrame([[np.nan, 1.0, 'A'], [2.0, np.nan, 3.0]])
        >>> df.style.format(na_rep='MISS', precision=3)
                0       1       2
        0    MISS   1.000       A
        1   2.000    MISS   3.000

        Using a ``formatter`` specification on consistent column dtypes

        >>> df.style.format('{:.2f}', na_rep='MISS', subset=[0,1])
                0      1          2
        0    MISS   1.00          A
        1    2.00   MISS   3.000000

        Using the default ``formatter`` for unspecified columns

        >>> df.style.format({0: '{:.2f}', 1: '£ {:.1f}'}, na_rep='MISS', precision=1)
                 0      1     2
        0    MISS   £ 1.0     A
        1    2.00    MISS   3.0

        Multiple ``na_rep`` or ``precision`` specifications under the default
        ``formatter``.

        >>> df.style.format(na_rep='MISS', precision=1, subset=[0])
        ...     .format(na_rep='PASS', precision=2, subset=[1, 2])
                0      1      2
        0    MISS   1.00      A
        1     2.0   PASS   3.00

        Using a callable ``formatter`` function.

        >>> func = lambda s: 'STRING' if isinstance(s, str) else 'FLOAT'
        >>> df.style.format({0: '{:.1f}', 2: func}, precision=4, na_rep='MISS')
                0        1        2
        0    MISS   1.0000   STRING
        1     2.0     MISS    FLOAT

        Using a ``formatter`` with HTML ``escape`` and ``na_rep``.

        >>> df = pd.DataFrame([['<div></div>', '"A&B"', None]])
        >>> s = df.style.format('<a href="a.com/{0}">{0}</a>', escape=True, na_rep="NA")
        >>> s.render()
        ...
        <td .. ><a href="a.com/&lt;div&gt;&lt;/div&gt;">&lt;div&gt;&lt;/div&gt;</a></td>
        <td .. ><a href="a.com/&#34;A&amp;B&#34;">&#34;A&amp;B&#34;</a></td>
        <td .. >NA</td>
        ...
        """
        if all(
            (
                formatter is None,
                subset is None,
                precision is None,
                na_rep is None,
                escape is False,
            )
        ):
            self._display_funcs.clear()
            return self  # clear the formatter / revert to default and avoid looping

        subset = slice(None) if subset is None else subset
        subset = non_reducing_slice(subset)
        data = self.data.loc[subset]

        columns = data.columns
        if not isinstance(formatter, dict):
            formatter = {col: formatter for col in columns}

        for col in columns:
            try:
                format_func = formatter[col]
            except KeyError:
                format_func = None
            format_func = _maybe_wrap_formatter(
                format_func, na_rep=na_rep, precision=precision, escape=escape
            )

            for row, value in data[[col]].itertuples():
                i, j = self.index.get_loc(row), self.columns.get_loc(col)
                self._display_funcs[(i, j)] = format_func

        return self


def _get_level_lengths(index, hidden_elements=None):
    """
    Given an index, find the level length for each element.

    Optional argument is a list of index positions which
    should not be visible.

    Result is a dictionary of (level, initial_position): span
    """
    if isinstance(index, MultiIndex):
        levels = index.format(sparsify=lib.no_default, adjoin=False)
    else:
        levels = index.format()

    if hidden_elements is None:
        hidden_elements = []

    lengths = {}
    if index.nlevels == 1:
        for i, value in enumerate(levels):
            if i not in hidden_elements:
                lengths[(0, i)] = 1
        return lengths

    for i, lvl in enumerate(levels):
        for j, row in enumerate(lvl):
            if not get_option("display.multi_sparse"):
                lengths[(i, j)] = 1
            elif (row is not lib.no_default) and (j not in hidden_elements):
                last_label = j
                lengths[(i, last_label)] = 1
            elif row is not lib.no_default:
                # even if its hidden, keep track of it in case
                # length >1 and later elements are visible
                last_label = j
                lengths[(i, last_label)] = 0
            elif j not in hidden_elements:
                lengths[(i, last_label)] += 1

    non_zero_lengths = {
        element: length for element, length in lengths.items() if length >= 1
    }

    return non_zero_lengths


def _is_visible(idx_row, idx_col, lengths) -> bool:
    """
    Index -> {(idx_row, idx_col): bool}).
    """
    return (idx_col, idx_row) in lengths


def _format_table_styles(styles: CSSStyles) -> CSSStyles:
    """
    looks for multiple CSS selectors and separates them:
    [{'selector': 'td, th', 'props': 'a:v;'}]
        ---> [{'selector': 'td', 'props': 'a:v;'},
              {'selector': 'th', 'props': 'a:v;'}]
    """
    return [
        item
        for sublist in [
            [  # this is a CSSDict when TypedDict is available to avoid cast.
                {"selector": x, "props": style["props"]}
                for x in cast(str, style["selector"]).split(",")
            ]
            for style in styles
        ]
        for item in sublist
    ]


def _default_formatter(x: Any, precision: int) -> Any:
    """
    Format the display of a value

    Parameters
    ----------
    x : Any
        Input variable to be formatted
    precision : Int
        Floating point precision used if ``x`` is float or complex.

    Returns
    -------
    value : Any
        Matches input type, or string if input is float or complex.
    """
    if isinstance(x, (float, complex)):
        return f"{x:.{precision}f}"
    return x


def _maybe_wrap_formatter(
    formatter: BaseFormatter | None = None,
    na_rep: str | None = None,
    precision: int | None = None,
    escape: bool = False,
) -> Callable:
    """
    Allows formatters to be expressed as str, callable or None, where None returns
    a default formatting function. wraps with na_rep, and precision where they are
    available.
    """
    if isinstance(formatter, str):
        formatter_func = lambda x: formatter.format(x)
    elif callable(formatter):
        formatter_func = formatter
    elif formatter is None:
        precision = get_option("display.precision") if precision is None else precision
        formatter_func = partial(_default_formatter, precision=precision)
    else:
        raise TypeError(f"'formatter' expected str or callable, got {type(formatter)}")

    def _str_escape(x, escape: bool):
        """if escaping: only use on str, else return input"""
        if escape and isinstance(x, str):
            return escape_func(x)
        else:
            return x

    display_func = lambda x: formatter_func(partial(_str_escape, escape=escape)(x))

    if na_rep is None:
        return display_func
    else:
        return lambda x: na_rep if isna(x) else display_func(x)


def non_reducing_slice(slice_):
    """
    Ensure that a slice doesn't reduce to a Series or Scalar.

    Any user-passed `subset` should have this called on it
    to make sure we're always working with DataFrames.
    """
    # default to column slice, like DataFrame
    # ['A', 'B'] -> IndexSlices[:, ['A', 'B']]
    kinds = (ABCSeries, np.ndarray, Index, list, str)
    if isinstance(slice_, kinds):
        slice_ = IndexSlice[:, slice_]

    def pred(part) -> bool:
        """
        Returns
        -------
        bool
            True if slice does *not* reduce,
            False if `part` is a tuple.
        """
        # true when slice does *not* reduce, False when part is a tuple,
        # i.e. MultiIndex slice
        if isinstance(part, tuple):
            # GH#39421 check for sub-slice:
            return any((isinstance(s, slice) or is_list_like(s)) for s in part)
        else:
            return isinstance(part, slice) or is_list_like(part)

    if not is_list_like(slice_):
        if not isinstance(slice_, slice):
            # a 1-d slice, like df.loc[1]
            slice_ = [[slice_]]
        else:
            # slice(a, b, c)
            slice_ = [slice_]  # to tuplize later
    else:
        slice_ = [part if pred(part) else [part] for part in slice_]
    return tuple(slice_)


def maybe_convert_css_to_tuples(style: CSSProperties) -> CSSList:
    """
    Convert css-string to sequence of tuples format if needed.
    'color:red; border:1px solid black;' -> [('color', 'red'),
                                             ('border','1px solid red')]
    """
    if isinstance(style, str):
        s = style.split(";")
        try:
            return [
                (x.split(":")[0].strip(), x.split(":")[1].strip())
                for x in s
                if x.strip() != ""
            ]
        except IndexError:
            raise ValueError(
                "Styles supplied as string must follow CSS rule formats, "
                f"for example 'attr: val;'. '{style}' was given."
            )
    return style


class Tooltips:
    """
    An extension to ``Styler`` that allows for and manipulates tooltips on hover
    of ``<td>`` cells in the HTML result.

    Parameters
    ----------
    css_name: str, default "pd-t"
        Name of the CSS class that controls visualisation of tooltips.
    css_props: list-like, default; see Notes
        List of (attr, value) tuples defining properties of the CSS class.
    tooltips: DataFrame, default empty
        DataFrame of strings aligned with underlying Styler data for tooltip
        display.

    Notes
    -----
    The default properties for the tooltip CSS class are:

        - visibility: hidden
        - position: absolute
        - z-index: 1
        - background-color: black
        - color: white
        - transform: translate(-20px, -20px)

    Hidden visibility is a key prerequisite to the hover functionality, and should
    always be included in any manual properties specification.
    """

    def __init__(
        self,
        css_props: CSSProperties = [
            ("visibility", "hidden"),
            ("position", "absolute"),
            ("z-index", 1),
            ("background-color", "black"),
            ("color", "white"),
            ("transform", "translate(-20px, -20px)"),
        ],
        css_name: str = "pd-t",
        tooltips: DataFrame = DataFrame(),
    ):
        self.class_name = css_name
        self.class_properties = css_props
        self.tt_data = tooltips
        self.table_styles: CSSStyles = []

    @property
    def _class_styles(self):
        """
        Combine the ``_Tooltips`` CSS class name and CSS properties to the format
        required to extend the underlying ``Styler`` `table_styles` to allow
        tooltips to render in HTML.

        Returns
        -------
        styles : List
        """
        return [
            {
                "selector": f".{self.class_name}",
                "props": maybe_convert_css_to_tuples(self.class_properties),
            }
        ]

    def _pseudo_css(self, uuid: str, name: str, row: int, col: int, text: str):
        """
        For every table data-cell that has a valid tooltip (not None, NaN or
        empty string) must create two pseudo CSS entries for the specific
        <td> element id which are added to overall table styles:
        an on hover visibility change and a content change
        dependent upon the user's chosen display string.

        For example:
            [{"selector": "T__row1_col1:hover .pd-t",
             "props": [("visibility", "visible")]},
            {"selector": "T__row1_col1 .pd-t::after",
             "props": [("content", "Some Valid Text String")]}]

        Parameters
        ----------
        uuid: str
            The uuid of the Styler instance
        name: str
            The css-name of the class used for styling tooltips
        row : int
            The row index of the specified tooltip string data
        col : int
            The col index of the specified tooltip string data
        text : str
            The textual content of the tooltip to be displayed in HTML.

        Returns
        -------
        pseudo_css : List
        """
        selector_id = "#T_" + uuid + "row" + str(row) + "_col" + str(col)
        return [
            {
                "selector": selector_id + f":hover .{name}",
                "props": [("visibility", "visible")],
            },
            {
                "selector": selector_id + f" .{name}::after",
                "props": [("content", f'"{text}"')],
            },
        ]

    def _translate(self, styler_data: FrameOrSeriesUnion, uuid: str, d: dict):
        """
        Mutate the render dictionary to allow for tooltips:

        - Add ``<span>`` HTML element to each data cells ``display_value``. Ignores
          headers.
        - Add table level CSS styles to control pseudo classes.

        Parameters
        ----------
        styler_data : DataFrame
            Underlying ``Styler`` DataFrame used for reindexing.
        uuid : str
            The underlying ``Styler`` uuid for CSS id.
        d : dict
            The dictionary prior to final render

        Returns
        -------
        render_dict : Dict
        """
        self.tt_data = self.tt_data.reindex_like(styler_data)

        if self.tt_data.empty:
            return d

        name = self.class_name

        mask = (self.tt_data.isna()) | (self.tt_data.eq(""))  # empty string = no ttip
        self.table_styles = [
            style
            for sublist in [
                self._pseudo_css(uuid, name, i, j, str(self.tt_data.iloc[i, j]))
                for i in range(len(self.tt_data.index))
                for j in range(len(self.tt_data.columns))
                if not mask.iloc[i, j]
            ]
            for style in sublist
        ]

        if self.table_styles:
            # add span class to every cell only if at least 1 non-empty tooltip
            for row in d["body"]:
                for item in row:
                    if item["type"] == "td":
                        item["display_value"] = (
                            str(item["display_value"])
                            + f'<span class="{self.class_name}"></span>'
                        )
            d["table_styles"].extend(self._class_styles)
            d["table_styles"].extend(self.table_styles)

        return d
