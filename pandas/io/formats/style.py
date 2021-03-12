"""
Module for applying conditional formatting to DataFrames and Series.
"""
from __future__ import annotations

from collections import defaultdict
from contextlib import contextmanager
import copy
from functools import partial
from typing import (
    Any,
    Callable,
    DefaultDict,
    Dict,
    Hashable,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)
from uuid import uuid4
import warnings

import numpy as np

from pandas._config import get_option

from pandas._libs import lib
from pandas._typing import (
    Axis,
    FrameOrSeries,
    FrameOrSeriesUnion,
    IndexLabel,
)
from pandas.compat._optional import import_optional_dependency
from pandas.util._decorators import doc

from pandas.core.dtypes.generic import ABCSeries

import pandas as pd
from pandas.api.types import is_list_like
from pandas.core import generic
import pandas.core.common as com
from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame
from pandas.core.indexes.api import Index

jinja2 = import_optional_dependency("jinja2", extra="DataFrame.style requires jinja2.")

BaseFormatter = Union[str, Callable]
ExtFormatter = Union[BaseFormatter, Dict[Any, Optional[BaseFormatter]]]
CSSPair = Tuple[str, Union[str, int, float]]
CSSList = List[CSSPair]
CSSProperties = Union[str, CSSList]
CSSStyles = List[Dict[str, CSSProperties]]  # = List[CSSDict]
# class CSSDict(TypedDict):  # available when TypedDict is valid in pandas
#     selector: str
#     props: CSSProperties

try:
    from matplotlib import colors
    import matplotlib.pyplot as plt

    has_mpl = True
except ImportError:
    has_mpl = False
    no_mpl_message = "{0} requires matplotlib."


@contextmanager
def _mpl(func: Callable):
    if has_mpl:
        yield plt, colors
    else:
        raise ImportError(no_mpl_message.format(func.__name__))


class Styler:
    """
    Helps style a DataFrame or Series according to the data with HTML and CSS.

    Parameters
    ----------
    data : Series or DataFrame
        Data to be styled - either a Series or DataFrame.
    precision : int
        Precision to round floats to, defaults to pd.options.display.precision.
    table_styles : list-like, default None
        List of {selector: (attr, value)} dicts; see Notes.
    uuid : str, default None
        A unique identifier to avoid CSS collisions; generated automatically.
    caption : str, default None
        Caption to attach to the table.
    table_attributes : str, default None
        Items that show up in the opening ``<table>`` tag
        in addition to automatic (by default) id.
    cell_ids : bool, default True
        If True, each cell will have an ``id`` attribute in their HTML tag.
        The ``id`` takes the form ``T_<uuid>_row<num_row>_col<num_col>``
        where ``<uuid>`` is the unique identifier, ``<num_row>`` is the row
        number and ``<num_col>`` is the column number.
    na_rep : str, optional
        Representation for missing values.
        If ``na_rep`` is None, no special formatting is applied.

        .. versionadded:: 1.0.0

    uuid_len : int, default 5
        If ``uuid`` is not specified, the length of the ``uuid`` to randomly generate
        expressed in hex characters, in range [0, 32].

        .. versionadded:: 1.2.0

    Attributes
    ----------
    env : Jinja2 jinja2.Environment
    template : Jinja2 Template
    loader : Jinja2 Loader

    See Also
    --------
    DataFrame.style : Return a Styler object containing methods for building
        a styled HTML representation for the DataFrame.

    Notes
    -----
    Most styling will be done by passing style functions into
    ``Styler.apply`` or ``Styler.applymap``. Style functions should
    return values with strings containing CSS ``'attr: value'`` that will
    be applied to the indicated cells.

    If using in the Jupyter notebook, Styler has defined a ``_repr_html_``
    to automatically render itself. Otherwise call Styler.render to get
    the generated HTML.

    CSS classes are attached to the generated HTML

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
    """

    loader = jinja2.PackageLoader("pandas", "io/formats/templates")
    env = jinja2.Environment(loader=loader, trim_blocks=True)
    template = env.get_template("html.tpl")

    def __init__(
        self,
        data: FrameOrSeriesUnion,
        precision: Optional[int] = None,
        table_styles: Optional[CSSStyles] = None,
        uuid: Optional[str] = None,
        caption: Optional[str] = None,
        table_attributes: Optional[str] = None,
        cell_ids: bool = True,
        na_rep: Optional[str] = None,
        uuid_len: int = 5,
    ):
        # validate ordered args
        if isinstance(data, pd.Series):
            data = data.to_frame()
        if not isinstance(data, DataFrame):
            raise TypeError("``data`` must be a Series or DataFrame")
        if not data.index.is_unique or not data.columns.is_unique:
            raise ValueError("style is not supported for non-unique indices.")
        self.data: DataFrame = data
        self.index: pd.Index = data.index
        self.columns: pd.Index = data.columns
        self.table_styles = table_styles
        if not isinstance(uuid_len, int) or not uuid_len >= 0:
            raise TypeError("``uuid_len`` must be an integer in range [0, 32].")
        self.uuid_len = min(32, uuid_len)
        self.uuid = (uuid or uuid4().hex[: self.uuid_len]) + "_"
        self.caption = caption
        self.table_attributes = table_attributes
        self.cell_ids = cell_ids

        # assign additional default vars
        self.hidden_index: bool = False
        self.hidden_columns: Sequence[int] = []
        self.ctx: DefaultDict[Tuple[int, int], CSSList] = defaultdict(list)
        self.cell_context: Dict[str, Any] = {}
        self._todo: List[Tuple[Callable, Tuple, Dict]] = []
        self.tooltips: Optional[_Tooltips] = None
        self._display_funcs: DefaultDict[  # maps (row, col) -> formatting function
            Tuple[int, int], Callable[[Any], str]
        ] = defaultdict(lambda: partial(_default_formatter, precision=None))
        self.precision = precision  # can be removed on set_precision depr cycle
        self.na_rep = na_rep  # can be removed on set_na_rep depr cycle
        self.format(formatter=None, precision=precision, na_rep=na_rep)

    def _repr_html_(self) -> str:
        """
        Hooks into Jupyter notebook rich display system.
        """
        return self.render()

    def set_tooltips(
        self,
        ttips: DataFrame,
        props: Optional[CSSProperties] = None,
        css_class: Optional[str] = None,
    ) -> Styler:
        """
        Add string based tooltips that will appear in the `Styler` HTML result. These
        tooltips are applicable only to`<td>` elements.

        .. versionadded:: 1.3.0

        Parameters
        ----------
        ttips : DataFrame
            DataFrame containing strings that will be translated to tooltips, mapped
            by identical column and index values that must exist on the underlying
            `Styler` data. None, NaN values, and empty strings will be ignored and
            not affect the rendered HTML.
        props : list-like or str, optional
            List of (attr, value) tuples or a valid CSS string. If ``None`` adopts
            the internal default values described in notes.
        css_class : str, optional
            Name of the tooltip class used in CSS, should conform to HTML standards.
            Only useful if integrating tooltips with external CSS. If ``None`` uses the
            internal default value 'pd-t'.

        Returns
        -------
        self : Styler

        Notes
        -----
        Tooltips are created by adding `<span class="pd-t"></span>` to each data cell
        and then manipulating the table level CSS to attach pseudo hover and pseudo
        after selectors to produce the required the results.

        The default properties for the tooltip CSS class are:

        - visibility: hidden
        - position: absolute
        - z-index: 1
        - background-color: black
        - color: white
        - transform: translate(-20px, -20px)

        The property 'visibility: hidden;' is a key prerequisite to the hover
        functionality, and should always be included in any manual properties
        specification, using the ``props`` argument.

        Tooltips are not designed to be efficient, and can add large amounts of
        additional HTML for larger tables, since they also require that ``cell_ids``
        is forced to `True`.

        Examples
        --------
        Basic application

        >>> df = pd.DataFrame(data=[[0, 1], [2, 3]])
        >>> ttips = pd.DataFrame(
        ...    data=[["Min", ""], [np.nan, "Max"]], columns=df.columns, index=df.index
        ... )
        >>> s = df.style.set_tooltips(ttips).render()

        Optionally controlling the tooltip visual display

        >>> df.style.set_tooltips(ttips, css_class='tt-add', props=[
        ...     ('visibility', 'hidden'),
        ...     ('position', 'absolute'),
        ...     ('z-index', 1)])
        >>> df.style.set_tooltips(ttips, css_class='tt-add',
        ...     props='visibility:hidden; position:absolute; z-index:1;')
        """
        if not self.cell_ids:
            # tooltips not optimised for individual cell check. requires reasonable
            # redesign and more extensive code for a feature that might be rarely used.
            raise NotImplementedError(
                "Tooltips can only render with 'cell_ids' is True."
            )
        if self.tooltips is None:  # create a default instance if necessary
            self.tooltips = _Tooltips()
        self.tooltips.tt_data = ttips
        if props:
            self.tooltips.class_properties = props
        if css_class:
            self.tooltips.class_name = css_class

        return self

    @doc(
        NDFrame.to_excel,
        klass="Styler",
        storage_options=generic._shared_docs["storage_options"],
    )
    def to_excel(
        self,
        excel_writer,
        sheet_name: str = "Sheet1",
        na_rep: str = "",
        float_format: Optional[str] = None,
        columns: Optional[Sequence[Hashable]] = None,
        header: Union[Sequence[Hashable], bool] = True,
        index: bool = True,
        index_label: Optional[IndexLabel] = None,
        startrow: int = 0,
        startcol: int = 0,
        engine: Optional[str] = None,
        merge_cells: bool = True,
        encoding: Optional[str] = None,
        inf_rep: str = "inf",
        verbose: bool = True,
        freeze_panes: Optional[Tuple[int, int]] = None,
    ) -> None:

        from pandas.io.formats.excel import ExcelFormatter

        formatter = ExcelFormatter(
            self,
            na_rep=na_rep,
            cols=columns,
            header=header,
            float_format=float_format,
            index=index,
            index_label=index_label,
            merge_cells=merge_cells,
            inf_rep=inf_rep,
        )
        formatter.write(
            excel_writer,
            sheet_name=sheet_name,
            startrow=startrow,
            startcol=startcol,
            freeze_panes=freeze_panes,
            engine=engine,
        )

    def _translate(self):
        """
        Convert the DataFrame in `self.data` and the attrs from `_build_styles`
        into a dictionary of {head, body, uuid, cellstyle}.
        """
        table_styles = self.table_styles or []
        caption = self.caption
        ctx = self.ctx
        hidden_index = self.hidden_index
        hidden_columns = self.hidden_columns
        uuid = self.uuid
        ROW_HEADING_CLASS = "row_heading"
        COL_HEADING_CLASS = "col_heading"
        INDEX_NAME_CLASS = "index_name"

        DATA_CLASS = "data"
        BLANK_CLASS = "blank"
        BLANK_VALUE = ""

        # for sparsifying a MultiIndex
        idx_lengths = _get_level_lengths(self.index)
        col_lengths = _get_level_lengths(self.columns, hidden_columns)

        cell_context = self.cell_context

        n_rlvls = self.data.index.nlevels
        n_clvls = self.data.columns.nlevels
        rlabels = self.data.index.tolist()
        clabels = self.data.columns.tolist()

        if n_rlvls == 1:
            rlabels = [[x] for x in rlabels]
        if n_clvls == 1:
            clabels = [[x] for x in clabels]
        clabels = list(zip(*clabels))

        cellstyle_map: DefaultDict[Tuple[CSSPair, ...], List[str]] = defaultdict(list)

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
                    cs = [
                        COL_HEADING_CLASS,
                        f"level{r}",
                        f"col{c}",
                    ]
                    cs.extend(
                        cell_context.get("col_headings", {}).get(r, {}).get(c, [])
                    )
                    es = {
                        "type": "th",
                        "value": value,
                        "display_value": value,
                        "class": " ".join(cs),
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
                cs = [DATA_CLASS, f"row{r}", f"col{c}"]
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
                    row_dict["id"] = "_".join(cs[1:])
                    props.extend(ctx[r, c])

                # add custom classes from cell context
                cs.extend(cell_context.get("data", {}).get(r, {}).get(c, []))
                row_dict["class"] = " ".join(cs)

                row_es.append(row_dict)
                if props:  # (), [] won't be in cellstyle_map, cellstyle respectively
                    cellstyle_map[tuple(props)].append(f"row{r}_col{c}")
            body.append(row_es)

        cellstyle: List[Dict[str, Union[CSSList, List[str]]]] = [
            {"props": list(props), "selectors": selectors}
            for props, selectors in cellstyle_map.items()
        ]

        table_attr = self.table_attributes
        use_mathjax = get_option("display.html.use_mathjax")
        if not use_mathjax:
            table_attr = table_attr or ""
            if 'class="' in table_attr:
                table_attr = table_attr.replace('class="', 'class="tex2jax_ignore ')
            else:
                table_attr += ' class="tex2jax_ignore"'

        d = {
            "head": head,
            "cellstyle": cellstyle,
            "body": body,
            "uuid": uuid,
            "table_styles": _format_table_styles(table_styles),
            "caption": caption,
            "table_attributes": table_attr,
        }
        if self.tooltips:
            d = self.tooltips._translate(self.data, self.uuid, d)

        return d

    def format(
        self,
        formatter: Optional[ExtFormatter] = None,
        subset: Optional[Union[slice, Sequence[Any]]] = None,
        na_rep: Optional[str] = None,
        precision: Optional[int] = None,
    ) -> Styler:
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

        Using a format specification on consistent column dtypes

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

        Using a callable formatting function

        >>> func = lambda s: 'STRING' if isinstance(s, str) else 'FLOAT'
        >>> df.style.format({0: '{:.1f}', 2: func}, precision=4, na_rep='MISS')
                0        1        2
        0    MISS   1.0000   STRING
        1     2.0     MISS    FLOAT
        """
        if all((formatter is None, subset is None, precision is None, na_rep is None)):
            self._display_funcs.clear()
            return self  # clear the formatter / revert to default and avoid looping

        subset = slice(None) if subset is None else subset
        subset = _non_reducing_slice(subset)
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
                format_func, na_rep=na_rep, precision=precision
            )

            for row, value in data[[col]].itertuples():
                i, j = self.index.get_loc(row), self.columns.get_loc(col)
                self._display_funcs[(i, j)] = format_func

        return self

    def set_td_classes(self, classes: DataFrame) -> Styler:
        """
        Add string based CSS class names to data cells that will appear within the
        `Styler` HTML result. These classes are added within specified `<td>` elements.

        Parameters
        ----------
        classes : DataFrame
            DataFrame containing strings that will be translated to CSS classes,
            mapped by identical column and index values that must exist on the
            underlying `Styler` data. None, NaN values, and empty strings will
            be ignored and not affect the rendered HTML.

        Returns
        -------
        self : Styler

        Examples
        --------
        >>> df = pd.DataFrame(data=[[1, 2, 3], [4, 5, 6]], columns=["A", "B", "C"])
        >>> classes = pd.DataFrame([
        ...     ["min-val red", "", "blue"],
        ...     ["red", None, "blue max-val"]
        ... ], index=df.index, columns=df.columns)
        >>> df.style.set_td_classes(classes)

        Using `MultiIndex` columns and a `classes` `DataFrame` as a subset of the
        underlying,

        >>> df = pd.DataFrame([[1,2],[3,4]], index=["a", "b"],
        ...     columns=[["level0", "level0"], ["level1a", "level1b"]])
        >>> classes = pd.DataFrame(["min-val"], index=["a"],
        ...     columns=[["level0"],["level1a"]])
        >>> df.style.set_td_classes(classes)

        Form of the output with new additional css classes,

        >>> df = pd.DataFrame([[1]])
        >>> css = pd.DataFrame(["other-class"])
        >>> s = Styler(df, uuid="_", cell_ids=False).set_td_classes(css)
        >>> s.hide_index().render()
        '<style  type="text/css" ></style>'
        '<table id="T__" >'
        '  <thead>'
        '    <tr><th class="col_heading level0 col0" >0</th></tr>'
        '  </thead>'
        '  <tbody>'
        '    <tr><td  class="data row0 col0 other-class" >1</td></tr>'
        '  </tbody>'
        '</table>'
        """
        classes = classes.reindex_like(self.data)

        mask = (classes.isna()) | (classes.eq(""))
        self.cell_context["data"] = {
            r: {
                c: [str(classes.iloc[r, c])]
                for c, cn in enumerate(classes.columns)
                if not mask.iloc[r, c]
            }
            for r, rn in enumerate(classes.index)
        }

        return self

    def render(self, **kwargs) -> str:
        """
        Render the built up styles to HTML.

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
        ``Styler`` objects have defined the ``_repr_html_`` method
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

    def _update_ctx(self, attrs: DataFrame) -> None:
        """
        Update the state of the Styler for data cells.

        Collects a mapping of {index_label: [('<property>', '<value>'), ..]}.

        Parameters
        ----------
        attrs : DataFrame
            should contain strings of '<property>: <value>;<prop2>: <val2>'
            Whitespace shouldn't matter and the final trailing ';' shouldn't
            matter.
        """
        for cn in attrs.columns:
            for rn, c in attrs[[cn]].itertuples():
                if not c:
                    continue
                css_list = _maybe_convert_css_to_tuples(c)
                i, j = self.index.get_loc(rn), self.columns.get_loc(cn)
                self.ctx[(i, j)].extend(css_list)

    def _copy(self, deepcopy: bool = False) -> Styler:
        styler = Styler(
            self.data,
            precision=self.precision,
            caption=self.caption,
            table_attributes=self.table_attributes,
            cell_ids=self.cell_ids,
            na_rep=self.na_rep,
        )

        styler.uuid = self.uuid
        styler.hidden_index = self.hidden_index

        if deepcopy:
            styler.ctx = copy.deepcopy(self.ctx)
            styler._todo = copy.deepcopy(self._todo)
            styler.table_styles = copy.deepcopy(self.table_styles)
            styler.hidden_columns = copy.copy(self.hidden_columns)
            styler.cell_context = copy.deepcopy(self.cell_context)
            styler.tooltips = copy.deepcopy(self.tooltips)
        else:
            styler.ctx = self.ctx
            styler._todo = self._todo
            styler.table_styles = self.table_styles
            styler.hidden_columns = self.hidden_columns
            styler.cell_context = self.cell_context
            styler.tooltips = self.tooltips

        return styler

    def __copy__(self) -> Styler:
        """
        Deep copy by default.
        """
        return self._copy(deepcopy=False)

    def __deepcopy__(self, memo) -> Styler:
        return self._copy(deepcopy=True)

    def clear(self) -> None:
        """
        Reset the styler, removing any previously applied styles.

        Returns None.
        """
        self.ctx.clear()
        self.tooltips = None
        self.cell_context = {}
        self._todo = []

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

    def _apply(
        self,
        func: Callable[..., Styler],
        axis: Optional[Axis] = 0,
        subset=None,
        **kwargs,
    ) -> Styler:
        subset = slice(None) if subset is None else subset
        subset = _non_reducing_slice(subset)
        data = self.data.loc[subset]
        if axis is not None:
            result = data.apply(func, axis=axis, result_type="expand", **kwargs)
            result.columns = data.columns
        else:
            result = func(data, **kwargs)
            if not isinstance(result, pd.DataFrame):
                if not isinstance(result, np.ndarray):
                    raise TypeError(
                        f"Function {repr(func)} must return a DataFrame or ndarray "
                        f"when passed to `Styler.apply` with axis=None"
                    )
                if not (data.shape == result.shape):
                    raise ValueError(
                        f"Function {repr(func)} returned ndarray with wrong shape.\n"
                        f"Result has shape: {result.shape}\n"
                        f"Expected shape: {data.shape}"
                    )
                result = DataFrame(result, index=data.index, columns=data.columns)
            elif not (
                result.index.equals(data.index) and result.columns.equals(data.columns)
            ):
                raise ValueError(
                    f"Result of {repr(func)} must have identical "
                    f"index and columns as the input"
                )

        if result.shape != data.shape:
            raise ValueError(
                f"Function {repr(func)} returned the wrong shape.\n"
                f"Result has shape: {result.shape}\n"
                f"Expected shape:   {data.shape}"
            )
        self._update_ctx(result)
        return self

    def apply(
        self,
        func: Callable[..., Styler],
        axis: Optional[Axis] = 0,
        subset=None,
        **kwargs,
    ) -> Styler:
        """
        Apply a CSS-styling function column-wise, row-wise, or table-wise.

        Updates the HTML representation with the result.

        Parameters
        ----------
        func : function
            ``func`` should take a Series or DataFrame (depending
            on ``axis``), and return an object with the same shape.
            Must return a DataFrame with identical index and
            column labels or an ndarray with same shape as input when ``axis=None``.

            .. versionchanged:: 1.3.0

        axis : {0 or 'index', 1 or 'columns', None}, default 0
            Apply to each column (``axis=0`` or ``'index'``), to each row
            (``axis=1`` or ``'columns'``), or to the entire DataFrame at once
            with ``axis=None``.
        subset : IndexSlice
            A valid indexer to limit ``data`` to *before* applying the
            function. Consider using a pandas.IndexSlice.
        **kwargs : dict
            Pass along to ``func``.

        Returns
        -------
        self : Styler

        Notes
        -----
        The output of ``func`` should be elements having CSS style as string or,
        if nothing is to be applied to that element, an empty string or ``None``.
        The output shape must match the input, i.e. if
        ``x`` is the input row, column, or table (depending on ``axis``),
        then ``func(x).shape == x.shape`` should be ``True``.

        This is similar to ``DataFrame.apply``, except that ``axis=None``
        applies the function to the entire DataFrame at once,
        rather than column-wise or row-wise.

        Examples
        --------
        >>> def highlight_max(x, color):
        ...     return np.where(x == np.nanmax(x.to_numpy()), f"color: {color};", None)
        >>> df = pd.DataFrame(np.random.randn(5, 2))
        >>> df.style.apply(highlight_max, color='red')
        >>> df.style.apply(highlight_max, color='blue', axis=1)
        >>> df.style.apply(highlight_max, color='green', axis=None)
        """
        self._todo.append(
            (lambda instance: getattr(instance, "_apply"), (func, axis, subset), kwargs)
        )
        return self

    def _applymap(self, func: Callable, subset=None, **kwargs) -> Styler:
        func = partial(func, **kwargs)  # applymap doesn't take kwargs?
        if subset is None:
            subset = pd.IndexSlice[:]
        subset = _non_reducing_slice(subset)
        result = self.data.loc[subset].applymap(func)
        self._update_ctx(result)
        return self

    def applymap(self, func: Callable, subset=None, **kwargs) -> Styler:
        """
        Apply a CSS-styling function elementwise.

        Updates the HTML representation with the result.

        Parameters
        ----------
        func : function
            ``func`` should take a scalar and return a scalar.
        subset : IndexSlice
            A valid indexer to limit ``data`` to *before* applying the
            function. Consider using a pandas.IndexSlice.
        **kwargs : dict
            Pass along to ``func``.

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.where: Updates the HTML representation with a style which is
            selected in accordance with the return value of a function.

        Notes
        -----
        The output of ``func`` should be a CSS style as string or, if nothing is to be
        applied, an empty string or ``None``.

        Examples
        --------
        >>> def color_negative(v, color):
        ...     return f"color: {color};" if v < 0 else None
        >>> df = pd.DataFrame(np.random.randn(5, 2))
        >>> df.style.applymap(color_negative, color='red')
        """
        self._todo.append(
            (lambda instance: getattr(instance, "_applymap"), (func, subset), kwargs)
        )
        return self

    def where(
        self,
        cond: Callable,
        value: str,
        other: Optional[str] = None,
        subset=None,
        **kwargs,
    ) -> Styler:
        """
        Apply a function elementwise.

        Updates the HTML representation with a style which is
        selected in accordance with the return value of a function.

        Parameters
        ----------
        cond : callable
            ``cond`` should take a scalar and return a boolean.
        value : str
            Applied when ``cond`` returns true.
        other : str
            Applied when ``cond`` returns false.
        subset : IndexSlice
            A valid indexer to limit ``data`` to *before* applying the
            function. Consider using a pandas.IndexSlice.
        **kwargs : dict
            Pass along to ``cond``.

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.applymap: Updates the HTML representation with the result.
        """
        if other is None:
            other = ""

        return self.applymap(
            lambda val: value if cond(val) else other, subset=subset, **kwargs
        )

    def set_precision(self, precision: int) -> Styler:
        """
        Set the precision used to display values.

        .. deprecated:: 1.3.0

        Parameters
        ----------
        precision : int

        Returns
        -------
        self : Styler

        Notes
        -----
        This method is deprecated see `Styler.format`.
        """
        warnings.warn(
            "this method is deprecated in favour of `Styler.format(precision=..)`",
            FutureWarning,
            stacklevel=2,
        )
        self.precision = precision
        return self.format(precision=precision, na_rep=self.na_rep)

    def set_table_attributes(self, attributes: str) -> Styler:
        """
        Set the table attributes.

        These are the items that show up in the opening ``<table>`` tag
        in addition to automatic (by default) id.

        Parameters
        ----------
        attributes : str

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

    def export(self) -> List[Tuple[Callable, Tuple, Dict]]:
        """
        Export the styles to applied to the current Styler.

        Can be applied to a second style with ``Styler.use``.

        Returns
        -------
        styles : list

        See Also
        --------
        Styler.use: Set the styles on the current Styler.
        """
        return self._todo

    def use(self, styles: List[Tuple[Callable, Tuple, Dict]]) -> Styler:
        """
        Set the styles on the current Styler.

        Possibly uses styles from ``Styler.export``.

        Parameters
        ----------
        styles : list
            List of style functions.

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.export : Export the styles to applied to the current Styler.
        """
        self._todo.extend(styles)
        return self

    def set_uuid(self, uuid: str) -> Styler:
        """
        Set the uuid for a Styler.

        Parameters
        ----------
        uuid : str

        Returns
        -------
        self : Styler
        """
        self.uuid = uuid
        return self

    def set_caption(self, caption: str) -> Styler:
        """
        Set the caption on a Styler.

        Parameters
        ----------
        caption : str

        Returns
        -------
        self : Styler
        """
        self.caption = caption
        return self

    def set_table_styles(
        self,
        table_styles: Union[Dict[Any, CSSStyles], CSSStyles],
        axis: int = 0,
        overwrite: bool = True,
    ) -> Styler:
        """
        Set the table styles on a Styler.

        These are placed in a ``<style>`` tag before the generated HTML table.

        This function can be used to style the entire table, columns, rows or
        specific HTML selectors.

        Parameters
        ----------
        table_styles : list or dict
            If supplying a list, each individual table_style should be a
            dictionary with ``selector`` and ``props`` keys. ``selector``
            should be a CSS selector that the style will be applied to
            (automatically prefixed by the table's UUID) and ``props``
            should be a list of tuples with ``(attribute, value)``.
            If supplying a dict, the dict keys should correspond to
            column names or index values, depending upon the specified
            `axis` argument. These will be mapped to row or col CSS
            selectors. MultiIndex values as dict keys should be
            in their respective tuple form. The dict values should be
            a list as specified in the form with CSS selectors and
            props that will be applied to the specified row or column.

            .. versionchanged:: 1.2.0

        axis : {0 or 'index', 1 or 'columns', None}, default 0
            Apply to each column (``axis=0`` or ``'index'``), to each row
            (``axis=1`` or ``'columns'``). Only used if `table_styles` is
            dict.

            .. versionadded:: 1.2.0

        overwrite : boolean, default True
            Styles are replaced if `True`, or extended if `False`. CSS
            rules are preserved so most recent styles set will dominate
            if selectors intersect.

            .. versionadded:: 1.2.0

        Returns
        -------
        self : Styler

        Examples
        --------
        >>> df = pd.DataFrame(np.random.randn(10, 4),
        ...                   columns=['A', 'B', 'C', 'D'])
        >>> df.style.set_table_styles(
        ...     [{'selector': 'tr:hover',
        ...       'props': [('background-color', 'yellow')]}]
        ... )

        Or with CSS strings

        >>> df.style.set_table_styles(
        ...     [{'selector': 'tr:hover',
        ...       'props': 'background-color: yellow; font-size: 1em;']}]
        ... )

        Adding column styling by name

        >>> df.style.set_table_styles({
        ...     'A': [{'selector': '',
        ...            'props': [('color', 'red')]}],
        ...     'B': [{'selector': 'td',
        ...            'props': 'color: blue;']}]
        ... }, overwrite=False)

        Adding row styling

        >>> df.style.set_table_styles({
        ...     0: [{'selector': 'td:hover',
        ...          'props': [('font-size', '25px')]}]
        ... }, axis=1, overwrite=False)
        """
        if isinstance(table_styles, dict):
            if axis in [0, "index"]:
                obj, idf = self.data.columns, ".col"
            else:
                obj, idf = self.data.index, ".row"

            table_styles = [
                {
                    "selector": str(s["selector"]) + idf + str(obj.get_loc(key)),
                    "props": _maybe_convert_css_to_tuples(s["props"]),
                }
                for key, styles in table_styles.items()
                for s in styles
            ]
        else:
            table_styles = [
                {
                    "selector": s["selector"],
                    "props": _maybe_convert_css_to_tuples(s["props"]),
                }
                for s in table_styles
            ]

        if not overwrite and self.table_styles is not None:
            self.table_styles.extend(table_styles)
        else:
            self.table_styles = table_styles
        return self

    def set_na_rep(self, na_rep: str) -> Styler:
        """
        Set the missing data representation on a Styler.

        .. versionadded:: 1.0.0

        .. deprecated:: 1.3.0

        Parameters
        ----------
        na_rep : str

        Returns
        -------
        self : Styler

        Notes
        -----
        This method is deprecated. See `Styler.format()`
        """
        warnings.warn(
            "this method is deprecated in favour of `Styler.format(na_rep=..)`",
            FutureWarning,
            stacklevel=2,
        )
        self.na_rep = na_rep
        return self.format(na_rep=na_rep, precision=self.precision)

    def hide_index(self) -> Styler:
        """
        Hide any indices from rendering.

        Returns
        -------
        self : Styler
        """
        self.hidden_index = True
        return self

    def hide_columns(self, subset) -> Styler:
        """
        Hide columns from rendering.

        Parameters
        ----------
        subset : IndexSlice
            An argument to ``DataFrame.loc`` that identifies which columns
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

    def background_gradient(
        self,
        cmap="PuBu",
        low: float = 0,
        high: float = 0,
        axis: Optional[Axis] = 0,
        subset=None,
        text_color_threshold: float = 0.408,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
    ) -> Styler:
        """
        Color the background in a gradient style.

        The background color is determined according
        to the data in each column (optionally row). Requires matplotlib.

        Parameters
        ----------
        cmap : str or colormap
            Matplotlib colormap.
        low : float
            Compress the range by the low.
        high : float
            Compress the range by the high.
        axis : {0 or 'index', 1 or 'columns', None}, default 0
            Apply to each column (``axis=0`` or ``'index'``), to each row
            (``axis=1`` or ``'columns'``), or to the entire DataFrame at once
            with ``axis=None``.
        subset : IndexSlice
            A valid slice for ``data`` to limit the style application to.
        text_color_threshold : float or int
            Luminance threshold for determining text color. Facilitates text
            visibility across varying background colors. From 0 to 1.
            0 = all text is dark colored, 1 = all text is light colored.

            .. versionadded:: 0.24.0

        vmin : float, optional
            Minimum data value that corresponds to colormap minimum value.
            When None (default): the minimum value of the data will be used.

            .. versionadded:: 1.0.0

        vmax : float, optional
            Maximum data value that corresponds to colormap maximum value.
            When None (default): the maximum value of the data will be used.

            .. versionadded:: 1.0.0

        Returns
        -------
        self : Styler

        Raises
        ------
        ValueError
            If ``text_color_threshold`` is not a value from 0 to 1.

        Notes
        -----
        Set ``text_color_threshold`` or tune ``low`` and ``high`` to keep the
        text legible by not using the entire range of the color map. The range
        of the data is extended by ``low * (x.max() - x.min())`` and ``high *
        (x.max() - x.min())`` before normalizing.
        """
        if subset is None:
            subset = self.data.select_dtypes(include=np.number).columns

        self.apply(
            self._background_gradient,
            cmap=cmap,
            subset=subset,
            axis=axis,
            low=low,
            high=high,
            text_color_threshold=text_color_threshold,
            vmin=vmin,
            vmax=vmax,
        )
        return self

    @staticmethod
    def _background_gradient(
        s,
        cmap="PuBu",
        low: float = 0,
        high: float = 0,
        text_color_threshold: float = 0.408,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
    ):
        """
        Color background in a range according to the data.
        """
        if (
            not isinstance(text_color_threshold, (float, int))
            or not 0 <= text_color_threshold <= 1
        ):
            msg = "`text_color_threshold` must be a value from 0 to 1."
            raise ValueError(msg)

        with _mpl(Styler.background_gradient) as (plt, colors):
            smin = np.nanmin(s.to_numpy()) if vmin is None else vmin
            smax = np.nanmax(s.to_numpy()) if vmax is None else vmax
            rng = smax - smin
            # extend lower / upper bounds, compresses color range
            norm = colors.Normalize(smin - (rng * low), smax + (rng * high))
            # matplotlib colors.Normalize modifies inplace?
            # https://github.com/matplotlib/matplotlib/issues/5427
            rgbas = plt.cm.get_cmap(cmap)(norm(s.to_numpy(dtype=float)))

            def relative_luminance(rgba) -> float:
                """
                Calculate relative luminance of a color.

                The calculation adheres to the W3C standards
                (https://www.w3.org/WAI/GL/wiki/Relative_luminance)

                Parameters
                ----------
                color : rgb or rgba tuple

                Returns
                -------
                float
                    The relative luminance as a value from 0 to 1
                """
                r, g, b = (
                    x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4
                    for x in rgba[:3]
                )
                return 0.2126 * r + 0.7152 * g + 0.0722 * b

            def css(rgba) -> str:
                dark = relative_luminance(rgba) < text_color_threshold
                text_color = "#f1f1f1" if dark else "#000000"
                return f"background-color: {colors.rgb2hex(rgba)};color: {text_color};"

            if s.ndim == 1:
                return [css(rgba) for rgba in rgbas]
            else:
                return pd.DataFrame(
                    [[css(rgba) for rgba in row] for row in rgbas],
                    index=s.index,
                    columns=s.columns,
                )

    def set_properties(self, subset=None, **kwargs) -> Styler:
        """
        Method to set one or more non-data dependent properties or each cell.

        Parameters
        ----------
        subset : IndexSlice
            A valid slice for ``data`` to limit the style application to.
        **kwargs : dict
            A dictionary of property, value pairs to be set for each cell.

        Returns
        -------
        self : Styler

        Examples
        --------
        >>> df = pd.DataFrame(np.random.randn(10, 4))
        >>> df.style.set_properties(color="white", align="right")
        >>> df.style.set_properties(**{'background-color': 'yellow'})
        """
        values = "".join(f"{p}: {v};" for p, v in kwargs.items())
        return self.applymap(lambda x: values, subset=subset)

    @staticmethod
    def _bar(
        s,
        align: str,
        colors: List[str],
        width: float = 100,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
    ):
        """
        Draw bar chart in dataframe cells.
        """
        # Get input value range.
        smin = np.nanmin(s.to_numpy()) if vmin is None else vmin
        smax = np.nanmax(s.to_numpy()) if vmax is None else vmax
        if align == "mid":
            smin = min(0, smin)
            smax = max(0, smax)
        elif align == "zero":
            # For "zero" mode, we want the range to be symmetrical around zero.
            smax = max(abs(smin), abs(smax))
            smin = -smax
        # Transform to percent-range of linear-gradient
        normed = width * (s.to_numpy(dtype=float) - smin) / (smax - smin + 1e-12)
        zero = -width * smin / (smax - smin + 1e-12)

        def css_bar(start: float, end: float, color: str) -> str:
            """
            Generate CSS code to draw a bar from start to end.
            """
            css = "width: 10em; height: 80%;"
            if end > start:
                css += "background: linear-gradient(90deg,"
                if start > 0:
                    css += f" transparent {start:.1f}%, {color} {start:.1f}%, "
                e = min(end, width)
                css += f"{color} {e:.1f}%, transparent {e:.1f}%)"
            return css

        def css(x):
            if pd.isna(x):
                return ""

            # avoid deprecated indexing `colors[x > zero]`
            color = colors[1] if x > zero else colors[0]

            if align == "left":
                return css_bar(0, x, color)
            else:
                return css_bar(min(x, zero), max(x, zero), color)

        if s.ndim == 1:
            return [css(x) for x in normed]
        else:
            return pd.DataFrame(
                [[css(x) for x in row] for row in normed],
                index=s.index,
                columns=s.columns,
            )

    def bar(
        self,
        subset=None,
        axis: Optional[Axis] = 0,
        color="#d65f5f",
        width: float = 100,
        align: str = "left",
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
    ) -> Styler:
        """
        Draw bar chart in the cell backgrounds.

        Parameters
        ----------
        subset : IndexSlice, optional
            A valid slice for `data` to limit the style application to.
        axis : {0 or 'index', 1 or 'columns', None}, default 0
            Apply to each column (``axis=0`` or ``'index'``), to each row
            (``axis=1`` or ``'columns'``), or to the entire DataFrame at once
            with ``axis=None``.
        color : str or 2-tuple/list
            If a str is passed, the color is the same for both
            negative and positive numbers. If 2-tuple/list is used, the
            first element is the color_negative and the second is the
            color_positive (eg: ['#d65f5f', '#5fba7d']).
        width : float, default 100
            A number between 0 or 100. The largest value will cover `width`
            percent of the cell's width.
        align : {'left', 'zero',' mid'}, default 'left'
            How to align the bars with the cells.

            - 'left' : the min value starts at the left of the cell.
            - 'zero' : a value of zero is located at the center of the cell.
            - 'mid' : the center of the cell is at (max-min)/2, or
              if values are all negative (positive) the zero is aligned
              at the right (left) of the cell.
        vmin : float, optional
            Minimum bar value, defining the left hand limit
            of the bar drawing range, lower values are clipped to `vmin`.
            When None (default): the minimum value of the data will be used.

            .. versionadded:: 0.24.0

        vmax : float, optional
            Maximum bar value, defining the right hand limit
            of the bar drawing range, higher values are clipped to `vmax`.
            When None (default): the maximum value of the data will be used.

            .. versionadded:: 0.24.0

        Returns
        -------
        self : Styler
        """
        if align not in ("left", "zero", "mid"):
            raise ValueError("`align` must be one of {'left', 'zero',' mid'}")

        if not (is_list_like(color)):
            color = [color, color]
        elif len(color) == 1:
            color = [color[0], color[0]]
        elif len(color) > 2:
            raise ValueError(
                "`color` must be string or a list-like "
                "of length 2: [`color_neg`, `color_pos`] "
                "(eg: color=['#d65f5f', '#5fba7d'])"
            )

        if subset is None:
            subset = self.data.select_dtypes(include=np.number).columns

        self.apply(
            self._bar,
            subset=subset,
            axis=axis,
            align=align,
            colors=color,
            width=width,
            vmin=vmin,
            vmax=vmax,
        )

        return self

    def highlight_null(
        self,
        null_color: str = "red",
        subset: Optional[IndexLabel] = None,
        props: Optional[str] = None,
    ) -> Styler:
        """
        Highlight missing values with a style.

        Parameters
        ----------
        null_color : str, default 'red'
        subset : label or list of labels, default None
            A valid slice for ``data`` to limit the style application to.

            .. versionadded:: 1.1.0

        props : str, default None
            CSS properties to use for highlighting. If ``props`` is given, ``color``
            is not used.

            .. versionadded:: 1.3.0

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.highlight_max: Highlight the maximum with a style.
        Styler.highlight_min: Highlight the minimum with a style.
        """

        def f(data: DataFrame, props: str) -> np.ndarray:
            return np.where(pd.isna(data).to_numpy(), props, "")

        if props is None:
            props = f"background-color: {null_color};"
        # error: Argument 1 to "apply" of "Styler" has incompatible type
        # "Callable[[DataFrame, str], ndarray]"; expected "Callable[..., Styler]"
        return self.apply(
            f, axis=None, subset=subset, props=props  # type: ignore[arg-type]
        )

    def highlight_max(
        self,
        subset: Optional[IndexLabel] = None,
        color: str = "yellow",
        axis: Optional[Axis] = 0,
        props: Optional[str] = None,
    ) -> Styler:
        """
        Highlight the maximum with a style.

        Parameters
        ----------
        subset : IndexSlice, default None
            A valid slice for ``data`` to limit the style application to.
        color : str, default 'yellow'
            Background color to use for highlighting.
        axis : {0 or 'index', 1 or 'columns', None}, default 0
            Apply to each column (``axis=0`` or ``'index'``), to each row
            (``axis=1`` or ``'columns'``), or to the entire DataFrame at once
            with ``axis=None``.
        props : str, default None
            CSS properties to use for highlighting. If ``props`` is given, ``color``
            is not used.

            .. versionadded:: 1.3.0

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.highlight_null: Highlight missing values with a style.
        Styler.highlight_min: Highlight the minimum with a style.
        """

        def f(data: FrameOrSeries, props: str) -> np.ndarray:
            return np.where(data == np.nanmax(data.to_numpy()), props, "")

        if props is None:
            props = f"background-color: {color};"
        # error: Argument 1 to "apply" of "Styler" has incompatible type
        # "Callable[[FrameOrSeries, str], ndarray]"; expected "Callable[..., Styler]"
        return self.apply(
            f, axis=axis, subset=subset, props=props  # type: ignore[arg-type]
        )

    def highlight_min(
        self,
        subset: Optional[IndexLabel] = None,
        color: str = "yellow",
        axis: Optional[Axis] = 0,
        props: Optional[str] = None,
    ) -> Styler:
        """
        Highlight the minimum with a style.

        Parameters
        ----------
        subset : IndexSlice, default None
            A valid slice for ``data`` to limit the style application to.
        color : str, default 'yellow'
            Background color to use for highlighting.
        axis : {0 or 'index', 1 or 'columns', None}, default 0
            Apply to each column (``axis=0`` or ``'index'``), to each row
            (``axis=1`` or ``'columns'``), or to the entire DataFrame at once
            with ``axis=None``.
        props : str, default None
            CSS properties to use for highlighting. If ``props`` is given, ``color``
            is not used.

            .. versionadded:: 1.3.0

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.highlight_null: Highlight missing values with a style.
        Styler.highlight_max: Highlight the maximum with a style.
        """

        def f(data: FrameOrSeries, props: str) -> np.ndarray:
            return np.where(data == np.nanmin(data.to_numpy()), props, "")

        if props is None:
            props = f"background-color: {color};"
        # error: Argument 1 to "apply" of "Styler" has incompatible type
        # "Callable[[FrameOrSeries, str], ndarray]"; expected "Callable[..., Styler]"
        return self.apply(
            f, axis=axis, subset=subset, props=props  # type: ignore[arg-type]
        )

    @classmethod
    def from_custom_template(cls, searchpath, name):
        """
        Factory function for creating a subclass of ``Styler``.

        Uses a custom template and Jinja environment.

        Parameters
        ----------
        searchpath : str or list
            Path or paths of directories containing the templates.
        name : str
            Name of your custom template to use for rendering.

        Returns
        -------
        MyStyler : subclass of Styler
            Has the correct ``env`` and ``template`` class attributes set.
        """
        loader = jinja2.ChoiceLoader([jinja2.FileSystemLoader(searchpath), cls.loader])

        # mypy doesn't like dynamically-defined classes
        # error: Variable "cls" is not valid as a type
        # error: Invalid base class "cls"
        class MyStyler(cls):  # type:ignore[valid-type,misc]
            env = jinja2.Environment(loader=loader)
            template = env.get_template(name)

        return MyStyler

    def pipe(self, func: Callable, *args, **kwargs):
        """
        Apply ``func(self, *args, **kwargs)``, and return the result.

        .. versionadded:: 0.24.0

        Parameters
        ----------
        func : function
            Function to apply to the Styler.  Alternatively, a
            ``(callable, keyword)`` tuple where ``keyword`` is a string
            indicating the keyword of ``callable`` that expects the Styler.
        *args : optional
            Arguments passed to `func`.
        **kwargs : optional
            A dictionary of keyword arguments passed into ``func``.

        Returns
        -------
        object :
            The value returned by ``func``.

        See Also
        --------
        DataFrame.pipe : Analogous method for DataFrame.
        Styler.apply : Apply a function row-wise, column-wise, or table-wise to
            modify the dataframe's styling.

        Notes
        -----
        Like :meth:`DataFrame.pipe`, this method can simplify the
        application of several user-defined functions to a styler.  Instead
        of writing:

        .. code-block:: python

            f(g(df.style.set_precision(3), arg1=a), arg2=b, arg3=c)

        users can write:

        .. code-block:: python

            (df.style.set_precision(3)
               .pipe(g, arg1=a)
               .pipe(f, arg2=b, arg3=c))

        In particular, this allows users to define functions that take a
        styler object, along with other parameters, and return the styler after
        making styling changes (such as calling :meth:`Styler.apply` or
        :meth:`Styler.set_properties`).  Using ``.pipe``, these user-defined
        style "transformations" can be interleaved with calls to the built-in
        Styler interface.

        Examples
        --------
        >>> def format_conversion(styler):
        ...     return (styler.set_properties(**{'text-align': 'right'})
        ...                   .format({'conversion': '{:.1%}'}))

        The user-defined ``format_conversion`` function above can be called
        within a sequence of other style modifications:

        >>> df = pd.DataFrame({'trial': list(range(5)),
        ...                    'conversion': [0.75, 0.85, np.nan, 0.7, 0.72]})
        >>> (df.style
        ...    .highlight_min(subset=['conversion'], color='yellow')
        ...    .pipe(format_conversion)
        ...    .set_caption("Results with minimum conversion highlighted."))
        """
        return com.pipe(self, func, *args, **kwargs)


class _Tooltips:
    """
    An extension to ``Styler`` that allows for and manipulates tooltips on hover
    of table data-cells in the HTML result.

    Parameters
    ----------
    css_name: str, default "pd-t"
        Name of the CSS class that controls visualisation of tooltips.
    css_props: list-like, default; see Notes
        List of (attr, value) tuples defining properties of the CSS class.
    tooltips: DataFrame, default empty
        DataFrame of strings aligned with underlying ``Styler`` data for tooltip
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
                "props": _maybe_convert_css_to_tuples(self.class_properties),
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

    def _translate(self, styler_data: FrameOrSeriesUnion, uuid: str, d: Dict):
        """
        Mutate the render dictionary to allow for tooltips:

        - Add `<span>` HTML element to each data cells `display_value`. Ignores
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


def _is_visible(idx_row, idx_col, lengths) -> bool:
    """
    Index -> {(idx_row, idx_col): bool}).
    """
    return (idx_col, idx_row) in lengths


def _get_level_lengths(index, hidden_elements=None):
    """
    Given an index, find the level length for each element.

    Optional argument is a list of index positions which
    should not be visible.

    Result is a dictionary of (level, initial_position): span
    """
    if isinstance(index, pd.MultiIndex):
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


def _default_formatter(x: Any, precision: Optional[int] = None) -> Any:
    """
    Format the display of a value

    Parameters
    ----------
    x : Any
        Input variable to be formatted
    precision : Int, optional
        Floating point precision used if ``x`` is float or complex.

    Returns
    -------
    value : Any
        Matches input type, or string if input is float or complex.
    """
    if precision is None:
        precision = get_option("display.precision")
    if isinstance(x, (float, complex)):
        return f"{x:.{precision}f}"
    return x


def _maybe_wrap_formatter(
    formatter: Optional[BaseFormatter] = None,
    na_rep: Optional[str] = None,
    precision: Optional[int] = None,
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
        formatter_func = partial(_default_formatter, precision=precision)
    else:
        raise TypeError(f"'formatter' expected str or callable, got {type(formatter)}")

    if na_rep is None:
        return formatter_func
    else:
        return lambda x: na_rep if pd.isna(x) else formatter_func(x)


def _maybe_convert_css_to_tuples(style: CSSProperties) -> CSSList:
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


def _non_reducing_slice(slice_):
    """
    Ensure that a slice doesn't reduce to a Series or Scalar.

    Any user-passed `subset` should have this called on it
    to make sure we're always working with DataFrames.
    """
    # default to column slice, like DataFrame
    # ['A', 'B'] -> IndexSlices[:, ['A', 'B']]
    kinds = (ABCSeries, np.ndarray, Index, list, str)
    if isinstance(slice_, kinds):
        slice_ = pd.IndexSlice[:, slice_]

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
