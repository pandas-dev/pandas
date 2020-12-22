"""
Module for applying conditional formatting to DataFrames and Series.
"""
from collections import defaultdict
from contextlib import contextmanager
import copy
from functools import partial
from itertools import product
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
)
from uuid import uuid4

import numpy as np

from pandas._config import get_option

from pandas._libs import lib
from pandas._typing import Axis, FrameOrSeries, FrameOrSeriesUnion, Label
from pandas.compat._optional import import_optional_dependency
from pandas.util._decorators import doc

from pandas.core.dtypes.common import is_float

import pandas as pd
from pandas.api.types import is_dict_like, is_list_like
from pandas.core import generic
import pandas.core.common as com
from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame
from pandas.core.indexing import maybe_numeric_slice, non_reducing_slice

jinja2 = import_optional_dependency("jinja2", extra="DataFrame.style requires jinja2.")


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
        table_styles: Optional[List[Dict[str, List[Tuple[str, str]]]]] = None,
        uuid: Optional[str] = None,
        caption: Optional[str] = None,
        table_attributes: Optional[str] = None,
        cell_ids: bool = True,
        na_rep: Optional[str] = None,
        uuid_len: int = 5,
    ):
        self.ctx: DefaultDict[Tuple[int, int], List[str]] = defaultdict(list)
        self._todo: List[Tuple[Callable, Tuple, Dict]] = []

        if not isinstance(data, (pd.Series, pd.DataFrame)):
            raise TypeError("``data`` must be a Series or DataFrame")
        if data.ndim == 1:
            data = data.to_frame()
        if not data.index.is_unique or not data.columns.is_unique:
            raise ValueError("style is not supported for non-unique indices.")

        self.data = data
        self.index = data.index
        self.columns = data.columns

        if not isinstance(uuid_len, int) or not uuid_len >= 0:
            raise TypeError("``uuid_len`` must be an integer in range [0, 32].")
        self.uuid_len = min(32, uuid_len)
        self.uuid = (uuid or uuid4().hex[: self.uuid_len]) + "_"
        self.table_styles = table_styles
        self.caption = caption
        if precision is None:
            precision = get_option("display.precision")
        self.precision = precision
        self.table_attributes = table_attributes
        self.hidden_index = False
        self.hidden_columns: Sequence[int] = []
        self.cell_ids = cell_ids
        self.na_rep = na_rep

        self.cell_context: Dict[str, Any] = {}

        # display_funcs maps (row, col) -> formatting function

        def default_display_func(x):
            if self.na_rep is not None and pd.isna(x):
                return self.na_rep
            elif is_float(x):
                display_format = f"{x:.{self.precision}f}"
                return display_format
            else:
                return x

        self._display_funcs: DefaultDict[
            Tuple[int, int], Callable[[Any], str]
        ] = defaultdict(lambda: default_display_func)

    def _repr_html_(self) -> str:
        """
        Hooks into Jupyter notebook rich display system.
        """
        return self.render()

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
        columns: Optional[Sequence[Label]] = None,
        header: Union[Sequence[Label], bool] = True,
        index: bool = True,
        index_label: Optional[Union[Label, Sequence[Label]]] = None,
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
        precision = self.precision
        hidden_index = self.hidden_index
        hidden_columns = self.hidden_columns
        uuid = self.uuid
        ROW_HEADING_CLASS = "row_heading"
        COL_HEADING_CLASS = "col_heading"
        INDEX_NAME_CLASS = "index_name"

        DATA_CLASS = "data"
        BLANK_CLASS = "blank"
        BLANK_VALUE = ""

        def format_attr(pair):
            return f"{pair['key']}={pair['value']}"

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

        cellstyle_map = defaultdict(list)
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
                        es["attributes"] = [
                            format_attr({"key": "colspan", "value": f'"{colspan}"'})
                        ]
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
                [{"type": "th", "value": BLANK_VALUE, "class": " ".join([BLANK_CLASS])}]
                * (len(clabels[0]) - len(hidden_columns))
            )

            head.append(index_header_row)

        body = []
        for r, idx in enumerate(self.data.index):
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
                    es["attributes"] = [
                        format_attr({"key": "rowspan", "value": f'"{rowspan}"'})
                    ]
                row_es.append(es)

            for c, col in enumerate(self.data.columns):
                cs = [DATA_CLASS, f"row{r}", f"col{c}"]
                cs.extend(cell_context.get("data", {}).get(r, {}).get(c, []))
                formatter = self._display_funcs[(r, c)]
                value = self.data.iloc[r, c]
                row_dict = {
                    "type": "td",
                    "value": value,
                    "class": " ".join(cs),
                    "display_value": formatter(value),
                    "is_visible": (c not in hidden_columns),
                }
                # only add an id if the cell has a style
                props = []
                if self.cell_ids or (r, c) in ctx:
                    row_dict["id"] = "_".join(cs[1:])
                    for x in ctx[r, c]:
                        # have to handle empty styles like ['']
                        if x.count(":"):
                            props.append(tuple(x.split(":")))
                        else:
                            props.append(("", ""))
                row_es.append(row_dict)
                cellstyle_map[tuple(props)].append(f"row{r}_col{c}")
            body.append(row_es)

        cellstyle = [
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

        return {
            "head": head,
            "cellstyle": cellstyle,
            "body": body,
            "uuid": uuid,
            "precision": precision,
            "table_styles": table_styles,
            "caption": caption,
            "table_attributes": table_attr,
        }

    def format(self, formatter, subset=None, na_rep: Optional[str] = None) -> "Styler":
        """
        Format the text display value of cells.

        Parameters
        ----------
        formatter : str, callable, dict or None
            If ``formatter`` is None, the default formatter is used.
        subset : IndexSlice
            An argument to ``DataFrame.loc`` that restricts which elements
            ``formatter`` is applied to.
        na_rep : str, optional
            Representation for missing values.
            If ``na_rep`` is None, no special formatting is applied.

            .. versionadded:: 1.0.0

        Returns
        -------
        self : Styler

        Notes
        -----
        ``formatter`` is either an ``a`` or a dict ``{column name: a}`` where
        ``a`` is one of

        - str: this will be wrapped in: ``a.format(x)``
        - callable: called with the value of an individual cell

        The default display value for numeric values is the "general" (``g``)
        format with ``pd.options.display.precision`` precision.

        Examples
        --------
        >>> df = pd.DataFrame(np.random.randn(4, 2), columns=['a', 'b'])
        >>> df.style.format("{:.2%}")
        >>> df['c'] = ['a', 'b', 'c', 'd']
        >>> df.style.format({'c': str.upper})
        """
        if formatter is None:
            assert self._display_funcs.default_factory is not None
            formatter = self._display_funcs.default_factory()

        if subset is None:
            row_locs = range(len(self.data))
            col_locs = range(len(self.data.columns))
        else:
            subset = non_reducing_slice(subset)
            if len(subset) == 1:
                subset = subset, self.data.columns

            sub_df = self.data.loc[subset]
            row_locs = self.data.index.get_indexer_for(sub_df.index)
            col_locs = self.data.columns.get_indexer_for(sub_df.columns)

        if is_dict_like(formatter):
            for col, col_formatter in formatter.items():
                # formatter must be callable, so '{}' are converted to lambdas
                col_formatter = _maybe_wrap_formatter(col_formatter, na_rep)
                col_num = self.data.columns.get_indexer_for([col])[0]

                for row_num in row_locs:
                    self._display_funcs[(row_num, col_num)] = col_formatter
        else:
            # single scalar to format all cells with
            formatter = _maybe_wrap_formatter(formatter, na_rep)
            locs = product(*(row_locs, col_locs))
            for i, j in locs:
                self._display_funcs[(i, j)] = formatter
        return self

    def set_td_classes(self, classes: DataFrame) -> "Styler":
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
            r: {c: [str(classes.iloc[r, c])]}
            for r, rn in enumerate(classes.index)
            for c, cn in enumerate(classes.columns)
            if not mask.iloc[r, c]
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
        * precision
        * table_styles
        * caption
        * table_attributes
        """
        self._compute()
        # TODO: namespace all the pandas keys
        d = self._translate()
        # filter out empty styles, every cell will have a class
        # but the list of props may just be [['', '']].
        # so we have the nested anys below
        trimmed = [x for x in d["cellstyle"] if any(any(y) for y in x["props"])]
        d["cellstyle"] = trimmed
        d.update(kwargs)
        return self.template.render(**d)

    def _update_ctx(self, attrs: DataFrame) -> None:
        """
        Update the state of the Styler.

        Collects a mapping of {index_label: ['<property>: <value>']}.

        Parameters
        ----------
        attrs : DataFrame
            should contain strings of '<property>: <value>;<prop2>: <val2>'
            Whitespace shouldn't matter and the final trailing ';' shouldn't
            matter.
        """
        coli = {k: i for i, k in enumerate(self.columns)}
        rowi = {k: i for i, k in enumerate(self.index)}
        for jj in range(len(attrs.columns)):
            cn = attrs.columns[jj]
            j = coli[cn]
            for rn, c in attrs[[cn]].itertuples():
                if not c:
                    continue
                c = c.rstrip(";")
                if not c:
                    continue
                i = rowi[rn]
                for pair in c.split(";"):
                    self.ctx[(i, j)].append(pair)

    def _copy(self, deepcopy: bool = False) -> "Styler":
        styler = Styler(
            self.data,
            precision=self.precision,
            caption=self.caption,
            uuid=self.uuid,
            table_styles=self.table_styles,
            na_rep=self.na_rep,
        )
        if deepcopy:
            styler.ctx = copy.deepcopy(self.ctx)
            styler._todo = copy.deepcopy(self._todo)
        else:
            styler.ctx = self.ctx
            styler._todo = self._todo
        return styler

    def __copy__(self) -> "Styler":
        """
        Deep copy by default.
        """
        return self._copy(deepcopy=False)

    def __deepcopy__(self, memo) -> "Styler":
        return self._copy(deepcopy=True)

    def clear(self) -> None:
        """
        Reset the styler, removing any previously applied styles.

        Returns None.
        """
        self.ctx.clear()
        self.cell_context = {}
        self._todo = []

    def _compute(self):
        """
        Execute the style functions built up in `self._todo`.

        Relies on the conventions that all style functions go through
        .apply or .applymap. The append styles to apply as tuples of

        (application method, *args, **kwargs)
        """
        r = self
        for func, args, kwargs in self._todo:
            r = func(self)(*args, **kwargs)
        return r

    def _apply(
        self,
        func: Callable[..., "Styler"],
        axis: Optional[Axis] = 0,
        subset=None,
        **kwargs,
    ) -> "Styler":
        subset = slice(None) if subset is None else subset
        subset = non_reducing_slice(subset)
        data = self.data.loc[subset]
        if axis is not None:
            result = data.apply(func, axis=axis, result_type="expand", **kwargs)
            result.columns = data.columns
        else:
            result = func(data, **kwargs)
            if not isinstance(result, pd.DataFrame):
                raise TypeError(
                    f"Function {repr(func)} must return a DataFrame when "
                    f"passed to `Styler.apply` with axis=None"
                )
            if not (
                result.index.equals(data.index) and result.columns.equals(data.columns)
            ):
                raise ValueError(
                    f"Result of {repr(func)} must have identical "
                    f"index and columns as the input"
                )

        result_shape = result.shape
        expected_shape = self.data.loc[subset].shape
        if result_shape != expected_shape:
            raise ValueError(
                f"Function {repr(func)} returned the wrong shape.\n"
                f"Result has shape: {result.shape}\n"
                f"Expected shape:   {expected_shape}"
            )
        self._update_ctx(result)
        return self

    def apply(
        self,
        func: Callable[..., "Styler"],
        axis: Optional[Axis] = 0,
        subset=None,
        **kwargs,
    ) -> "Styler":
        """
        Apply a function column-wise, row-wise, or table-wise.

        Updates the HTML representation with the result.

        Parameters
        ----------
        func : function
            ``func`` should take a Series or DataFrame (depending
            on ``axis``), and return an object with the same shape.
            Must return a DataFrame with identical index and
            column labels when ``axis=None``.
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
        The output shape of ``func`` should match the input, i.e. if
        ``x`` is the input row, column, or table (depending on ``axis``),
        then ``func(x).shape == x.shape`` should be true.

        This is similar to ``DataFrame.apply``, except that ``axis=None``
        applies the function to the entire DataFrame at once,
        rather than column-wise or row-wise.

        Examples
        --------
        >>> def highlight_max(x):
        ...     return ['background-color: yellow' if v == x.max() else ''
                        for v in x]
        ...
        >>> df = pd.DataFrame(np.random.randn(5, 2))
        >>> df.style.apply(highlight_max)
        """
        self._todo.append(
            (lambda instance: getattr(instance, "_apply"), (func, axis, subset), kwargs)
        )
        return self

    def _applymap(self, func: Callable, subset=None, **kwargs) -> "Styler":
        func = partial(func, **kwargs)  # applymap doesn't take kwargs?
        if subset is None:
            subset = pd.IndexSlice[:]
        subset = non_reducing_slice(subset)
        result = self.data.loc[subset].applymap(func)
        self._update_ctx(result)
        return self

    def applymap(self, func: Callable, subset=None, **kwargs) -> "Styler":
        """
        Apply a function elementwise.

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
    ) -> "Styler":
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

    def set_precision(self, precision: int) -> "Styler":
        """
        Set the precision used to render.

        Parameters
        ----------
        precision : int

        Returns
        -------
        self : Styler
        """
        self.precision = precision
        return self

    def set_table_attributes(self, attributes: str) -> "Styler":
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

    def use(self, styles: List[Tuple[Callable, Tuple, Dict]]) -> "Styler":
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

    def set_uuid(self, uuid: str) -> "Styler":
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

    def set_caption(self, caption: str) -> "Styler":
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

    def set_table_styles(self, table_styles, axis=0, overwrite=True) -> "Styler":
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

        Adding column styling by name

        >>> df.style.set_table_styles({
        ...     'A': [{'selector': '',
        ...            'props': [('color', 'red')]}],
        ...     'B': [{'selector': 'td',
        ...            'props': [('color', 'blue')]}]
        ... }, overwrite=False)

        Adding row styling

        >>> df.style.set_table_styles({
        ...     0: [{'selector': 'td:hover',
        ...          'props': [('font-size', '25px')]}]
        ... }, axis=1, overwrite=False)
        """
        if is_dict_like(table_styles):
            if axis in [0, "index"]:
                obj, idf = self.data.columns, ".col"
            else:
                obj, idf = self.data.index, ".row"

            table_styles = [
                {
                    "selector": s["selector"] + idf + str(obj.get_loc(key)),
                    "props": s["props"],
                }
                for key, styles in table_styles.items()
                for s in styles
            ]

        if not overwrite and self.table_styles is not None:
            self.table_styles.extend(table_styles)
        else:
            self.table_styles = table_styles
        return self

    def set_na_rep(self, na_rep: str) -> "Styler":
        """
        Set the missing data representation on a Styler.

        .. versionadded:: 1.0.0

        Parameters
        ----------
        na_rep : str

        Returns
        -------
        self : Styler
        """
        self.na_rep = na_rep
        return self

    def hide_index(self) -> "Styler":
        """
        Hide any indices from rendering.

        Returns
        -------
        self : Styler
        """
        self.hidden_index = True
        return self

    def hide_columns(self, subset) -> "Styler":
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
        subset = non_reducing_slice(subset)
        hidden_df = self.data.loc[subset]
        self.hidden_columns = self.columns.get_indexer_for(hidden_df.columns)
        return self

    # -----------------------------------------------------------------------
    # A collection of "builtin" styles
    # -----------------------------------------------------------------------

    @staticmethod
    def _highlight_null(v, null_color: str) -> str:
        return f"background-color: {null_color}" if pd.isna(v) else ""

    def highlight_null(
        self,
        null_color: str = "red",
        subset: Optional[Union[Label, Sequence[Label]]] = None,
    ) -> "Styler":
        """
        Shade the background ``null_color`` for missing values.

        Parameters
        ----------
        null_color : str, default 'red'
        subset : label or list of labels, default None
            A valid slice for ``data`` to limit the style application to.

            .. versionadded:: 1.1.0

        Returns
        -------
        self : Styler
        """
        self.applymap(self._highlight_null, null_color=null_color, subset=subset)
        return self

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
    ) -> "Styler":
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
        subset = maybe_numeric_slice(self.data, subset)
        subset = non_reducing_slice(subset)
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
                    x / 12.92 if x <= 0.03928 else ((x + 0.055) / 1.055 ** 2.4)
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

    def set_properties(self, subset=None, **kwargs) -> "Styler":
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
        values = ";".join(f"{p}: {v}" for p, v in kwargs.items())
        f = lambda x: values
        return self.applymap(f, subset=subset)

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
    ) -> "Styler":
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

        subset = maybe_numeric_slice(self.data, subset)
        subset = non_reducing_slice(subset)
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

    def highlight_max(
        self, subset=None, color: str = "yellow", axis: Optional[Axis] = 0
    ) -> "Styler":
        """
        Highlight the maximum by shading the background.

        Parameters
        ----------
        subset : IndexSlice, default None
            A valid slice for ``data`` to limit the style application to.
        color : str, default 'yellow'
        axis : {0 or 'index', 1 or 'columns', None}, default 0
            Apply to each column (``axis=0`` or ``'index'``), to each row
            (``axis=1`` or ``'columns'``), or to the entire DataFrame at once
            with ``axis=None``.

        Returns
        -------
        self : Styler
        """
        return self._highlight_handler(subset=subset, color=color, axis=axis, max_=True)

    def highlight_min(
        self, subset=None, color: str = "yellow", axis: Optional[Axis] = 0
    ) -> "Styler":
        """
        Highlight the minimum by shading the background.

        Parameters
        ----------
        subset : IndexSlice, default None
            A valid slice for ``data`` to limit the style application to.
        color : str, default 'yellow'
        axis : {0 or 'index', 1 or 'columns', None}, default 0
            Apply to each column (``axis=0`` or ``'index'``), to each row
            (``axis=1`` or ``'columns'``), or to the entire DataFrame at once
            with ``axis=None``.

        Returns
        -------
        self : Styler
        """
        return self._highlight_handler(
            subset=subset, color=color, axis=axis, max_=False
        )

    def _highlight_handler(
        self,
        subset=None,
        color: str = "yellow",
        axis: Optional[Axis] = None,
        max_: bool = True,
    ) -> "Styler":
        subset = non_reducing_slice(maybe_numeric_slice(self.data, subset))
        self.apply(
            self._highlight_extrema, color=color, axis=axis, subset=subset, max_=max_
        )
        return self

    @staticmethod
    def _highlight_extrema(
        data: FrameOrSeries, color: str = "yellow", max_: bool = True
    ):
        """
        Highlight the min or max in a Series or DataFrame.
        """
        attr = f"background-color: {color}"

        if max_:
            extrema = data == np.nanmax(data.to_numpy())
        else:
            extrema = data == np.nanmin(data.to_numpy())

        if data.ndim == 1:  # Series from .apply
            return [attr if v else "" for v in extrema]
        else:  # DataFrame from .tee
            return pd.DataFrame(
                np.where(extrema, attr, ""), index=data.index, columns=data.columns
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

        # mypy doesnt like dynamically-defined class
        # error: Variable "cls" is not valid as a type  [valid-type]
        # error: Invalid base class "cls"  [misc]
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


def _maybe_wrap_formatter(
    formatter: Union[Callable, str], na_rep: Optional[str]
) -> Callable:
    if isinstance(formatter, str):
        formatter_func = lambda x: formatter.format(x)
    elif callable(formatter):
        formatter_func = formatter
    else:
        msg = f"Expected a template string or callable, got {formatter} instead"
        raise TypeError(msg)

    if na_rep is None:
        return formatter_func
    elif isinstance(na_rep, str):
        return lambda x: na_rep if pd.isna(x) else formatter_func(x)
    else:
        msg = f"Expected a string, got {na_rep} instead"
        raise TypeError(msg)
