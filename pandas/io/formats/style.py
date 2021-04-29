"""
Module for applying conditional formatting to DataFrames and Series.
"""
from __future__ import annotations

from contextlib import contextmanager
import copy
from functools import partial
import operator
from typing import (
    Any,
    Callable,
    Hashable,
    Sequence,
)
import warnings

import numpy as np

from pandas._typing import (
    Axis,
    FrameOrSeries,
    FrameOrSeriesUnion,
    IndexLabel,
    Scalar,
)
from pandas.compat._optional import import_optional_dependency
from pandas.util._decorators import doc

import pandas as pd
from pandas.api.types import is_list_like
from pandas.core import generic
import pandas.core.common as com
from pandas.core.frame import (
    DataFrame,
    Series,
)
from pandas.core.generic import NDFrame

jinja2 = import_optional_dependency("jinja2", extra="DataFrame.style requires jinja2.")

from pandas.io.formats.style_render import (
    CSSProperties,
    CSSStyles,
    StylerRenderer,
    Tooltips,
    maybe_convert_css_to_tuples,
    non_reducing_slice,
)

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


class Styler(StylerRenderer):
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

    decimal : str, default "."
        Character used as decimal separator for floats, complex and integers

        .. versionadded:: 1.3.0

    thousands : str, optional, default None
        Character used as thousands separator for floats, complex and integers

        .. versionadded:: 1.3.0

    escape : bool, default False
        Replace the characters ``&``, ``<``, ``>``, ``'``, and ``"`` in cell display
        strings with HTML-safe sequences.

        ... versionadded:: 1.3.0

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

    def __init__(
        self,
        data: FrameOrSeriesUnion,
        precision: int | None = None,
        table_styles: CSSStyles | None = None,
        uuid: str | None = None,
        caption: str | None = None,
        table_attributes: str | None = None,
        cell_ids: bool = True,
        na_rep: str | None = None,
        uuid_len: int = 5,
        decimal: str = ".",
        thousands: str | None = None,
        escape: bool = False,
    ):
        super().__init__(
            data=data,
            uuid=uuid,
            uuid_len=uuid_len,
            table_styles=table_styles,
            table_attributes=table_attributes,
            caption=caption,
            cell_ids=cell_ids,
        )

        # validate ordered args
        self.precision = precision  # can be removed on set_precision depr cycle
        self.na_rep = na_rep  # can be removed on set_na_rep depr cycle
        self.format(
            formatter=None,
            precision=precision,
            na_rep=na_rep,
            escape=escape,
            decimal=decimal,
            thousands=thousands,
        )

    def _repr_html_(self) -> str:
        """
        Hooks into Jupyter notebook rich display system.
        """
        return self._render_html()

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
        return self._render_html(**kwargs)

    def set_tooltips(
        self,
        ttips: DataFrame,
        props: CSSProperties | None = None,
        css_class: str | None = None,
    ) -> Styler:
        """
        Set the DataFrame of strings on ``Styler`` generating ``:hover`` tooltips.

        These string based tooltips are only applicable to ``<td>`` HTML elements,
        and cannot be used for column or index headers.

        .. versionadded:: 1.3.0

        Parameters
        ----------
        ttips : DataFrame
            DataFrame containing strings that will be translated to tooltips, mapped
            by identical column and index values that must exist on the underlying
            Styler data. None, NaN values, and empty strings will be ignored and
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
            self.tooltips = Tooltips()
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
        float_format: str | None = None,
        columns: Sequence[Hashable] | None = None,
        header: Sequence[Hashable] | bool = True,
        index: bool = True,
        index_label: IndexLabel | None = None,
        startrow: int = 0,
        startcol: int = 0,
        engine: str | None = None,
        merge_cells: bool = True,
        encoding: str | None = None,
        inf_rep: str = "inf",
        verbose: bool = True,
        freeze_panes: tuple[int, int] | None = None,
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

    def set_td_classes(self, classes: DataFrame) -> Styler:
        """
        Set the DataFrame of strings added to the ``class`` attribute of ``<td>``
        HTML elements.

        Parameters
        ----------
        classes : DataFrame
            DataFrame containing strings that will be translated to CSS classes,
            mapped by identical column and index key values that must exist on the
            underlying Styler data. None, NaN values, and empty strings will
            be ignored and not affect the rendered HTML.

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.set_table_styles: Set the table styles included within the ``<style>``
            HTML element.
        Styler.set_table_attributes: Set the table attributes added to the ``<table>``
            HTML element.

        Notes
        -----
        Can be used in combination with ``Styler.set_table_styles`` to define an
        internal CSS solution without reference to external CSS files.

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
        >>> css = pd.DataFrame([["other-class"]])
        >>> s = Styler(df, uuid="_", cell_ids=False).set_td_classes(css)
        >>> s.hide_index().render()
        '<style type="text/css"></style>'
        '<table id="T__">'
        '  <thead>'
        '    <tr><th class="col_heading level0 col0" >0</th></tr>'
        '  </thead>'
        '  <tbody>'
        '    <tr><td class="data row0 col0 other-class" >1</td></tr>'
        '  </tbody>'
        '</table>'
        """
        classes = classes.reindex_like(self.data)

        for r, row_tup in enumerate(classes.itertuples()):
            for c, value in enumerate(row_tup[1:]):
                if not (pd.isna(value) or value == ""):
                    self.cell_context[(r, c)] = str(value)

        return self

    def _update_ctx(self, attrs: DataFrame) -> None:
        """
        Update the state of the ``Styler`` for data cells.

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
                css_list = maybe_convert_css_to_tuples(c)
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
        Reset the ``Styler``, removing any previously applied styles.

        Returns None.
        """
        self.ctx.clear()
        self.tooltips = None
        self.cell_context.clear()
        self._todo.clear()

        self.hidden_index = False
        self.hidden_columns = []
        # self.format and self.table_styles may be dependent on user
        # input in self.__init__()

    def _apply(
        self,
        func: Callable[..., Styler],
        axis: Axis | None = 0,
        subset=None,
        **kwargs,
    ) -> Styler:
        subset = slice(None) if subset is None else subset
        subset = non_reducing_slice(subset)
        data = self.data.loc[subset]
        if axis is not None:
            result = data.apply(func, axis=axis, result_type="expand", **kwargs)
            result.columns = data.columns
        else:
            result = func(data, **kwargs)
            if not isinstance(result, DataFrame):
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
        axis: Axis | None = 0,
        subset=None,
        **kwargs,
    ) -> Styler:
        """
        Apply a CSS-styling function column-wise, row-wise, or table-wise.

        Updates the HTML representation with the result.

        Parameters
        ----------
        func : function
            ``func`` should take a Series if ``axis`` in [0,1] and return an object
            of same length, also with identical index if the object is a Series.
            ``func`` should take a DataFrame if ``axis`` is ``None`` and return either
            an ndarray with the same shape or a DataFrame with identical columns and
            index.

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

        See Also
        --------
        Styler.applymap: Apply a CSS-styling function elementwise.

        Notes
        -----
        The elements of the output of ``func`` should be CSS styles as strings, in the
        format 'attribute: value; attribute2: value2; ...' or,
        if nothing is to be applied to that element, an empty string or ``None``.

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
        subset = non_reducing_slice(subset)
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
        Styler.apply: Apply a CSS-styling function column-wise, row-wise, or table-wise.

        Notes
        -----
        The elements of the output of ``func`` should be CSS styles as strings, in the
        format 'attribute: value; attribute2: value2; ...' or,
        if nothing is to be applied to that element, an empty string or ``None``.

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
        other: str | None = None,
        subset=None,
        **kwargs,
    ) -> Styler:
        """
        Apply CSS-styles based on a conditional function elementwise.

        .. deprecated:: 1.3.0

        Updates the HTML representation with a style which is
        selected in accordance with the return value of a function.

        Parameters
        ----------
        cond : callable
            ``cond`` should take a scalar, and optional keyword arguments, and return
            a boolean.
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
        Styler.applymap: Apply a CSS-styling function elementwise.
        Styler.apply: Apply a CSS-styling function column-wise, row-wise, or table-wise.

        Notes
        -----
        This method is deprecated.

        This method is a convenience wrapper for :meth:`Styler.applymap`, which we
        recommend using instead.

        The example:
        >>> df = pd.DataFrame([[1, 2], [3, 4]])
        >>> def cond(v, limit=4):
        ...     return v > 1 and v != limit
        >>> df.style.where(cond, value='color:green;', other='color:red;')

        should be refactored to:
        >>> def style_func(v, value, other, limit=4):
        ...     cond = v > 1 and v != limit
        ...     return value if cond else other
        >>> df.style.applymap(style_func, value='color:green;', other='color:red;')
        """
        warnings.warn(
            "this method is deprecated in favour of `Styler.applymap()`",
            FutureWarning,
            stacklevel=2,
        )

        if other is None:
            other = ""

        return self.applymap(
            lambda val: value if cond(val, **kwargs) else other,
            subset=subset,
        )

    def set_precision(self, precision: int) -> StylerRenderer:
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
        Set the table attributes added to the ``<table>`` HTML element.

        These are items in addition to automatic (by default) ``id`` attribute.

        Parameters
        ----------
        attributes : str

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.set_table_styles: Set the table styles included within the ``<style>``
            HTML element.
        Styler.set_td_classes: Set the DataFrame of strings added to the ``class``
            attribute of ``<td>`` HTML elements.

        Examples
        --------
        >>> df = pd.DataFrame(np.random.randn(10, 4))
        >>> df.style.set_table_attributes('class="pure-table"')
        # ... <table class="pure-table"> ...
        """
        self.table_attributes = attributes
        return self

    def export(self) -> list[tuple[Callable, tuple, dict]]:
        """
        Export the styles applied to the current ``Styler``.

        Can be applied to a second Styler with ``Styler.use``.

        Returns
        -------
        styles : list

        See Also
        --------
        Styler.use: Set the styles on the current ``Styler``.
        """
        return self._todo

    def use(self, styles: list[tuple[Callable, tuple, dict]]) -> Styler:
        """
        Set the styles on the current ``Styler``.

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
        Styler.export : Export the styles to applied to the current ``Styler``.
        """
        self._todo.extend(styles)
        return self

    def set_uuid(self, uuid: str) -> Styler:
        """
        Set the uuid applied to ``id`` attributes of HTML elements.

        Parameters
        ----------
        uuid : str

        Returns
        -------
        self : Styler

        Notes
        -----
        Almost all HTML elements within the table, and including the ``<table>`` element
        are assigned ``id`` attributes. The format is ``T_uuid_<extra>`` where
        ``<extra>`` is typically a more specific identifier, such as ``row1_col2``.
        """
        self.uuid = uuid
        return self

    def set_caption(self, caption: str) -> Styler:
        """
        Set the text added to a ``<caption>`` HTML element.

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
        table_styles: dict[Any, CSSStyles] | CSSStyles,
        axis: int = 0,
        overwrite: bool = True,
    ) -> Styler:
        """
        Set the table styles included within the ``<style>`` HTML element.

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

        overwrite : bool, default True
            Styles are replaced if `True`, or extended if `False`. CSS
            rules are preserved so most recent styles set will dominate
            if selectors intersect.

            .. versionadded:: 1.2.0

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.set_td_classes: Set the DataFrame of strings added to the ``class``
            attribute of ``<td>`` HTML elements.
        Styler.set_table_attributes: Set the table attributes added to the ``<table>``
            HTML element.

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
                    "props": maybe_convert_css_to_tuples(s["props"]),
                }
                for key, styles in table_styles.items()
                for s in styles
            ]
        else:
            table_styles = [
                {
                    "selector": s["selector"],
                    "props": maybe_convert_css_to_tuples(s["props"]),
                }
                for s in table_styles
            ]

        if not overwrite and self.table_styles is not None:
            self.table_styles.extend(table_styles)
        else:
            self.table_styles = table_styles
        return self

    def set_na_rep(self, na_rep: str) -> StylerRenderer:
        """
        Set the missing data representation on a ``Styler``.

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
        subset = non_reducing_slice(subset)
        hidden_df = self.data.loc[subset]
        hcols = self.columns.get_indexer_for(hidden_df.columns)
        # error: Incompatible types in assignment (expression has type
        # "ndarray", variable has type "Sequence[int]")
        self.hidden_columns = hcols  # type: ignore[assignment]
        return self

    # -----------------------------------------------------------------------
    # A collection of "builtin" styles
    # -----------------------------------------------------------------------

    def background_gradient(
        self,
        cmap="PuBu",
        low: float = 0,
        high: float = 0,
        axis: Axis | None = 0,
        subset=None,
        text_color_threshold: float = 0.408,
        vmin: float | None = None,
        vmax: float | None = None,
        gmap: Sequence | None = None,
    ) -> Styler:
        """
        Color the background in a gradient style.

        The background color is determined according
        to the data in each column, row or frame, or by a given
        gradient map. Requires matplotlib.

        Parameters
        ----------
        cmap : str or colormap
            Matplotlib colormap.
        low : float
            Compress the color range at the low end. This is a multiple of the data
            range to extend below the minimum; good values usually in [0, 1],
            defaults to 0.
        high : float
            Compress the color range at the high end. This is a multiple of the data
            range to extend above the maximum; good values usually in [0, 1],
            defaults to 0.
        axis : {0 or 'index', 1 or 'columns', None}, default 0
            Apply to each column (``axis=0`` or ``'index'``), to each row
            (``axis=1`` or ``'columns'``), or to the entire DataFrame at once
            with ``axis=None``.
        subset : IndexSlice
            A valid slice for ``data`` to limit the style application to.
        text_color_threshold : float or int
            Luminance threshold for determining text color in [0, 1]. Facilitates text
            visibility across varying background colors. All text is dark if 0, and
            light if 1, defaults to 0.408.

            .. versionadded:: 0.24.0

        vmin : float, optional
            Minimum data value that corresponds to colormap minimum value.
            If not specified the minimum value of the data (or gmap) will be used.

            .. versionadded:: 1.0.0

        vmax : float, optional
            Maximum data value that corresponds to colormap maximum value.
            If not specified the maximum value of the data (or gmap) will be used.

            .. versionadded:: 1.0.0

        gmap : array-like, optional
            Gradient map for determining the background colors. If not supplied
            will use the underlying data from rows, columns or frame. If given as an
            ndarray or list-like must be an identical shape to the underlying data
            considering ``axis`` and ``subset``. If given as DataFrame or Series must
            have same index and column labels considering ``axis`` and ``subset``.
            If supplied, ``vmin`` and ``vmax`` should be given relative to this
            gradient map.

            .. versionadded:: 1.3.0

        Returns
        -------
        self : Styler

        Notes
        -----
        When using ``low`` and ``high`` the range
        of the gradient, given by the data if ``gmap`` is not given or by ``gmap``,
        is extended at the low end effectively by
        `map.min - low * map.range` and at the high end by
        `map.max + high * map.range` before the colors are normalized and determined.

        If combining with ``vmin`` and ``vmax`` the `map.min`, `map.max` and
        `map.range` are replaced by values according to the values derived from
        ``vmin`` and ``vmax``.

        This method will preselect numeric columns and ignore non-numeric columns
        unless a ``gmap`` is supplied in which case no preselection occurs.

        Examples
        --------
        >>> df = pd.DataFrame({
        ...          'City': ['Stockholm', 'Oslo', 'Copenhagen'],
        ...          'Temp (c)': [21.6, 22.4, 24.5],
        ...          'Rain (mm)': [5.0, 13.3, 0.0],
        ...          'Wind (m/s)': [3.2, 3.1, 6.7]
        ... })

        Shading the values column-wise, with ``axis=0``, preselecting numeric columns

        >>> df.style.background_gradient(axis=0)

        .. figure:: ../../_static/style/bg_ax0.png

        Shading all values collectively using ``axis=None``

        >>> df.style.background_gradient(axis=None)

        .. figure:: ../../_static/style/bg_axNone.png

        Compress the color map from the both ``low`` and ``high`` ends

        >>> df.style.background_gradient(axis=None, low=0.75, high=1.0)

        .. figure:: ../../_static/style/bg_axNone_lowhigh.png

        Manually setting ``vmin`` and ``vmax`` gradient thresholds

        >>> df.style.background_gradient(axis=None, vmin=6.7, vmax=21.6)

        .. figure:: ../../_static/style/bg_axNone_vminvmax.png

        Setting a ``gmap`` and applying to all columns with another ``cmap``

        >>> df.style.background_gradient(axis=0, gmap=df['Temp (c)'], cmap='YlOrRd')

        .. figure:: ../../_static/style/bg_gmap.png

        Setting the gradient map for a dataframe (i.e. ``axis=None``), we need to
        explicitly state ``subset`` to match the ``gmap`` shape

        >>> gmap = np.array([[1,2,3], [2,3,4], [3,4,5]])
        >>> df.style.background_gradient(axis=None, gmap=gmap,
        ...     cmap='YlOrRd', subset=['Temp (c)', 'Rain (mm)', 'Wind (m/s)']
        ... )

        .. figure:: ../../_static/style/bg_axNone_gmap.png
        """
        if subset is None and gmap is None:
            subset = self.data.select_dtypes(include=np.number).columns

        self.apply(
            _background_gradient,
            cmap=cmap,
            subset=subset,
            axis=axis,
            low=low,
            high=high,
            text_color_threshold=text_color_threshold,
            vmin=vmin,
            vmax=vmax,
            gmap=gmap,
        )
        return self

    def set_properties(self, subset=None, **kwargs) -> Styler:
        """
        Set defined CSS-properties to each ``<td>`` HTML element within the given
        subset.

        Parameters
        ----------
        subset : IndexSlice
            A valid slice for ``data`` to limit the style application to.
        **kwargs : dict
            A dictionary of property, value pairs to be set for each cell.

        Returns
        -------
        self : Styler

        Notes
        -----
        This is a convenience methods which wraps the :meth:`Styler.applymap` calling a
        function returning the CSS-properties independently of the data.

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
        colors: list[str],
        width: float = 100,
        vmin: float | None = None,
        vmax: float | None = None,
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
            return DataFrame(
                [[css(x) for x in row] for row in normed],
                index=s.index,
                columns=s.columns,
            )

    def bar(
        self,
        subset=None,
        axis: Axis | None = 0,
        color="#d65f5f",
        width: float = 100,
        align: str = "left",
        vmin: float | None = None,
        vmax: float | None = None,
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
        subset: IndexLabel | None = None,
        props: str | None = None,
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
        Styler.highlight_between: Highlight a defined range with a style.
        Styler.highlight_quantile: Highlight values defined by a quantile with a style.
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
        subset: IndexLabel | None = None,
        color: str = "yellow",
        axis: Axis | None = 0,
        props: str | None = None,
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
        Styler.highlight_between: Highlight a defined range with a style.
        Styler.highlight_quantile: Highlight values defined by a quantile with a style.
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
        subset: IndexLabel | None = None,
        color: str = "yellow",
        axis: Axis | None = 0,
        props: str | None = None,
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
        Styler.highlight_between: Highlight a defined range with a style.
        Styler.highlight_quantile: Highlight values defined by a quantile with a style.
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

    def highlight_between(
        self,
        subset: IndexLabel | None = None,
        color: str = "yellow",
        axis: Axis | None = 0,
        left: Scalar | Sequence | None = None,
        right: Scalar | Sequence | None = None,
        inclusive: str = "both",
        props: str | None = None,
    ) -> Styler:
        """
        Highlight a defined range with a style.

        .. versionadded:: 1.3.0

        Parameters
        ----------
        subset : IndexSlice, default None
            A valid slice for ``data`` to limit the style application to.
        color : str, default 'yellow'
            Background color to use for highlighting.
        axis : {0 or 'index', 1 or 'columns', None}, default 0
            If ``left`` or ``right`` given as sequence, axis along which to apply those
            boundaries. See examples.
        left : scalar or datetime-like, or sequence or array-like, default None
            Left bound for defining the range.
        right : scalar or datetime-like, or sequence or array-like, default None
            Right bound for defining the range.
        inclusive : {'both', 'neither', 'left', 'right'}
            Identify whether bounds are closed or open.
        props : str, default None
            CSS properties to use for highlighting. If ``props`` is given, ``color``
            is not used.

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.highlight_null: Highlight missing values with a style.
        Styler.highlight_max: Highlight the maximum with a style.
        Styler.highlight_min: Highlight the minimum with a style.
        Styler.highlight_quantile: Highlight values defined by a quantile with a style.

        Notes
        -----
        If ``left`` is ``None`` only the right bound is applied.
        If ``right`` is ``None`` only the left bound is applied. If both are ``None``
        all values are highlighted.

        ``axis`` is only needed if ``left`` or ``right`` are provided as a sequence or
        an array-like object for aligning the shapes. If ``left`` and ``right`` are
        both scalars then all ``axis`` inputs will give the same result.

        This function only works with compatible ``dtypes``. For example a datetime-like
        region can only use equivalent datetime-like ``left`` and ``right`` arguments.
        Use ``subset`` to control regions which have multiple ``dtypes``.

        Examples
        --------
        Basic usage

        >>> df = pd.DataFrame({
        ...     'One': [1.2, 1.6, 1.5],
        ...     'Two': [2.9, 2.1, 2.5],
        ...     'Three': [3.1, 3.2, 3.8],
        ... })
        >>> df.style.highlight_between(left=2.1, right=2.9)

        .. figure:: ../../_static/style/hbetw_basic.png

        Using a range input sequnce along an ``axis``, in this case setting a ``left``
        and ``right`` for each column individually

        >>> df.style.highlight_between(left=[1.4, 2.4, 3.4], right=[1.6, 2.6, 3.6],
        ...     axis=1, color="#fffd75")

        .. figure:: ../../_static/style/hbetw_seq.png

        Using ``axis=None`` and providing the ``left`` argument as an array that
        matches the input DataFrame, with a constant ``right``

        >>> df.style.highlight_between(left=[[2,2,3],[2,2,3],[3,3,3]], right=3.5,
        ...     axis=None, color="#fffd75")

        .. figure:: ../../_static/style/hbetw_axNone.png

        Using ``props`` instead of default background coloring

        >>> df.style.highlight_between(left=1.5, right=3.5,
        ...     props='font-weight:bold;color:#e83e8c')

        .. figure:: ../../_static/style/hbetw_props.png
        """
        if props is None:
            props = f"background-color: {color};"
        return self.apply(
            _highlight_between,  # type: ignore[arg-type]
            axis=axis,
            subset=subset,
            props=props,
            left=left,
            right=right,
            inclusive=inclusive,
        )

    def highlight_quantile(
        self,
        subset: IndexLabel | None = None,
        color: str = "yellow",
        axis: Axis | None = 0,
        q_left: float = 0.0,
        q_right: float = 1.0,
        interpolation: str = "linear",
        inclusive: str = "both",
        props: str | None = None,
    ) -> Styler:
        """
        Highlight values defined by a quantile with a style.

        .. versionadded:: 1.3.0

        Parameters
        ----------
        subset : IndexSlice, default None
            A valid slice for ``data`` to limit the style application to.
        color : str, default 'yellow'
            Background color to use for highlighting
        axis : {0 or 'index', 1 or 'columns', None}, default 0
            Axis along which to determine and highlight quantiles. If ``None`` quantiles
            are measured over the entire DataFrame. See examples.
        q_left : float, default 0
            Left bound, in [0, q_right), for the target quantile range.
        q_right : float, default 1
            Right bound, in (q_left, 1], for the target quantile range.
        interpolation : {‘linear’, ‘lower’, ‘higher’, ‘midpoint’, ‘nearest’}
            Argument passed to ``Series.quantile`` or ``DataFrame.quantile`` for
            quantile estimation.
        inclusive : {'both', 'neither', 'left', 'right'}
            Identify whether quantile bounds are closed or open.
        props : str, default None
            CSS properties to use for highlighting. If ``props`` is given, ``color``
            is not used.

        Returns
        -------
        self : Styler

        See Also
        --------
        Styler.highlight_null: Highlight missing values with a style.
        Styler.highlight_max: Highlight the maximum with a style.
        Styler.highlight_min: Highlight the minimum with a style.
        Styler.highlight_between: Highlight a defined range with a style.

        Notes
        -----
        This function does not work with ``str`` dtypes.

        Examples
        --------
        Using ``axis=None`` and apply a quantile to all collective data

        >>> df = pd.DataFrame(np.arange(10).reshape(2,5) + 1)
        >>> df.style.highlight_quantile(axis=None, q_left=0.8, color="#fffd75")

        .. figure:: ../../_static/style/hq_axNone.png

        Or highlight quantiles row-wise or column-wise, in this case by row-wise

        >>> df.style.highlight_quantile(axis=1, q_left=0.8, color="#fffd75")

        .. figure:: ../../_static/style/hq_ax1.png

        Use ``props`` instead of default background coloring

        >>> df.style.highlight_quantile(axis=None, q_left=0.2, q_right=0.8,
        ...     props='font-weight:bold;color:#e83e8c')

        .. figure:: ../../_static/style/hq_props.png
        """
        subset_ = slice(None) if subset is None else subset
        subset_ = non_reducing_slice(subset_)
        data = self.data.loc[subset_]

        # after quantile is found along axis, e.g. along rows,
        # applying the calculated quantile to alternate axis, e.g. to each column
        kwargs = {"q": [q_left, q_right], "interpolation": interpolation}
        if axis in [0, "index"]:
            q = data.quantile(axis=axis, numeric_only=False, **kwargs)
            axis_apply: int | None = 1
        elif axis in [1, "columns"]:
            q = data.quantile(axis=axis, numeric_only=False, **kwargs)
            axis_apply = 0
        else:  # axis is None
            q = Series(data.to_numpy().ravel()).quantile(**kwargs)
            axis_apply = None

        if props is None:
            props = f"background-color: {color};"
        return self.apply(
            _highlight_between,  # type: ignore[arg-type]
            axis=axis_apply,
            subset=subset,
            props=props,
            left=q.iloc[0],
            right=q.iloc[1],
            inclusive=inclusive,
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
            template_html = env.get_template(name)

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
        Styler.apply : Apply a CSS-styling function column-wise, row-wise, or
            table-wise.

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


def _validate_apply_axis_arg(
    arg: FrameOrSeries | Sequence | np.ndarray,
    arg_name: str,
    dtype: Any | None,
    data: FrameOrSeries,
) -> np.ndarray:
    """
    For the apply-type methods, ``axis=None`` creates ``data`` as DataFrame, and for
    ``axis=[1,0]`` it creates a Series. Where ``arg`` is expected as an element
    of some operator with ``data`` we must make sure that the two are compatible shapes,
    or raise.

    Parameters
    ----------
    arg : sequence, Series or DataFrame
        the user input arg
    arg_name : string
        name of the arg for use in error messages
    dtype : numpy dtype, optional
        forced numpy dtype if given
    data : Series or DataFrame
        underling subset of Styler data on which operations are performed

    Returns
    -------
    ndarray
    """
    dtype = {"dtype": dtype} if dtype else {}
    # raise if input is wrong for axis:
    if isinstance(arg, Series) and isinstance(data, DataFrame):
        raise ValueError(
            f"'{arg_name}' is a Series but underlying data for operations "
            f"is a DataFrame since 'axis=None'"
        )
    elif isinstance(arg, DataFrame) and isinstance(data, Series):
        raise ValueError(
            f"'{arg_name}' is a DataFrame but underlying data for "
            f"operations is a Series with 'axis in [0,1]'"
        )
    elif isinstance(arg, (Series, DataFrame)):  # align indx / cols to data
        arg = arg.reindex_like(data, method=None).to_numpy(**dtype)
    else:
        arg = np.asarray(arg, **dtype)
        assert isinstance(arg, np.ndarray)  # mypy requirement
        if arg.shape != data.shape:  # check valid input
            raise ValueError(
                f"supplied '{arg_name}' is not correct shape for data over "
                f"selected 'axis': got {arg.shape}, "
                f"expected {data.shape}"
            )
    return arg


def _background_gradient(
    data,
    cmap="PuBu",
    low: float = 0,
    high: float = 0,
    text_color_threshold: float = 0.408,
    vmin: float | None = None,
    vmax: float | None = None,
    gmap: Sequence | np.ndarray | FrameOrSeries | None = None,
):
    """
    Color background in a range according to the data or a gradient map
    """
    if gmap is None:  # the data is used the gmap
        gmap = data.to_numpy(dtype=float)
    else:  # else validate gmap against the underlying data
        gmap = _validate_apply_axis_arg(gmap, "gmap", float, data)

    with _mpl(Styler.background_gradient) as (plt, colors):
        smin = np.nanmin(gmap) if vmin is None else vmin
        smax = np.nanmax(gmap) if vmax is None else vmax
        rng = smax - smin
        # extend lower / upper bounds, compresses color range
        norm = colors.Normalize(smin - (rng * low), smax + (rng * high))
        rgbas = plt.cm.get_cmap(cmap)(norm(gmap))

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

        if data.ndim == 1:
            return [css(rgba) for rgba in rgbas]
        else:
            return DataFrame(
                [[css(rgba) for rgba in row] for row in rgbas],
                index=data.index,
                columns=data.columns,
            )


def _highlight_between(
    data: FrameOrSeries,
    props: str,
    left: Scalar | Sequence | np.ndarray | FrameOrSeries | None = None,
    right: Scalar | Sequence | np.ndarray | FrameOrSeries | None = None,
    inclusive: bool | str = True,
) -> np.ndarray:
    """
    Return an array of css props based on condition of data values within given range.
    """
    if np.iterable(left) and not isinstance(left, str):
        left = _validate_apply_axis_arg(
            left, "left", None, data  # type: ignore[arg-type]
        )

    if np.iterable(right) and not isinstance(right, str):
        right = _validate_apply_axis_arg(
            right, "right", None, data  # type: ignore[arg-type]
        )

    # get ops with correct boundary attribution
    if inclusive == "both":
        ops = (operator.ge, operator.le)
    elif inclusive == "neither":
        ops = (operator.gt, operator.lt)
    elif inclusive == "left":
        ops = (operator.ge, operator.lt)
    elif inclusive == "right":
        ops = (operator.gt, operator.le)
    else:
        raise ValueError(
            f"'inclusive' values can be 'both', 'left', 'right', or 'neither' "
            f"got {inclusive}"
        )

    g_left = (
        ops[0](data, left)
        if left is not None
        else np.full(data.shape, True, dtype=bool)
    )
    l_right = (
        ops[1](data, right)
        if right is not None
        else np.full(data.shape, True, dtype=bool)
    )
    return np.where(g_left & l_right, props, "")
