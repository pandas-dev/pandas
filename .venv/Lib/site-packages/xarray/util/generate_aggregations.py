"""Generate module and stub file for arithmetic operators of various xarray classes.

For internal xarray development use only.

Usage:
    python xarray/util/generate_aggregations.py
    pytest --doctest-modules xarray/core/_aggregations.py --accept || true
    pytest --doctest-modules xarray/core/_aggregations.py

This requires [pytest-accept](https://github.com/max-sixty/pytest-accept).
The second run of pytest is deliberate, since the first will return an error
while replacing the doctests.

"""

import collections
import textwrap
from dataclasses import dataclass, field

MODULE_PREAMBLE = '''\
"""Mixin classes with reduction operations."""

# This file was generated using xarray.util.generate_aggregations. Do not edit manually.

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Callable

from xarray.core import duck_array_ops
from xarray.core.options import OPTIONS
from xarray.core.types import Dims, Self
from xarray.core.utils import contains_only_chunked_or_numpy, module_available

if TYPE_CHECKING:
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

flox_available = module_available("flox")
'''

NAMED_ARRAY_MODULE_PREAMBLE = '''\
"""Mixin classes with reduction operations."""
# This file was generated using xarray.util.generate_aggregations. Do not edit manually.

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Callable

from xarray.core import duck_array_ops
from xarray.core.types import Dims, Self
'''

AGGREGATIONS_PREAMBLE = """

class {obj}{cls}Aggregations:
    __slots__ = ()

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> Self:
        raise NotImplementedError()"""

NAMED_ARRAY_AGGREGATIONS_PREAMBLE = """

class {obj}{cls}Aggregations:
    __slots__ = ()

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> Self:
        raise NotImplementedError()"""


GROUPBY_PREAMBLE = """

class {obj}{cls}Aggregations:
    _obj: {obj}

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> {obj}:
        raise NotImplementedError()

    def _flox_reduce(
        self,
        dim: Dims,
        **kwargs: Any,
    ) -> {obj}:
        raise NotImplementedError()"""

RESAMPLE_PREAMBLE = """

class {obj}{cls}Aggregations:
    _obj: {obj}

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> {obj}:
        raise NotImplementedError()

    def _flox_reduce(
        self,
        dim: Dims,
        **kwargs: Any,
    ) -> {obj}:
        raise NotImplementedError()"""

TEMPLATE_REDUCTION_SIGNATURE = '''
    def {method}(
        self,
        dim: Dims = None,{kw_only}{extra_kwargs}{keep_attrs}
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this {obj}'s data by applying ``{method}`` along some dimension(s).

        Parameters
        ----------'''

TEMPLATE_REDUCTION_SIGNATURE_GROUPBY = '''
    def {method}(
        self,
        dim: Dims = None,
        *,{extra_kwargs}
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> {obj}:
        """
        Reduce this {obj}'s data by applying ``{method}`` along some dimension(s).

        Parameters
        ----------'''

TEMPLATE_RETURNS = """
        Returns
        -------
        reduced : {obj}
            New {obj} with ``{method}`` applied to its data and the
            indicated dimension(s) removed"""

TEMPLATE_SEE_ALSO = """
        See Also
        --------
{see_also_methods}
        :ref:`{docref}`
            User guide on {docref_description}."""

TEMPLATE_NOTES = """
        Notes
        -----
{notes}"""

_DIM_DOCSTRING = """dim : str, Iterable of Hashable, "..." or None, default: None
    Name of dimension[s] along which to apply ``{method}``. For e.g. ``dim="x"``
    or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions."""

_DIM_DOCSTRING_GROUPBY = """dim : str, Iterable of Hashable, "..." or None, default: None
    Name of dimension[s] along which to apply ``{method}``. For e.g. ``dim="x"``
    or ``dim=["x", "y"]``. If None, will reduce over the {cls} dimensions.
    If "...", will reduce over all dimensions."""

_SKIPNA_DOCSTRING = """skipna : bool or None, optional
    If True, skip missing values (as marked by NaN). By default, only
    skips missing values for float dtypes; other dtypes either do not
    have a sentinel missing value (int) or ``skipna=True`` has not been
    implemented (object, datetime64 or timedelta64)."""

_MINCOUNT_DOCSTRING = """min_count : int or None, optional
    The required number of valid values to perform the operation. If
    fewer than min_count non-NA values are present the result will be
    NA. Only used if skipna is set to True or defaults to True for the
    array's dtype. Changed in version 0.17.0: if specified on an integer
    array and skipna=True, the result will be a float array."""

_DDOF_DOCSTRING = """ddof : int, default: 0
    “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
    where ``N`` represents the number of elements."""

_KEEP_ATTRS_DOCSTRING = """keep_attrs : bool or None, optional
    If True, ``attrs`` will be copied from the original
    object to the new one.  If False, the new object will be
    returned without attributes."""

_KWARGS_DOCSTRING = """**kwargs : Any
    Additional keyword arguments passed on to the appropriate array
    function for calculating ``{method}`` on this object's data.
    These could include dask-specific kwargs like ``split_every``."""

_NUMERIC_ONLY_NOTES = "Non-numeric variables will be removed prior to reducing."

_FLOX_NOTES_TEMPLATE = """Use the ``flox`` package to significantly speed up {kind} computations,
especially with dask arrays. Xarray will use flox by default if installed.
Pass flox-specific keyword arguments in ``**kwargs``.
See the `flox documentation <https://flox.readthedocs.io>`_ for more."""
_FLOX_GROUPBY_NOTES = _FLOX_NOTES_TEMPLATE.format(kind="groupby")
_FLOX_RESAMPLE_NOTES = _FLOX_NOTES_TEMPLATE.format(kind="resampling")

ExtraKwarg = collections.namedtuple("ExtraKwarg", "docs kwarg call example")
skipna = ExtraKwarg(
    docs=_SKIPNA_DOCSTRING,
    kwarg="skipna: bool | None = None,",
    call="skipna=skipna,",
    example="""\n
        Use ``skipna`` to control whether NaNs are ignored.

        >>> {calculation}(skipna=False)""",
)
min_count = ExtraKwarg(
    docs=_MINCOUNT_DOCSTRING,
    kwarg="min_count: int | None = None,",
    call="min_count=min_count,",
    example="""\n
        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> {calculation}(skipna=True, min_count=2)""",
)
ddof = ExtraKwarg(
    docs=_DDOF_DOCSTRING,
    kwarg="ddof: int = 0,",
    call="ddof=ddof,",
    example="""\n
        Specify ``ddof=1`` for an unbiased estimate.

        >>> {calculation}(skipna=True, ddof=1)""",
)


@dataclass
class DataStructure:
    name: str
    create_example: str
    example_var_name: str
    numeric_only: bool = False
    see_also_modules: tuple[str] = tuple


class Method:
    def __init__(
        self,
        name,
        bool_reduce=False,
        extra_kwargs=tuple(),
        numeric_only=False,
        see_also_modules=("numpy", "dask.array"),
        min_flox_version=None,
    ):
        self.name = name
        self.extra_kwargs = extra_kwargs
        self.numeric_only = numeric_only
        self.see_also_modules = see_also_modules
        self.min_flox_version = min_flox_version
        if bool_reduce:
            self.array_method = f"array_{name}"
            self.np_example_array = """
        ...     np.array([True, True, True, True, True, False], dtype=bool)"""

        else:
            self.array_method = name
            self.np_example_array = """
        ...     np.array([1, 2, 3, 0, 2, np.nan])"""


@dataclass
class AggregationGenerator:
    _dim_docstring = _DIM_DOCSTRING
    _template_signature = TEMPLATE_REDUCTION_SIGNATURE

    cls: str
    datastructure: DataStructure
    methods: tuple[Method, ...]
    docref: str
    docref_description: str
    example_call_preamble: str
    definition_preamble: str
    has_keep_attrs: bool = True
    notes: str = ""
    preamble: str = field(init=False)

    def __post_init__(self):
        self.preamble = self.definition_preamble.format(
            obj=self.datastructure.name, cls=self.cls
        )

    def generate_methods(self):
        yield [self.preamble]
        for method in self.methods:
            yield self.generate_method(method)

    def generate_method(self, method):
        has_kw_only = method.extra_kwargs or self.has_keep_attrs

        template_kwargs = dict(
            obj=self.datastructure.name,
            method=method.name,
            keep_attrs=(
                "\n        keep_attrs: bool | None = None,"
                if self.has_keep_attrs
                else ""
            ),
            kw_only="\n        *," if has_kw_only else "",
        )

        if method.extra_kwargs:
            extra_kwargs = "\n        " + "\n        ".join(
                [kwarg.kwarg for kwarg in method.extra_kwargs if kwarg.kwarg]
            )
        else:
            extra_kwargs = ""

        yield self._template_signature.format(
            **template_kwargs,
            extra_kwargs=extra_kwargs,
        )

        for text in [
            self._dim_docstring.format(method=method.name, cls=self.cls),
            *(kwarg.docs for kwarg in method.extra_kwargs if kwarg.docs),
            _KEEP_ATTRS_DOCSTRING if self.has_keep_attrs else None,
            _KWARGS_DOCSTRING.format(method=method.name),
        ]:
            if text:
                yield textwrap.indent(text, 8 * " ")

        yield TEMPLATE_RETURNS.format(**template_kwargs)

        # we want Dataset.count to refer to DataArray.count
        # but we also want DatasetGroupBy.count to refer to Dataset.count
        # The generic aggregations have self.cls == ''
        others = (
            self.datastructure.see_also_modules
            if self.cls == ""
            else (self.datastructure.name,)
        )
        see_also_methods = "\n".join(
            " " * 8 + f"{mod}.{method.name}"
            for mod in (method.see_also_modules + others)
        )
        # Fixes broken links mentioned in #8055
        yield TEMPLATE_SEE_ALSO.format(
            **template_kwargs,
            docref=self.docref,
            docref_description=self.docref_description,
            see_also_methods=see_also_methods,
        )

        notes = self.notes
        if method.numeric_only:
            if notes != "":
                notes += "\n\n"
            notes += _NUMERIC_ONLY_NOTES

        if notes != "":
            yield TEMPLATE_NOTES.format(notes=textwrap.indent(notes, 8 * " "))

        yield textwrap.indent(self.generate_example(method=method), "")
        yield '        """'

        yield self.generate_code(method, self.has_keep_attrs)

    def generate_example(self, method):
        created = self.datastructure.create_example.format(
            example_array=method.np_example_array
        )
        calculation = f"{self.datastructure.example_var_name}{self.example_call_preamble}.{method.name}"
        if method.extra_kwargs:
            extra_examples = "".join(
                kwarg.example for kwarg in method.extra_kwargs if kwarg.example
            ).format(calculation=calculation, method=method.name)
        else:
            extra_examples = ""

        return f"""
        Examples
        --------{created}
        >>> {self.datastructure.example_var_name}

        >>> {calculation}(){extra_examples}"""


class GroupByAggregationGenerator(AggregationGenerator):
    _dim_docstring = _DIM_DOCSTRING_GROUPBY
    _template_signature = TEMPLATE_REDUCTION_SIGNATURE_GROUPBY

    def generate_code(self, method, has_keep_attrs):
        extra_kwargs = [kwarg.call for kwarg in method.extra_kwargs if kwarg.call]

        if self.datastructure.numeric_only:
            extra_kwargs.append(f"numeric_only={method.numeric_only},")

        # median isn't enabled yet, because it would break if a single group was present in multiple
        # chunks. The non-flox code path will just rechunk every group to a single chunk and execute the median
        method_is_not_flox_supported = method.name in ("median", "cumsum", "cumprod")
        if method_is_not_flox_supported:
            indent = 12
        else:
            indent = 16

        if extra_kwargs:
            extra_kwargs = textwrap.indent("\n" + "\n".join(extra_kwargs), indent * " ")
        else:
            extra_kwargs = ""

        if method_is_not_flox_supported:
            return f"""\
        return self.reduce(
            duck_array_ops.{method.array_method},
            dim=dim,{extra_kwargs}
            keep_attrs=keep_attrs,
            **kwargs,
        )"""

        min_version_check = f"""
            and module_available("flox", minversion="{method.min_flox_version}")"""

        return (
            """\
        if (
            flox_available
            and OPTIONS["use_flox"]"""
            + (min_version_check if method.min_flox_version is not None else "")
            + f"""
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="{method.name}",
                dim=dim,{extra_kwargs}
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.{method.array_method},
                dim=dim,{extra_kwargs}
                keep_attrs=keep_attrs,
                **kwargs,
            )"""
        )


class GenericAggregationGenerator(AggregationGenerator):
    def generate_code(self, method, has_keep_attrs):
        extra_kwargs = [kwarg.call for kwarg in method.extra_kwargs if kwarg.call]

        if self.datastructure.numeric_only:
            extra_kwargs.append(f"numeric_only={method.numeric_only},")

        if extra_kwargs:
            extra_kwargs = textwrap.indent("\n" + "\n".join(extra_kwargs), 12 * " ")
        else:
            extra_kwargs = ""
        keep_attrs = (
            "\n" + 12 * " " + "keep_attrs=keep_attrs," if has_keep_attrs else ""
        )
        return f"""\
        return self.reduce(
            duck_array_ops.{method.array_method},
            dim=dim,{extra_kwargs}{keep_attrs}
            **kwargs,
        )"""


AGGREGATION_METHODS = (
    # Reductions:
    Method("count", see_also_modules=("pandas.DataFrame", "dask.dataframe.DataFrame")),
    Method("all", bool_reduce=True),
    Method("any", bool_reduce=True),
    Method("max", extra_kwargs=(skipna,)),
    Method("min", extra_kwargs=(skipna,)),
    Method("mean", extra_kwargs=(skipna,), numeric_only=True),
    Method("prod", extra_kwargs=(skipna, min_count), numeric_only=True),
    Method("sum", extra_kwargs=(skipna, min_count), numeric_only=True),
    Method("std", extra_kwargs=(skipna, ddof), numeric_only=True),
    Method("var", extra_kwargs=(skipna, ddof), numeric_only=True),
    Method(
        "median", extra_kwargs=(skipna,), numeric_only=True, min_flox_version="0.9.2"
    ),
    # Cumulatives:
    Method("cumsum", extra_kwargs=(skipna,), numeric_only=True),
    Method("cumprod", extra_kwargs=(skipna,), numeric_only=True),
)


DATASET_OBJECT = DataStructure(
    name="Dataset",
    create_example="""
        >>> da = xr.DataArray({example_array},
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))""",
    example_var_name="ds",
    numeric_only=True,
    see_also_modules=("DataArray",),
)
DATAARRAY_OBJECT = DataStructure(
    name="DataArray",
    create_example="""
        >>> da = xr.DataArray({example_array},
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )""",
    example_var_name="da",
    numeric_only=False,
    see_also_modules=("Dataset",),
)
DATASET_GENERATOR = GenericAggregationGenerator(
    cls="",
    datastructure=DATASET_OBJECT,
    methods=AGGREGATION_METHODS,
    docref="agg",
    docref_description="reduction or aggregation operations",
    example_call_preamble="",
    definition_preamble=AGGREGATIONS_PREAMBLE,
)
DATAARRAY_GENERATOR = GenericAggregationGenerator(
    cls="",
    datastructure=DATAARRAY_OBJECT,
    methods=AGGREGATION_METHODS,
    docref="agg",
    docref_description="reduction or aggregation operations",
    example_call_preamble="",
    definition_preamble=AGGREGATIONS_PREAMBLE,
)
DATAARRAY_GROUPBY_GENERATOR = GroupByAggregationGenerator(
    cls="GroupBy",
    datastructure=DATAARRAY_OBJECT,
    methods=AGGREGATION_METHODS,
    docref="groupby",
    docref_description="groupby operations",
    example_call_preamble='.groupby("labels")',
    definition_preamble=GROUPBY_PREAMBLE,
    notes=_FLOX_GROUPBY_NOTES,
)
DATAARRAY_RESAMPLE_GENERATOR = GroupByAggregationGenerator(
    cls="Resample",
    datastructure=DATAARRAY_OBJECT,
    methods=AGGREGATION_METHODS,
    docref="resampling",
    docref_description="resampling operations",
    example_call_preamble='.resample(time="3ME")',
    definition_preamble=RESAMPLE_PREAMBLE,
    notes=_FLOX_RESAMPLE_NOTES,
)
DATASET_GROUPBY_GENERATOR = GroupByAggregationGenerator(
    cls="GroupBy",
    datastructure=DATASET_OBJECT,
    methods=AGGREGATION_METHODS,
    docref="groupby",
    docref_description="groupby operations",
    example_call_preamble='.groupby("labels")',
    definition_preamble=GROUPBY_PREAMBLE,
    notes=_FLOX_GROUPBY_NOTES,
)
DATASET_RESAMPLE_GENERATOR = GroupByAggregationGenerator(
    cls="Resample",
    datastructure=DATASET_OBJECT,
    methods=AGGREGATION_METHODS,
    docref="resampling",
    docref_description="resampling operations",
    example_call_preamble='.resample(time="3ME")',
    definition_preamble=RESAMPLE_PREAMBLE,
    notes=_FLOX_RESAMPLE_NOTES,
)

NAMED_ARRAY_OBJECT = DataStructure(
    name="NamedArray",
    create_example="""
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray(
        ...     "x",{example_array},
        ... )""",
    example_var_name="na",
    numeric_only=False,
    see_also_modules=("Dataset", "DataArray"),
)

NAMED_ARRAY_GENERATOR = GenericAggregationGenerator(
    cls="",
    datastructure=NAMED_ARRAY_OBJECT,
    methods=AGGREGATION_METHODS,
    docref="agg",
    docref_description="reduction or aggregation operations",
    example_call_preamble="",
    definition_preamble=NAMED_ARRAY_AGGREGATIONS_PREAMBLE,
    has_keep_attrs=False,
)


def write_methods(filepath, generators, preamble):
    with open(filepath, mode="w", encoding="utf-8") as f:
        f.write(preamble)
        for gen in generators:
            for lines in gen.generate_methods():
                for line in lines:
                    f.write(line + "\n")


if __name__ == "__main__":
    import os
    from pathlib import Path

    p = Path(os.getcwd())
    write_methods(
        filepath=p.parent / "xarray" / "xarray" / "core" / "_aggregations.py",
        generators=[
            DATASET_GENERATOR,
            DATAARRAY_GENERATOR,
            DATASET_GROUPBY_GENERATOR,
            DATASET_RESAMPLE_GENERATOR,
            DATAARRAY_GROUPBY_GENERATOR,
            DATAARRAY_RESAMPLE_GENERATOR,
        ],
        preamble=MODULE_PREAMBLE,
    )
    write_methods(
        filepath=p.parent / "xarray" / "xarray" / "namedarray" / "_aggregations.py",
        generators=[NAMED_ARRAY_GENERATOR],
        preamble=NAMED_ARRAY_MODULE_PREAMBLE,
    )
    # filepath = p.parent / "core" / "_aggregations.py"  # Run from script location
