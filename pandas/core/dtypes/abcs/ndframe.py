from __future__ import annotations

from typing import TYPE_CHECKING

from pandas._libs import lib
from pandas.util._decorators import rewrite_axis_style_signature

from pandas.core import indexing
from pandas.core.base import PandasObject, SelectionMixin

if TYPE_CHECKING:
    from abc import ABC
    from datetime import timedelta
    import pickle
    from typing import (
        Any,
        Callable,
        Dict,
        FrozenSet,
        Hashable,
        List,
        Mapping,
        Optional,
        Sequence,
        Set,
        Tuple,
        Type,
        Union,
    )

    import numpy as np

    from pandas._libs.tslibs import BaseOffset
    from pandas._typing import (
        Axis,
        CompressionOptions,
        FilePathOrBuffer,
        FrameOrSeries,
        IndexKeyFunc,
        IndexLabel,
        JSONSerializable,
        Label,
        Level,
        Renamer,
        StorageOptions,
        TimedeltaConvertibleTypes,
        TimestampConvertibleTypes,
        ValueKeyFunc,
    )
    from pandas.errors import AbstractMethodError

    from pandas.core.dtypes.generic import ABCSeries

    from pandas.core.flags import Flags
    from pandas.core.indexes.api import Index, MultiIndex
    from pandas.core.internals import BlockManager
    from pandas.core.resample import Resampler
    from pandas.core.series import Series
    from pandas.core.window import Expanding, ExponentialMovingWindow
    from pandas.core.window.indexers import BaseIndexer

bool_t = bool  # Need alias because NDFrame has def bool:


class ABCNDFrame(ABC, PandasObject, SelectionMixin, indexing.IndexingMixin):
    _internal_names: List[str] = [
        "_mgr",
        "_cacher",
        "_item_cache",
        "_cache",
        "_is_copy",
        "_subtyp",
        "_name",
        "_index",
        "_default_kind",
        "_default_fill_value",
        "_metadata",
        "__array_struct__",
        "__array_interface__",
        "_flags",
    ]
    _internal_names_set: Set[str] = set(_internal_names)
    _accessors: Set[str] = set()
    _hidden_attrs: FrozenSet[str] = frozenset(["get_values", "tshift"])
    _metadata: List[str] = []
    _is_copy = None
    _mgr: BlockManager
    _attrs: Dict[Optional[Hashable], Any]
    _typ: str

    # ----------------------------------------------------------------------
    # Constructors

    @classmethod
    def _init_mgr(cls, mgr, axes, dtype=None, copy: bool = False) -> BlockManager:
        raise NotImplementedError

    # ----------------------------------------------------------------------
    # attrs and flags

    @property
    def attrs(self) -> Dict[Optional[Hashable], Any]:
        raise NotImplementedError

    @attrs.setter
    def attrs(self, value: Mapping[Optional[Hashable], Any]) -> None:
        raise NotImplementedError

    @property
    def flags(self) -> Flags:
        raise NotImplementedError

    def set_flags(
        self: FrameOrSeries,
        *,
        copy: bool = False,
        allows_duplicate_labels: Optional[bool] = None,
    ):
        raise NotImplementedError

    @classmethod
    def _validate_dtype(cls, dtype):
        raise NotImplementedError

    # ----------------------------------------------------------------------
    # Construction

    @property
    def _constructor(self: FrameOrSeries) -> Type[FrameOrSeries]:
        raise AbstractMethodError(self)

    @property
    def _constructor_sliced(self):
        raise AbstractMethodError(self)

    @property
    def _constructor_expanddim(self):
        raise NotImplementedError

    # ----------------------------------------------------------------------
    # Internals

    @property
    def _data(self):
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Axis
    _stat_axis_number = 0
    _stat_axis_name = "index"
    _ix = None
    _AXIS_ORDERS: List[str]
    _AXIS_TO_AXIS_NUMBER: Dict[Axis, int] = {0: 0, "index": 0, "rows": 0}
    _AXIS_REVERSED: bool
    _info_axis_number: int
    _info_axis_name: str
    _AXIS_LEN: int

    @property
    def _AXIS_NUMBERS(self) -> Dict[str, int]:
        raise AbstractMethodError(self)

    @property
    def _AXIS_NAMES(self) -> Dict[int, str]:
        raise AbstractMethodError(self)

    def _construct_axes_dict(self, axes=None, **kwargs):
        raise AbstractMethodError(self)

    @classmethod
    def _construct_axes_from_arguments(
        cls, args, kwargs, require_all: bool = False, sentinel=None
    ):
        raise AbstractMethodError(cls)

    @classmethod
    def _get_axis_number(cls, axis: Axis) -> int:
        raise AbstractMethodError(cls)

    @classmethod
    def _get_axis_name(cls, axis: Axis) -> str:
        raise AbstractMethodError(cls)

    def _get_axis(self, axis: Axis) -> Index:
        raise AbstractMethodError(self)

    @classmethod
    def _get_block_manager_axis(cls, axis: Axis) -> int:
        raise AbstractMethodError(cls)

    def _get_axis_resolvers(self, axis: str) -> Dict[str, Union[Series, MultiIndex]]:
        raise AbstractMethodError(self)

    def _get_index_resolvers(self) -> Dict[str, Union[Series, MultiIndex]]:
        raise AbstractMethodError(self)

    def _get_cleaned_column_resolvers(self) -> Dict[str, ABCSeries]:
        raise AbstractMethodError(self)

    @property
    def _info_axis(self) -> Index:
        raise AbstractMethodError(self)

    @property
    def _stat_axis(self) -> Index:
        raise AbstractMethodError(self)

    @property
    def shape(self) -> Tuple[int, ...]:
        raise AbstractMethodError(self)

    @property
    def axes(self) -> List[Index]:
        raise AbstractMethodError(self)

    @property
    def ndim(self) -> int:
        raise AbstractMethodError(self)

    @property
    def size(self) -> int:
        raise AbstractMethodError(self)

    @property
    def _selected_obj(self: FrameOrSeries) -> FrameOrSeries:
        raise AbstractMethodError(self)

    @property
    def _obj_with_exclusions(self: FrameOrSeries) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def set_axis(self, labels, axis: Axis = 0, inplace: bool = False):
        raise AbstractMethodError(self)

    def _set_axis_nocheck(self, labels, axis: Axis, inplace: bool):
        raise AbstractMethodError(self)

    def _set_axis(self, axis: int, labels: Index) -> None:
        raise AbstractMethodError(self)

    def swapaxes(self: FrameOrSeries, axis1, axis2, copy=True) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def droplevel(self: FrameOrSeries, level, axis=0) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def pop(self, item: Label) -> Union[Series, Any]:
        raise AbstractMethodError(self)

    def squeeze(self, axis=None):
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Rename

    def rename(
        self: FrameOrSeries,
        mapper: Optional[Renamer] = None,
        *,
        index: Optional[Renamer] = None,
        columns: Optional[Renamer] = None,
        axis: Optional[Axis] = None,
        copy: bool = True,
        inplace: bool = False,
        level: Optional[Level] = None,
        errors: str = "ignore",
    ) -> Optional[FrameOrSeries]:
        raise AbstractMethodError(self)

    @rewrite_axis_style_signature("mapper", [("copy", True), ("inplace", False)])
    def rename_axis(self, mapper=lib.no_default, **kwargs):
        raise AbstractMethodError(self)

    def _set_axis_name(self, name, axis=0, inplace=False):
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Comparison Methods

    def _indexed_same(self, other) -> bool:
        raise AbstractMethodError(self)

    def equals(self, other: object) -> bool:
        raise AbstractMethodError(self)

    # -------------------------------------------------------------------------
    # Unary Methods

    def __neg__(self):
        raise AbstractMethodError(self)

    def __pos__(self):
        raise AbstractMethodError(self)

    def __invert__(self):
        raise AbstractMethodError(self)

    def __nonzero__(self):
        raise AbstractMethodError(self)

    __bool__ = __nonzero__

    def bool(self):
        raise AbstractMethodError(self)

    def __abs__(self: FrameOrSeries) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def __round__(self: FrameOrSeries, decimals: int = 0) -> FrameOrSeries:
        raise AbstractMethodError(self)

    # -------------------------------------------------------------------------
    # Label or Level Combination Helpers
    #
    # A collection of helper methods for DataFrame/Series operations that
    # accept a combination of column/index labels and levels.  All such
    # operations should utilize/extend these methods when possible so that we
    # have consistent precedence and validation logic throughout the library.

    def _is_level_reference(self, key, axis=0):
        raise AbstractMethodError(self)

    def _is_label_reference(self, key, axis=0) -> bool_t:
        raise AbstractMethodError(self)

    def _is_label_or_level_reference(self, key: str, axis: int = 0) -> bool_t:
        raise AbstractMethodError(self)

    def _check_label_or_level_ambiguity(self, key, axis: int = 0) -> None:
        raise AbstractMethodError(self)

    def _get_label_or_level_values(self, key: str, axis: int = 0) -> np.ndarray:
        raise AbstractMethodError(self)

    def _drop_labels_or_levels(self, keys, axis: int = 0):
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Iteration

    def __hash__(self) -> int:
        raise AbstractMethodError(self)

    def __iter__(self):
        raise AbstractMethodError(self)

    # can we get a better explanation of this?
    def keys(self):
        raise AbstractMethodError(self)

    def items(self):
        raise AbstractMethodError(self)

    def iteritems(self):
        raise AbstractMethodError(self)

    def __len__(self) -> int:
        raise AbstractMethodError(self)

    def __contains__(self, key) -> bool_t:
        raise AbstractMethodError(self)

    @property
    def empty(self) -> bool_t:
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Array Interface

    # This is also set in IndexOpsMixin
    # GH#23114 Ensure ndarray.__op__(DataFrame) returns NotImplemented
    __array_priority__ = 1000

    def __array__(self, dtype=None) -> np.ndarray:
        raise AbstractMethodError(self)

    def __array_wrap__(
        self,
        result: np.ndarray,
        context: Optional[Tuple[Callable, Tuple[Any, ...], int]] = None,
    ):
        raise AbstractMethodError(self)

    # ideally we would define this to avoid the getattr checks, but
    # is slower
    # @property
    # def __array_interface__(self):
    #    """ provide numpy array interface method """
    #    values = self.values
    #    return dict(typestr=values.dtype.str,shape=values.shape,data=values)

    # ----------------------------------------------------------------------
    # Picklability

    def __getstate__(self) -> Dict[str, Any]:
        raise AbstractMethodError(self)

    def __setstate__(self, state):
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Rendering Methods

    def __repr__(self) -> str:
        raise AbstractMethodError(self)

    def _repr_latex_(self):
        raise AbstractMethodError(self)

    def _repr_data_resource_(self):
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # I/O Methods

    def to_excel(
        self,
        excel_writer,
        sheet_name="Sheet1",
        na_rep="",
        float_format=None,
        columns=None,
        header=True,
        index=True,
        index_label=None,
        startrow=0,
        startcol=0,
        engine=None,
        merge_cells=True,
        encoding=None,
        inf_rep="inf",
        verbose=True,
        freeze_panes=None,
    ) -> None:
        raise AbstractMethodError(self)

    def to_json(
        self,
        path_or_buf: Optional[FilePathOrBuffer] = None,
        orient: Optional[str] = None,
        date_format: Optional[str] = None,
        double_precision: int = 10,
        force_ascii: bool_t = True,
        date_unit: str = "ms",
        default_handler: Optional[Callable[[Any], JSONSerializable]] = None,
        lines: bool_t = False,
        compression: CompressionOptions = "infer",
        index: bool_t = True,
        indent: Optional[int] = None,
        storage_options: StorageOptions = None,
    ) -> Optional[str]:
        raise AbstractMethodError(self)

    def to_hdf(
        self,
        path_or_buf,
        key: str,
        mode: str = "a",
        complevel: Optional[int] = None,
        complib: Optional[str] = None,
        append: bool_t = False,
        format: Optional[str] = None,
        index: bool_t = True,
        min_itemsize: Optional[Union[int, Dict[str, int]]] = None,
        nan_rep=None,
        dropna: Optional[bool_t] = None,
        data_columns: Optional[Union[bool_t, List[str]]] = None,
        errors: str = "strict",
        encoding: str = "UTF-8",
    ) -> None:
        raise AbstractMethodError(self)

    def to_sql(
        self,
        name: str,
        con,
        schema=None,
        if_exists: str = "fail",
        index: bool_t = True,
        index_label=None,
        chunksize=None,
        dtype=None,
        method=None,
    ) -> None:
        raise AbstractMethodError(self)

    def to_pickle(
        self,
        path,
        compression: CompressionOptions = "infer",
        protocol: int = pickle.HIGHEST_PROTOCOL,
        storage_options: StorageOptions = None,
    ) -> None:
        raise AbstractMethodError(self)

    def to_clipboard(
        self, excel: bool_t = True, sep: Optional[str] = None, **kwargs
    ) -> None:
        raise AbstractMethodError(self)

    def to_xarray(self):
        raise AbstractMethodError(self)

    def to_latex(
        self,
        buf=None,
        columns=None,
        col_space=None,
        header=True,
        index=True,
        na_rep="NaN",
        formatters=None,
        float_format=None,
        sparsify=None,
        index_names=True,
        bold_rows=False,
        column_format=None,
        longtable=None,
        escape=None,
        encoding=None,
        decimal=".",
        multicolumn=None,
        multicolumn_format=None,
        multirow=None,
        caption=None,
        label=None,
        position=None,
    ):
        raise AbstractMethodError(self)

    def to_csv(
        self,
        path_or_buf: Optional[FilePathOrBuffer] = None,
        sep: str = ",",
        na_rep: str = "",
        float_format: Optional[str] = None,
        columns: Optional[Sequence[Label]] = None,
        header: Union[bool_t, List[str]] = True,
        index: bool_t = True,
        index_label: Optional[IndexLabel] = None,
        mode: str = "w",
        encoding: Optional[str] = None,
        compression: CompressionOptions = "infer",
        quoting: Optional[int] = None,
        quotechar: str = '"',
        line_terminator: Optional[str] = None,
        chunksize: Optional[int] = None,
        date_format: Optional[str] = None,
        doublequote: bool_t = True,
        escapechar: Optional[str] = None,
        decimal: str = ".",
        errors: str = "strict",
        storage_options: StorageOptions = None,
    ) -> Optional[str]:
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Lookup Caching

    def _set_as_cached(self, item, cacher) -> None:
        raise AbstractMethodError(self)

    def _reset_cacher(self) -> None:
        raise AbstractMethodError(self)

    def _maybe_cache_changed(self, item, value) -> None:
        raise AbstractMethodError(self)

    @property
    def _is_cached(self) -> bool_t:
        raise AbstractMethodError(self)

    def _get_cacher(self):
        raise AbstractMethodError(self)

    def _maybe_update_cacher(
        self, clear: bool_t = False, verify_is_copy: bool_t = True
    ) -> None:
        raise AbstractMethodError(self)

    def _clear_item_cache(self) -> None:
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Indexing Methods

    def take(
        self: FrameOrSeries, indices, axis=0, is_copy: Optional[bool_t] = None, **kwargs
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def _take_with_is_copy(self: FrameOrSeries, indices, axis=0) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def xs(self, key, axis=0, level=None, drop_level: bool_t = True):
        raise AbstractMethodError(self)

    def __getitem__(self, item):
        raise AbstractMethodError(self)

    def _get_item_cache(self, item):
        raise AbstractMethodError(self)

    def _slice(self: FrameOrSeries, slobj: slice, axis=0) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def _iset_item(self, loc: int, value) -> None:
        raise AbstractMethodError(self)

    def _set_item(self, key, value) -> None:
        raise AbstractMethodError(self)

    def _set_is_copy(self, ref, copy: bool_t = True) -> None:
        raise AbstractMethodError(self)

    def _check_is_chained_assignment_possible(self) -> bool_t:
        raise AbstractMethodError(self)

    def _check_setitem_copy(self, stacklevel=4, t="setting", force=False):
        raise AbstractMethodError(self)

    def __delitem__(self, key) -> None:
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Unsorted

    def _check_inplace_and_allows_duplicate_labels(self, inplace):
        raise AbstractMethodError(self)

    def get(self, key, default=None):
        raise AbstractMethodError(self)

    @property
    def _is_view(self) -> bool_t:
        raise AbstractMethodError(self)

    def reindex_like(
        self: FrameOrSeries,
        other,
        method: Optional[str] = None,
        copy: bool_t = True,
        limit=None,
        tolerance=None,
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def drop(
        self,
        labels=None,
        axis=0,
        index=None,
        columns=None,
        level=None,
        inplace: bool_t = False,
        errors: str = "raise",
    ):
        raise AbstractMethodError(self)

    def _drop_axis(
        self: FrameOrSeries, labels, axis, level=None, errors: str = "raise"
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def _update_inplace(self, result, verify_is_copy: bool_t = True) -> None:
        raise AbstractMethodError(self)

    def add_prefix(self: FrameOrSeries, prefix: str) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def add_suffix(self: FrameOrSeries, suffix: str) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def sort_values(
        self,
        axis=0,
        ascending=True,
        inplace: bool_t = False,
        kind: str = "quicksort",
        na_position: str = "last",
        ignore_index: bool_t = False,
        key: ValueKeyFunc = None,
    ):
        raise AbstractMethodError(self)

    def sort_index(
        self,
        axis=0,
        level=None,
        ascending: bool_t = True,
        inplace: bool_t = False,
        kind: str = "quicksort",
        na_position: str = "last",
        sort_remaining: bool_t = True,
        ignore_index: bool_t = False,
        key: IndexKeyFunc = None,
    ):
        raise AbstractMethodError(self)

    def reindex(self: FrameOrSeries, *args, **kwargs) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def _reindex_axes(
        self: FrameOrSeries, axes, level, limit, tolerance, method, fill_value, copy
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def _needs_reindex_multi(self, axes, method, level) -> bool_t:
        raise AbstractMethodError(self)

    def _reindex_multi(self, axes, copy, fill_value):
        raise AbstractMethodError(self)

    def _reindex_with_indexers(
        self: FrameOrSeries,
        reindexers,
        fill_value=None,
        copy: bool_t = False,
        allow_dups: bool_t = False,
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def filter(
        self: FrameOrSeries,
        items=None,
        like: Optional[str] = None,
        regex: Optional[str] = None,
        axis=None,
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def head(self: FrameOrSeries, n: int = 5) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def tail(self: FrameOrSeries, n: int = 5) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def sample(
        self: FrameOrSeries,
        n=None,
        frac=None,
        replace=False,
        weights=None,
        random_state=None,
        axis=None,
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def pipe(self, func, *args, **kwargs):
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Attribute access

    def __finalize__(
        self: FrameOrSeries, other, method: Optional[str] = None, **kwargs
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def __getattr__(self, name: str):
        raise AbstractMethodError(self)

    def __setattr__(self, name: str, value) -> None:
        raise AbstractMethodError(self)

    def _dir_additions(self) -> Set[str]:
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Consolidation of internals

    def _protect_consolidate(self, f):
        raise AbstractMethodError(self)

    def _consolidate_inplace(self) -> None:
        raise AbstractMethodError(self)

    def _consolidate(self):
        raise AbstractMethodError(self)

    @property
    def _is_mixed_type(self) -> bool_t:
        raise AbstractMethodError(self)

    def _check_inplace_setting(self, value) -> bool_t:
        raise AbstractMethodError(self)

    def _get_numeric_data(self):
        raise AbstractMethodError(self)

    def _get_bool_data(self):
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Internal Interface Methods

    @property
    def values(self) -> np.ndarray:
        raise AbstractMethodError(self)

    @property
    def _values(self) -> np.ndarray:
        raise AbstractMethodError(self)

    @property
    def dtypes(self):
        raise AbstractMethodError(self)

    def _to_dict_of_blocks(self, copy: bool_t = True):
        raise AbstractMethodError(self)

    def astype(
        self: FrameOrSeries, dtype, copy: bool_t = True, errors: str = "raise"
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def copy(self: FrameOrSeries, deep: bool_t = True) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def __copy__(self: FrameOrSeries, deep: bool_t = True) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def __deepcopy__(self: FrameOrSeries, memo=None) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def _convert(
        self: FrameOrSeries,
        datetime: bool_t = False,
        numeric: bool_t = False,
        timedelta: bool_t = False,
        coerce: bool_t = False,
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def infer_objects(self: FrameOrSeries) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def convert_dtypes(
        self: FrameOrSeries,
        infer_objects: bool_t = True,
        convert_string: bool_t = True,
        convert_integer: bool_t = True,
        convert_boolean: bool_t = True,
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Filling NA's

    def fillna(
        self: FrameOrSeries,
        value=None,
        method=None,
        axis=None,
        inplace: bool_t = False,
        limit=None,
        downcast=None,
    ) -> Optional[FrameOrSeries]:
        raise AbstractMethodError(self)

    def ffill(
        self: FrameOrSeries,
        axis=None,
        inplace: bool_t = False,
        limit=None,
        downcast=None,
    ) -> Optional[FrameOrSeries]:
        raise AbstractMethodError(self)

    pad = ffill

    def bfill(
        self: FrameOrSeries,
        axis=None,
        inplace: bool_t = False,
        limit=None,
        downcast=None,
    ) -> Optional[FrameOrSeries]:
        raise AbstractMethodError(self)

    backfill = bfill

    def replace(
        self,
        to_replace=None,
        value=None,
        inplace: bool_t = False,
        limit: Optional[int] = None,
        regex=False,
        method="pad",
    ):
        raise AbstractMethodError(self)

    def interpolate(
        self: FrameOrSeries,
        method: str = "linear",
        axis: Axis = 0,
        limit: Optional[int] = None,
        inplace: bool_t = False,
        limit_direction: Optional[str] = None,
        limit_area: Optional[str] = None,
        downcast: Optional[str] = None,
        **kwargs,
    ) -> Optional[FrameOrSeries]:
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Timeseries methods Methods

    def asof(self, where, subset=None):
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Action Methods

    def isna(self: FrameOrSeries) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def isnull(self: FrameOrSeries) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def notna(self: FrameOrSeries) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def notnull(self: FrameOrSeries) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def _clip_with_scalar(self, lower, upper, inplace: bool_t = False):
        raise AbstractMethodError(self)

    def _clip_with_one_bound(self, threshold, method, axis, inplace):
        raise AbstractMethodError(self)

    def clip(
        self: FrameOrSeries,
        lower=None,
        upper=None,
        axis=None,
        inplace: bool_t = False,
        *args,
        **kwargs,
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def asfreq(
        self: FrameOrSeries,
        freq,
        method=None,
        how: Optional[str] = None,
        normalize: bool_t = False,
        fill_value=None,
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def at_time(
        self: FrameOrSeries, time, asof: bool_t = False, axis=None
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def between_time(
        self: FrameOrSeries,
        start_time,
        end_time,
        include_start: bool_t = True,
        include_end: bool_t = True,
        axis=None,
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def resample(
        self,
        rule,
        axis=0,
        closed: Optional[str] = None,
        label: Optional[str] = None,
        convention: str = "start",
        kind: Optional[str] = None,
        loffset=None,
        base: Optional[int] = None,
        on=None,
        level=None,
        origin: Union[str, TimestampConvertibleTypes] = "start_day",
        offset: Optional[TimedeltaConvertibleTypes] = None,
    ) -> Resampler:
        raise AbstractMethodError(self)

    def first(self: FrameOrSeries, offset) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def last(self: FrameOrSeries, offset) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def rank(
        self: FrameOrSeries,
        axis=0,
        method: str = "average",
        numeric_only: Optional[bool_t] = None,
        na_option: str = "keep",
        ascending: bool_t = True,
        pct: bool_t = False,
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def compare(
        self,
        other,
        align_axis: Axis = 1,
        keep_shape: bool_t = False,
        keep_equal: bool_t = False,
    ):
        raise AbstractMethodError(self)

    def align(
        self,
        other,
        join="outer",
        axis=None,
        level=None,
        copy=True,
        fill_value=None,
        method=None,
        limit=None,
        fill_axis=0,
        broadcast_axis=None,
    ):
        raise AbstractMethodError(self)

    def _align_frame(
        self,
        other,
        join="outer",
        axis=None,
        level=None,
        copy: bool_t = True,
        fill_value=None,
        method=None,
        limit=None,
        fill_axis=0,
    ):
        raise AbstractMethodError(self)

    def _align_series(
        self,
        other,
        join="outer",
        axis=None,
        level=None,
        copy: bool_t = True,
        fill_value=None,
        method=None,
        limit=None,
        fill_axis=0,
    ):
        raise AbstractMethodError(self)

    def _where(
        self,
        cond,
        other=np.nan,
        inplace=False,
        axis=None,
        level=None,
        errors="raise",
        try_cast=False,
    ):
        raise AbstractMethodError(self)

    def where(
        self,
        cond,
        other=np.nan,
        inplace=False,
        axis=None,
        level=None,
        errors="raise",
        try_cast=False,
    ):
        raise AbstractMethodError(self)

    def mask(
        self,
        cond,
        other=np.nan,
        inplace=False,
        axis=None,
        level=None,
        errors="raise",
        try_cast=False,
    ):
        raise AbstractMethodError(self)

    def shift(
        self: FrameOrSeries, periods=1, freq=None, axis=0, fill_value=None
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def slice_shift(self: FrameOrSeries, periods: int = 1, axis=0) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def tshift(
        self: FrameOrSeries, periods: int = 1, freq=None, axis: Axis = 0
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def truncate(
        self: FrameOrSeries, before=None, after=None, axis=None, copy: bool_t = True
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def tz_convert(
        self: FrameOrSeries, tz, axis=0, level=None, copy: bool_t = True
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def tz_localize(
        self: FrameOrSeries,
        tz,
        axis=0,
        level=None,
        copy: bool_t = True,
        ambiguous="raise",
        nonexistent: str = "raise",
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Numeric Methods

    def abs(self: FrameOrSeries) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def describe(
        self: FrameOrSeries,
        percentiles=None,
        include=None,
        exclude=None,
        datetime_is_numeric=False,
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def pct_change(
        self: FrameOrSeries,
        periods=1,
        fill_method="pad",
        limit=None,
        freq=None,
        **kwargs,
    ) -> FrameOrSeries:
        raise AbstractMethodError(self)

    def _agg_by_level(self, name, axis=0, level=0, skipna=True, **kwargs):
        raise AbstractMethodError(self)

    def _logical_func(
        self, name: str, func, axis=0, bool_only=None, skipna=True, level=None, **kwargs
    ):
        raise AbstractMethodError(self)

    def any(self, axis=0, bool_only=None, skipna=True, level=None, **kwargs):
        raise AbstractMethodError(self)

    def all(self, axis=0, bool_only=None, skipna=True, level=None, **kwargs):
        raise AbstractMethodError(self)

    def _accum_func(self, name: str, func, axis=None, skipna=True, *args, **kwargs):
        raise AbstractMethodError(self)

    def cummax(self, axis=None, skipna=True, *args, **kwargs):
        raise AbstractMethodError(self)

    def cummin(self, axis=None, skipna=True, *args, **kwargs):
        raise AbstractMethodError(self)

    def cumsum(self, axis=None, skipna=True, *args, **kwargs):
        raise AbstractMethodError(self)

    def cumprod(self, axis=None, skipna=True, *args, **kwargs):
        raise AbstractMethodError(self)

    def _stat_function_ddof(
        self,
        name: str,
        func,
        axis=None,
        skipna=None,
        level=None,
        ddof=1,
        numeric_only=None,
        **kwargs,
    ):
        raise AbstractMethodError(self)

    def sem(
        self, axis=None, skipna=None, level=None, ddof=1, numeric_only=None, **kwargs
    ):
        raise AbstractMethodError(self)

    def var(
        self, axis=None, skipna=None, level=None, ddof=1, numeric_only=None, **kwargs
    ):
        raise AbstractMethodError(self)

    def std(
        self, axis=None, skipna=None, level=None, ddof=1, numeric_only=None, **kwargs
    ):
        raise AbstractMethodError(self)

    def _stat_function(
        self,
        name: str,
        func,
        axis=None,
        skipna=None,
        level=None,
        numeric_only=None,
        **kwargs,
    ):
        raise AbstractMethodError(self)

    def min(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs):
        raise AbstractMethodError(self)

    def max(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs):
        raise AbstractMethodError(self)

    def mean(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs):
        raise AbstractMethodError(self)

    def median(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs):
        raise AbstractMethodError(self)

    def skew(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs):
        raise AbstractMethodError(self)

    def kurt(self, axis=None, skipna=None, level=None, numeric_only=None, **kwargs):
        raise AbstractMethodError(self)

    kurtosis = kurt

    def _min_count_stat_function(
        self,
        name: str,
        func,
        axis=None,
        skipna=None,
        level=None,
        numeric_only=None,
        min_count=0,
        **kwargs,
    ):
        raise AbstractMethodError(self)

    def sum(
        self,
        axis=None,
        skipna=None,
        level=None,
        numeric_only=None,
        min_count=0,
        **kwargs,
    ):
        raise AbstractMethodError(self)

    def prod(
        self,
        axis=None,
        skipna=None,
        level=None,
        numeric_only=None,
        min_count=0,
        **kwargs,
    ):
        raise AbstractMethodError(self)

    product = prod

    def mad(self, axis=None, skipna=None, level=None):
        raise AbstractMethodError(self)

    @classmethod
    def _add_numeric_operations(cls):
        raise AbstractMethodError(cls)

    def rolling(
        self,
        window: Union[int, timedelta, BaseOffset, BaseIndexer],
        min_periods: Optional[int] = None,
        center: bool_t = False,
        win_type: Optional[str] = None,
        on: Optional[str] = None,
        axis: Axis = 0,
        closed: Optional[str] = None,
    ):
        raise AbstractMethodError(self)

    def expanding(
        self, min_periods: int = 1, center: Optional[bool_t] = None, axis: Axis = 0
    ) -> Expanding:
        raise AbstractMethodError(self)

    def ewm(
        self,
        com: Optional[float] = None,
        span: Optional[float] = None,
        halflife: Optional[Union[float, TimedeltaConvertibleTypes]] = None,
        alpha: Optional[float] = None,
        min_periods: int = 0,
        adjust: bool_t = True,
        ignore_na: bool_t = False,
        axis: Axis = 0,
        times: Optional[Union[str, np.ndarray, FrameOrSeries]] = None,
    ) -> ExponentialMovingWindow:
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Arithmetic Methods

    def _inplace_method(self, other, op):
        raise AbstractMethodError(self)

    def __iadd__(self, other):
        raise AbstractMethodError(self)

    def __isub__(self, other):
        raise AbstractMethodError(self)

    def __imul__(self, other):
        raise AbstractMethodError(self)

    def __itruediv__(self, other):
        raise AbstractMethodError(self)

    def __ifloordiv__(self, other):
        raise AbstractMethodError(self)

    def __imod__(self, other):
        raise AbstractMethodError(self)

    def __ipow__(self, other):
        raise AbstractMethodError(self)

    def __iand__(self, other):
        raise AbstractMethodError(self)

    def __ior__(self, other):
        raise AbstractMethodError(self)

    def __ixor__(self, other):
        raise AbstractMethodError(self)

    # ----------------------------------------------------------------------
    # Misc methods

    def _find_valid_index(self, how: str):
        raise AbstractMethodError(self)

    def first_valid_index(self):
        raise AbstractMethodError(self)

    def last_valid_index(self):
        raise AbstractMethodError(self)
