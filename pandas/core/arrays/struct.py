from __future__ import annotations

from copy import deepcopy
import operator
from collections import namedtuple, OrderedDict
from typing import TYPE_CHECKING, Type, Union, TypeVar, Tuple, List, NamedTuple, Iterable, Sequence, Dict, Any

from numpy import dtype

import re

import numpy as np
import pandas as pd
from pandas.api.extensions import (
    register_extension_dtype,
    ExtensionDtype,
    ExtensionArray,
)

# from pandas.core.construction import extract_array

try:
    import pyarrow as pa
except ImportError:
    pa = None

import logging

log = logging.getLogger(__name__)

NoneType = type(None)

StructArrayT = TypeVar("StructArrayT", bound="StructArray")


__all__ = [
    "StructDtype",
    "StructArray",
    # "ArrowStructType",
]

def _ZERO_VALUE_FOR_TYPE(dtype):
    # convert to pandas dtype
    dtype = pd.api.types.pandas_dtype(dtype)
    if hasattr(dtype, "zero_value"):
        return dtype.zero_value
    elif (
            pd.api.types.is_numeric_dtype(dtype) or
            pd.api.types.is_datetime64_dtype(dtype) or
            pd.api.types.is_timedelta64_dtype(dtype) or
            pd.api.types.is_categorical_dtype(dtype)
    ):
        return 0
    elif pd.api.types.is_bool_dtype(dtype):
        return False
    elif pd.api.types.is_string_dtype(dtype):
        return ""
    else:
        return ValueError(f"Unknown type: '{dtype}'")


def _EMPTY_ARRAY(size, dtype):
    zero_value = _ZERO_VALUE_FOR_TYPE(dtype)
    return pd.array([zero_value], dtype=dtype).repeat(size)


def _INFER_DTYPE_FROM_SCALAR(value):
    # this does not work
    # return pd.api.types.pandas_dtype(pd.api.types.infer_dtype([v]))

    # ugly hack abusing pd.array:
    return pd.array([value]).dtype


# def _ARRAY_CLASS_FROM_TYPE(dtype):
#     if pd.api.types.is_extension_type(dtype):
#         cls = dtype.construct_array_type()
#
#     return cls

def _CONCAT_SAME_TYPE_FN(dtype):
    if hasattr(dtype, "construct_array_type"):
        cls = dtype.construct_array_type()
        if hasattr(cls, "_concat_same_type"):
            return cls._concat_same_type

    # workaround via factorization
    return lambda values: pd.array(np.concatenate(values), dtype=dtype)


@register_extension_dtype
class StructDtype(ExtensionDtype):
    name = "Struct"

    kind = "O"
    str = "|O08"
    base = np.dtype("O")
    num = 103
    _metadata = ("fields", "nullable")
    # e.g. Struct([(x: Integer), (y: Integer)], True)
    _match = re.compile(r"[S|s]truct\(\[(?P<fields>.+)\](,(?P<nullable>.*))?\)")
    _cache: Dict[str, ExtensionDtype] = {}

    _fields: OrderedDict[str, ExtensionDtype]
    _type: NamedTuple

    def __new__(cls, fields: Union[List[Tuple[str, Type]], Dict[str, Type], str], nullable=None):
        # from pandas.core.dtypes.common import is_string_dtype, pandas_dtype

        if isinstance(fields, StructDtype):
            return fields
        elif pa is not None and isinstance(fields, pa.StructType):
            pa_structtype: pa.StructType = fields
            fields = OrderedDict()
            for idx in range(pa_structtype.num_fields):
                field = pa_structtype[idx]
                field_name = field.name
                field_dtype = field.type.to_pandas_dtype()
                fields[field_name] = field_dtype
        elif fields is None:
            raise ValueError("Missing field types")
        elif isinstance(fields, str):
            m = cls._match.search(fields)
            if m is not None:
                # parse field types
                fields = []
                for field in m.group("fields").split(","):
                    s = field.split(":")
                    if len(s) != 2:
                        raise ValueError("Wrong field description: %s" % field)
                    field_name: str = s[0].strip()
                    field_type = pd.api.types.pandas_dtype(s[1].strip())

                    fields.append((field_name, field_type))
                # get nullability of struct
                nullable = m.group(nullable)
                if nullable is not None:
                    nullable = bool(nullable.strip())

        # make sure that every field type is actually a Pandas type
        fields = OrderedDict(fields)
        pd_type_fields = OrderedDict()
        for f_name, f_type in fields.items():
            f_type = pd.api.types.pandas_dtype(f_type)
            pd_type_fields[f_name] = f_type
        fields = pd_type_fields

        # default for struct is nullable=True
        if nullable is None:
            nullable = True

        # deduplication of dtype
        cache_key = str(fields) + str(nullable)
        try:
            return cls._cache[cache_key]
        except KeyError:
            u = object.__new__(cls)
            # set dtype properties
            u._fields = fields
            u._nullable = nullable
            u._type = None

            cls._cache[cache_key] = u
            return u

    @property
    def type(self) -> Type[NamedTuple]:
        if self._type is None:
            # create scalar struct type
            self._type: NamedTuple = namedtuple("Struct", self.fields.keys())
            self._type.dtype = self

        return self._type

    @property
    def fields(self) -> OrderedDict[str, ExtensionDtype]:
        """
        The fields contained in this struct.
        """
        return self._fields

    @property
    def field_names(self):
        return self._fields.keys()

    @property
    def field_dtypes(self):
        return self._fields.values()

    @property
    def nullable(self):
        return self._nullable

    @property
    def zero_value(self):
        return self.type(*[_ZERO_VALUE_FOR_TYPE(x) for x in self.field_dtypes])

    @classmethod
    def construct_array_type(cls) -> Type["StructArray"]:
        """
        Return the array type associated with this dtype.
        Returns
        -------
        type
        """
        return StructArray

    def __str__(self) -> str:
        field_str_reprs = [
            f"{f_name}: {str(f_type)}" for f_name, f_type in self.fields.items()
        ]
        return f'''struct([{", ".join(field_str_reprs)}], nullable={self.nullable}'''

    def __hash__(self) -> int:
        # make myself hashable
        return hash(str(self))

    def __repr__(self) -> str:
        return f"dtype('{str(self)}'"

    def __from_arrow__(
            self, array: Union["pa.Array", "pa.ChunkedArray"]
    ) -> "StructArray":
        """
        Construct VariantArray from pyarrow Array/ChunkedArray.
        """
        import pyarrow as pa

        if isinstance(array, pa.Array):
            chunks = [array]
        else:
            # pyarrow.ChunkedArray
            chunks = array.chunks

        results = []
        for arr in chunks:
            data = OrderedDict()
            for f_name, f_type in self.fields:
                pd_array = pd.array(arr.storage.field(f_name), dtype=f_type)
                data[f_name] = pd_array

            if self.nullable:
                # set missing values correctly
                mask = ~ np.asarray(arr.is_valid())
            else:
                mask = None

            iarr = StructArray(data, mask, dtype=self)

            results.append(iarr)

        return StructArray._concat_same_type(results)


# # x = StructDtype({"a": int, "b": int})
# class ArrowStructType(pa.ExtensionType):
#     def __init__(self, storage_type: pa.StructType):
#         pa.ExtensionType.__init__(self, storage_type, "struct")
#
#     def __arrow_ext_serialize__(self):
#         # metadata = {"subtype": str(self.subtype), "closed": self.closed}
#         # return json.dumps(metadata).encode()
#         return b'{}'
#
#     @classmethod
#     def __arrow_ext_deserialize__(cls, storage_type: pa.StructType, serialized):
#         # metadata = json.loads(serialized.decode())
#         # subtype = pyarrow.type_for_alias(metadata["subtype"])
#         # closed = metadata["closed"]
#         return ArrowStructType(storage_type)
#
#     def __eq__(self, other):
#         if isinstance(other, pa.BaseExtensionType):
#             return type(self) == type(other)
#         else:
#             return NotImplemented
#
#     def __hash__(self):
#         return hash(str(self))
#
#     def to_pandas_dtype(self):
#         return StructDtype(self.storage_type)
#
#
# import pyarrow as pa
# # register the type with a dummy instance
# _arrow_variant_type = ArrowVariantType()
# pa.register_extension_type(_arrow_variant_type)

class StructArray(ExtensionArray):
    """
    TODO
    Array representation of a list of Variant objects.

    Internally, stores the variants in a column-based representation.

    Attributes:
        - chrom: array of chromosome names
        - start: array of variant starts (0-based, inclusive)
        - end: array of variant ends (1-based, exclusive)
        - ref: array of reference sequences
        - alt: array of alternative sequences

    In addition, can be converted from/to VCF-formatted (chrom, pos, ref, alt) representations.
    """
    _fields: OrderedDict[str, ExtensionArray]
    "fields that this StructArray holds"
    # _field_names: List[str]
    # "ordered list of field names that this StructArray holds"
    _mask: Union[NoneType, np.ndarray]
    "boolean array that denotes missing values"
    _dtype: StructDtype
    "data type of this StructArray"

    na_value = pd.NA

    def __init__(self, data, mask=None, dtype=None, copy=False, validate=True):
        if isinstance(data, StructArray):
            if dtype is not None and dtype != data.dtype:
                raise ValueError(f"Passed StructArray but dtypes do not match: '{dtype}' != '{data.dtype}'")
            else:
                dtype = data.dtype

            fields = data.fields
            if dtype.nullable and mask is None:
                mask = data._mask

            if copy:
                fields = deepcopy(fields)
                mask = deepcopy(mask)
        elif isinstance(data, OrderedDict) and not copy:
            fields = data
            if dtype is None:
                field_types = self._infer_field_dtypes_from_arrays(fields.values, names=fields.keys())
                nullable = (mask is not None)
                dtype = StructDtype(field_types, nullable=nullable)

            if copy:
                fields = deepcopy(fields)
                mask = deepcopy(mask)
        elif isinstance(data, dict) or isinstance(data, pd.DataFrame):
            fields = OrderedDict([(k, pd.array(v, copy=copy)) for k, v in data.items()])
            if dtype is None:
                field_types = self._infer_field_dtypes_from_arrays(fields.values, names=fields.keys())
                nullable = (mask is not None)
                dtype = StructDtype(field_types, nullable=nullable)

            if copy:
                # no need to copy fields as this will already be done by pd.array if necessary
                mask = deepcopy(mask)
        elif self._is_valid_scalar(data):
            # if dtype is None:
            #     field_types = self._infer_field_dtypes_from_scalar(data)
            #     nullable = (mask is not None)
            #     dtype = StructDtype(field_types, nullable=nullable)
            # fields = self._parse_scalar(data, dtype=dtype)
            raise ValueError(f"Cannot pass scalar '{data}' to '{self}'.")
        else:
            # assuming list of scalars
            if dtype is None:
                dtype = self._infer_dtype_from_list_of_scalars(data)

            data = self._parse_listlike(data, dtype)
            fields = data._fields
            # no need to copy 'fields' since it is a copy anyways

            if mask is None:
                mask = data._mask
            else:
                if copy:
                    # need to copy the specified mask
                    mask = deepcopy(mask)

        # TODO: Copy-on-write when copy=False?

        self._fields = fields
        self._mask = mask
        self._dtype = dtype

        if validate:
            self._validate()

    def _validate(self):
        validated_fields = OrderedDict()
        for f, v in self._fields.items():
            validated_fields[f] = pd.array(v, dtype=self.dtype.fields[f], copy=False)
        self._fields = validated_fields

    def __getattr__(self, item):
        if item in self.fields:
            return self.fields[item]
        else:
            raise AttributeError

    @property
    def fields(self):
        return self._fields

    @property
    def field_names(self):
        return self._fields.keys()

    @property
    def field_values(self):
        return self._fields.values()

    @property
    def field_types(self):
        return self._dtype.fields

    @classmethod
    def _validate_dtypes(cls, data: Iterable[Tuple[str, any]], dtype: Iterable[Tuple[str, Type]]):
        """
        ensure correct array types
        """
        return [
            (k, pd.array(v, dtype=v_dtype)) for ((k, v), (k2, v_dtype)) in zip(data, dtype)
        ]

    def _filter(self, cond, inplace=False):
        # handle scalar case
        if pd.api.types.is_scalar(cond):
            if not inplace:
                if np.all(pd.isna(self)):
                    return self.na_value
                else:
                    return self.dtype.type(*[v[cond] for v in self.field_values])
            else:
                cond = [cond]

        if inplace:
            retval = self
        else:
            retval = self.copy()

        for f_name, f_data in retval._fields.items():
            retval._fields[f_name] = f_data[cond]
        if retval.dtype.nullable:
            retval._mask = retval._mask[cond]

        return retval

    @classmethod
    def _parse_listlike(cls, values: Iterable, dtype: StructDtype, validate=True) -> StructArrayT:
        # list-like of structs
        if not isinstance(values, list):
            values = list(values)
        length = len(values)

        casted_values = [cls._parse_scalar(v, dtype) for v in values]

        if dtype.nullable:
            mask = np.array([pd.isna(v) for v in values])
        else:
            # validate type
            if validate and np.any([pd.isna(v) for v in values]):
                raise ValueError("Data type '{dtype}' is not nullable but passed values contain NA!")

            mask = None

        # transpose List[Struct] to Tuple[Array]
        fields_list = tuple(_EMPTY_ARRAY(size=length, dtype=d) for f, d in dtype.fields.items())
        for idx, scalar in enumerate(casted_values):
            for field_idx, scalar_field in enumerate(scalar):
                fields_list[field_idx][idx] = scalar_field

        data = OrderedDict()
        for f, v in zip(dtype.field_names, fields_list):
            data[f] = v

        return cls(data, mask=mask, dtype=dtype, copy=False)

    @classmethod
    def _parse_scalar(cls, value, dtype: StructDtype) -> NamedTuple:
        """

        Args:
            value:
            dtype:

        Returns: dtype.type objects or dtype.zero_value if element is NA

        """
        if isinstance(value, tuple):
            retval = dtype.type(*value)
        elif isinstance(value, dict):
            retval = dtype.type(**value)
        elif pd.isna(value):
            retval = dtype.zero_value
        else:
            raise TypeError(
                "can only parse Struct-like objects"
            )

        # TODO: validate types

        return retval

    # @classmethod
    # def _parse_fill_value(cls, value) -> Tuple:
    #     return cls._parse_scalar(value)

    @classmethod
    def _parse_setitem_value(cls, value, dtype: StructDtype):
        if cls._is_valid_scalar(value):
            return cls._parse_scalar(value, dtype)
        else:
            return cls._parse_listlike(value, dtype)

    @staticmethod
    def _is_valid_scalar(value):
        return (
                isinstance(value, tuple) or
                isinstance(value, dict)
        )

    @classmethod
    def _infer_field_dtypes_from_scalar(cls, scalar) -> OrderedDict[str, Any[dtype, ExtensionDtype]]:
        if hasattr(scalar, "dtype"):
            return scalar.dtype

        if isinstance(scalar, NamedTuple):
            field_names = scalar._fields
            values = scalar
        elif isinstance(scalar, tuple):
            # simply enumerate fields
            field_names = [f"f_{idx}" for idx in range(len(scalar))]
            values = scalar
        elif isinstance(scalar, dict):
            field_names = scalar.keys()
            values = scalar.values()
        else:
            raise ValueError(f"Unable to infer field dtypes from '{scalar}'")

        field_dtypes = OrderedDict()
        for f, v in zip(field_names, values):
            field_dtypes[str(f)] = _INFER_DTYPE_FROM_SCALAR(v)

        return field_dtypes

    @classmethod
    def _infer_field_dtypes_from_arrays(cls, arrays, names=None) -> OrderedDict[str, Any[dtype, ExtensionDtype]]:
        field_types = OrderedDict()
        for field_idx, array in enumerate(arrays):
            # use pandas Series to infer dtype
            if not isinstance(array, pd.Series):
                array = pd.Series(array)

            if names is not None:
                field_name = names[field_idx]
            elif array.name is not None:
                # if array is already a series, this also allows to keep its name as field name
                field_name = array.name
            else:
                # just use field index as name
                field_name = f"f_{field_idx}"

            # make sure that field_name is a string
            field_name = str(field_name)
            # get dtype from series
            field_dtype = array.dtype

            if field_name in field_types:
                raise ValueError(f"Duplicate field name: '{field_name}'!")

            field_types[field_name] = field_dtype

        return field_types

    @classmethod
    def _infer_dtype_from_list_of_scalars(cls, values: List[Tuple]) -> StructDtype:
        values = np.asarray(values)
        mask = np.array([pd.isna(v) for v in values])
        nullable = np.any(mask)

        nonnull_values = values[~ mask]
        df = pd.DataFrame.from_records(nonnull_values)

        return StructDtype(
            [(f"f_{field_name}", df.dtypes[field_name]) for field_name in df.columns],
            nullable=nullable
        )

        # for scalar in values:
        #     if pd.isna(scalar):
        #         continue
        #
        #     field_dtypes = cls._infer_field_dtypes_from_scalar(scalar)
        #
        #     return StructDtype(field_dtypes, nullable=nullable)
        #
        # raise ValueError(f"Unable to infer field dtypes!")

    @classmethod
    def from_arrays(
            cls: type[StructArrayT],
            list_of_fields,
            dtype: StructDtype = None,
            mask=None,
            verify=True,
            copy: bool = False,
    ) -> StructArrayT:
        if dtype is None:
            field_types = cls._infer_field_dtypes_from_arrays(list_of_fields)
            nullable = (mask is not None)

            dtype = StructDtype(field_types, nullable=nullable)

        data = OrderedDict()
        for f, v in zip(list_of_fields, dtype.field_names):
            data[f] = v

        return cls(data=data, mask=mask, dtype=dtype, copy=copy)

    @classmethod
    def from_tuples(
            cls: type[StructArrayT],
            data: Iterable[Tuple],
            dtype: StructDtype = None,
    ) -> StructArrayT:
        if dtype is None:
            # infer dtype
            if not isinstance(data, list):
                data = list(data)
            if len(data) == 0:
                raise ValueError("Empty list specified!")

            # if cls._is_valid_scalar(data):
            #     field_types = cls._infer_field_dtypes_from_scalar(data)
            #     return StructDtype(field_types, nullable=True)

            dtype = cls._infer_dtype_from_list_of_scalars(data)

        return cls._parse_listlike(data, dtype=dtype)

    @classmethod
    def from_df(cls, df: pd.DataFrame, is_null_column="is_null", dtype: StructDtype = None) -> StructArrayT:
        """
        Creates a new StructArray from a DataFrame.

        :param df: Pandas dataframe
        :param dtype: optional StructType datatype
        :return: StructArray
        """
        cols = df.columns

        if dtype is not None:
            if not isinstance(dtype, StructDtype):
                raise ValueError("Invalid dtype: '%s' is no instance of StructType" % dtype)
            nullable = dtype.nullable
        else:
            nullable = is_null_column in cols
            dtype = StructDtype(
                [(field_name, df.dtypes[field_name]) for field_name in cols if field_name != is_null_column]
            )

        if nullable:
            # validate is_null column
            if is_null_column not in cols:
                raise ValueError("dtype is nullable but column '%s' is missing!" % is_null_column)
            elif not pd.api.types.is_bool_dtype(df.dtypes[is_null_column]):
                raise ValueError("dtype of '%s' column is not boolean" % is_null_column)

            mask = df[is_null_column]
        else:
            mask = None

        return StructArray(
            data={k: df[k].astype(v) for k, v in dtype.fields.items()},
            mask=mask,
            dtype=dtype,
            copy=False
        )

    def as_frame(self, is_null_column="is_null"):
        retval = pd.DataFrame({k: v for k, v in self.fields.items()})
        if self.dtype.nullable:
            if is_null_column in self.field_names:
                raise ValueError(
                    f"is_null_column='{is_null_column}' cannot be used since there is an equally-named struct field!"
                )
            retval[is_null_column] = self.isna()
        return retval

    # ---------------------------------------------------------------------
    # ExtensionArray interface

    @classmethod
    def _concat_same_type(
            cls: Type[StructArrayT], to_concat: Sequence[StructArrayT]
    ) -> StructArrayT:
        """
        Concatenate multiple arrays of this dtype.

        Parameters
        ----------
        to_concat: sequence of this type

        Returns
        -------
        ExtensionArray
        """
        dtype = to_concat[0].dtype
        data = OrderedDict()
        for fname, ftype in dtype.fields.items():
            # get array type of field
            concat_fn = _CONCAT_SAME_TYPE_FN(ftype)

            # concatenate same types
            data[fname] = concat_fn([v.fields[fname] for v in to_concat])

        if dtype.nullable:
            mask = np.concatenate([pd.isna(v) for v in to_concat])
        else:
            mask = None

        return cls(data, mask=mask, dtype=dtype, copy=False)

    @classmethod
    def _from_sequence(
            cls: type[StructArrayT],
            scalars,
            *,
            dtype=None,
            copy: bool = False,
    ) -> StructArrayT:
        """
        Construct a new ExtensionArray from a sequence of scalars.

        Parameters
        ----------
        scalars: Sequence

            Each element will be an instance of the scalar type for this array, ``cls.dtype.type``.
        dtype: dtype, optional

            Construct for this particular dtype.This should be a Dtype compatible with the ExtensionArray.
        copy: boolean, default False

            If True, copy the underlying data.

        Returns
        -------
        ExtensionArray
        """
        if isinstance(scalars, cls) and scalars.dtype == dtype:
            return cls(scalars, dtype=dtype, copy=copy)

        if dtype is None:
            dtype = cls._infer_dtype_from_list_of_scalars(scalars)
        return cls._parse_listlike(scalars, dtype=dtype)

    def _values_for_factorize(self):
        # type: () -> Tuple[np.ndarray, Any]
        """
        Return an array and missing value suitable for factorization.

        Returns
        -------
        values : ndarray

            An array suitable for factorization. This should maintain order
            and be a supported dtype (Float64, Int64, UInt64, String, Object).
            By default, the extension array is cast to object dtype.
        na_value : object
            The value in `values` to consider missing. This will be treated
            as NA in the factorization routines, so it will be coded as
            `na_sentinel` and not included in `uniques`. By default,
            ``np.nan`` is used.

        Notes
        -----
        The values returned by this method are also used in
        :func:`pandas.util.hash_pandas_object`.
        """
        return np.asarray(self), self.na_value

    @classmethod
    def _from_factorized(
            cls: type[StructArrayT], values: np.ndarray, original: StructArrayT
    ) -> StructArrayT:
        """
        Reconstruct an ExtensionArray after factorization.

        Parameters
        ----------
        values : ndarray
            An integer ndarray with the factorized values.
        original : ExtensionArray
            The original ExtensionArray that factorize was called on.

        See Also
        --------
        factorize : Top-level factorize method that dispatches here.
        ExtensionArray.factorize : Encode the extension array as an enumerated type.
        """
        return cls._from_sequence(values, dtype=original.dtype)

    @classmethod
    def _from_sequence_of_strings(cls, strings, *, dtype=None, copy=False):
        return cls.from_str(strings)

    def __array__(self, dtype=None) -> np.ndarray:
        """
        Return the StructArray's data as a numpy array of Variant objects (with dtype='object')
        """
        mask = self.isna()

        result = np.empty(len(self), dtype=object)
        for i in range(len(self)):
            if mask[i]:
                result[i] = pd.NA
            else:
                result[i] = self.take(i)

        if dtype != None:
            result = result.astype(dtype)
        return result

    def __arrow_array__(self, type=None):
        """
        Convert myself into a pyarrow Array.
        """
        import pyarrow as pa

        storage_array = pa.StructArray.from_arrays(
            arrays=[pa.array(f, from_pandas=True) for f in self.field_values],
            names=list(self.field_names),
        )

        mask = self.isna()
        if mask.any():
            # if there are missing values, set validity bitmap also on the array level
            null_bitmap = pa.array(~mask).buffers()[1]
            storage_array = pa.StructArray.from_buffers(
                storage_array.type,
                len(storage_array),
                [null_bitmap],
                children=[storage_array.field(i) for i in range(len(self.fields))],
            )

        if type is not None:
            if type.equals(storage_array.type):
                return storage_array
            else:
                raise TypeError(
                    f"Not supported to convert StructArray to '{type}' type"
                )

        return storage_array

    # ---------------------------------------------------------------------
    # Descriptive
    def copy(self: StructArrayT) -> StructArrayT:
        """
        Return a copy of the array.
        """
        return StructArray(self, copy=True)

    @property
    def dtype(self) -> StructDtype:
        return self._dtype

    @property
    def nbytes(self) -> int:
        size = sum(f.nbytes for f in self.fields.values())
        if self.dtype.nullable:
            size += self._mask.nbytes
        return size

    @property
    def size(self) -> int:
        # Avoid materializing self.values
        return self._fields.values()[0].size

    def __iter__(self):
        return iter(np.asarray(self))

    def __len__(self) -> int:
        return len(next(iter(self.fields.values())))

    def __getitem__(self, key):
        key = pd.api.indexers.check_array_indexer(self, key)
        return self._filter(key)

    def __setitem__(self, key, value):
        key = pd.api.indexers.check_array_indexer(self, key)
        value_array = self._parse_setitem_value(value, self.dtype)

        for f, v in value_array.fields.items():
            self._fields[f][key] = v

    def take(self, indices, allow_fill=False, fill_value=None):
        from pandas.api.extensions import take

        if allow_fill:
            if fill_value is None or pd.isna(fill_value):
                fill_value = None
            elif self._is_valid_scalar(fill_value):
                fill_value = self._parse_scalar(fill_value, dtype=self.dtype)
            else:
                raise TypeError(f"'{type(fill_value)}' is not a valid fill value!")

        # scalar case
        if pd.api.types.is_scalar(indices):
            if self._mask[indices]:  # missing value at 'indices'
                if allow_fill and fill_value is not None:
                    return fill_value
                else:
                    return self.na_value
            else:
                # return Struct of this type
                return self.dtype.type(*[v[indices] for v in self.field_values])

        fields = OrderedDict()
        for f, v in self.fields.items():
            field_fill_value = fill_value[f] if fill_value is not None else None
            fields[f] = take(v, indices, allow_fill=allow_fill, fill_value=field_fill_value)

        if self.dtype.nullable:
            mask = take(self._mask, indices)
        else:
            mask = None

        return type(self)(fields, mask=mask, dtype=self.dtype, copy=False)

    def isna(self):
        """
        Returns boolean NumPy array indicating if eachvalue is missing
        """
        return self._mask

    def unique(self):
        """
        Compute the ExtensionArray of unique values.

        Returns
        -------
        uniques: ExtensionArray
        """
        # TODO: test alternatives, e.g.:
        # uniques = StructArray(self.as_frame().drop_duplicates())

        factors, uniques = pd.factorize(self)
        if np.any(factors < 0):
            uniques = self._concat_same_type([uniques, StructArray._from_sequence([self.na_value], dtype=self.dtype)])

        return uniques

    def astype(self, dtype, copy=True):
        dtype = pd.api.types.pandas_dtype(dtype)

        if pd.api.types.is_string_dtype(dtype):
            return pd.array(self.to_str(), dtype=dtype)
        elif dtype == self.dtype:
            if copy:
                return self.copy()
            else:
                return self
        else:
            raise TypeError(f"Cannot cast '{self.dtype}' to '{dtype}'")
            # return super().astype(dtype, copy)

    def min(self, axis=None, skipna: bool = True, **kwargs) -> NamedTuple:
        raise NotImplementedError()

    def max(self, axis=None, skipna: bool = True, **kwargs) -> NamedTuple:
        raise NotImplementedError()

    def value_counts(self, dropna=False):
        return pd.value_counts(np.asarray(self), dropna=dropna).astype("Int64")

    def _cmp_method(self, other, op, fail_on_missing=True):
        # ensure pandas array for list-like and eliminate non-compatible scalars
        if self._is_valid_scalar(other):
            other = self._parse_scalar(other, dtype=self.dtype)
            # directly use 'other' as iterable tuple
            other_field_values = other
        else:
            # ensure same dtype
            other = pd.array(other, dtype=self.dtype)

            if len(self) != len(other):
                raise ValueError("Lengths must match to compare")

            other_field_values = other.field_values

        # if other.dtype != self.dtype:
        #     raise ValueError(f"Unknown type of other: {other.dtype}")

        if op is operator.eq or op is operator.ne:
            is_eq = np.array(True)
            for this_f, other_f in zip(self.field_values, other_field_values):
                # element is equal if all fields are equal
                # Has to be casted to numpy boolean array; otherwise we get
                #   `TypeError: unsupported operand type(s) for &: 'PandasArray' and 'PandasArray'`
                # TODO: switch to Pandas boolean array, once these support all binary ops
                is_eq = np.asarray(this_f == other_f) & is_eq

            if self.dtype.nullable:
                is_eq = is_eq | (pd.isna(self) & pd.isna(other))  # return true if both values are NA

            if op is operator.ne:
                return ~is_eq
            else:
                return is_eq
        elif op in {operator.gt, operator.ge, operator.le, operator.lt}:
            prev_is_eq = np.array(True)
            self_is_gt_other = np.array(False)
            for this_f, other_f in zip(self.field_values, other_field_values):
                # TODO: switch to Pandas boolean array, once these support all binary ops
                cur_is_eq = np.asarray(this_f == other_f)
                cur_is_gt = np.asarray(this_f > other_f)

                # self > other if the first non-equal field is greater
                self_is_gt_other = (cur_is_gt & prev_is_eq) | self_is_gt_other
                prev_is_eq = cur_is_eq
            if op is operator.gt:
                retval = self_is_gt_other
            elif op is operator.ge:
                retval = self_is_gt_other | prev_is_eq
            elif op is operator.lt:
                retval = ~ (self_is_gt_other | prev_is_eq)
            else:  # op is operator.le
                retval = ~self_is_gt_other

            if self.dtype.nullable:
                any_of_both_missing = (pd.isna(self) | pd.isna(other))
                if fail_on_missing and np.any(any_of_both_missing):
                        raise TypeError("boolean value of NA is ambiguous")

                # return false if any of the values is NA
                return retval & ~any_of_both_missing
            else:
                return retval
        else:
            raise ValueError(f"Unknown op {op}")

    _arith_method = _cmp_method

    # @unpack_zerodim_and_defer("__eq__")
    def __eq__(self, other):
        return self._cmp_method(other, operator.eq)

    # @unpack_zerodim_and_defer("__ne__")
    def __ne__(self, other):
        return self._cmp_method(other, operator.ne)

    # @unpack_zerodim_and_defer("__gt__")
    def __gt__(self, other):
        return self._cmp_method(other, operator.gt)

    # @unpack_zerodim_and_defer("__ge__")
    def __ge__(self, other):
        return self._cmp_method(other, operator.ge)

    # @unpack_zerodim_and_defer("__lt__")
    def __lt__(self, other):
        return self._cmp_method(other, operator.lt)

    # @unpack_zerodim_and_defer("__le__")
    def __le__(self, other):
        return self._cmp_method(other, operator.le)

    # ---------------------------------------------------------------------
    # Rendering Methods

    def _format_data(self) -> str:
        n = len(self)

        max_seq_items = min((pd.get_option("display.max_seq_items") or n) // 10, 10)

        formatter = str

        if n == 0:
            summary = "[]"
        elif n == 1:
            first = formatter(self[0])
            summary = f"[{first}]"
        elif n == 2:
            first = formatter(self[0])
            last = formatter(self[-1])
            summary = f"[{first}, {last}]"
        else:
            if n > max_seq_items:
                n = min(max_seq_items // 2, 10)
                head = [formatter(x) for x in self[:n]]
                tail = [formatter(x) for x in self[-n:]]
                head_str = ", ".join(head)
                tail_str = ", ".join(tail)
                summary = f"[{head_str} ... {tail_str}]"
            else:
                tail = [formatter(x) for x in self]
                tail_str = ", ".join(tail)
                summary = f"[{tail_str}]"

        return summary

    def __repr__(self) -> str:
        # the short repr has no trailing newline, while the truncated
        # repr does. So we include a newline in our template, and strip
        # any trailing newlines from format_object_summary
        data = self._format_data()
        class_name = f"<{type(self).__name__}>\n"

        template = f"{class_name}{data}\nLength: {len(self)}, dtype: {self.dtype}"
        return template

    def _format_space(self) -> str:
        space = " " * (len(type(self).__name__) + 1)

        return f"\n{space}"

    def __str__(self):
        return self.__repr__()


