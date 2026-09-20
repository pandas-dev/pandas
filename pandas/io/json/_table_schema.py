"""
Table Schema builders

https://specs.frictionlessdata.io/table-schema/
"""

from __future__ import annotations

from collections import Counter
from datetime import timezone
from typing import (
    TYPE_CHECKING,
    Any,
    cast,
)
import warnings

from pandas._config import option_context

from pandas._libs import lib
from pandas._libs._ujson import ujson_loads
from pandas._libs.tslibs import timezones
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.base import _registry as registry
from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_integer_dtype,
    is_numeric_dtype,
    is_string_dtype,
)
from pandas.core.dtypes.dtypes import (
    CategoricalDtype,
    DatetimeTZDtype,
    ExtensionDtype,
    PeriodDtype,
)

from pandas import (
    DataFrame,
    Index,
)
import pandas.core.common as com

from pandas.tseries.frequencies import to_offset

if TYPE_CHECKING:
    from collections.abc import (
        Hashable,
        Sequence,
    )

    from pandas._typing import (
        DtypeObj,
        JSONSerializable,
    )

    from pandas import Series
    from pandas.core.indexes.multi import MultiIndex


TABLE_SCHEMA_VERSION = "1.4.0"


def as_json_table_type(x: DtypeObj) -> str:
    """
    Convert a NumPy / pandas type to its corresponding json_table.

    Parameters
    ----------
    x : np.dtype or ExtensionDtype

    Returns
    -------
    str
        the Table Schema data types

    Notes
    -----
    This table shows the relationship between NumPy / pandas dtypes,
    and Table Schema dtypes.

    ==============  =================
    Pandas type     Table Schema type
    ==============  =================
    int64           integer
    float64         number
    bool            boolean
    datetime64[ns]  datetime
    timedelta64[ns] duration
    object          str
    categorical     any
    =============== =================
    """
    if is_integer_dtype(x):
        return "integer"
    elif is_bool_dtype(x):
        return "boolean"
    elif is_numeric_dtype(x):
        return "number"
    elif lib.is_np_dtype(x, "M") or isinstance(x, (DatetimeTZDtype, PeriodDtype)):
        return "datetime"
    elif lib.is_np_dtype(x, "m"):
        return "duration"
    elif is_string_dtype(x):
        return "string"
    else:
        return "any"


def set_default_names(data):
    """Sets index names to 'index' for regular, or 'level_x' for Multi"""
    if com.all_not_none(*data.index.names):
        nms = data.index.names
        if len(nms) == 1 and data.index.name == "index":
            warnings.warn(
                "Index name of 'index' is not round-trippable.",
                stacklevel=find_stack_level(),
            )
        elif len(nms) > 1 and any(
            isinstance(x, str) and x.startswith("level_") for x in nms
        ):
            warnings.warn(
                "Index names beginning with 'level_' are not round-trippable.",
                stacklevel=find_stack_level(),
            )
        return data

    data = data.copy(deep=False)
    if data.index.nlevels > 1:
        data.index.names = com.fill_missing_names(data.index.names)
    else:
        data.index.name = data.index.name or "index"
    return data


def convert_pandas_type_to_json_field(arr) -> dict[str, JSONSerializable]:
    dtype = arr.dtype
    name: JSONSerializable
    if arr.name is None:
        name = "values"
    else:
        name = arr.name
    field: dict[str, JSONSerializable] = {
        "name": name,
        "type": as_json_table_type(dtype),
    }

    if isinstance(dtype, CategoricalDtype):
        cats = dtype.categories
        ordered = dtype.ordered

        field["constraints"] = {"enum": list(cats)}
        field["ordered"] = ordered
    elif isinstance(dtype, PeriodDtype):
        field["freq"] = dtype.freq.freqstr
    elif isinstance(dtype, DatetimeTZDtype):
        if timezones.is_utc(dtype.tz):
            field["tz"] = "UTC"
        else:
            zone = timezones.get_timezone(dtype.tz)
            if isinstance(zone, str):
                field["tz"] = zone
            elif isinstance(zone, timezone):
                # fixed-offset stdlib timezone, e.g. timezone(timedelta(hours=1));
                # str gives a round-trippable "UTC+HH:MM" (GH#39537)
                field["tz"] = str(zone)
    elif isinstance(dtype, ExtensionDtype):
        field["extDtype"] = dtype.name
    return field


def convert_json_field_to_pandas_type(field) -> str | CategoricalDtype:
    """
    Converts a JSON field descriptor into its corresponding NumPy / pandas type

    Parameters
    ----------
    field
        A JSON field descriptor

    Returns
    -------
    dtype

    Raises
    ------
    ValueError
        If the type of the provided field is unknown or currently unsupported

    Examples
    --------
    >>> convert_json_field_to_pandas_type({"name": "an_int", "type": "integer"})
    'int64'

    >>> convert_json_field_to_pandas_type(
    ...     {
    ...         "name": "a_categorical",
    ...         "type": "any",
    ...         "constraints": {"enum": ["a", "b", "c"]},
    ...         "ordered": True,
    ...     }
    ... )
    CategoricalDtype(categories=['a', 'b', 'c'], ordered=True, categories_dtype=str)

    >>> convert_json_field_to_pandas_type({"name": "a_datetime", "type": "datetime"})
    'datetime64[ns]'

    >>> convert_json_field_to_pandas_type(
    ...     {"name": "a_datetime_with_tz", "type": "datetime", "tz": "US/Central"}
    ... )
    'datetime64[ns, US/Central]'
    """
    typ = field["type"]
    if typ == "string":
        return field.get("extDtype", None)
    elif typ == "integer":
        return field.get("extDtype", "int64")
    elif typ == "number":
        return field.get("extDtype", "float64")
    elif typ == "boolean":
        return field.get("extDtype", "bool")
    elif typ == "duration":
        return "timedelta64"
    elif typ == "datetime":
        if field.get("tz"):
            return f"datetime64[ns, {field['tz']}]"
        elif field.get("freq"):
            # GH#9586 rename frequency M to ME for offsets
            offset = to_offset(field["freq"])
            freq = PeriodDtype(offset)._freqstr
            # GH#47747 using datetime over period to minimize the change surface
            return f"period[{freq}]"
        else:
            return "datetime64[ns]"
    elif typ == "any":
        if "constraints" in field and "ordered" in field:
            return CategoricalDtype(
                categories=field["constraints"]["enum"], ordered=field["ordered"]
            )
        elif "extDtype" in field:
            return registry.find(field["extDtype"])
        else:
            return "object"

    raise ValueError(f"Unsupported or invalid field type: {typ}")


def build_table_schema(
    data: DataFrame | Series,
    index: bool = True,
    primary_key: bool | None = None,
    version: bool = True,
) -> dict[str, JSONSerializable]:
    """
    Create a Table schema from ``data``.

    This method is a utility to generate a JSON-serializable schema
    representation of a pandas Series or DataFrame, compatible with the
    Table Schema specification. It enables structured data to be shared
    and validated in various applications, ensuring consistency and
    interoperability.

    Parameters
    ----------
    data : Series or DataFrame
        The input data for which the table schema is to be created.
    index : bool, default True
        Whether to include ``data.index`` in the schema.
    primary_key : bool or None, default True
        Column names to designate as the primary key.
        The default `None` will set `'primaryKey'` to the index
        level or levels if the index is unique.
    version : bool, default True
        Whether to include a field `pandas_version` with the version
        of pandas that last revised the table schema. This version
        can be different from the installed pandas version.

    Returns
    -------
    dict
        A dictionary representing the Table schema.

    See Also
    --------
    DataFrame.to_json : Convert the object to a JSON string.
    read_json : Convert a JSON string to pandas object.

    Notes
    -----
    See `Table Schema
    <https://pandas.pydata.org/docs/user_guide/io.html#table-schema>`__ for
    conversion types.
    Timedeltas as converted to ISO8601 duration format with
    9 decimal places after the seconds field for nanosecond precision.

    Categoricals are converted to the `any` dtype, and use the `enum` field
    constraint to list the allowed values. The `ordered` attribute is included
    in an `ordered` field.

    Examples
    --------
    >>> from pandas.io.json._table_schema import build_table_schema
    >>> df = pd.DataFrame(
    ...     {'A': [1, 2, 3],
    ...      'B': ['a', 'b', 'c'],
    ...      'C': pd.date_range('2016-01-01', freq='D', periods=3),
    ...      }, index=pd.Index(range(3), name='idx'))
    >>> build_table_schema(df)
    {'fields': \
[{'name': 'idx', 'type': 'integer'}, \
{'name': 'A', 'type': 'integer'}, \
{'name': 'B', 'type': 'string', 'extDtype': 'str'}, \
{'name': 'C', 'type': 'datetime'}], \
'primaryKey': ['idx'], \
'pandas_version': '1.4.0'}
    """
    if index is True:
        data = set_default_names(data)

    schema: dict[str, Any] = {}
    fields = []

    if index:
        if data.index.nlevels > 1:
            data.index = cast("MultiIndex", data.index)
            for level, name in zip(data.index.levels, data.index.names, strict=True):
                new_field = convert_pandas_type_to_json_field(level)
                new_field["name"] = name
                fields.append(new_field)
        else:
            fields.append(convert_pandas_type_to_json_field(data.index))

    if data.ndim > 1:
        for column, s in data.items():
            fields.append(convert_pandas_type_to_json_field(s))
    else:
        fields.append(convert_pandas_type_to_json_field(data))

    schema["fields"] = fields
    if index and data.index.is_unique and primary_key is None:
        if data.index.nlevels == 1:
            schema["primaryKey"] = [data.index.name]
        else:
            schema["primaryKey"] = data.index.names
    elif primary_key is not None:
        schema["primaryKey"] = primary_key

    if version:
        schema["pandas_version"] = TABLE_SCHEMA_VERSION
    return schema


def _unmatched_names(
    names: Sequence[Hashable], records: Sequence[Any]
) -> list[Hashable]:
    """
    Non-string field names with no key holding their values in "data".
    """
    # JSON object keys are always strings, so a non-string field name is
    #  keyed by its string form in "data" (GH#19129).
    unmatched = []
    for name in names:
        if isinstance(name, str):
            continue
        key = str(name)
        if not any(isinstance(record, dict) and key in record for record in records):
            unmatched.append(name)
    return unmatched


def _float_names_from_keys(
    perturbed: Sequence[Hashable], records: Sequence[Any]
) -> dict[float, float]:
    """
    Map each perturbed float name to the exact label its "data" key spells.
    """
    # a key that parses to the same double as the field name spells that
    #  label exactly, so float() on the key recovers it (GH#19129)
    wanted = set(perturbed)
    matches: dict[float, list[float]] = {}
    record_keys = {
        key for record in records if isinstance(record, dict) for key in record
    }
    for key in record_keys:
        try:
            value = ujson_loads(key, precise_float=False)
        except ValueError:
            # not a JSON literal, so not the spelling of a float label
            continue
        if isinstance(value, float) and value in wanted:
            matches.setdefault(value, []).append(float(key))
    # a name two keys could spell is left alone, for the unmatched check to catch
    return {value: keys[0] for value, keys in matches.items() if len(keys) == 1}


def _field_position(col_order: Sequence[Hashable], key: Hashable) -> int | None:
    """
    Index of the field `key` names, or None if it names no field.
    """
    # matched on type as well as value, since pandas conflates 1, 1.0 and True
    #  but "primaryKey" names exactly one of them (GH#19129)
    for pos, col in enumerate(col_order):
        if type(col) is type(key) and col == key:
            return pos
    return None


def parse_table_schema(json, precise_float: bool) -> DataFrame:
    """
    Builds a DataFrame from a given schema

    Parameters
    ----------
    json :
        A JSON table schema
    precise_float : bool
        Flag controlling precision when decoding string to double values, as
        dictated by ``read_json``

    Returns
    -------
    df : DataFrame

    Raises
    ------
    NotImplementedError
        If the JSON table schema contains either timezone or timedelta data
    ValueError
        If a field name cannot be matched to the data it labels, or if
        "primaryKey" names a field the schema does not declare

    Notes
    -----
        Because :func:`DataFrame.to_json` uses the string 'index' to denote a
        name-less :class:`Index`, this function sets the name of the returned
        :class:`DataFrame` to ``None`` when said string is encountered with a
        normal :class:`Index`. For a :class:`MultiIndex`, the same limitation
        applies to any strings beginning with 'level_'. Therefore, an
        :class:`Index` name of 'index'  and :class:`MultiIndex` names starting
        with 'level_' are not supported.

    See Also
    --------
    build_table_schema : Inverse function.
    pandas.read_json
    """
    table = ujson_loads(json, precise_float=precise_float)
    schema = table["schema"]
    records = table["data"]
    names = [field["name"] for field in schema["fields"]]
    # only an object record is keyed by label; with none of them there is no
    #  key to match a label against, so DataFrame takes the labels as they are
    rows = records if isinstance(records, list) else []
    keyed = any(isinstance(record, dict) for record in rows)
    positional = isinstance(records, list) and not keyed
    # a float label often does not match the "data" key holding its values:
    #  the fast parser perturbs the schema literal (GH#19129). precise_float
    #  governs the data values, not the labels
    recovered: dict[float, float] = {}
    unmatched_floats = _unmatched_names(
        [name for name in names if isinstance(name, float)], rows
    )
    if unmatched_floats and keyed:
        # the key spells the label, so parsing it back recovers the label --
        #  including a sign the writer dropped, as it does for -0.0
        recovered = _float_names_from_keys(unmatched_floats, rows)
        unmatched_floats = [name for name in unmatched_floats if name not in recovered]
    if unmatched_floats and not precise_float:
        # the key did not spell the label, so the schema literal itself has to
        #  be parsed exactly, which costs a second decode of the whole document
        try:
            schema = ujson_loads(json, precise_float=True)["schema"]
        except ValueError:
            # some literal in the document is out of range for the exact
            #  parser, so the fast parser's labels are all there is
            pass
        else:
            names = [field["name"] for field in schema["fields"]]
    if recovered:
        names = [
            recovered.get(name, name) if isinstance(name, float) else name
            for name in names
        ]
    fields = schema["fields"]
    # a keyed record is looked up by the label's string form, so the frame is
    #  built under those keys and the labels restored at the end
    col_order: list[Hashable] = (
        list(names) if positional else [str(name) for name in names]
    )
    if keyed:
        # two labels sharing that string form are indistinguishable in "data"
        counts = Counter(col_order)
        collisions = [
            name for name, col in zip(names, col_order, strict=True) if counts[col] > 1
        ]
        if collisions:
            raise ValueError(
                f"Field names {collisions} share a string form, which is how "
                "'data' keys its values, so they cannot be read back"
            )
        unmatched = _unmatched_names(names, rows)
        # a label absent from every record is ambiguous: its values may all be
        #  missing, or a key may spell it differently. A key that no field
        #  claims is the evidence for the second, and without one the column
        #  reads as all-missing, which is what a string label already does
        claimed = {str(name) for name in names}
        orphans = sorted(
            key
            for record in rows
            if isinstance(record, dict)
            for key in record
            if key not in claimed
        )
        if unmatched and orphans:
            msg = f"Field names {unmatched} have no matching key in 'data'"
            if any(isinstance(name, float) for name in unmatched):
                msg += (
                    "; to_json writes a float label at 'double_precision' "
                    "digits, 15 at most, but keys its values by the full repr"
                )
            raise ValueError(f"{msg}; {orphans} match no field name")
    df = DataFrame(records, columns=col_order)
    # address the frame by position from here on: pandas conflates 1, 1.0 and
    #  True and reads None back as NaN, so a label is not a safe key
    df.columns = range(len(names))

    dtypes = {
        pos: convert_json_field_to_pandas_type(field)
        for pos, field in enumerate(fields)
    }

    # No ISO constructor for Timedelta as of yet, so need to raise
    if "timedelta64" in dtypes.values():
        raise NotImplementedError(
            'table="orient" can not yet read ISO-formatted Timedelta data'
        )

    with option_context("future.distinguish_nan_and_na", False):
        df = df.astype(dtypes)

    if "primaryKey" in schema:
        pkey = schema["primaryKey"]
        if not isinstance(pkey, list):
            # the spec allows a bare field name as well as an array of them,
            #  and build_table_schema(primary_key=0) writes a non-string one
            pkey = [pkey]
        # primaryKey repeats the field name, so it was perturbed the same way
        pkey = [
            recovered.get(key, key) if isinstance(key, float) else key for key in pkey
        ]
        primary_key = pkey if positional else [str(key) for key in pkey]
        key_positions = []
        for raw, key in zip(pkey, primary_key, strict=True):
            pos = _field_position(col_order, key)
            if pos is None:
                raise ValueError(f"'primaryKey' names {raw!r}, which is not a field")
            key_positions.append(pos)
        df = df.set_index(key_positions)
        df.index.names = [names[pos] for pos in key_positions]
        if len(df.index.names) == 1:
            if df.index.name == "index":
                df.index.name = None
        else:
            df.index.names = [
                None if isinstance(name, str) and name.startswith("level_") else name
                for name in df.index.names
            ]

    # undo the stringification; rebuilding from a plain list also lets a
    #  uniform label type infer its own dtype rather than staying object
    if df.columns.empty:
        # nothing to infer from, so borrow the field names' own dtype
        df.columns = Index(names)[:0]
    else:
        df.columns = [names[pos] for pos in df.columns]

    return df
