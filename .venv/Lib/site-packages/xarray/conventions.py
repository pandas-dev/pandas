from __future__ import annotations

import itertools
import warnings
from collections import defaultdict
from collections.abc import Hashable, Iterable, Mapping, MutableMapping
from typing import TYPE_CHECKING, Any, Literal, TypeVar, Union

import numpy as np

from xarray.coders import CFDatetimeCoder, CFTimedeltaCoder
from xarray.coding import strings, variables
from xarray.coding.variables import SerializationWarning, pop_to
from xarray.core import indexing
from xarray.core.common import (
    _contains_datetime_like_objects,
    contains_cftime_datetimes,
)
from xarray.core.utils import emit_user_level_warning
from xarray.core.variable import IndexVariable, Variable
from xarray.namedarray.utils import is_duck_array

CF_RELATED_DATA = (
    "bounds",
    "grid_mapping",
    "climatology",
    "geometry",
    "node_coordinates",
    "node_count",
    "part_node_count",
    "interior_ring",
    "cell_measures",
    "formula_terms",
)
CF_RELATED_DATA_NEEDS_PARSING = (
    "grid_mapping",
    "cell_measures",
    "formula_terms",
)


if TYPE_CHECKING:
    from xarray.backends.common import AbstractDataStore
    from xarray.core.dataset import Dataset

    T_VarTuple = tuple[tuple[Hashable, ...], Any, dict, dict]
    T_Name = Union[Hashable, None]
    T_Variables = Mapping[Any, Variable]
    T_Attrs = MutableMapping[Any, Any]
    T_DropVariables = Union[str, Iterable[Hashable], None]
    T_DatasetOrAbstractstore = Union[Dataset, AbstractDataStore]


def ensure_not_multiindex(var: Variable, name: T_Name = None) -> None:
    # only the pandas multi-index dimension coordinate cannot be serialized (tuple values)
    if isinstance(var._data, indexing.PandasMultiIndexingAdapter):
        if name is None and isinstance(var, IndexVariable):
            name = var.name
        if var.dims == (name,):
            raise NotImplementedError(
                f"variable {name!r} is a MultiIndex, which cannot yet be "
                "serialized. Instead, either use reset_index() "
                "to convert MultiIndex levels into coordinate variables instead "
                "or use https://cf-xarray.readthedocs.io/en/latest/coding.html."
            )


def encode_cf_variable(
    var: Variable, needs_copy: bool = True, name: T_Name = None
) -> Variable:
    """
    Converts a Variable into a Variable which follows some
    of the CF conventions:

        - Nans are masked using _FillValue (or the deprecated missing_value)
        - Rescaling via: scale_factor and add_offset
        - datetimes are converted to the CF 'units since time' format
        - dtype encodings are enforced.

    Parameters
    ----------
    var : Variable
        A variable holding un-encoded data.

    Returns
    -------
    out : Variable
        A variable which has been encoded as described above.
    """
    ensure_not_multiindex(var, name=name)

    for coder in [
        CFDatetimeCoder(),
        CFTimedeltaCoder(),
        variables.CFScaleOffsetCoder(),
        variables.CFMaskCoder(),
        variables.NativeEnumCoder(),
        variables.NonStringCoder(),
        variables.DefaultFillvalueCoder(),
        variables.BooleanCoder(),
    ]:
        var = coder.encode(var, name=name)

    for attr_name in CF_RELATED_DATA:
        pop_to(var.encoding, var.attrs, attr_name)
    return var


def decode_cf_variable(
    name: Hashable,
    var: Variable,
    concat_characters: bool = True,
    mask_and_scale: bool = True,
    decode_times: bool | CFDatetimeCoder = True,
    decode_endianness: bool = True,
    stack_char_dim: bool = True,
    use_cftime: bool | None = None,
    decode_timedelta: bool | CFTimedeltaCoder | None = None,
) -> Variable:
    """
    Decodes a variable which may hold CF encoded information.

    This includes variables that have been masked and scaled, which
    hold CF style time variables (this is almost always the case if
    the dataset has been serialized) and which have strings encoded
    as character arrays.

    Parameters
    ----------
    name : str
        Name of the variable. Used for better error messages.
    var : Variable
        A variable holding potentially CF encoded information.
    concat_characters : bool
        Should character arrays be concatenated to strings, for
        example: ["h", "e", "l", "l", "o"] -> "hello"
    mask_and_scale : bool
        Lazily scale (using scale_factor and add_offset) and mask
        (using _FillValue). If the _Unsigned attribute is present
        treat integer arrays as unsigned.
    decode_times : bool or CFDatetimeCoder
        Decode cf times ("hours since 2000-01-01") to np.datetime64.
    decode_endianness : bool
        Decode arrays from non-native to native endianness.
    stack_char_dim : bool
        Whether to stack characters into bytes along the last dimension of this
        array. Passed as an argument because we need to look at the full
        dataset to figure out if this is appropriate.
    use_cftime : bool, optional
        Only relevant if encoded dates come from a standard calendar
        (e.g. "gregorian", "proleptic_gregorian", "standard", or not
        specified).  If None (default), attempt to decode times to
        ``np.datetime64[ns]`` objects; if this is not possible, decode times to
        ``cftime.datetime`` objects. If True, always decode times to
        ``cftime.datetime`` objects, regardless of whether or not they can be
        represented using ``np.datetime64[ns]`` objects.  If False, always
        decode times to ``np.datetime64[ns]`` objects; if this is not possible
        raise an error.

        .. deprecated:: 2025.01.1
           Please pass a :py:class:`coders.CFDatetimeCoder` instance initialized with ``use_cftime`` to the ``decode_times`` kwarg instead.
    decode_timedelta : None, bool, or CFTimedeltaCoder
        Decode cf timedeltas ("hours") to np.timedelta64.

    Returns
    -------
    out : Variable
        A variable holding the decoded equivalent of var.
    """
    # Ensure datetime-like Variables are passed through unmodified (GH 6453)
    if _contains_datetime_like_objects(var):
        return var

    original_dtype = var.dtype

    decode_timedelta_was_none = decode_timedelta is None
    if decode_timedelta is None:
        if isinstance(decode_times, CFDatetimeCoder):
            decode_timedelta = CFTimedeltaCoder(time_unit=decode_times.time_unit)
        else:
            decode_timedelta = bool(decode_times)

    if concat_characters:
        if stack_char_dim:
            var = strings.CharacterArrayCoder().decode(var, name=name)
        var = strings.EncodedStringCoder().decode(var)

    if original_dtype.kind == "O":
        var = variables.ObjectVLenStringCoder().decode(var)
        original_dtype = var.dtype

    if original_dtype.kind == "T":
        var = variables.Numpy2StringDTypeCoder().decode(var)

    if mask_and_scale:
        for coder in [
            variables.CFMaskCoder(
                decode_times=decode_times, decode_timedelta=decode_timedelta
            ),
            variables.CFScaleOffsetCoder(
                decode_times=decode_times, decode_timedelta=decode_timedelta
            ),
        ]:
            var = coder.decode(var, name=name)

    if decode_timedelta:
        if isinstance(decode_timedelta, bool):
            decode_timedelta = CFTimedeltaCoder(
                decode_via_units=decode_timedelta, decode_via_dtype=decode_timedelta
            )
        decode_timedelta._emit_decode_timedelta_future_warning = (
            decode_timedelta_was_none
        )
        var = decode_timedelta.decode(var, name=name)
    if decode_times:
        # remove checks after end of deprecation cycle
        if not isinstance(decode_times, CFDatetimeCoder):
            if use_cftime is not None:
                emit_user_level_warning(
                    "Usage of 'use_cftime' as a kwarg is deprecated. "
                    "Please pass a 'CFDatetimeCoder' instance initialized "
                    "with 'use_cftime' to the 'decode_times' kwarg instead.\n"
                    "Example usage:\n"
                    "    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)\n"
                    "    ds = xr.open_dataset(decode_times=time_coder)\n",
                    DeprecationWarning,
                )
            decode_times = CFDatetimeCoder(use_cftime=use_cftime)
        elif use_cftime is not None:
            raise TypeError(
                "Usage of 'use_cftime' as a kwarg is not allowed "
                "if a 'CFDatetimeCoder' instance is passed to "
                "'decode_times'. Please set 'use_cftime' "
                "when initializing 'CFDatetimeCoder' instead.\n"
                "Example usage:\n"
                "    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)\n"
                "    ds = xr.open_dataset(decode_times=time_coder)\n",
            )
        var = decode_times.decode(var, name=name)

    if decode_endianness and not var.dtype.isnative:
        var = variables.EndianCoder().decode(var)
        original_dtype = var.dtype

    var = variables.BooleanCoder().decode(var)

    dimensions, data, attributes, encoding = variables.unpack_for_decoding(var)

    encoding.setdefault("dtype", original_dtype)

    if (
        # we don't need to lazily index duck arrays
        not is_duck_array(data)
        # These arrays already support lazy indexing
        # OR for IndexingAdapters, it makes no sense to wrap them
        and not isinstance(data, indexing.ExplicitlyIndexedNDArrayMixin)
    ):
        # this path applies to bare BackendArray objects.
        # It is not hit for any internal Xarray backend
        data = indexing.LazilyIndexedArray(data)

    return Variable(dimensions, data, attributes, encoding=encoding, fastpath=True)


def _update_bounds_attributes(variables: T_Variables) -> None:
    """Adds time attributes to time bounds variables.

    Variables handling time bounds ("Cell boundaries" in the CF
    conventions) do not necessarily carry the necessary attributes to be
    decoded. This copies the attributes from the time variable to the
    associated boundaries.

    See Also:

    http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/
         cf-conventions.html#cell-boundaries

    https://github.com/pydata/xarray/issues/2565
    """

    # For all time variables with bounds
    for v in variables.values():
        attrs = v.attrs
        units = attrs.get("units")
        has_date_units = isinstance(units, str) and "since" in units
        if has_date_units and "bounds" in attrs and attrs["bounds"] in variables:
            bounds_attrs = variables[attrs["bounds"]].attrs
            bounds_attrs.setdefault("units", attrs["units"])
            if "calendar" in attrs:
                bounds_attrs.setdefault("calendar", attrs["calendar"])


def _update_bounds_encoding(variables: T_Variables) -> None:
    """Adds time encoding to time bounds variables.

    Variables handling time bounds ("Cell boundaries" in the CF
    conventions) do not necessarily carry the necessary attributes to be
    decoded. This copies the encoding from the time variable to the
    associated bounds variable so that we write CF-compliant files.

    See Also:

    http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/
         cf-conventions.html#cell-boundaries

    https://github.com/pydata/xarray/issues/2565
    """

    # For all time variables with bounds
    for name, v in variables.items():
        attrs = v.attrs
        encoding = v.encoding
        has_date_units = "units" in encoding and "since" in encoding["units"]
        is_datetime_type = np.issubdtype(
            v.dtype, np.datetime64
        ) or contains_cftime_datetimes(v)

        if (
            is_datetime_type
            and not has_date_units
            and "bounds" in attrs
            and attrs["bounds"] in variables
        ):
            emit_user_level_warning(
                f"Variable {name} has datetime type and a "
                f"bounds variable but {name}.encoding does not have "
                f"units specified. The units encodings for {name} "
                f"and {attrs['bounds']} will be determined independently "
                "and may not be equal, counter to CF-conventions. "
                "If this is a concern, specify a units encoding for "
                f"{name} before writing to a file.",
            )

        if has_date_units and "bounds" in attrs and attrs["bounds"] in variables:
            bounds_encoding = variables[attrs["bounds"]].encoding
            bounds_encoding.setdefault("units", encoding["units"])
            if "calendar" in encoding:
                bounds_encoding.setdefault("calendar", encoding["calendar"])


T = TypeVar("T")
U = TypeVar("U")


def _item_or_default(obj: Mapping[Any, T | U] | T, key: Hashable, default: T) -> T | U:
    """
    Return item by key if obj is mapping and key is present, else return default value.
    """
    return obj.get(key, default) if isinstance(obj, Mapping) else obj


def decode_cf_variables(
    variables: T_Variables,
    attributes: T_Attrs,
    concat_characters: bool | Mapping[str, bool] = True,
    mask_and_scale: bool | Mapping[str, bool] = True,
    decode_times: bool | CFDatetimeCoder | Mapping[str, bool | CFDatetimeCoder] = True,
    decode_coords: bool | Literal["coordinates", "all"] = True,
    drop_variables: T_DropVariables = None,
    use_cftime: bool | Mapping[str, bool] | None = None,
    decode_timedelta: bool
    | CFTimedeltaCoder
    | Mapping[str, bool | CFTimedeltaCoder]
    | None = None,
) -> tuple[T_Variables, T_Attrs, set[Hashable]]:
    """
    Decode several CF encoded variables.

    See: decode_cf_variable
    """
    # Only emit one instance of the decode_timedelta default change
    # FutureWarning. This can be removed once this change is made.
    warnings.filterwarnings("once", "decode_timedelta", FutureWarning)

    dimensions_used_by = defaultdict(list)
    for v in variables.values():
        for d in v.dims:
            dimensions_used_by[d].append(v)

    def stackable(dim: Hashable) -> bool:
        # figure out if a dimension can be concatenated over
        if dim in variables:
            return False
        for v in dimensions_used_by[dim]:
            if v.dtype.kind != "S" or dim != v.dims[-1]:
                return False
        return True

    coord_names = set()

    if isinstance(drop_variables, str):
        drop_variables = [drop_variables]
    elif drop_variables is None:
        drop_variables = []
    drop_variables = set(drop_variables)

    # Time bounds coordinates might miss the decoding attributes
    if decode_times:
        _update_bounds_attributes(variables)

    new_vars = {}
    for k, v in variables.items():
        if k in drop_variables:
            continue
        stack_char_dim = (
            _item_or_default(concat_characters, k, True)
            and v.dtype == "S1"
            and v.ndim > 0
            and stackable(v.dims[-1])
        )
        try:
            new_vars[k] = decode_cf_variable(
                k,
                v,
                concat_characters=_item_or_default(concat_characters, k, True),
                mask_and_scale=_item_or_default(mask_and_scale, k, True),
                decode_times=_item_or_default(decode_times, k, True),
                stack_char_dim=stack_char_dim,
                use_cftime=_item_or_default(use_cftime, k, None),
                decode_timedelta=_item_or_default(decode_timedelta, k, None),
            )
        except Exception as e:
            raise type(e)(f"Failed to decode variable {k!r}: {e}") from e
        if decode_coords in [True, "coordinates", "all"]:
            var_attrs = new_vars[k].attrs
            if "coordinates" in var_attrs:
                var_coord_names = [
                    c for c in var_attrs["coordinates"].split() if c in variables
                ]
                # propagate as is
                new_vars[k].encoding["coordinates"] = var_attrs["coordinates"]
                del var_attrs["coordinates"]
                # but only use as coordinate if existing
                if var_coord_names:
                    coord_names.update(var_coord_names)

        if decode_coords == "all":
            for attr_name in CF_RELATED_DATA:
                if attr_name in var_attrs:
                    # fixes stray colon
                    attr_val = var_attrs[attr_name].replace(" :", ":")
                    var_names = attr_val.split()
                    # if grid_mapping is a single string, do not enter here
                    if (
                        attr_name in CF_RELATED_DATA_NEEDS_PARSING
                        and len(var_names) > 1
                    ):
                        # map the keys to list of strings
                        # "A: b c d E: f g" returns
                        # {"A": ["b", "c", "d"], "E": ["f", "g"]}
                        roles_and_names = defaultdict(list)
                        key = None
                        for vname in var_names:
                            if ":" in vname:
                                key = vname.strip(":")
                            else:
                                if key is None:
                                    raise ValueError(
                                        f"First element {vname!r} of [{attr_val!r}] misses ':', "
                                        f"cannot decode {attr_name!r}."
                                    )
                                roles_and_names[key].append(vname)
                        # for grid_mapping keys are var_names
                        if attr_name == "grid_mapping":
                            var_names = list(roles_and_names.keys())
                        else:
                            # for cell_measures and formula_terms values are var names
                            var_names = list(itertools.chain(*roles_and_names.values()))
                            # consistency check (one element per key)
                            if len(var_names) != len(roles_and_names.keys()):
                                emit_user_level_warning(
                                    f"Attribute {attr_name!r} has malformed content [{attr_val!r}], "
                                    f"decoding {var_names!r} to coordinates."
                                )
                    if all(var_name in variables for var_name in var_names):
                        new_vars[k].encoding[attr_name] = attr_val
                        coord_names.update(var_names)
                    else:
                        referenced_vars_not_in_variables = [
                            proj_name
                            for proj_name in var_names
                            if proj_name not in variables
                        ]
                        emit_user_level_warning(
                            f"Variable(s) referenced in {attr_name} not in variables: {referenced_vars_not_in_variables}",
                        )
                    del var_attrs[attr_name]

    if decode_coords and isinstance(attributes.get("coordinates", None), str):
        attributes = dict(attributes)
        crds = attributes.pop("coordinates")
        coord_names.update(crds.split())

    return new_vars, attributes, coord_names


def decode_cf(
    obj: T_DatasetOrAbstractstore,
    concat_characters: bool = True,
    mask_and_scale: bool = True,
    decode_times: bool | CFDatetimeCoder | Mapping[str, bool | CFDatetimeCoder] = True,
    decode_coords: bool | Literal["coordinates", "all"] = True,
    drop_variables: T_DropVariables = None,
    use_cftime: bool | None = None,
    decode_timedelta: bool
    | CFTimedeltaCoder
    | Mapping[str, bool | CFTimedeltaCoder]
    | None = None,
) -> Dataset:
    """Decode the given Dataset or Datastore according to CF conventions into
    a new Dataset.

    Parameters
    ----------
    obj : Dataset or DataStore
        Object to decode.
    concat_characters : bool, optional
        Should character arrays be concatenated to strings, for
        example: ["h", "e", "l", "l", "o"] -> "hello"
    mask_and_scale : bool, optional
        Lazily scale (using scale_factor and add_offset) and mask
        (using _FillValue).
    decode_times : bool | CFDatetimeCoder | Mapping[str, bool | CFDatetimeCoder], optional
        Decode cf times (e.g., integers since "hours since 2000-01-01") to
        np.datetime64.
    decode_coords : bool or {"coordinates", "all"}, optional
        Controls which variables are set as coordinate variables:

        - "coordinates" or True: Set variables referred to in the
          ``'coordinates'`` attribute of the datasets or individual variables
          as coordinate variables.
        - "all": Set variables referred to in  ``'grid_mapping'``, ``'bounds'`` and
          other attributes as coordinate variables.
    drop_variables : str or iterable, optional
        A variable or list of variables to exclude from being parsed from the
        dataset. This may be useful to drop variables with problems or
        inconsistent values.
    use_cftime : bool, optional
        Only relevant if encoded dates come from a standard calendar
        (e.g. "gregorian", "proleptic_gregorian", "standard", or not
        specified).  If None (default), attempt to decode times to
        ``np.datetime64[ns]`` objects; if this is not possible, decode times to
        ``cftime.datetime`` objects. If True, always decode times to
        ``cftime.datetime`` objects, regardless of whether or not they can be
        represented using ``np.datetime64[ns]`` objects.  If False, always
        decode times to ``np.datetime64[ns]`` objects; if this is not possible
        raise an error.

        .. deprecated:: 2025.01.1
           Please pass a :py:class:`coders.CFDatetimeCoder` instance initialized with ``use_cftime`` to the ``decode_times`` kwarg instead.

    decode_timedelta : bool | CFTimedeltaCoder | Mapping[str, bool | CFTimedeltaCoder], optional
        If True or :py:class:`CFTimedeltaCoder`, decode variables and
        coordinates with time units in
        {"days", "hours", "minutes", "seconds", "milliseconds", "microseconds"}
        into timedelta objects. If False, leave them encoded as numbers.
        If None (default), assume the same behavior as decode_times. The
        resolution of the decoded timedeltas can be configured with the
        ``time_unit`` argument in the :py:class:`CFTimedeltaCoder` passed.

    Returns
    -------
    decoded : Dataset
    """
    from xarray.backends.common import AbstractDataStore
    from xarray.core.dataset import Dataset

    vars: T_Variables
    attrs: T_Attrs
    if isinstance(obj, Dataset):
        vars = obj._variables
        attrs = obj.attrs
        extra_coords = set(obj.coords)
        close = obj._close
        encoding = obj.encoding
    elif isinstance(obj, AbstractDataStore):
        vars, attrs = obj.load()
        extra_coords = set()
        close = obj.close
        encoding = obj.get_encoding()
    else:
        raise TypeError("can only decode Dataset or DataStore objects")

    vars, attrs, coord_names = decode_cf_variables(
        vars,
        attrs,
        concat_characters,
        mask_and_scale,
        decode_times,
        decode_coords,
        drop_variables=drop_variables,
        use_cftime=use_cftime,
        decode_timedelta=decode_timedelta,
    )
    ds = Dataset(vars, attrs=attrs)
    ds = ds.set_coords(coord_names.union(extra_coords).intersection(vars))
    ds.set_close(close)
    ds.encoding = encoding

    return ds


def cf_decoder(
    variables: T_Variables,
    attributes: T_Attrs,
    concat_characters: bool = True,
    mask_and_scale: bool = True,
    decode_times: bool | CFDatetimeCoder | Mapping[str, bool | CFDatetimeCoder] = True,
) -> tuple[T_Variables, T_Attrs]:
    """
    Decode a set of CF encoded variables and attributes.

    Parameters
    ----------
    variables : dict
        A dictionary mapping from variable name to xarray.Variable
    attributes : dict
        A dictionary mapping from attribute name to value
    concat_characters : bool
        Should character arrays be concatenated to strings, for
        example: ["h", "e", "l", "l", "o"] -> "hello"
    mask_and_scale : bool
        Lazily scale (using scale_factor and add_offset) and mask
        (using _FillValue).
    decode_times : bool | CFDatetimeCoder | Mapping[str, bool | CFDatetimeCoder]
        Decode cf times ("hours since 2000-01-01") to np.datetime64.

    Returns
    -------
    decoded_variables : dict
        A dictionary mapping from variable name to xarray.Variable objects.
    decoded_attributes : dict
        A dictionary mapping from attribute name to values.

    See Also
    --------
    decode_cf_variable
    """
    variables, attributes, _ = decode_cf_variables(
        variables,
        attributes,
        concat_characters,
        mask_and_scale,
        decode_times,
    )
    return variables, attributes


def _encode_coordinates(
    variables: T_Variables, attributes: T_Attrs, non_dim_coord_names
):
    # calculate global and variable specific coordinates
    non_dim_coord_names = set(non_dim_coord_names)

    for name in list(non_dim_coord_names):
        if isinstance(name, str) and " " in name:
            emit_user_level_warning(
                f"coordinate {name!r} has a space in its name, which means it "
                "cannot be marked as a coordinate on disk and will be "
                "saved as a data variable instead",
                category=SerializationWarning,
            )
            non_dim_coord_names.discard(name)

    global_coordinates = non_dim_coord_names.copy()
    variable_coordinates = defaultdict(set)
    not_technically_coordinates = set()
    for coord_name in non_dim_coord_names:
        target_dims = variables[coord_name].dims
        for k, v in variables.items():
            if (
                k not in non_dim_coord_names
                and k not in v.dims
                and set(target_dims) <= set(v.dims)
            ):
                variable_coordinates[k].add(coord_name)

            if any(
                coord_name in v.encoding.get(attr_name, tuple())
                for attr_name in CF_RELATED_DATA
            ):
                not_technically_coordinates.add(coord_name)
                global_coordinates.discard(coord_name)

    variables = {k: v.copy(deep=False) for k, v in variables.items()}

    # keep track of variable names written to file under the "coordinates" attributes
    written_coords = set()
    for name, var in variables.items():
        encoding = var.encoding
        attrs = var.attrs
        if "coordinates" in attrs and "coordinates" in encoding:
            raise ValueError(
                f"'coordinates' found in both attrs and encoding for variable {name!r}."
            )

        # if coordinates set to None, don't write coordinates attribute
        if ("coordinates" in attrs and attrs.get("coordinates") is None) or (
            "coordinates" in encoding and encoding.get("coordinates") is None
        ):
            # make sure "coordinates" is removed from attrs/encoding
            attrs.pop("coordinates", None)
            encoding.pop("coordinates", None)
            continue

        # this will copy coordinates from encoding to attrs if "coordinates" in attrs
        # after the next line, "coordinates" is never in encoding
        # we get support for attrs["coordinates"] for free.
        coords_str = pop_to(encoding, attrs, "coordinates") or attrs.get("coordinates")
        if not coords_str and variable_coordinates[name]:
            coordinates_text = " ".join(
                str(coord_name)
                for coord_name in sorted(variable_coordinates[name])
                if coord_name not in not_technically_coordinates
            )
            if coordinates_text:
                attrs["coordinates"] = coordinates_text
        if "coordinates" in attrs:
            written_coords.update(attrs["coordinates"].split())

    # These coordinates are not associated with any particular variables, so we
    # save them under a global 'coordinates' attribute so xarray can roundtrip
    # the dataset faithfully. Because this serialization goes beyond CF
    # conventions, only do it if necessary.
    # Reference discussion:
    # https://cfconventions.org/mailing-list-archive/Data/7400.html
    global_coordinates.difference_update(written_coords)
    if global_coordinates:
        attributes = dict(attributes)
        if "coordinates" in attributes:
            emit_user_level_warning(
                f"cannot serialize global coordinates {global_coordinates!r} because the global "
                f"attribute 'coordinates' already exists. This may prevent faithful roundtripping"
                f"of xarray datasets",
                category=SerializationWarning,
            )
        else:
            attributes["coordinates"] = " ".join(sorted(map(str, global_coordinates)))

    return variables, attributes


def encode_dataset_coordinates(dataset: Dataset):
    """Encode coordinates on the given dataset object into variable specific
    and global attributes.

    When possible, this is done according to CF conventions.

    Parameters
    ----------
    dataset : Dataset
        Object to encode.

    Returns
    -------
    variables : dict
    attrs : dict
    """
    non_dim_coord_names = set(dataset.coords) - set(dataset.dims)
    return _encode_coordinates(
        dataset._variables, dataset.attrs, non_dim_coord_names=non_dim_coord_names
    )


def cf_encoder(variables: T_Variables, attributes: T_Attrs):
    """
    Encode a set of CF encoded variables and attributes.
    Takes a dicts of variables and attributes and encodes them
    to conform to CF conventions as much as possible.
    This includes masking, scaling, character array handling,
    and CF-time encoding.

    Parameters
    ----------
    variables : dict
        A dictionary mapping from variable name to xarray.Variable
    attributes : dict
        A dictionary mapping from attribute name to value

    Returns
    -------
    encoded_variables : dict
        A dictionary mapping from variable name to xarray.Variable,
    encoded_attributes : dict
        A dictionary mapping from attribute name to value

    See Also
    --------
    decode_cf_variable, encode_cf_variable
    """

    # add encoding for time bounds variables if present.
    _update_bounds_encoding(variables)

    new_vars = {}
    for k, v in variables.items():
        try:
            new_vars[k] = encode_cf_variable(v, name=k)
        except Exception as e:
            e.add_note(f"Raised while encoding variable {k!r} with value {v!r}")
            raise

    # Remove attrs from bounds variables (issue #2921)
    for var in new_vars.values():
        bounds = var.attrs.get("bounds")
        if bounds and bounds in new_vars:
            # see http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
            for attr in [
                "units",
                "standard_name",
                "axis",
                "positive",
                "calendar",
                "long_name",
                "leap_month",
                "leap_year",
                "month_lengths",
            ]:
                if (
                    attr in new_vars[bounds].attrs
                    and attr in var.attrs
                    and new_vars[bounds].attrs[attr] == var.attrs[attr]
                ):
                    new_vars[bounds].attrs.pop(attr)

    return new_vars, attributes
