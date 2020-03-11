# ---------------------------------------------------------------------
# JSON normalization routines

from collections import defaultdict
import copy
from typing import Any, DefaultDict, Dict, Iterable, List, Optional, Union

import numpy as np

from pandas._libs.writers import convert_json_to_lines
from pandas._typing import Scalar
from pandas.util._decorators import deprecate

import pandas as pd
from pandas import DataFrame


def convert_to_line_delimits(s):
    """
    Helper function that converts JSON lists to line delimited JSON.
    """

    # Determine we have a JSON list to turn to lines otherwise just return the
    # json object, only lists can
    if not s[0] == "[" and s[-1] == "]":
        return s
    s = s[1:-1]

    return convert_json_to_lines(s)


def nested_to_record(
    ds,
    prefix: str = "",
    sep: str = ".",
    level: int = 0,
    max_level: Optional[int] = None,
):
    """
    A simplified json_normalize

    Converts a nested dict into a flat dict ("record"), unlike json_normalize,
    it does not attempt to extract a subset of the data.

    Parameters
    ----------
    ds : dict or list of dicts
    prefix: the prefix, optional, default: ""
    sep : str, default '.'
        Nested records will generate names separated by sep,
        e.g., for sep='.', { 'foo' : { 'bar' : 0 } } -> foo.bar
    level: int, optional, default: 0
        The number of levels in the json string.

    max_level: int, optional, default: None
        The max depth to normalize.

        .. versionadded:: 0.25.0

    Returns
    -------
    d - dict or list of dicts, matching `ds`

    Examples
    --------

    IN[52]: nested_to_record(dict(flat1=1,dict1=dict(c=1,d=2),
                                  nested=dict(e=dict(c=1,d=2),d=2)))
    Out[52]:
    {'dict1.c': 1,
     'dict1.d': 2,
     'flat1': 1,
     'nested.d': 2,
     'nested.e.c': 1,
     'nested.e.d': 2}
    """
    singleton = False
    if isinstance(ds, dict):
        ds = [ds]
        singleton = True
    new_ds = []
    for d in ds:
        new_d = copy.deepcopy(d)
        for k, v in d.items():
            # each key gets renamed with prefix
            if not isinstance(k, str):
                k = str(k)
            if level == 0:
                newkey = k
            else:
                newkey = prefix + sep + k

            # flatten if type is dict and
            # current dict level  < maximum level provided and
            # only dicts gets recurse-flattened
            # only at level>1 do we rename the rest of the keys
            if not isinstance(v, dict) or (
                max_level is not None and level >= max_level
            ):
                if level != 0:  # so we skip copying for top level, common case
                    v = new_d.pop(k)
                    new_d[newkey] = v
                continue
            else:
                v = new_d.pop(k)
                new_d.update(nested_to_record(v, newkey, sep, level + 1, max_level))
        new_ds.append(new_d)

    if singleton:
        return new_ds[0]
    return new_ds


def _json_normalize(
    data: Union[Dict, List[Dict]],
    record_path: Optional[Union[str, List]] = None,
    meta: Optional[Union[str, List[Union[str, List[str]]]]] = None,
    meta_prefix: Optional[str] = None,
    record_prefix: Optional[str] = None,
    errors: Optional[str] = "raise",
    sep: str = ".",
    max_level: Optional[int] = None,
) -> "DataFrame":
    """
    Normalize semi-structured JSON data into a flat table.

    Parameters
    ----------
    data : dict or list of dicts
        Unserialized JSON objects.
    record_path : str or list of str, default None
        Path in each object to list of records. If not passed, data will be
        assumed to be an array of records.
    meta : list of paths (str or list of str), default None
        Fields to use as metadata for each record in resulting table.
    meta_prefix : str, default None
        If True, prefix records with dotted (?) path, e.g. foo.bar.field if
        meta is ['foo', 'bar'].
    record_prefix : str, default None
        If True, prefix records with dotted (?) path, e.g. foo.bar.field if
        path to records is ['foo', 'bar'].
    errors : {'raise', 'ignore'}, default 'raise'
        Configures error handling.

        * 'ignore' : will ignore KeyError if keys listed in meta are not
          always present.
        * 'raise' : will raise KeyError if keys listed in meta are not
          always present.
    sep : str, default '.'
        Nested records will generate names separated by sep.
        e.g., for sep='.', {'foo': {'bar': 0}} -> foo.bar.
    max_level : int, default None
        Max number of levels(depth of dict) to normalize.
        if None, normalizes all levels.

        .. versionadded:: 0.25.0

    Returns
    -------
    frame : DataFrame
    Normalize semi-structured JSON data into a flat table.

    Examples
    --------

    >>> from pandas.io.json import json_normalize
    >>> data = [{'id': 1, 'name': {'first': 'Coleen', 'last': 'Volk'}},
    ...         {'name': {'given': 'Mose', 'family': 'Regner'}},
    ...         {'id': 2, 'name': 'Faye Raker'}]
    >>> json_normalize(data)
        id        name name.family name.first name.given name.last
    0  1.0         NaN         NaN     Coleen        NaN      Volk
    1  NaN         NaN      Regner        NaN       Mose       NaN
    2  2.0  Faye Raker         NaN        NaN        NaN       NaN

    >>> data = [{'id': 1,
    ...          'name': "Cole Volk",
    ...          'fitness': {'height': 130, 'weight': 60}},
    ...         {'name': "Mose Reg",
    ...          'fitness': {'height': 130, 'weight': 60}},
    ...         {'id': 2, 'name': 'Faye Raker',
    ...          'fitness': {'height': 130, 'weight': 60}}]
    >>> json_normalize(data, max_level=0)
                fitness                 id        name
    0   {'height': 130, 'weight': 60}  1.0   Cole Volk
    1   {'height': 130, 'weight': 60}  NaN    Mose Reg
    2   {'height': 130, 'weight': 60}  2.0  Faye Raker

    Normalizes nested data up to level 1.

    >>> data = [{'id': 1,
    ...          'name': "Cole Volk",
    ...          'fitness': {'height': 130, 'weight': 60}},
    ...         {'name': "Mose Reg",
    ...          'fitness': {'height': 130, 'weight': 60}},
    ...         {'id': 2, 'name': 'Faye Raker',
    ...          'fitness': {'height': 130, 'weight': 60}}]
    >>> json_normalize(data, max_level=1)
      fitness.height  fitness.weight   id    name
    0   130              60          1.0    Cole Volk
    1   130              60          NaN    Mose Reg
    2   130              60          2.0    Faye Raker

    >>> data = [{'state': 'Florida',
    ...          'shortname': 'FL',
    ...          'info': {'governor': 'Rick Scott'},
    ...          'counties': [{'name': 'Dade', 'population': 12345},
    ...                       {'name': 'Broward', 'population': 40000},
    ...                       {'name': 'Palm Beach', 'population': 60000}]},
    ...         {'state': 'Ohio',
    ...          'shortname': 'OH',
    ...          'info': {'governor': 'John Kasich'},
    ...          'counties': [{'name': 'Summit', 'population': 1234},
    ...                       {'name': 'Cuyahoga', 'population': 1337}]}]
    >>> result = json_normalize(data, 'counties', ['state', 'shortname',
    ...                                            ['info', 'governor']])
    >>> result
             name  population    state shortname info.governor
    0        Dade       12345   Florida    FL    Rick Scott
    1     Broward       40000   Florida    FL    Rick Scott
    2  Palm Beach       60000   Florida    FL    Rick Scott
    3      Summit        1234   Ohio       OH    John Kasich
    4    Cuyahoga        1337   Ohio       OH    John Kasich

    >>> data = {'A': [1, 2]}
    >>> json_normalize(data, 'A', record_prefix='Prefix.')
        Prefix.0
    0          1
    1          2

    Returns normalized data with columns prefixed with the given string.
    """

    def _pull_field(
        js: Dict[str, Any], spec: Union[List, str]
    ) -> Union[Scalar, Iterable]:
        """Internal function to pull field"""
        result = js  # type: ignore
        if isinstance(spec, list):
            for field in spec:
                result = result[field]
        else:
            result = result[spec]
        return result

    def _pull_records(js: Dict[str, Any], spec: Union[List, str]) -> Iterable:
        """
        Interal function to pull field for records, and similar to
        _pull_field, but require to return Iterable. And will raise error
        if has non iterable value.
        """
        result = _pull_field(js, spec)

        # GH 31507 GH 30145, if result is not Iterable, raise TypeError if not
        # null, otherwise return an empty list
        if not isinstance(result, Iterable):
            if pd.isnull(result):
                result = []  # type: ignore
            else:
                raise TypeError(
                    f"{js} has non iterable value {result} for path {spec}. "
                    "Must be iterable or null."
                )
        return result

    if isinstance(data, list) and not data:
        return DataFrame()

    # A bit of a hackjob
    if isinstance(data, dict):
        data = [data]

    if record_path is None:
        if any([isinstance(x, dict) for x in y.values()] for y in data):
            # naive normalization, this is idempotent for flat records
            # and potentially will inflate the data considerably for
            # deeply nested structures:
            #  {VeryLong: { b: 1,c:2}} -> {VeryLong.b:1 ,VeryLong.c:@}
            #
            # TODO: handle record value which are lists, at least error
            #       reasonably
            data = nested_to_record(data, sep=sep, max_level=max_level)
        return DataFrame(data)
    elif not isinstance(record_path, list):
        record_path = [record_path]

    if meta is None:
        meta = []
    elif not isinstance(meta, list):
        meta = [meta]

    _meta = [m if isinstance(m, list) else [m] for m in meta]

    # Disastrously inefficient for now
    records: List = []
    lengths = []

    meta_vals: DefaultDict = defaultdict(list)
    meta_keys = [sep.join(val) for val in _meta]

    def _recursive_extract(data, path, seen_meta, level=0):
        if isinstance(data, dict):
            data = [data]
        if len(path) > 1:
            for obj in data:
                for val, key in zip(_meta, meta_keys):
                    if level + 1 == len(val):
                        seen_meta[key] = _pull_field(obj, val[-1])

                _recursive_extract(obj[path[0]], path[1:], seen_meta, level=level + 1)
        else:
            for obj in data:
                recs = _pull_records(obj, path[0])
                recs = [
                    nested_to_record(r, sep=sep, max_level=max_level)
                    if isinstance(r, dict)
                    else r
                    for r in recs
                ]

                # For repeating the metadata later
                lengths.append(len(recs))
                for val, key in zip(_meta, meta_keys):
                    if level + 1 > len(val):
                        meta_val = seen_meta[key]
                    else:
                        try:
                            meta_val = _pull_field(obj, val[level:])
                        except KeyError as e:
                            if errors == "ignore":
                                meta_val = np.nan
                            else:
                                raise KeyError(
                                    "Try running with "
                                    "errors='ignore' as key "
                                    f"{e} is not always present"
                                )
                    meta_vals[key].append(meta_val)
                records.extend(recs)

    _recursive_extract(data, record_path, {}, level=0)

    result = DataFrame(records)

    if record_prefix is not None:
        result = result.rename(columns=lambda x: f"{record_prefix}{x}")

    # Data types, a problem
    for k, v in meta_vals.items():
        if meta_prefix is not None:
            k = meta_prefix + k

        if k in result:
            raise ValueError(
                f"Conflicting metadata name {k}, need distinguishing prefix "
            )
        result[k] = np.array(v, dtype=object).repeat(lengths)
    return result


json_normalize = deprecate(
    "pandas.io.json.json_normalize", _json_normalize, "1.0.0", "pandas.json_normalize"
)
