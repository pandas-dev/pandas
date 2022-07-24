# ---------------------------------------------------------------------
# JSON normalization routines
from __future__ import annotations

from collections import (
    abc,
    defaultdict,
)
import copy
from typing import (
    Any,
    DefaultDict,
    Iterable,
)

import numpy as np

from pandas._libs.writers import convert_json_to_lines
from pandas._typing import * # Other imports raise errors
from pandas.util._decorators import deprecate

import pandas as pd
from pandas import DataFrame

import numpy as np
# import json

# Helper Functions
def get_dict_path(kh: dict, ignore: dict, ig_key: str)-> dict[str, Any]:
   """
   Helper Function to get all parent keys up to the point
   of deletion in a dictionary. Reorginises current disorganised path dict (kh)
   and checks for ignore values on the same and consecutive levels.
   """
   new_keys = []
   new_lvls = []
   ignore_count = 0

   set_lvls = set(list(kh.values()))
   for lvl in set_lvls:
      idxs = [(idx_in_dict, key) for idx_in_dict, (key, val) in enumerate(list(kh.items())) if val == lvl]
      idx_ig_key = [i for i, key in idxs if key == ig_key]
      idx_ig_key = idx_ig_key[0] if len(idx_ig_key) >= 1 else idx_ig_key
      temp_kh = {k: v for k, v in list(kh.items()) if v == lvl}
      
      # Finds num ignore cols in one level
      ig_sl_count = 0
      for key in list(temp_kh):
         if key in ignore["cols"]:
            ig_sl_count += 1
            if ig_sl_count > 1:
               break
      
      if not ig_sl_count > 1:
         # If cols to ignore exist on consecutive lvls, keep the first col
         new_keys.append(list(kh)[idxs[-1][0]])
         new_lvls.append(list(kh.values())[idxs[-1][0]])
         if list(kh)[idxs[-1][0]] in ignore["cols"]:
            ignore_count += 1
         if ignore_count > 1:
            # Prioritises first ignore seen for consecutive lvls
            new_lvls = new_lvls[:-1] 
            new_keys = new_keys[:-1]
      else:
         new_keys.append(list(kh)[idx_ig_key])
         new_lvls.append(list(kh.values())[idx_ig_key])

   # Resetting history
   kh = {} 
   kh = {k: v for k, v in zip(new_keys, new_lvls)}
   return kh

def get_locs(kh: dict)-> list:
   """
   Function to get strings of all parent keys up to del point
   """
   return [[key for key in list(kh)],\
            [lvl for lvl in list(kh.values())]] 

def get_multi_dict(splits: list, val: Any)-> dict[str, Any]:
   """
   Function to recursively create multilevel dicitonary
   from a list of key locations. Takes desired value to be put at that location.
   """
   if len(splits) == 0:
      return val
   elif isinstance(splits, list):
      first_value = splits[0] 
      splits_dict = {first_value : get_multi_dict(splits[1:], val)}
      return splits_dict

# Helper Functions to separte each deleted val into a separted dictionary, 
# with user chosen lvls for names. These operate like a pivot table (in sql).
def locs_to_val(d_dict: dict, loc_arr: list, lvl: int = 0)-> Any:
   """
   Uses recursion to grab value of dictionary given a list, its the reverse of
   get_multi_dict(), and gets a value rather than creating a dictionary.
   """
   if isinstance(d_dict, dict) and lvl != len(loc_arr)-1:
      for key, val in d_dict.items():
         if isinstance(val, dict):
            if lvl < len(loc_arr):
               if key == loc_arr[lvl]:
                  new_val = val
                  return locs_to_val(new_val, loc_arr, lvl + 1)
   elif isinstance(d_dict, dict) and lvl == len(loc_arr)-1:
      return d_dict[loc_arr[-1]]
   return d_dict[loc_arr[-1]]

def update_dict(del_dict: dict, temp_dict: dict, replace: bool = True, add_val: bool = False,\
            sep: str = "_", suffix: str = "0", prefix: str = "")-> dict[str, Any]:
   """
   Function that updates multilevel dictionary without resetting values. Much like an
   implementation of the list.append() function but for dictionaries. User can
   specify whether to replace, ignore or append certain additions.
   """
   if isinstance(temp_dict, dict):
      for k, v in temp_dict.items():
         if not isinstance(k, str):
            k = str(k)
         if isinstance(v, dict) and isinstance(del_dict.get(k, {}), dict):
            del_dict[k] = {**del_dict.get(k, {}), **update_dict(del_dict.get(k, {}), v)}
         elif isinstance(v, tuple):
            del_dict[k] = del_dict.get(k, ()) + v
         elif isinstance(v, list):
            if isinstance(del_dict.get(k, ()), tuple):
               del_dict[k] = [*list(del_dict.get(k, ())), *v]
            elif isinstance(del_dict.get(k, []), list):
               del_dict[k] = [*del_dict.get(k, []), *v]
         else:
            if k in list(del_dict):
               if replace and not add_val:
                  del_dict[k] = v
               elif not replace and add_val:
                  del_dict[prefix + k + sep + suffix] = v
               elif not replace and not add_val:
                  pass
            else:
               del_dict[k] = v
   return del_dict

def get_del(ignore_col: str, del_dict: dict, dels: list, max_level: int or None = None,\
   ignore_loc: bool = False)-> dict[str, Any]:
   """
   Helper Function to create a "pivot_table" (see earlier) of the delete values.
   It uses recursion in locs_to_val() to search for a value to put in this new
   ordered dictionary. The result is appended to the main dicitonary containing all the pivots.
   It can either be used with a list of location lists or a sinlge location list (to reduce
   time and soze complexities for large dictionaries).
   """
   # max_level always starts from 0
   dict_cols = {}
   dels = np.array(dels)

   if len(dels.shape) == 2: # If a list of paths is passed
      for idx in range(len(dels)):
         if ignore_col in dels[idx]:
            if len(dels[idx]) >= (max_level):
               dict_cols[str(dels[idx][max_level])] = locs_to_val(del_dict, dels[idx])
            elif max_level == None:
               # If max_level is none, use lowest level available
               dict_cols[str(dels[idx][0])] = locs_to_val(del_dict, dels[idx])
            else:
               if not ignore_loc:
                  # If max_level is greater than len(loc_lvls) use level at maximum idx
                  dict_cols[str(dels[idx][-1])] = locs_to_val(del_dict, dels[idx])
               else:
                  # Ignore a column if incorrect max_level supplied
                  pass
   elif len(dels.shape) == 1: # If a single path is passed
      if ignore_col in dels:
         if len(dels) >= (max_level):
            dict_cols[str(dels[max_level])] = locs_to_val(del_dict, dels)
         elif max_level == None:
            # If max_level is none, use lowest level available
            dict_cols[str(dels[0])] = locs_to_val(del_dict, dels)
         else:
            if not ignore_loc:
               # If max_level is greater than len(loc_lvls) use level at maximum idx
               dict_cols[str(dels[-1])] = locs_to_val(del_dict, dels)
            else:
               # Ignore a column if incorrect max_level supplied
               pass
   return dict_cols
   
def convert_to_line_delimits(s: str) -> str:
    """
    Helper function that converts JSON lists to line delimited JSON.
    """
    # Determine we have a JSON list to turn to lines otherwise just return the
    # json object, only lists can
    if not s[0] == "[" and s[-1] == "]":
        return s
    s = s[1:-1]

    return convert_json_to_lines(s)

# Returns normalised_json, (unordered_dict_of_deleted_vals, locations_of_deleted_vals,
# pivot_tables_of_deleted_vals, path_idx)
def nested_ignore_cols_to_record(
   ds,
   prefix: str = "",
   sep: str = ".",
   level: int = 0,
   max_level: int or None = None,
   keys_hist: dict = {},
   del_dict: dict = {},
   dels: list = [],
   pivot_dels: dict = {},
   return_dels: bool = False,
   path_idx: int = 0,
   first_update: bool = True,
   ignore: dict = {"cols": None,\
                   "name_lvls": None}
):
   """
    A more commplex version of nested_to_record(), user can pass in columns to
    ignore in the input dictionary (before flattening).

    Can also return the deleted values so they can be used for other cases.

    If no column to be ignored is specified, function just performs the 
    un-edited nested_to_record().

    Only truly deletes values if no max_level is specified (max_level == None).

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
    keys_hist: dict, optional, default: {}
        This saves the path to a value to be deleted, gets reset
        after the value gets deleted.
    del_dict: dict, optional, default: {}
        This saves the values deleted from the original data in their
        own dictionary. THIS DICTIONARY IS UNORDERED.
    dels: list, optional, default: []
        This is an extension of del_dict, saves the locations of each deleted
        value in an array as well as their levels leading up to the deleted value.
    pivot_dels: dict, optional, default: {}
        This is the ordered version of del_dict. The top level keys are the deleted value
        keys/columns and the values are dictionaries with user-specified "name_lvl" keys
        in ignore (see below).
    return_dels: bool, optional, default: True
        If true, returns the deleted values, as well as the flattened input dictionary,
        in the format: flattened_dict, (del_dict, dels, pivot_dels, path_idx)
    path_idx: int, default: 0
        Internal variable used to find the updated running idx of values in dels.
    ignore: dict, optional, default: {"cols": None, "name_lvls": None}
        This contains the columns to be deleted and the desired levels (in the array of the
        locations of deleted values) to be used as keys in the level below the top level of
        the ordered pivot_dels dictionary (see above).

    .. versionadded:: x.xx.x --- Please fill in

    Returns
    -------
    1. d - dict or list of dicts, matching `ds`
    or
    2. (1.) and del_tuple - tuple of (del_dict, dels, pivot_dels, path_idx) (see above)

    Example (1.)
    --------
    nested_ignore_cols_to_record(
      dict(flat1=1, dict1=dict(c=1, d=2), nested=dict(e=dict(c=1, d=2), d=2))

            {
        'flat1': 1, 
        'dict1.c': 1, 
        'dict1.d': 2, 
        'nested.e.c': 1, 
        'nested.e.d': 2, 
        'nested.d': 2
        }

    Example (2.)
    --------
    nested_ignore_cols_to_record(
         dict(flat1=1, dict1=dict(c=1, d=2), nested=dict(e=dict(c=1, d=2), d=2),
                ignore = {"cols": ["e"], "name_lvls": None},
                return_dels = True)

       returns: NORMALISED DICT WITH 'e' col ignored:
               {'flat1': 1, 'dict1.c': 1, 'dict1.d': 2, 'nested.d': 2},

                FULL_DELS_TUPLE:
               ({'nested': {'e': {'c': 1, 'd': 2}}}, [['nested', 'e']], {'e': {'nested': {'c': 1, 'd': 2}}}, 1)
               
               FULL_DELS_TUPLE[0] = del_dict
               FULL_DELS_TUPLE[1] = dels
               FULL_DELS_TUPLE[2] = pivot_dels
               FULL_DELS_TUPLE[3] = path_idx
   """
   singleton = False
   if isinstance(ds, dict):
      ds = [ds]
      singleton = True
   new_ds = []
   for d in ds:
      new_d = copy.deepcopy(d)
      same_key = 0
      path_arr = []
      for k, v in d.items():
         # each key gets renamed with prefix
         if not isinstance(k, str):
            k = str(k)

         if len(list(keys_hist)) > 0:
            # Checking for same consecutive keys
            if k == list(keys_hist)[-1]:
               same_key += 1
               keys_hist[k + "_" + str(same_key)] = level
            else:
               keys_hist[k] = level

         if ignore["cols"] is None:
            if level == 0:
               newkey = k
               keys_hist[k] = level
            else:
               newkey = prefix + sep + k
         elif k not in ignore["cols"]:
            if level == 0:
               newkey = k
               keys_hist[k] = level
            else:
               newkey = prefix + sep + k
         elif k in ignore["cols"]:
            if first_update == True:
               del_dict = {}
               pivot_dels = {}
               dels = []
               path_arr = []

               keys_hist = get_dict_path(keys_hist, ignore, k)
               locs, _ = get_locs(keys_hist)
               dels.append(locs)
               path_arr.append((locs, path_idx))
               path_idx += 1
               first_update = False
            else:
               keys_hist = get_dict_path(keys_hist, ignore, k)
               locs, _ = get_locs(keys_hist)
               dels.append(locs)
               path_arr.append((locs, path_idx))
               path_idx += 1
            continue

         # flatten if type is dict and
         # current dict level < maximum level provided and
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
            keys_hist = get_dict_path(keys_hist, ignore, k)

            if return_dels == True:
               nest_dict, d_tuple = \
               nested_ignore_cols_to_record(v, newkey, sep, level + 1, max_level,\
                           keys_hist, del_dict, dels, pivot_dels, return_dels, path_idx, first_update = first_update, ignore = ignore)
               del_dict = d_tuple[0]
               dels = d_tuple[1]
               pivot_dels = d_tuple[2]
               path_idx = d_tuple[3]
               
               # Updating first_update
               for item in d_tuple:
                  if isinstance(item, dict):
                     if len(item.items()) > 0:
                        first_update = False
                  elif isinstance(item, list):
                     if len(item) > 0:
                        first_update = False
                  if isinstance(item, int) or isinstance(item, float):
                     if item > 0:
                        first_update = False
                  break
            else:
               nest_dict = \
               nested_ignore_cols_to_record(v, newkey, sep, level + 1, max_level,\
                           keys_hist, del_dict, dels, pivot_dels, return_dels, path_idx, first_update, ignore = ignore)
            new_d.update(nest_dict)
      new_ds.append(new_d)
   
   if singleton:
      if (ignore["cols"] is not None and len(path_arr) > 0
          and ignore["name_lvls"] is not None
          ):
         # Resetting all dels_tup for repeated implementations
         del_tups = [t[0] for t in path_arr]
         temp_d = {}
         for del_tup in del_tups:
            # Assuming del_val is last in each del_tup array
            last_val = del_tup[-1] 
            if last_val in list(new_ds[0]) and last_val in ignore["cols"]:
               ig_ml_tup = [(ig, ml) for ig, ml in zip(ignore["cols"], ignore["name_lvls"])\
                           if ig == last_val]
               ignore_col = ig_ml_tup[0][0]
               ml = ig_ml_tup[0][1]
               val = new_ds[0].pop(last_val)
               temp = get_multi_dict(del_tup, val)
               del_dict = update_dict(del_dict, temp)

               if ignore["name_lvls"] is not None:
                  temp_d[ignore_col] = get_del(ignore_col, del_dict, del_tup, ml)
               else:
                  temp_d[ignore_col] = get_del(ignore_col, del_dict, del_tup, 0)

         pivot_dels = update_dict(pivot_dels, temp_d)

      if return_dels == True: # User defined
         return new_ds[0], (del_dict, dels, pivot_dels, path_idx)
      return new_ds[0]

   if (ignore["cols"] is not None and len(path_arr) > 0
       and ignore["name_lvls"] is not None
          ):
      del_tups = [t[0] for t in path_arr]
      temp_d = {}
      for del_tup in del_tups:
         # Assuming del_val is last in each del_tup array
         last_val = del_tup[-1] 
         if last_val in list(new_ds) and last_val in ignore["cols"]:
            ig_ml_tup = [(ig, ml) for ig, ml in zip(ignore["cols"], ignore["name_lvls"])\
                           if ig == last_val]
            ignore_col = ig_ml_tup[0][0]
            ml = ig_ml_tup[0][1]
            val = new_ds.pop(last_val)
            temp = get_multi_dict(del_tup, val)
            del_dict = update_dict(del_dict, temp)

            if ignore["name_lvls"] is not None:
               temp_d[ignore_col] = get_del(ignore_col, del_dict, del_tup, ml)
            else:
               temp_d[ignore_col] = get_del(ignore_col, del_dict, del_tup, 0)

      pivot_dels = update_dict(pivot_dels, temp_d)

   if return_dels == True: # User defined
      return new_ds, (del_dict, dels, pivot_dels, path_idx)
   return new_ds

def nested_to_record(
    ds,
    prefix: str = "",
    sep: str = ".",
    level: int = 0,
    max_level: int | None = None,
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
    >>> nested_to_record(
    ...     dict(flat1=1, dict1=dict(c=1, d=2), nested=dict(e=dict(c=1, d=2), d=2))
    ... )
    {\
'flat1': 1, \
'dict1.c': 1, \
'dict1.d': 2, \
'nested.e.c': 1, \
'nested.e.d': 2, \
'nested.d': 2\
}
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


def _normalise_json(
    data: Any,
    key_string: str,
    normalized_dict: dict[str, Any],
    separator: str,
) -> dict[str, Any]:
    """
    Main recursive function
    Designed for the most basic use case of pd.json_normalize(data)
    intended as a performance improvement, see #15621

    Parameters
    ----------
    data : Any
        Type dependent on types contained within nested Json
    key_string : str
        New key (with separator(s) in) for data
    normalized_dict : dict
        The new normalized/flattened Json dict
    separator : str, default '.'
        Nested records will generate names separated by sep,
        e.g., for sep='.', { 'foo' : { 'bar' : 0 } } -> foo.bar
    """
    if isinstance(data, dict):
        for key, value in data.items():
            new_key = f"{key_string}{separator}{key}"
            _normalise_json(
                data=value,
                # to avoid adding the separator to the start of every key
                # GH#43831 avoid adding key if key_string blank
                key_string=new_key
                if new_key[: len(separator)] != separator
                else new_key[len(separator) :],
                normalized_dict=normalized_dict,
                separator=separator,
            )
    else:
        normalized_dict[key_string] = data
    return normalized_dict


def _normalise_json_ordered(data: dict[str, Any], separator: str) -> dict[str, Any]:
    """
    Order the top level keys and then recursively go to depth

    Parameters
    ----------
    data : dict or list of dicts
    separator : str, default '.'
        Nested records will generate names separated by sep,
        e.g., for sep='.', { 'foo' : { 'bar' : 0 } } -> foo.bar

    Returns
    -------
    dict or list of dicts, matching `normalised_json_object`
    """
    top_dict_ = {k: v for k, v in data.items() if not isinstance(v, dict)}
    nested_dict_ = _normalise_json(
        data={k: v for k, v in data.items() if isinstance(v, dict)},
        key_string="",
        normalized_dict={},
        separator=separator,
    )
    return {**top_dict_, **nested_dict_}


def _simple_json_normalize(
    ds: dict | list[dict],
    sep: str = ".",
) -> dict | list[dict] | Any:
    """
    A optimized basic json_normalize

    Converts a nested dict into a flat dict ("record"), unlike
    json_normalize and nested_to_record it doesn't do anything clever.
    But for the most basic use cases it enhances performance.
    E.g. pd.json_normalize(data)

    Parameters
    ----------
    ds : dict or list of dicts
    sep : str, default '.'
        Nested records will generate names separated by sep,
        e.g., for sep='.', { 'foo' : { 'bar' : 0 } } -> foo.bar

    Returns
    -------
    frame : DataFrame
    d - dict or list of dicts, matching `normalised_json_object`

    Examples
    --------
    >>> _simple_json_normalize(
    ...     {
    ...         "flat1": 1,
    ...         "dict1": {"c": 1, "d": 2},
    ...         "nested": {"e": {"c": 1, "d": 2}, "d": 2},
    ...     }
    ... )
    {\
'flat1': 1, \
'dict1.c': 1, \
'dict1.d': 2, \
'nested.e.c': 1, \
'nested.e.d': 2, \
'nested.d': 2\
}

    """
    normalised_json_object = {}
    # expect a dictionary, as most jsons are. However, lists are perfectly valid
    if isinstance(ds, dict):
        normalised_json_object = _normalise_json_ordered(data=ds, separator=sep)
    elif isinstance(ds, list):
        normalised_json_list = [_simple_json_normalize(row, sep=sep) for row in ds]
        return normalised_json_list
    return normalised_json_object


def _json_normalize(
    data: dict | list[dict],
    record_path: str | list | None = None,
    meta: str | list[str | list[str]] | None = None,
    meta_prefix: str | None = None,
    record_prefix: str | None = None,
    errors: IgnoreRaise = "raise",
    sep: str = ".",
    max_level: int | None = None,
) -> DataFrame:
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
    >>> data = [
    ...     {"id": 1, "name": {"first": "Coleen", "last": "Volk"}},
    ...     {"name": {"given": "Mark", "family": "Regner"}},
    ...     {"id": 2, "name": "Faye Raker"},
    ... ]
    >>> pd.json_normalize(data)
        id name.first name.last name.given name.family        name
    0  1.0     Coleen      Volk        NaN         NaN         NaN
    1  NaN        NaN       NaN       Mark      Regner         NaN
    2  2.0        NaN       NaN        NaN         NaN  Faye Raker

    >>> data = [
    ...     {
    ...         "id": 1,
    ...         "name": "Cole Volk",
    ...         "fitness": {"height": 130, "weight": 60},
    ...     },
    ...     {"name": "Mark Reg", "fitness": {"height": 130, "weight": 60}},
    ...     {
    ...         "id": 2,
    ...         "name": "Faye Raker",
    ...         "fitness": {"height": 130, "weight": 60},
    ...     },
    ... ]
    >>> pd.json_normalize(data, max_level=0)
        id        name                        fitness
    0  1.0   Cole Volk  {'height': 130, 'weight': 60}
    1  NaN    Mark Reg  {'height': 130, 'weight': 60}
    2  2.0  Faye Raker  {'height': 130, 'weight': 60}

    Normalizes nested data up to level 1.

    >>> data = [
    ...     {
    ...         "id": 1,
    ...         "name": "Cole Volk",
    ...         "fitness": {"height": 130, "weight": 60},
    ...     },
    ...     {"name": "Mark Reg", "fitness": {"height": 130, "weight": 60}},
    ...     {
    ...         "id": 2,
    ...         "name": "Faye Raker",
    ...         "fitness": {"height": 130, "weight": 60},
    ...     },
    ... ]
    >>> pd.json_normalize(data, max_level=1)
        id        name  fitness.height  fitness.weight
    0  1.0   Cole Volk             130              60
    1  NaN    Mark Reg             130              60
    2  2.0  Faye Raker             130              60

    >>> data = [
    ...     {
    ...         "state": "Florida",
    ...         "shortname": "FL",
    ...         "info": {"governor": "Rick Scott"},
    ...         "counties": [
    ...             {"name": "Dade", "population": 12345},
    ...             {"name": "Broward", "population": 40000},
    ...             {"name": "Palm Beach", "population": 60000},
    ...         ],
    ...     },
    ...     {
    ...         "state": "Ohio",
    ...         "shortname": "OH",
    ...         "info": {"governor": "John Kasich"},
    ...         "counties": [
    ...             {"name": "Summit", "population": 1234},
    ...             {"name": "Cuyahoga", "population": 1337},
    ...         ],
    ...     },
    ... ]
    >>> result = pd.json_normalize(
    ...     data, "counties", ["state", "shortname", ["info", "governor"]]
    ... )
    >>> result
             name  population    state shortname info.governor
    0        Dade       12345   Florida    FL    Rick Scott
    1     Broward       40000   Florida    FL    Rick Scott
    2  Palm Beach       60000   Florida    FL    Rick Scott
    3      Summit        1234   Ohio       OH    John Kasich
    4    Cuyahoga        1337   Ohio       OH    John Kasich

    >>> data = {"A": [1, 2]}
    >>> pd.json_normalize(data, "A", record_prefix="Prefix.")
        Prefix.0
    0          1
    1          2

    Returns normalized data with columns prefixed with the given string.
    """

    def _pull_field(
        js: dict[str, Any], spec: list | str, extract_record: bool = False
    ) -> Scalar | Iterable:
        """Internal function to pull field"""
        result = js
        try:
            if isinstance(spec, list):
                for field in spec:
                    if result is None:
                        raise KeyError(field)
                    result = result[field]
            else:
                result = result[spec]
        except KeyError as e:
            if extract_record:
                raise KeyError(
                    f"Key {e} not found. If specifying a record_path, all elements of "
                    f"data should have the path."
                ) from e
            elif errors == "ignore":
                return np.nan
            else:
                raise KeyError(
                    f"Key {e} not found. To replace missing values of {e} with "
                    f"np.nan, pass in errors='ignore'"
                ) from e

        return result

    def _pull_records(js: dict[str, Any], spec: list | str) -> list:
        """
        Internal function to pull field for records, and similar to
        _pull_field, but require to return list. And will raise error
        if has non iterable value.
        """
        result = _pull_field(js, spec, extract_record=True)

        # GH 31507 GH 30145, GH 26284 if result is not list, raise TypeError if not
        # null, otherwise return an empty list
        if not isinstance(result, list):
            if pd.isnull(result):
                result = []
            else:
                raise TypeError(
                    f"{js} has non list value {result} for path {spec}. "
                    "Must be list or null."
                )
        return result

    if isinstance(data, list) and not data:
        return DataFrame()
    elif isinstance(data, dict):
        # A bit of a hackjob
        data = [data]
    elif isinstance(data, abc.Iterable) and not isinstance(data, str):
        # GH35923 Fix pd.json_normalize to not skip the first element of a
        # generator input
        data = list(data)
    else:
        raise NotImplementedError

    # check to see if a simple recursive function is possible to
    # improve performance (see #15621) but only for cases such
    # as pd.Dataframe(data) or pd.Dataframe(data, sep)
    if (
        record_path is None
        and meta is None
        and meta_prefix is None
        and record_prefix is None
        and max_level is None
    ):
        return DataFrame(_simple_json_normalize(data, sep=sep))

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
    records: list = []
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
                        meta_val = _pull_field(obj, val[level:])
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
