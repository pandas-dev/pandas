# #############################################################################
# Copyright 2018 Hoffmann-La Roche
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# #############################################################################

from collections.abc import Callable, Iterator
import multiprocessing as mp
from itertools import chain
from os import PathLike
from typing import TYPE_CHECKING, Any, Literal, TypeAlias, overload, Protocol #, Concatenate: see later

import narwhals.stable.v2 as nw

from ._readstat_parser import parser_entry_point, PyreadstatError
from ._readstat_writer import writer_entry_point
from .worker import worker
from .pyclasses import metadata_container, MissingRange
from .pyfunctions import set_value_labels, set_catalog_to_sas

# Typing interface

if TYPE_CHECKING:
    from typing import Concatenate  # TODO: move back to top-level import when dropping Python 3.10 support
    # Setup type aliases for the public interface.
    # These are not executed at runtime, but they help type checkers understand
    # the expected types of the public functions and classes.

    # Since pyreadstat can work with both pandas and polars, we define a DataFrame type that can be either.
    try:
        from pandas import DataFrame as PandasDataFrame  # type: ignore
    except ImportError:
        # Define a dummy DataFrame class to avoid accepting any type as PandasDataFrame when pandas is not installed
        class PandasDataFrame:
            pass

    try:
        from polars import DataFrame as PolarsDataFrame  # type: ignore
    except ImportError:
        # Define a dummy DataFrame class to avoid accepting any type as PolarsDataFrame when polars is not installed
        class PolarsDataFrame:
            pass

DataFrame: TypeAlias = "PandasDataFrame | PolarsDataFrame"  # Define type at runtime for introspection

class FileLike(Protocol):
    """Protocol for file-like objects accepted by pyreadstat"""

    # Should work with any file-like object that has read and seek methods, such as those returned by open() or io.BytesIO
    def read(self, size: int | None = -1, /) -> bytes: ...
    def seek(self, pos: int, whence: int = 0, /) -> int: ...


FilePathLike: TypeAlias = str | bytes | PathLike[str] | PathLike[bytes]
FilePathorBuffer: TypeAlias = FilePathLike | FileLike

DictOutput: TypeAlias = dict[str, list[Any]]

# TODO: when dropping Python 3.10 support, remove the string quotes and move Concatenate back to the top-level import:
#   PyreadstatReadFunction: TypeAlias = Callable[Concatenate[FilePathorBuffer, ...], tuple[DataFrame | DictOutput, metadata_container]]
PyreadstatReadFunction: TypeAlias = "Callable[Concatenate[FilePathorBuffer, ...], tuple[DataFrame | DictOutput, metadata_container]]"


# Public interface

# Parsing functions


@overload
def read_sas7bdat(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    catalog_file: FilePathorBuffer | None = ...,
    formats_as_category: bool = ...,
    formats_as_ordered_category: bool = ...,
    encoding: str | None = ...,
    usecols: list[str] | None = ...,
    user_missing: bool = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["pandas"] | None = ...,
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> "tuple[PandasDataFrame, metadata_container]": ...
@overload
def read_sas7bdat(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    catalog_file: FilePathorBuffer | None = ...,
    formats_as_category: bool = ...,
    formats_as_ordered_category: bool = ...,
    encoding: str | None = ...,
    usecols: list[str] | None = ...,
    user_missing: bool = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["polars"] = "polars",
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> "tuple[PolarsDataFrame, metadata_container]": ...
@overload
def read_sas7bdat(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    catalog_file: FilePathorBuffer | None = ...,
    formats_as_category: bool = ...,
    formats_as_ordered_category: bool = ...,
    encoding: str | None = ...,
    usecols: list[str] | None = ...,
    user_missing: bool = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["dict"] = "dict",
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> tuple[DictOutput, metadata_container]: ...
def read_sas7bdat(
    filename_path: FilePathorBuffer,
    metadataonly: bool = False,
    dates_as_pandas_datetime: bool = False,
    catalog_file: FilePathorBuffer | None = None,
    formats_as_category: bool = True,
    formats_as_ordered_category: bool = False,
    encoding: str | None = None,
    usecols: list[str] | None = None,
    user_missing: bool = False,
    disable_datetime_conversion: bool = False,
    row_limit: int = 0,
    row_offset: int = 0,
    output_format: Literal["pandas", "polars", "dict"] | None = None,
    extra_datetime_formats: list[str] | None = None,
    extra_date_formats: list[str] | None = None,
    extra_time_formats: list[str] | None = None,
) -> "tuple[DataFrame | DictOutput, metadata_container]":
    r"""
    Read a SAS sas7bdat file.
    It accepts the path to a sas7bcat.

    Parameters
    ----------
        filename_path : str, bytes, Path-like object or file-like object
            path to the file or file-like object. In python 2.7 the string is assumed to be utf-8 encoded.
        metadataonly : bool, optional
            by default False. IF true, no data will be read but only metadata, so that you can get all elements in the
            metadata object. The data frame will be set with the correct column names but no data.
        dates_as_pandas_datetime : bool, optional
            by default False. If true dates will be transformed to pandas datetime64 instead of date, effective only for pandas.
        catalog_file : str, bytes, Path-like object or file-like object, optional
            path to a sas7bcat file or file-like object. By default is None. If not None, will parse the catalog file and replace the values
            by the formats in the catalog, if any appropiate is found. If this is not the behavior you are looking for,
            Use read_sas7bcat to parse the catalog independently
            of the sas7bdat and set_catalog_to_sas to apply the resulting format into sas7bdat files.
        formats_as_category : bool, optional
            Will take effect only if the catalog_file was specified. If True the variables whose values were replaced
            by the formats will be transformed into categories.
        formats_as_ordered_category : bool, optional
            defaults to False. If True the variables having formats will be transformed into ordered categories/enums.
            it has precedence over formats_as_category, meaning if this is True, it will take effect irrespective of
            the value of formats_as_category.
        encoding : str, optional
            Defaults to None. If set, the system will use the defined encoding instead of guessing it. It has to be an
            iconv-compatible name
        usecols : list, optional
            a list with column names to read from the file. Only those columns will be imported. Case sensitive!
        user_missing : bool, optional
            by default False, in this case user defined missing values are delivered as nan. If true, the missing values
            will be deliver as is, and an extra piece of information will be set in the metadata (missing_user_values)
            to be able to interpret those values as missing.
        disable_datetime_conversion : bool, optional
            if True pyreadstat will not attempt to convert dates, datetimes and times to python objects but those columns
            will remain as numbers. In order to convert them later to an appropiate python object, the user can use the
            information about the original variable format stored in the metadata object in original_variable_types.
            Disabling datetime conversion speeds up reading files. In addition it helps to overcome situations where
            there are datetimes that are beyond the limits of python datetime (which is limited to year 10,000, dates
            beyond that will rise an Overflow error in pyreadstat).
        row_limit : int, optional
            maximum number of rows to read. The default is 0 meaning unlimited.
        row_offset : int, optional
            start reading rows after this offset. By default 0, meaning start with the first row not skipping anything.
        output_format : str, optional
            one of 'pandas' (default), 'polars' or 'dict'. If 'dict' a dictionary with numpy arrays as values will be returned, the
            user can then convert it to her preferred data format. Using dict is faster as the other types as the conversion to a
            dataframe is avoided.
        extra_datetime_formats: list of str, optional
            formats to be parsed as python datetime objects
        extra_date_formats: list of str, optional
            formats to be parsed as python date objects
        extra_time_formats: list of str, optional
            formats to be parsed as python time objects


    Returns
    -------
        data_frame : dataframe or dict
            a dataframe or dict with the data.
        metadata
            object with metadata. The members variables_value_labels will be empty unless a valid catalog file is
            supplied.
            Look at the documentation for more information.
    """
    parser_format = "sas7bdat"
    data_frame, metadata = parser_entry_point(
        filename_path,
        parser_format,
        metadataonly=metadataonly,
        dates_as_pandas_datetime=dates_as_pandas_datetime,
        formats_as_category=formats_as_category,
        formats_as_ordered_category=formats_as_ordered_category,
        encoding=encoding,
        usecols=usecols,
        user_missing=user_missing,
        disable_datetime_conversion=disable_datetime_conversion,
        row_limit=row_limit,
        row_offset=row_offset,
        output_format=output_format,
        extra_datetime_formats=extra_datetime_formats,
        extra_date_formats=extra_date_formats,
        extra_time_formats=extra_time_formats,
    )

    metadata.file_format = parser_format

    if catalog_file:
        _, catalog = read_sas7bcat(catalog_file, encoding=encoding)
        data_frame, metadata = set_catalog_to_sas(
            data_frame,
            metadata,
            catalog,
            formats_as_category=formats_as_category,
            formats_as_ordered_category=formats_as_ordered_category,
        )

    return data_frame, metadata


@overload
def read_xport(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    encoding: str | None = ...,
    usecols: list[str] | None = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["pandas"] | None = ...,
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> "tuple[PandasDataFrame, metadata_container]": ...
@overload
def read_xport(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    encoding: str | None = ...,
    usecols: list[str] | None = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["polars"] = "polars",
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> "tuple[PolarsDataFrame, metadata_container]": ...
@overload
def read_xport(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    encoding: str | None = ...,
    usecols: list[str] | None = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["dict"] = "dict",
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> tuple[DictOutput, metadata_container]: ...
def read_xport(
    filename_path: FilePathorBuffer,
    metadataonly: bool = False,
    dates_as_pandas_datetime: bool = False,
    encoding: str | None = None,
    usecols: list[str] | None = None,
    disable_datetime_conversion: bool = False,
    row_limit: int = 0,
    row_offset: int = 0,
    output_format: Literal["pandas", "polars", "dict"] | None = None,
    extra_datetime_formats: list[str] | None = None,
    extra_date_formats: list[str] | None = None,
    extra_time_formats: list[str] | None = None,
) -> "tuple[DataFrame | DictOutput, metadata_container]":
    r"""
    Read a SAS xport file.

    Parameters
    ----------
        filename_path : str, bytes, Path-like object or file-like object
            path to the file or file-like object. In python 2.7 the string is assumed to be utf-8 encoded.
        metadataonly : bool, optional
            by default False. IF true, no data will be read but only metadata, so that you can get all elements in the
            metadata object. The data frame will be set with the correct column names but no data.
        dates_as_pandas_datetime : bool, optional
            by default False. If true dates will be transformed to pandas datetime64 instead of date, effective only for pandas.
        encoding : str, optional
            Defaults to None. If set, the system will use the defined encoding instead of guessing it. It has to be an
            iconv-compatible name
        usecols : list, optional
            a list with column names to read from the file. Only those columns will be imported. Case sensitive!
        disable_datetime_conversion : bool, optional
            if True pyreadstat will not attempt to convert dates, datetimes and times to python objects but those columns
            will remain as numbers. In order to convert them later to an appropiate python object, the user can use the
            information about the original variable format stored in the metadata object in original_variable_types.
            Disabling datetime conversion speeds up reading files. In addition it helps to overcome situations where
            there are datetimes that are beyond the limits of python datetime (which is limited to year 10,000, dates
            beyond that will rise an Overflow error in pyreadstat).
        row_limit : int, optional
            maximum number of rows to read. The default is 0 meaning unlimited.
        row_offset : int, optional
            start reading rows after this offset. By default 0, meaning start with the first row not skipping anything.
        output_format : str, optional
            one of 'pandas' (default), 'polars' or 'dict'. If 'dict' a dictionary with numpy arrays as values will be returned, the
            user can then convert it to her preferred data format. Using dict is faster as the other types as the conversion to a
            dataframe is avoided.
        extra_datetime_formats: list of str, optional
            formats to be parsed as python datetime objects
        extra_date_formats: list of str, optional
            formats to be parsed as python date objects
        extra_time_formats: list of str, optional
            formats to be parsed as python time objects

    Returns
    -------
        data_frame : dataframe or dict
            a dataframe or dict with the data.
        metadata :
            object with metadata. Look at the documentation for more information.
    """
    parser_format = "xport"
    data_frame, metadata = parser_entry_point(
        filename_path,
        parser_format,
        metadataonly=metadataonly,
        dates_as_pandas_datetime=dates_as_pandas_datetime,
        encoding=encoding,
        usecols=usecols,
        disable_datetime_conversion=disable_datetime_conversion,
        row_limit=row_limit,
        row_offset=row_offset,
        output_format=output_format,
        extra_datetime_formats=extra_datetime_formats,
        extra_date_formats=extra_date_formats,
        extra_time_formats=extra_time_formats,
    )

    metadata.file_format = parser_format

    return data_frame, metadata


@overload
def read_dta(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    apply_value_formats: bool = ...,
    formats_as_category: bool = ...,
    formats_as_ordered_category: bool = ...,
    encoding: str | None = ...,
    usecols: list[str] | None = ...,
    user_missing: bool = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["pandas"] | None = ...,
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> "tuple[PandasDataFrame, metadata_container]": ...
@overload
def read_dta(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    apply_value_formats: bool = ...,
    formats_as_category: bool = ...,
    formats_as_ordered_category: bool = ...,
    encoding: str | None = ...,
    usecols: list[str] | None = ...,
    user_missing: bool = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["polars"] = "polars",
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> "tuple[PolarsDataFrame, metadata_container]": ...
@overload
def read_dta(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    apply_value_formats: bool = ...,
    formats_as_category: bool = ...,
    formats_as_ordered_category: bool = ...,
    encoding: str | None = ...,
    usecols: list[str] | None = ...,
    user_missing: bool = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["dict"] = "dict",
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> tuple[DictOutput, metadata_container]: ...
def read_dta(
    filename_path: FilePathorBuffer,
    metadataonly: bool = False,
    dates_as_pandas_datetime: bool = False,
    apply_value_formats: bool = False,
    formats_as_category: bool = True,
    formats_as_ordered_category: bool = False,
    encoding: str | None = None,
    usecols: list[str] | None = None,
    user_missing: bool = False,
    disable_datetime_conversion: bool = False,
    row_limit: int = 0,
    row_offset: int = 0,
    output_format: Literal["pandas", "polars", "dict"] | None = None,
    extra_datetime_formats: list[str] | None = None,
    extra_date_formats: list[str] | None = None,
    extra_time_formats: list[str] | None = None,
) -> "tuple[DataFrame | DictOutput, metadata_container]":
    r"""
    Read a STATA dta file

    Parameters
    ----------
        filename_path : str, bytes, Path-like object or file-like object
            path to the file or file-like object. In python 2.7 the string is assumed to be utf-8 encoded.
        metadataonly : bool, optional
            by default False. IF true, no data will be read but only metadata, so that you can get all elements in the
            metadata object. The data frame will be set with the correct column names but no data.
        dates_as_pandas_datetime : bool, optional
            by default False. If true dates will be transformed to pandas datetime64 instead of date, effective only for pandas.
        apply_value_formats : bool, optional
            by default False. If true it will change values in the dataframe for they value labels in the metadata,
            if any appropiate are found.
        formats_as_category : bool, optional
            by default True. Takes effect only if apply_value_formats is True. If True, variables with values changed
            for their formatted version will be transformed into categories.
        formats_as_ordered_category : bool, optional
            defaults to False. If True the variables having formats will be transformed into ordered categories/enum.
            it has precedence over formats_as_category, meaning if this is True, it will take effect irrespective of
            the value of formats_as_category.
        encoding : str, optional
            Defaults to None. If set, the system will use the defined encoding instead of guessing it. It has to be an
            iconv-compatible name
        usecols : list, optional
            a list with column names to read from the file. Only those columns will be imported. Case sensitive!
        user_missing : bool, optional
            by default False, in this case user defined missing values are delivered as nan. If true, the missing values
            will be deliver as is, and an extra piece of information will be set in the metadata (missing_user_values)
            to be able to interpret those values as missing.
        disable_datetime_conversion : bool, optional
            if True pyreadstat will not attempt to convert dates, datetimes and times to python objects but those columns
            will remain as numbers. In order to convert them later to an appropiate python object, the user can use the
            information about the original variable format stored in the metadata object in original_variable_types.
            Disabling datetime conversion speeds up reading files. In addition it helps to overcome situations where
            there are datetimes that are beyond the limits of python datetime (which is limited to year 10,000, dates
            beyond that will rise an Overflow error in pyreadstat).
        row_limit : int, optional
            maximum number of rows to read. The default is 0 meaning unlimited.
        row_offset : int, optional
            start reading rows after this offset. By default 0, meaning start with the first row not skipping anything.
        output_format : str, optional
            one of 'pandas' (default), 'polars' or 'dict'. If 'dict' a dictionary with numpy arrays as values will be returned, the
            user can then convert it to her preferred data format. Using dict is faster as the other types as the conversion to a
            dataframe is avoided.
        extra_datetime_formats: list of str, optional
            formats to be parsed as python datetime objects
        extra_date_formats: list of str, optional
            formats to be parsed as python date objects
        extra_time_formats: list of str, optional
            formats to be parsed as python time objects

    Returns
    -------
        data_frame : dataframe or dict
            a dataframe or dict with the data.
        metadata :
            object with metadata. Look at the documentation for more information.
    """
    parser_format = "dta"
    data_frame, metadata = parser_entry_point(
        filename_path,
        parser_format=parser_format,
        metadataonly=metadataonly,
        dates_as_pandas_datetime=dates_as_pandas_datetime,
        formats_as_category=formats_as_category,
        formats_as_ordered_category=formats_as_ordered_category,
        encoding=encoding,
        usecols=usecols,
        user_missing=user_missing,
        disable_datetime_conversion=disable_datetime_conversion,
        row_limit=row_limit,
        row_offset=row_offset,
        output_format=output_format,
        extra_datetime_formats=extra_datetime_formats,
        extra_date_formats=extra_date_formats,
        extra_time_formats=extra_time_formats,
    )

    metadata.file_format = parser_format

    if apply_value_formats:
        data_frame = set_value_labels(
            data_frame,
            metadata,
            formats_as_category=formats_as_category,
            formats_as_ordered_category=formats_as_ordered_category,
        )

    return data_frame, metadata


@overload
def read_sav(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    apply_value_formats: bool = ...,
    formats_as_category: bool = ...,
    formats_as_ordered_category: bool = ...,
    encoding: str | None = ...,
    usecols: list[str] | None = ...,
    user_missing: bool = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["pandas"] | None = ...,
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> "tuple[PandasDataFrame, metadata_container]": ...
@overload
def read_sav(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    apply_value_formats: bool = ...,
    formats_as_category: bool = ...,
    formats_as_ordered_category: bool = ...,
    encoding: str | None = ...,
    usecols: list[str] | None = ...,
    user_missing: bool = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["polars"] = "polars",
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> "tuple[PolarsDataFrame, metadata_container]": ...
@overload
def read_sav(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    apply_value_formats: bool = ...,
    formats_as_category: bool = ...,
    formats_as_ordered_category: bool = ...,
    encoding: str | None = ...,
    usecols: list[str] | None = ...,
    user_missing: bool = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["dict"] = "dict",
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> tuple[DictOutput, metadata_container]: ...
def read_sav(
    filename_path: FilePathorBuffer,
    metadataonly: bool = False,
    dates_as_pandas_datetime: bool = False,
    apply_value_formats: bool = False,
    formats_as_category: bool = True,
    formats_as_ordered_category: bool = False,
    encoding: str | None = None,
    usecols: list[str] | None = None,
    user_missing: bool = False,
    disable_datetime_conversion: bool = False,
    row_limit: int = 0,
    row_offset: int = 0,
    output_format: Literal["pandas", "polars", "dict"] | None = None,
    extra_datetime_formats: list[str] | None = None,
    extra_date_formats: list[str] | None = None,
    extra_time_formats: list[str] | None = None,
) -> "tuple[DataFrame | DictOutput, metadata_container]":
    r"""
    Read a SPSS sav or zsav (compressed) files

    Parameters
    ----------
        filename_path : str, bytes, Path-like object or file-like object
            path to the file or file-like object. In python 2.7 the string is assumed to be utf-8 encoded.
        metadataonly : bool, optional
            by default False. IF true, no data will be read but only metadata, so that you can get all elements in the
            metadata object. The data frame will be set with the correct column names but no data.
        dates_as_pandas_datetime : bool, optional
            by default False. If true dates will be transformed to pandas datetime64 instead of date, effective only for pandas.
        apply_value_formats : bool, optional
            by default False. If true it will change values in the dataframe for they value labels in the metadata,
            if any appropiate are found.
        formats_as_category : bool, optional
            by default True. Takes effect only if apply_value_formats is True. If True, variables with values changed
            for their formatted version will be transformed into categories.
        formats_as_ordered_category : bool, optional
            defaults to False. If True the variables having formats will be transformed into ordered categories/enum.
            it has precedence over formats_as_category, meaning if this is True, it will take effect irrespective of
            the value of formats_as_category.
        encoding : str, optional
            Defaults to None. If set, the system will use the defined encoding instead of guessing it. It has to be an
            iconv-compatible name
        usecols : list, optional
            a list with column names to read from the file. Only those columns will be imported. Case sensitive!
        user_missing : bool, optional
            by default False, in this case user defined missing values are delivered as nan. If true, the missing values
            will be deliver as is, and an extra piece of information will be set in the metadata (missing_ranges)
            to be able to interpret those values as missing.
        disable_datetime_conversion : bool, optional
            if True pyreadstat will not attempt to convert dates, datetimes and times to python objects but those columns
            will remain as numbers. In order to convert them later to an appropiate python object, the user can use the
            information about the original variable format stored in the metadata object in original_variable_types.
            Disabling datetime conversion speeds up reading files. In addition it helps to overcome situations where
            there are datetimes that are beyond the limits of python datetime (which is limited to year 10,000, dates
            beyond that will rise an Overflow error in pyreadstat).
        row_limit : int, optional
            maximum number of rows to read. The default is 0 meaning unlimited.
        row_offset : int, optional
            start reading rows after this offset. By default 0, meaning start with the first row not skipping anything.
        output_format : str, optional
            one of 'pandas' (default), 'polars' or 'dict'. If 'dict' a dictionary with numpy arrays as values will be returned, the
            user can then convert it to her preferred data format. Using dict is faster as the other types as the conversion to a
            dataframe is avoided.
        extra_datetime_formats: list of str, optional
            formats to be parsed as python datetime objects
        extra_date_formats: list of str, optional
            formats to be parsed as python date objects
        extra_time_formats: list of str, optional
            formats to be parsed as python time objects

    Returns
    -------
        data_frame : dataframe or dict
            a dataframe or dict with the data.
        metadata :
            object with metadata. Look at the documentation for more information.
    """
    parser_format = "sav/zsav"

    data_frame, metadata = parser_entry_point(
        filename_path,
        parser_format=parser_format,
        metadataonly=metadataonly,
        dates_as_pandas_datetime=dates_as_pandas_datetime,
        formats_as_category=formats_as_category,
        formats_as_ordered_category=formats_as_ordered_category,
        encoding=encoding,
        usecols=usecols,
        user_missing=user_missing,
        disable_datetime_conversion=disable_datetime_conversion,
        row_limit=row_limit,
        row_offset=row_offset,
        output_format=output_format,
        extra_datetime_formats=extra_datetime_formats,
        extra_date_formats=extra_date_formats,
        extra_time_formats=extra_time_formats,
    )

    metadata.file_format = parser_format

    if apply_value_formats:
        data_frame = set_value_labels(
            data_frame,
            metadata,
            formats_as_category=formats_as_category,
            formats_as_ordered_category=formats_as_ordered_category,
        )

    return data_frame, metadata


@overload
def read_por(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    apply_value_formats: bool = ...,
    formats_as_category: bool = ...,
    formats_as_ordered_category: bool = ...,
    usecols: list[str] | None = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["pandas"] | None = ...,
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> "tuple[PandasDataFrame, metadata_container]": ...
@overload
def read_por(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    apply_value_formats: bool = ...,
    formats_as_category: bool = ...,
    formats_as_ordered_category: bool = ...,
    usecols: list[str] | None = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["polars"] = "polars",
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> "tuple[PolarsDataFrame, metadata_container]": ...
@overload
def read_por(
    filename_path: FilePathorBuffer,
    metadataonly: bool = ...,
    dates_as_pandas_datetime: bool = ...,
    apply_value_formats: bool = ...,
    formats_as_category: bool = ...,
    formats_as_ordered_category: bool = ...,
    usecols: list[str] | None = ...,
    disable_datetime_conversion: bool = ...,
    row_limit: int = ...,
    row_offset: int = ...,
    output_format: Literal["dict"] = "dict",
    extra_datetime_formats: list[str] | None = ...,
    extra_date_formats: list[str] | None = ...,
    extra_time_formats: list[str] | None = ...,
) -> tuple[DictOutput, metadata_container]: ...
def read_por(
    filename_path: FilePathorBuffer,
    metadataonly: bool = False,
    dates_as_pandas_datetime: bool = False,
    apply_value_formats: bool = False,
    formats_as_category: bool = True,
    formats_as_ordered_category: bool = False,
    usecols: list[str] | None = None,
    disable_datetime_conversion: bool = False,
    row_limit: int = 0,
    row_offset: int = 0,
    output_format: Literal["pandas", "polars", "dict"] | None = None,
    extra_datetime_formats: list[str] | None = None,
    extra_date_formats: list[str] | None = None,
    extra_time_formats: list[str] | None = None,
) -> "tuple[DataFrame | DictOutput, metadata_container]":
    r"""
    Read a SPSS por file. Files are assumed to be UTF-8 encoded, the encoding cannot be set to other.

    Parameters
    ----------
        filename_path : str, bytes, Path-like object or file-like object
            path to the file or file-like object. In python 2.7 the string is assumed to be utf-8 encoded.
        metadataonly : bool, optional
            by default False. IF true, no data will be read but only metadata, so that you can get all elements in the
            metadata object. The data frame will be set with the correct column names but no data.
            Notice that number_rows will be None as por files do not have the number of rows recorded in the file metadata.
        dates_as_pandas_datetime : bool, optional
            by default False. If true dates will be transformed to pandas datetime64 instead of date, effective only for pandas.
        apply_value_formats : bool, optional
            by default False. If true it will change values in the dataframe for they value labels in the metadata,
            if any appropiate are found.
        formats_as_category : bool, optional
            by default True. Takes effect only if apply_value_formats is True. If True, variables with values changed
            for their formatted version will be transformed into categories.
        formats_as_ordered_category : bool, optional
            defaults to False. If True the variables having formats will be transformed into ordered categories/enum.
            it has precedence over formats_as_category, meaning if this is True, it will take effect irrespective of
            the value of formats_as_category.
        usecols : list, optional
            a list with column names to read from the file. Only those columns will be imported. Case sensitive!
        disable_datetime_conversion : bool, optional
            if True pyreadstat will not attempt to convert dates, datetimes and times to python objects but those columns
            will remain as numbers. In order to convert them later to an appropiate python object, the user can use the
            information about the original variable format stored in the metadata object in original_variable_types.
            Disabling datetime conversion speeds up reading files. In addition it helps to overcome situations where
            there are datetimes that are beyond the limits of python datetime (which is limited to year 10,000, dates
            beyond that will rise an Overflow error in pyreadstat).
        row_limit : int, optional
            maximum number of rows to read. The default is 0 meaning unlimited.
        row_offset : int, optional
            start reading rows after this offset. By default 0, meaning start with the first row not skipping anything.
        output_format : str, optional
            one of 'pandas' (default), 'polars' or 'dict'. If 'dict' a dictionary with numpy arrays as values will be returned, the
            user can then convert it to her preferred data format. Using dict is faster as the other types as the conversion to a
            dataframe is avoided.
        extra_datetime_formats: list of str, optional
            formats to be parsed as python datetime objects
        extra_date_formats: list of str, optional
            formats to be parsed as python date objects
        extra_time_formats: list of str, optional
            formats to be parsed as python time objects

    Returns
    -------
        data_frame : dataframe or dict
            a dataframe or dict with the data.
        metadata :
            object with metadata. Look at the documentation for more information.
    """
    parser_format = "por"
    data_frame, metadata = parser_entry_point(
        filename_path,
        parser_format=parser_format,
        metadataonly=metadataonly,
        dates_as_pandas_datetime=dates_as_pandas_datetime,
        formats_as_category=formats_as_category,
        formats_as_ordered_category=formats_as_ordered_category,
        usecols=usecols,
        disable_datetime_conversion=disable_datetime_conversion,
        row_limit=row_limit,
        row_offset=row_offset,
        output_format=output_format,
        extra_datetime_formats=extra_datetime_formats,
        extra_date_formats=extra_date_formats,
        extra_time_formats=extra_time_formats,
    )

    metadata.file_format = parser_format

    if apply_value_formats:
        data_frame = set_value_labels(data_frame, metadata, formats_as_category=formats_as_category)

    return data_frame, metadata


@overload
def read_sas7bcat(
    filename_path: FilePathorBuffer,
    encoding: str | None = ...,
    output_format: Literal["pandas"] | None = ...,
) -> "tuple[PandasDataFrame, metadata_container]": ...
@overload
def read_sas7bcat(
    filename_path: FilePathorBuffer,
    encoding: str | None = ...,
    output_format: Literal["polars"] = "polars",
) -> "tuple[PolarsDataFrame, metadata_container]": ...
@overload
def read_sas7bcat(
    filename_path: FilePathorBuffer,
    encoding: str | None = ...,
    output_format: Literal["dict"] = "dict",
) -> tuple[DictOutput, metadata_container]: ...
def read_sas7bcat(
    filename_path: FilePathorBuffer,
    encoding: str | None = None,
    output_format: Literal["pandas", "polars", "dict"] | None = None,
) -> "tuple[DataFrame | DictOutput, metadata_container]":
    r"""
    Read a SAS sas7bcat file. The returning dataframe will be empty. The metadata object will contain a dictionary
    value_labels that contains the formats. When parsing the sas7bdat file, in the metadata, the dictionary
    variable_to_label contains a map from variable name to the formats.
    In order to apply the catalog to the sas7bdat file use set_catalog_to_sas or pass the catalog file as an argument
    to read_sas7bdat directly.
    SAS catalog files are difficult ones, some of them can be read only in specific SAS version, may contain strange
    encodings etc. Therefore it may be that many catalog files are not readable from this application.

    Parameters
    ----------
        filename_path : str, bytes, Path-like object or file-like object
            path to the file or file-like object. In python 2.7 the string is assumed to be utf-8 encoded.
        encoding : str, optional
            Defaults to None. If set, the system will use the defined encoding instead of guessing it. It has to be an
            iconv-compatible name
        output_format : str, optional
            one of 'pandas' (default), 'polars' or 'dict'. If 'dict' a dictionary with numpy arrays as values will be returned.
            Notice that for this function the resulting object is always empty, this is done for consistency with other functions
            but has no impact on performance.

    Returns
    -------
        data_frame : dataframe or dict
            a dataframe with the data (no data in this case, so will be always empty).
        metadata :
            object with metadata. The member value_labels is the one that contains the formats.
            Look at the documentation for more information.
    """
    parser_format = "sas7bcat"
    data_frame, metadata = parser_entry_point(
        filename_path,
        parser_format=parser_format,
        encoding=encoding,
        output_format=output_format,
    )

    metadata.file_format = parser_format

    return data_frame, metadata


# convenience functions to read in chunks


@overload
def read_file_in_chunks(
    read_function: PyreadstatReadFunction,
    file_path: FilePathLike,
    chunksize: int = ...,
    offset: int = ...,
    limit: int = ...,
    multiprocess: bool = ...,
    num_processes: int = ...,
    num_rows: int | None = ...,
    *,
    output_format: Literal["pandas"] | None = ...,
    **kwargs: Any,
) -> "Iterator[tuple[PandasDataFrame, metadata_container]]": ...
@overload
def read_file_in_chunks(
    read_function: PyreadstatReadFunction,
    file_path: FilePathLike,
    chunksize: int = ...,
    offset: int = ...,
    limit: int = ...,
    multiprocess: bool = ...,
    num_processes: int = ...,
    num_rows: int | None = ...,
    *,
    output_format: Literal["polars"] = "polars",
    **kwargs: Any,
) -> "Iterator[tuple[PolarsDataFrame, metadata_container]]": ...
@overload
def read_file_in_chunks(
    read_function: PyreadstatReadFunction,
    file_path: FilePathLike,
    chunksize: int = ...,
    offset: int = ...,
    limit: int = ...,
    multiprocess: bool = ...,
    num_processes: int = ...,
    num_rows: int | None = ...,
    *,
    output_format: Literal["dict"] = "dict",
    **kwargs: Any,
) -> Iterator[tuple[DictOutput, metadata_container]]: ...
def read_file_in_chunks(
    read_function: PyreadstatReadFunction,
    file_path: FilePathLike,
    chunksize: int = 100000,
    offset: int = 0,
    limit: int = 0,
    multiprocess: bool = False,
    num_processes: int = 4,
    num_rows: int | None = None,
    **kwargs: Any,
) -> "Iterator[tuple[DataFrame | DictOutput, metadata_container]]":
    """
    Returns a generator that will allow to read a file in chunks.

    If using multiprocessing, for Xport, Por and some defective sav files where the number of rows in the dataset canot be obtained from the metadata,
    the parameter num_rows must be set to a number equal or larger than the number of rows in the dataset. That information must
    be obtained by the user before running this function.

    Parameters
    ----------
        read_function : pyreadstat function
            a pyreadstat reading function
        file_path : str, bytes or Path-like object
            path to the file to be read
        chunksize : integer, optional
            size of the chunks to read
        offset : integer, optional
            start reading the file after certain number of rows
        limit : integer, optional
            stop reading the file after certain number of rows, will be added to offset
        multiprocess: bool, optional
            use multiprocessing to read each chunk?
        num_processes: integer, optional
            in case multiprocess is true, how many workers/processes to spawn?
        num_rows: integer, optional
            number of rows in the dataset. If using multiprocessing it is obligatory for files where
            the number of rows cannot be obtained from the medatata, such as por and
            some defective xport and sav files. The user must obtain this value by reading the file without multiprocessing first or any other means. A number
            larger than the actual number of rows will work as well. Discarded if the number of rows can be obtained from the metadata or not using
            multiprocessing.
        kwargs : dict, optional
            any other keyword argument to pass to the read_function. row_limit and row_offset will be discarded if present.

    Yields
    -------
        data_frame : dataframe
            a dataframe with the data
        metadata :
            object with metadata.
            Look at the documentation for more information.

        it : generator
            A generator that reads the file in chunks.
    """

    if read_function == read_sas7bcat:
        raise Exception("read_sas7bcat not supported")

    if "row_offset" in kwargs:
        _ = kwargs.pop("row_offset")

    if "row_limit" in kwargs:
        _ = kwargs.pop("row_limit")

    if "num_processes" in kwargs:
        _ = kwargs.pop("num_processes")

    _, meta = read_function(file_path, metadataonly=True)
    numrows = meta.number_rows
    if numrows:
        if not limit:
            limit = numrows
        else:
            limit = min(offset + limit, numrows)
    else:
        if limit:
            limit = offset + limit
    df = [0]
    while len(df):
        if limit and (offset >= limit):
            break
        if multiprocess:
            df, meta = read_file_multiprocessing(
                read_function,
                file_path,
                num_processes=num_processes,
                row_offset=offset,
                row_limit=chunksize,
                num_rows=num_rows,
                **kwargs,
            )
        else:
            df, meta = read_function(file_path, row_offset=offset, row_limit=chunksize, **kwargs)
        if len(df):
            yield df, meta
            offset += chunksize


@overload
def read_file_multiprocessing(
    read_function: PyreadstatReadFunction,
    file_path: FilePathLike,
    num_processes: int | None = ...,
    num_rows: int | None = ...,
    *,
    output_format: Literal["pandas"] | None = ...,
    **kwargs: Any,
) -> "tuple[PandasDataFrame, metadata_container]": ...
@overload
def read_file_multiprocessing(
    read_function: PyreadstatReadFunction,
    file_path: FilePathLike,
    num_processes: int | None = ...,
    num_rows: int | None = ...,
    *,
    output_format: Literal["polars"] = "polars",
    **kwargs: Any,
) -> "tuple[PolarsDataFrame, metadata_container]": ...
@overload
def read_file_multiprocessing(
    read_function: PyreadstatReadFunction,
    file_path: FilePathLike,
    num_processes: int | None = ...,
    num_rows: int | None = ...,
    *,
    output_format: Literal["dict"] = "dict",
    **kwargs: Any,
) -> tuple[DictOutput, metadata_container]: ...
def read_file_multiprocessing(
    read_function: PyreadstatReadFunction,
    file_path: FilePathLike,
    num_processes: int | None = None,
    num_rows: int | None = None,
    **kwargs: Any,
) -> "tuple[DataFrame | DictOutput, metadata_container]":
    """
    Reads a file in parallel using multiprocessing.
    For Xport, Por and some defective sav files where the number of rows in the dataset canot be obtained from the metadata,
    the parameter num_rows must be set to a number equal or larger than the number of rows in the dataset. That information must
    be obtained by the user before running this function.

    Parameters
    ----------
        read_function : pyreadstat function
            a pyreadstat reading function
        file_path : str, bytes or Path-like object
            path to the file to be read
        num_processes : integer, optional
            number of processes to spawn, by default the min 4 and the max cores on the computer
        num_rows: integer, optional
            number of rows in the dataset. Obligatory for files where the number of rows cannot be obtained from the medatata, such as por and
            some defective xport and sav files. The user must obtain this value by reading the file without multiprocessing first or any other means. A number
            larger than the actual number of rows will work as well. Discarded if the number of rows can be obtained from the metadata.
        kwargs : dict, optional
            any other keyword argument to pass to the read_function.

    Returns
    -------
        data_frame : dataframe
            a dataframe with the data
        metadata :
            object with metadata. Look at the documentation for more information.
    """

    if read_function in (read_sas7bcat,):
        raise Exception("read_sas7bcat is not supported")

    if read_function == read_por and num_rows is None:
        raise Exception(
            "num_rows must be specified for read_por to be a number equal or larger than the number of rows in the dataset."
        )

    if not num_processes:
        # let's be conservative with the number of workers
        num_processes = min(mp.cpu_count(), 4)
    _ = kwargs.pop("metadataonly", None)
    row_offset = kwargs.pop("row_offset", 0)
    row_limit = kwargs.pop("row_limit", float("inf"))
    _, meta = read_function(file_path, metadataonly=True, **kwargs)
    numrows = meta.number_rows

    if numrows is None:
        if num_rows is None:
            raise Exception(
                "The number of rows of the file cannot be determined from the file's metadata. If you still want to proceed, please set num_rows to a number equal or larger than the number of rows of your data"
            )
        numrows = num_rows
    elif numrows == 0:
        final, meta = read_function(file_path, **kwargs)

    numrows = min(max(numrows - row_offset, 0), row_limit)
    divs = [numrows // num_processes + (1 if x < numrows % num_processes else 0) for x in range(num_processes)]
    offsets = list()
    prev_offset = row_offset
    prev_div = 0
    for indx, div in enumerate(divs):
        offset = prev_offset + prev_div
        prev_offset = offset
        prev_div = div
        offsets.append((offset, div))
    jobs = [(read_function, file_path, offset, chunksize, kwargs) for offset, chunksize in offsets]
    pool = mp.Pool(processes=num_processes)
    try:
        chunks = pool.map(worker, jobs)
    except:
        raise
    finally:
        pool.close()
    output_format = kwargs.get("output_format")
    if output_format == "dict":
        keys = chunks[0].keys()
        final = {key: list(chain.from_iterable(chunk[key] for chunk in chunks)) for key in keys}
    else:
        # final = pd.concat(chunks, axis=0, ignore_index=True)
        chunks = [nw.from_native(x) for x in chunks]
        final = nw.concat(chunks, how="vertical")
        ispandas = False
        if final.implementation.is_pandas():
            ispandas = True
        final = final.to_native()
        if ispandas:
            final = final.reset_index(drop=True)
    return final, meta


# Write API


def write_sav(
    df: "DataFrame",
    dst_path: FilePathLike,
    file_label: str = "",
    column_labels: list[str] | dict[str, str] | None = None,
    compress: bool = False,
    row_compress: bool = False,
    note: str | list[str] | None = None,
    variable_value_labels: dict[str, dict[int | float, str]] | None = None,
    missing_ranges: dict[str, list[int | float | str | MissingRange]] | None = None,
    variable_display_width: dict[str, int] | None = None,
    variable_measure: dict[str, str] | None = None,
    variable_format: dict[str, str] | None = None,
) -> None:
    """
    Writes a dataframe to a SPSS sav or zsav file.

    Parameters
    ----------
    df : dataframe
        dataframe to write to sav or zsav
    dst_path : str, bytes or Path-like object
        full path to the result sav or zsav file
    file_label : str, optional
        a label for the file
    column_labels : list or dict, optional
        labels for columns (variables), if list must be the same length as the number of columns. Variables with no
        labels must be represented by None. If dict values must be variable names and values variable labels.
        In such case there is no need to include all variables; labels for non existent
        variables will be ignored with no warning or error.
    compress : boolean, optional
        if true a zsav will be written, by default False, a sav is written
    row_compress : boolean, optional
        if true it applies row compression, by default False, compress and row_compress cannot be both true at the same time
    note : str or list of str, optional
        a note or list of notes to add to the file
    variable_value_labels : dict, optional
        value labels, a dictionary with key variable name and value a dictionary with key values and
        values labels. Variable names must match variable names in the dataframe otherwise will be
        ignored. Value types must match the type of the column in the dataframe.
    missing_ranges : dict, optional
        user defined missing values. Must be a dictionary with keys as variable names matching variable
        names in the dataframe. The values must be a list. Each element in that list can either be
        either a discrete numeric or string value (max 3 per variable) or a dictionary with keys 'hi' and 'lo' to
        indicate the upper and lower range for numeric values (max 1 range value + 1 discrete value per
        variable). hi and lo may also be the same value in which case it will be interpreted as a discrete
        missing value.
        For this to be effective, values in the dataframe must be the same as reported here and not NaN.
    variable_display_width : dict, optional
        set the display width for variables. Must be a dictonary with keys being variable names and
        values being integers.
    variable_measure: dict, optional
        sets the measure type for a variable. Must be a dictionary with keys being variable names and
        values being strings one of "nominal", "ordinal", "scale" or "unknown" (default).
    variable_format: dict, optional
        sets the format of a variable. Must be a dictionary with keys being the variable names and
        values being strings defining the format. See README, setting variable formats section,
        for more information.
    """
    writer_format = "sav"

    # formats
    formats_presets = {"restricted_integer": "N{var_width}", "integer": "F{var_width}.0"}
    if variable_format:
        for col_name, col_format in variable_format.items():
            if col_format in formats_presets.keys() and col_name in df.columns:
                var_width = str(len(str(max(df[col_name]))))
                variable_format[col_name] = formats_presets[col_format].format(var_width=var_width)

    writer_entry_point(
        df,
        dst_path,
        writer_format=writer_format,
        file_label=file_label,
        column_labels=column_labels,
        compress=compress,
        row_compress=row_compress,
        note=note,
        variable_value_labels=variable_value_labels,
        missing_ranges=missing_ranges,
        variable_display_width=variable_display_width,
        variable_measure=variable_measure,
        variable_format=variable_format,
    )


def write_dta(
    df: "DataFrame",
    dst_path: FilePathLike,
    file_label: str = "",
    column_labels: list[str] | dict[str, str] | None = None,
    version: int = 15,
    variable_value_labels: dict[str, dict[int | float, str]] | None = None,
    missing_user_values: dict[str, list[str]] | None = None,
    variable_format: dict[str, str] | None = None,
) -> None:
    """
    Writes a dataframe to a STATA dta file

    Parameters
    ----------
    df : dataframe
        dataframe to write to sav or zsav
    dst_path : str, bytes or Path-like object
        full path to the result dta file
    file_label : str, optional
        a label for the file
    column_labels : list or dict, optional
        labels for columns (variables), if list must be the same length as the number of columns. Variables with no
        labels must be represented by None. If dict values must be variable names and values variable labels.
        In such case there is no need to include all variables; labels for non existent
        variables will be ignored with no warning or error.
    version : int, optional
        dta file version, supported from 8 to 15, default is 15
    variable_value_labels : dict, optional
        value labels, a dictionary with key variable name and value a dictionary with key values and
        values labels. Variable names must match variable names in the dataframe otherwise will be
        ignored. Value types must match the type of the column in the dataframe.
    missing_user_values : dict, optional
        user defined missing values for numeric variables. Must be a dictionary with keys being variable
        names and values being a list of missing values. Missing values must be a single character
        between a and z.
    variable_format: dict, optional
        sets the format of a variable. Must be a dictionary with keys being the variable names and
        values being strings defining the format. See README, setting variable formats section,
        for more information.
    """

    writer_format = "dta"
    writer_entry_point(
        df,
        dst_path,
        writer_format=writer_format,
        file_label=file_label,
        column_labels=column_labels,
        version=version,
        variable_value_labels=variable_value_labels,
        missing_user_values=missing_user_values,
        variable_format=variable_format,
    )


def write_xport(
    df: "DataFrame",
    dst_path: FilePathLike,
    file_label: str = "",
    column_labels: list[str] | dict[str, str] | None = None,
    table_name: str | None = None,
    file_format_version: Literal[5, 8] = 8,
    variable_format: dict[str, str] | None = None,
) -> None:
    """
    Writes a dataframe to a SAS Xport (xpt) file.
    If no table_name is specified the dataset has by default the name DATASET (take it into account if
    reading the file from SAS.)
    Versions 5 and 8 are supported, default is 8.

    Parameters
    ----------
    df : dataframe
        dataframe to write to xport
    dst_path : str, bytes or Path-like object
        full path to the result xport file
    file_label : str, optional
        a label for the file
    column_labels : list or dict, optional
        labels for columns (variables), if list must be the same length as the number of columns. Variables with no
        labels must be represented by None. If dict values must be variable names and values variable labels.
        In such case there is no need to include all variables; labels for non existent
        variables will be ignored with no warning or error.
    table_name : str, optional
        name of the dataset, by default DATASET
    file_format_version : int, optional
        XPORT file version, either 8 or 5, default is 8
    variable_format: dict, optional
        sets the format of a variable. Must be a dictionary with keys being the variable names and
        values being strings defining the format. See README, setting variable formats section,
        for more information.
    """

    writer_format = "xport"
    writer_entry_point(
        df,
        dst_path,
        writer_format=writer_format,
        file_label=file_label,
        column_labels=column_labels,
        version=file_format_version,
        table_name=table_name,
        variable_format=variable_format,
    )


def write_por(
    df: "DataFrame",
    dst_path: FilePathLike,
    file_label: str = "",
    column_labels: list[str] | dict[str, str] | None = None,
    variable_format: dict[str, str] | None = None,
) -> None:
    """
    Writes a dataframe to a SPSS POR file.

    Parameters
    ----------
    df : dataframe
        data frame to write to por
    dst_path : str, bytes or Path-like object
        full path to the result por file
    file_label : str, optional
        a label for the file
    column_labels : list or dict, optional
        labels for columns (variables), if list must be the same length as the number of columns. Variables with no
        labels must be represented by None. If dict values must be variable names and values variable labels.
        In such case there is no need to include all variables; labels for non existent
        variables will be ignored with no warning or error.
    variable_format: dict, optional
        sets the format of a variable. Must be a dictionary with keys being the variable names and
        values being strings defining the format. See README, setting variable formats section,
        for more information.
    """

    writer_format = "por"
    writer_entry_point(
        df,
        dst_path,
        writer_format=writer_format,
        file_label=file_label,
        column_labels=column_labels,
        variable_format=variable_format,
    )
