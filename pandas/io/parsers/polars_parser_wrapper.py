"""
Polars parser wrapper for reading CSV files.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from pandas.compat._optional import import_optional_dependency
from pandas.errors import ParserError

from pandas.io.parsers.base_parser import ParserBase

if TYPE_CHECKING:
    from pandas._typing import ReadBuffer

    from pandas import DataFrame


class PolarsParserWrapper(ParserBase):
    """
    Wrapper for the polars engine for read_csv()
    """

    def __init__(self, src: ReadBuffer[bytes] | ReadBuffer[str], **kwds) -> None:
        super().__init__(kwds)
        self.kwds = kwds
        self.src = src

        self._parse_kwds()

    def _parse_kwds(self) -> None:
        """
        Validates keywords before passing to polars.
        """
        encoding: str | None = self.kwds.get("encoding")
        self.encoding = "utf-8" if encoding is None else encoding

        na_values = self.kwds["na_values"]
        if isinstance(na_values, dict):
            raise ValueError(
                "The polars engine doesn't support passing a dict for na_values"
            )
        self.na_values = list(self.kwds["na_values"])

    def _get_polars_options(self) -> dict:
        """
        Map pandas options to polars read_csv options.
        """
        # Import polars
        pl = import_optional_dependency("polars")
        
        polars_options = {}
        
        # Basic options mapping
        if self.kwds.get("sep") is not None:
            polars_options["separator"] = self.kwds["sep"]
        
        if self.kwds.get("header") is not None:
            header = self.kwds["header"]
            if header is None:
                polars_options["has_header"] = False
            elif header == 0:
                polars_options["has_header"] = True
            else:
                # For multi-line headers, skip rows and assume header
                polars_options["has_header"] = True
                polars_options["skip_rows"] = header
        
        if self.kwds.get("skiprows") is not None:
            skiprows = self.kwds["skiprows"]
            if isinstance(skiprows, int):
                polars_options["skip_rows"] = skiprows
        
        if self.kwds.get("na_values") is not None:
            na_vals = self.kwds["na_values"]
            if isinstance(na_vals, str):
                polars_options["null_values"] = [na_vals]
            elif hasattr(na_vals, '__iter__'):
                polars_options["null_values"] = list(na_vals)
        
        if self.kwds.get("quotechar") is not None:
            polars_options["quote_char"] = self.kwds["quotechar"]
        
        if self.kwds.get("comment") is not None:
            polars_options["comment_prefix"] = self.kwds["comment"]
        
        if self.kwds.get("encoding") is not None:
            polars_options["encoding"] = self.kwds["encoding"]
        
        # Handle usecols - only column names are supported 
        if self.kwds.get("usecols") is not None:
            usecols = self.kwds["usecols"]
            if callable(usecols):
                raise ValueError(
                    "The polars engine does not support callable usecols"
                )
            polars_options["columns"] = usecols

        # Handle nrows
        if self.kwds.get("nrows") is not None:
            polars_options["n_rows"] = self.kwds["nrows"]

        # Handle dtype mapping
        if self.kwds.get("dtype") is not None:
            dtype = self.kwds["dtype"]
            if isinstance(dtype, dict):
                # Convert pandas dtypes to polars dtypes
                polars_schema = {}
                for col, dt in dtype.items():
                    polars_schema[col] = self._convert_dtype_to_polars(dt)
                polars_options["schema"] = polars_schema
            # Single dtype for all columns will be handled after reading

        return polars_options

    def _convert_dtype_to_polars(self, pandas_dtype_str):
        """
        Convert pandas dtype string to polars dtype.
        """
        pl = import_optional_dependency("polars")
        
        # Map common pandas dtypes to polars dtypes
        dtype_mapping = {
            "object": pl.Utf8,
            "str": pl.Utf8,
            "string": pl.Utf8,
            "int64": pl.Int64,
            "int32": pl.Int32,
            "int16": pl.Int16,
            "int8": pl.Int8,
            "uint64": pl.UInt64,
            "uint32": pl.UInt32,
            "uint16": pl.UInt16,
            "uint8": pl.UInt8,
            "float64": pl.Float64,
            "float32": pl.Float32,
            "bool": pl.Boolean,
            "datetime64[ns]": pl.Datetime("ns"),
            "category": pl.Categorical,
        }
        
        # Handle string representation
        if isinstance(pandas_dtype_str, str):
            return dtype_mapping.get(pandas_dtype_str, pl.Utf8)
        else:
            # For actual dtype objects, convert to string first
            dtype_str = str(pandas_dtype_str)
            return dtype_mapping.get(dtype_str, pl.Utf8)

    def _adjust_column_names(self, df) -> bool:
        """
        Adjust column names if needed.
        """
        multi_index_named = True
        
        # Handle custom column names
        if self.names is not None:
            if len(self.names) != len(df.columns):
                raise ValueError(
                    f"Number of names ({len(self.names)}) does not match "
                    f"number of columns ({len(df.columns)})"
                )
            df = df.select([
                df[old_name].alias(new_name) 
                for old_name, new_name in zip(df.columns, self.names)
            ])
            
        return multi_index_named, df

    def _finalize_index(self, frame: DataFrame, multi_index_named: bool) -> DataFrame:
        """
        Set up the index if index_col is specified.
        """
        if self.index_col is not None:
            if isinstance(self.index_col, list):
                # MultiIndex case
                frame.set_index(self.index_col, drop=True, inplace=True)
            else:
                # Single index
                frame.set_index(self.index_col, drop=True, inplace=True)
                
            # Clear names if headerless and no name given
            if self.header is None and not multi_index_named:
                frame.index.names = [None] * len(frame.index.names)

        return frame

    def _finalize_dtype(self, frame: DataFrame) -> DataFrame:
        """
        Apply any remaining dtype conversions.
        """
        if self.dtype is not None and not isinstance(self.dtype, dict):
            # Single dtype for all columns
            try:
                for col in frame.columns:
                    if col not in (self.index_col or []):
                        frame[col] = frame[col].astype(self.dtype)
            except (TypeError, ValueError) as err:
                raise ValueError(f"Error converting dtypes: {err}") from err
                
        return frame

    def _apply_filtering(self, lazy_frame):
        """
        Apply column selection and row filtering using lazy operations.
        """
        # Column selection (usecols equivalent)
        if self.kwds.get("usecols") is not None:
            usecols = self.kwds["usecols"]
            if not callable(usecols):
                try:
                    lazy_frame = lazy_frame.select(usecols)
                except Exception as e:
                    # Fallback to pandas-style selection after collection
                    pass
        
        # Row filtering could be added here for predicate pushdown
        # For now, we'll handle skiprows and nrows in the scan_csv call
        
        return lazy_frame

    def read(self) -> DataFrame:
        """
        Reads the contents of a CSV file into a DataFrame using Polars
        and converts it to pandas.

        Returns
        -------
        DataFrame
            The DataFrame created from the CSV file.
        """
        pl = import_optional_dependency("polars")
        
        try:
            # Get polars options
            polars_options = self._get_polars_options()
            
            # For file-like objects, read content and use read_csv
            if hasattr(self.src, 'read'):
                # For file-like objects, we need to get the content
                content = self.src.read()
                if isinstance(content, bytes):
                    content = content.decode(self.encoding)
                
                # Use read_csv with string content
                from io import StringIO
                polars_df = pl.read_csv(StringIO(content), **polars_options)
            else:
                # For file paths, we can use scan_csv for lazy evaluation
                if isinstance(self.src, str):
                    # Use lazy reading for better performance
                    lazy_df = pl.scan_csv(self.src, **polars_options)
                    lazy_df = self._apply_filtering(lazy_df)
                    polars_df = lazy_df.collect()
                else:
                    # Fallback to read_csv for other cases
                    polars_df = pl.read_csv(self.src, **polars_options)
            
            # Convert to pandas DataFrame
            frame = polars_df.to_pandas()
            
        except Exception as e:
            if "polars" in str(e).lower() or "pl." in str(e):
                raise ParserError(f"Polars parsing error: {e}") from e
            else:
                raise ParserError(f"Error reading CSV with polars engine: {e}") from e

        # Adjust column names if needed
        multi_index_named, frame = self._adjust_column_names_pandas(frame)
        
        # Apply date conversions
        frame = self._do_date_conversions(frame.columns, frame)
        
        # Set up index
        frame = self._finalize_index(frame, multi_index_named)
        
        # Apply remaining dtype conversions
        frame = self._finalize_dtype(frame)
        
        return frame

    def _adjust_column_names_pandas(self, frame: DataFrame) -> tuple[bool, DataFrame]:
        """
        Adjust column names for pandas DataFrame after conversion from Polars.
        """
        multi_index_named = True
        
        # Handle custom column names
        if self.names is not None:
            if len(self.names) != len(frame.columns):
                raise ValueError(
                    f"Number of names ({len(self.names)}) does not match "
                    f"number of columns ({len(frame.columns)})"
                )
            frame.columns = self.names
            
        return multi_index_named, frame

    def close(self) -> None:
        """
        Close any open resources.
        """
        # Polars doesn't require explicit cleanup for most cases
        pass
