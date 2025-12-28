from __future__ import annotations

import functools
import warnings
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Literal,
)

import numpy as np

from pandas._config import using_copy_on_write

from pandas.util._decorators import doc

from pandas.core.dtypes.common import (
    is_datetime64_any_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_object_dtype,
    is_string_dtype,
)
from pandas.core.dtypes.dtypes import CategoricalDtype

import pandas as pd
from pandas.api.extensions import (
    register_dataframe_accessor,
    register_series_accessor,
)

if TYPE_CHECKING:
    from pandas import (
        DataFrame,
        Series,
    )


@register_dataframe_accessor("modern")
class ModernAccessor:
    """
    The 'Modern Pandas' Accessor.
    Access via ``df.modern...`` to use enhanced features.
    """

    def __init__(self, pandas_obj: DataFrame | Series):
        self._obj = pandas_obj
        self._validate_queue = []

    # ------------------------------------------------------------------
    # 1. OPTIMIZE: Memory Reduction
    # ------------------------------------------------------------------
    def optimize(self, verbose: bool = True) -> DataFrame:
        """
        Aggressively optimize memory usage by downcasting types.
        """
        df = self._obj.copy()
        start_mem = df.memory_usage(deep=True).sum()

        for col in df.columns:
            # 1. Float downcasting (float64 -> float32)
            if is_float_dtype(df[col]):
                df[col] = pd.to_numeric(df[col], downcast="float")

            # 2. Integer downcasting (int64 -> int8/16/32)
            elif is_integer_dtype(df[col]):
                df[col] = pd.to_numeric(df[col], downcast="integer")

            # 3. Object/String Optimization
            elif is_object_dtype(df[col]) or is_string_dtype(df[col]):
                num_unique = len(df[col].unique())
                num_total = len(df[col])
                
                # If low cardinality (<50%), convert to Category
                if num_unique / num_total < 0.5:
                    df[col] = df[col].astype("category")
                else:
                    # Otherwise verify it's pyarrow string if available (pandas 2+)
                    try:
                        df[col] = df[col].astype("string[pyarrow]")
                    except (ImportError, TypeError):
                        pass # Fallback to object

        end_mem = df.memory_usage(deep=True).sum()
        if verbose:
            saved = (start_mem - end_mem) / start_mem * 100
            print(f"Memory Optimized: {saved:.2f}% reduction detected.")
            print(f"Start: {start_mem / 1024**2:.2f} MB | End: {end_mem / 1024**2:.2f} MB")
        
        return df

    # ------------------------------------------------------------------
    # 2. CLEAN: Smart Cleaning
    # ------------------------------------------------------------------
    def clean_smart(self) -> DataFrame:
        """
        Heuristic-based data cleaning. 
        - Fixes common NULL strings
        - Strips whitespace
        - Standardizes column names
        """
        df = self._obj.copy()

        # 1. Standardize NULLs
        null_strings = ["na", "n/a", "null", "none", "-", "?", "missing", "nan"]
        df = df.replace(null_strings, np.nan, regex=False)

        # 2. Clean Strings in Object Columns
        for col in df.select_dtypes(include=["object", "string"]):
            # Strip whitespace
            df[col] = df[col].astype(str).str.strip()
            # Convert "nan" string back to real NaN after strip
            mask = df[col].str.lower() == "nan"
            df.loc[mask, col] = np.nan

        # 3. Clean Column Names (snake_case)
        df.columns = (
            df.columns.astype(str)
            .str.lower()
            .str.replace(r"\s+", "_", regex=True)
            .str.replace(r"[^a-z0-9_]", "", regex=True)
        )

        return df

    # ------------------------------------------------------------------
    # 3. EXPECT: Safety Net Validation
    # ------------------------------------------------------------------
    @property
    def expect(self) -> ModernAccessor:
        """Start a validation chain."""
        return self

    def no_nulls(self, subset: list[str] | str | None = None):
        """Queue a check for no nulls."""
        self._validate_queue.append(("no_nulls", subset))
        return self

    def unique(self, subset: list[str] | str):
        """Queue a check for uniqueness."""
        self._validate_queue.append(("unique", subset))
        return self
    
    def validate(self) -> DataFrame:
        """Execute all queued expectations."""
        df = self._obj
        errors = []

        for check, arg in self._validate_queue:
            if check == "no_nulls":
                cols = arg if arg else df.columns
                if df[cols].isnull().any().any():
                    errors.append(f"Validation Failed: Nulls found in {cols}")
            elif check == "unique":
                cols = [arg] if isinstance(arg, str) else arg
                for c in cols:
                    if not df[c].is_unique:
                        errors.append(f"Validation Failed: Column '{c}' is not unique.")
        
        self._validate_queue = [] # Clear queue
        
        if errors:
            raise ValueError("\n".join(errors))
        return df

    # ------------------------------------------------------------------
    # 4. FUZZY: Merge
    # ------------------------------------------------------------------
    def merge_fuzzy(self, right: DataFrame, on: str, threshold: int = 80) -> DataFrame:
        """
        Perform a fuzzy merge on a string column.
        Requires 'rapidfuzz'.
        """
        try:
            from rapidfuzz import process, fuzz
        except ImportError:
            raise ImportError("Please install 'rapidfuzz' to use df.modern.merge_fuzzy()")

        df_left = self._obj.copy()
        df_right = right.copy()

        # Get unique values to match
        left_keys = df_left[on].dropna().unique()
        right_keys = df_right[on].dropna().unique()

        # Find matches
        matches = {}
        for key in left_keys:
            # process.extractOne returns (match, score, index)
            res = process.extractOne(key, right_keys, scorer=fuzz.WRatio)
            if res and res[1] >= threshold:
                matches[key] = res[0]

        # Map back to dataframe
        match_col = f"{on}_matched"
        df_left[match_col] = df_left[on].map(matches)
        
        # Merge
        return df_left.merge(df_right, left_on=match_col, right_on=on, suffixes=("", "_fuzzy"))

    # ------------------------------------------------------------------
    # 5. BRIDGE: Interoperability
    # ------------------------------------------------------------------
    @property
    def bridge(self):
        return BridgeAccessor(self._obj)

    # ------------------------------------------------------------------
    # 6. AI: LLM Interface
    # ------------------------------------------------------------------
    @property
    def ai(self):
        return AIAccessor(self._obj)


class BridgeAccessor:
    def __init__(self, obj):
        self._obj = obj

    def to_polars(self):
        try:
            import polars as pl
            return pl.from_pandas(self._obj)
        except ImportError:
            raise ImportError("Please install 'polars' to use df.modern.bridge.to_polars()")

    def to_duckdb(self, table_name="pandas_df"):
        try:
            import duckdb
            con = duckdb.connect()
            con.register(table_name, self._obj)
            return con
        except ImportError:
            raise ImportError("Please install 'duckdb' to use df.modern.bridge.to_duckdb()")


class AIAccessor:
    def __init__(self, obj):
        self._obj = obj

    def ask(self, query: str, api_key: str | None = None):
        """
        Natural language query. 
        Note: In a real implementation, this would call OpenAI/Gemini APIs.
        For this prototype, we simulate the response generation logic.
        """
        # Simulation Logic
        print(f"🤖 AI Request: '{query}'")
        print("Analyzing Data Structure...")
        print(f"Columns: {list(self._obj.columns)}")
        
        if "plot" in query.lower():
            print("Action: Generating Plot Code...")
            print("--> df.plot(kind='bar') # Simulated")
        elif "max" in query.lower() or "mean" in query.lower():
            print("Action: Aggregation...")
            print("--> df.groupby('...').mean()")
        else:
            print("Action: General Query")
            
        print("Note: Set OpenAI API Key to execute real queries.")
        return "AI_Response_Object"

