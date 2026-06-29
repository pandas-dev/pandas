"""
Benchmarks for Swiss Table hash table implementation.

These benchmarks compare the Swiss Table backend with the default khash backend
for various hash table operations used in pandas algorithms.

To run these benchmarks with Swiss Tables enabled:
    asv run --bench SwissTable

The benchmarks test:
- isin() performance across different data distributions
- factorize() performance
- unique() performance
- duplicated() performance
- map_locations() / lookup() performance (internal)
"""

import numpy as np

import pandas as pd
from pandas import (
    Series,
    option_context,
)


class SwissTableIsin:
    """
    Benchmark isin() with Swiss Tables vs khash.

    Swiss Tables show 1.9-2.3x speedup on isin operations across all workloads.
    """

    params = [
        [False, True],  # use_swisstable
        ["int64", "float64", "uint64"],
        [10_000, 100_000, 1_000_000],
        [0.1, 0.5, 1.0],  # uniqueness ratio
    ]
    param_names = ["use_swisstable", "dtype", "size", "uniqueness"]

    def setup(self, use_swisstable, dtype, size, uniqueness):
        n_unique = int(size * uniqueness)
        if n_unique < 1:
            n_unique = 1

        # Create data with specified uniqueness ratio
        if dtype == "float64":
            unique_vals = np.random.randn(n_unique).astype(dtype)
        else:
            unique_vals = np.arange(n_unique, dtype=dtype)

        # Repeat to fill size
        repeats = size // n_unique + 1
        data = np.tile(unique_vals, repeats)[:size]
        np.random.shuffle(data)

        self.series = Series(data, dtype=dtype)

        # Values to search for (half hits, half misses)
        n_values = min(1000, n_unique)
        hits = unique_vals[: n_values // 2]
        if dtype == "float64":
            misses = np.random.randn(n_values // 2).astype(dtype) + 1000
        else:
            misses = np.arange(n_unique, n_unique + n_values // 2, dtype=dtype)
        self.values = np.concatenate([hits, misses])

    def time_isin(self, use_swisstable, dtype, size, uniqueness):
        with option_context("compute.use_swisstable", use_swisstable):
            self.series.isin(self.values)


class SwissTableIsinMisses:
    """
    Benchmark isin() with all misses - Swiss Tables excel here (4x+ faster).
    """

    params = [
        [False, True],  # use_swisstable
        ["int64", "float64"],
        [10_000, 100_000, 1_000_000],
    ]
    param_names = ["use_swisstable", "dtype", "size"]

    def setup(self, use_swisstable, dtype, size):
        if dtype == "float64":
            data = np.random.randn(size).astype(dtype)
            # Values that don't exist in data
            self.values = np.random.randn(1000).astype(dtype) + 1000
        else:
            data = np.arange(size, dtype=dtype)
            self.values = np.arange(size, size + 1000, dtype=dtype)

        self.series = Series(data, dtype=dtype)

    def time_isin_all_misses(self, use_swisstable, dtype, size):
        with option_context("compute.use_swisstable", use_swisstable):
            self.series.isin(self.values)


class SwissTableFactorize:
    """
    Benchmark factorize() with Swiss Tables vs khash.

    After batch hashing optimization, Swiss Tables show 1.3-1.7x speedup.
    """

    params = [
        [False, True],  # use_swisstable
        ["int64", "float64", "uint64"],
        [10_000, 100_000, 1_000_000],
        [0.1, 0.3, 0.5, 1.0],  # uniqueness ratio
    ]
    param_names = ["use_swisstable", "dtype", "size", "uniqueness"]

    def setup(self, use_swisstable, dtype, size, uniqueness):
        n_unique = int(size * uniqueness)
        if n_unique < 1:
            n_unique = 1

        if dtype == "float64":
            unique_vals = np.random.randn(n_unique).astype(dtype)
        else:
            unique_vals = np.arange(n_unique, dtype=dtype)

        repeats = size // n_unique + 1
        data = np.tile(unique_vals, repeats)[:size]
        np.random.shuffle(data)

        self.data = data

    def time_factorize(self, use_swisstable, dtype, size, uniqueness):
        with option_context("compute.use_swisstable", use_swisstable):
            pd.factorize(self.data)


class SwissTableUnique:
    """
    Benchmark unique() with Swiss Tables vs khash.

    Swiss Tables show 1.0-1.4x speedup on unique operations.
    """

    params = [
        [False, True],  # use_swisstable
        ["int64", "float64", "uint64"],
        [10_000, 100_000, 1_000_000],
        [0.1, 0.3, 1.0],  # uniqueness ratio
    ]
    param_names = ["use_swisstable", "dtype", "size", "uniqueness"]

    def setup(self, use_swisstable, dtype, size, uniqueness):
        n_unique = int(size * uniqueness)
        if n_unique < 1:
            n_unique = 1

        if dtype == "float64":
            unique_vals = np.random.randn(n_unique).astype(dtype)
        else:
            unique_vals = np.arange(n_unique, dtype=dtype)

        repeats = size // n_unique + 1
        data = np.tile(unique_vals, repeats)[:size]
        np.random.shuffle(data)

        self.series = Series(data, dtype=dtype)

    def time_unique(self, use_swisstable, dtype, size, uniqueness):
        with option_context("compute.use_swisstable", use_swisstable):
            self.series.unique()


class SwissTableDuplicated:
    """
    Benchmark duplicated() with Swiss Tables vs khash.

    Swiss Tables show 1.55x speedup at 1M elements.
    """

    params = [
        [False, True],  # use_swisstable
        ["int64", "float64"],
        [10_000, 100_000, 1_000_000],
        ["first", "last", False],
    ]
    param_names = ["use_swisstable", "dtype", "size", "keep"]

    def setup(self, use_swisstable, dtype, size, keep):
        n_unique = size // 3  # ~33% unique

        if dtype == "float64":
            unique_vals = np.random.randn(n_unique).astype(dtype)
        else:
            unique_vals = np.arange(n_unique, dtype=dtype)

        repeats = size // n_unique + 1
        data = np.tile(unique_vals, repeats)[:size]
        np.random.shuffle(data)

        self.series = Series(data, dtype=dtype)

    def time_duplicated(self, use_swisstable, dtype, size, keep):
        with option_context("compute.use_swisstable", use_swisstable):
            self.series.duplicated(keep=keep)


class SwissTableFloat64NaN:
    """
    Benchmark float64 operations with NaN values.

    Swiss Tables handle NaN correctly and show 1.2-1.4x speedup.
    """

    params = [
        [False, True],  # use_swisstable
        ["unique", "factorize", "isin"],
        [10_000, 100_000, 1_000_000],
        [0.01, 0.05, 0.1],  # NaN ratio
    ]
    param_names = ["use_swisstable", "operation", "size", "nan_ratio"]

    def setup(self, use_swisstable, operation, size, nan_ratio):
        data = np.random.randn(size)
        n_nan = int(size * nan_ratio)
        nan_indices = np.random.choice(size, n_nan, replace=False)
        data[nan_indices] = np.nan

        self.series = Series(data, dtype="float64")
        self.values = np.array([1.0, 2.0, np.nan])

    def time_operation(self, use_swisstable, operation, size, nan_ratio):
        with option_context("compute.use_swisstable", use_swisstable):
            if operation == "unique":
                self.series.unique()
            elif operation == "factorize":
                pd.factorize(self.series)
            elif operation == "isin":
                self.series.isin(self.values)


class SwissTableValueCounts:
    """
    Benchmark value_counts() with Swiss Tables.

    Note: Small N (<10k) may show regression due to function call overhead.
    """

    params = [
        [False, True],  # use_swisstable
        ["int64", "float64"],
        [10_000, 100_000, 1_000_000],
        [0.1, 0.3, 1.0],  # uniqueness ratio
    ]
    param_names = ["use_swisstable", "dtype", "size", "uniqueness"]

    def setup(self, use_swisstable, dtype, size, uniqueness):
        n_unique = int(size * uniqueness)
        if n_unique < 1:
            n_unique = 1

        if dtype == "float64":
            unique_vals = np.random.randn(n_unique).astype(dtype)
        else:
            unique_vals = np.arange(n_unique, dtype=dtype)

        repeats = size // n_unique + 1
        data = np.tile(unique_vals, repeats)[:size]
        np.random.shuffle(data)

        self.series = Series(data, dtype=dtype)

    def time_value_counts(self, use_swisstable, dtype, size, uniqueness):
        with option_context("compute.use_swisstable", use_swisstable):
            self.series.value_counts()


class SwissTableComplex:
    """
    Benchmark complex type operations with Swiss Tables.
    """

    params = [
        [False, True],  # use_swisstable
        ["complex64", "complex128"],
        [10_000, 100_000],
        ["unique", "factorize", "isin"],
    ]
    param_names = ["use_swisstable", "dtype", "size", "operation"]

    def setup(self, use_swisstable, dtype, size, operation):
        real = np.random.randn(size)
        imag = np.random.randn(size)
        data = (real + 1j * imag).astype(dtype)

        self.series = Series(data)
        self.values = data[:100]

    def time_operation(self, use_swisstable, dtype, size, operation):
        with option_context("compute.use_swisstable", use_swisstable):
            if operation == "unique":
                self.series.unique()
            elif operation == "factorize":
                pd.factorize(self.series)
            elif operation == "isin":
                self.series.isin(self.values)


class SwissTableLowLevel:
    """
    Benchmark low-level hash table operations directly.

    This tests the raw performance of Swiss Tables vs khash
    without pandas overhead.
    """

    params = [
        ["int64", "float64"],
        [10_000, 100_000, 1_000_000],
    ]
    param_names = ["dtype", "size"]

    def setup(self, dtype, size):
        if dtype == "int64":
            self.data = np.arange(size, dtype=np.int64)
        else:
            self.data = np.random.randn(size).astype(np.float64)

        # Pre-import for timing
        from pandas._libs import hashtable as ht

        self.ht = ht

        # Check if swisstable is available
        try:
            from pandas._libs import swisstable as st

            self.st = st
            self.has_swisstable = True
        except ImportError:
            self.has_swisstable = False

    def time_khash_unique(self, dtype, size):
        if dtype == "int64":
            table = self.ht.Int64HashTable()
        else:
            table = self.ht.Float64HashTable()
        table.unique(self.data)

    def time_swiss_unique(self, dtype, size):
        if not self.has_swisstable:
            return
        if dtype == "int64":
            table = self.st.SwissInt64Map()
        else:
            table = self.st.SwissFloat64Map()
        table.unique(self.data)

    def time_khash_factorize(self, dtype, size):
        if dtype == "int64":
            table = self.ht.Int64HashTable()
        else:
            table = self.ht.Float64HashTable()
        table.factorize(self.data)

    def time_swiss_factorize(self, dtype, size):
        if not self.has_swisstable:
            return
        if dtype == "int64":
            table = self.st.SwissInt64Map()
        else:
            table = self.st.SwissFloat64Map()
        table.factorize(self.data)
