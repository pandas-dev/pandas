from functools import wraps
import threading

import numpy as np

from pandas import (
    DataFrame,
    Index,
    Series,
    date_range,
    factorize,
    read_csv,
)
from pandas.core.algorithms import take_nd

try:
    from pandas import (
        rolling_kurt,
        rolling_max,
        rolling_mean,
        rolling_median,
        rolling_min,
        rolling_skew,
        rolling_std,
        rolling_var,
    )

    have_rolling_methods = True
except ImportError:
    have_rolling_methods = False
try:
    from pandas._libs import algos
except ImportError:
    from pandas import algos

from .pandas_vb_common import BaseIO  # isort:skip


def run_parallel(num_threads=2, kwargs_list=None):
    """
    Decorator to run the same function multiple times in parallel.

    Parameters
    ----------
    num_threads : int, optional
        The number of times the function is run in parallel.
    kwargs_list : list of dicts, optional
        The list of kwargs to update original
        function kwargs on different threads.

    Notes
    -----
    This decorator does not pass the return value of the decorated function.

    Original from scikit-image:

    https://github.com/scikit-image/scikit-image/pull/1519

    """
    assert num_threads > 0
    has_kwargs_list = kwargs_list is not None
    if has_kwargs_list:
        assert len(kwargs_list) == num_threads

    def wrapper(func):
        @wraps(func)
        def inner(*args, **kwargs):
            if has_kwargs_list:
                update_kwargs = lambda i: dict(kwargs, **kwargs_list[i])
            else:
                update_kwargs = lambda i: kwargs
            threads = []
            for i in range(num_threads):
                updated_kwargs = update_kwargs(i)
                thread = threading.Thread(target=func, args=args, kwargs=updated_kwargs)
                threads.append(thread)
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()

        return inner

    return wrapper


class ParallelGroupbyMethods:
    params = ([2, 4, 8], ["count", "last", "max", "mean", "min", "prod", "sum", "var"])
    param_names = ["threads", "method"]

    def setup(self, threads, method):
        N = 10**6
        ngroups = 10**3
        df = DataFrame(
            {"key": np.random.randint(0, ngroups, size=N), "data": np.random.randn(N)}
        )

        @run_parallel(num_threads=threads)
        def parallel():
            getattr(df.groupby("key")["data"], method)()

        self.parallel = parallel

        def loop():
            getattr(df.groupby("key")["data"], method)()

        self.loop = loop

    def time_parallel(self, threads, method):
        self.parallel()

    def time_loop(self, threads, method):
        for i in range(threads):
            self.loop()


class ParallelGroups:
    params = [2, 4, 8]
    param_names = ["threads"]

    def setup(self, threads):
        size = 2**22
        ngroups = 10**3
        data = Series(np.random.randint(0, ngroups, size=size))

        @run_parallel(num_threads=threads)
        def get_groups():
            data.groupby(data).groups

        self.get_groups = get_groups

    def time_get_groups(self, threads):
        self.get_groups()


class ParallelTake1D:
    params = ["int64", "float64"]
    param_names = ["dtype"]

    def setup(self, dtype):
        N = 10**6
        df = DataFrame({"col": np.arange(N, dtype=dtype)})
        indexer = np.arange(100, len(df) - 100)

        @run_parallel(num_threads=2)
        def parallel_take1d():
            take_nd(df["col"].values, indexer)

        self.parallel_take1d = parallel_take1d

    def time_take1d(self, dtype):
        self.parallel_take1d()


class ParallelKth:
    # This depends exclusively on code in _libs/, could go in libs.py

    number = 1
    repeat = 5

    def setup(self):
        N = 10**7
        k = 5 * 10**5
        kwargs_list = [{"arr": np.random.randn(N)}, {"arr": np.random.randn(N)}]

        @run_parallel(num_threads=2, kwargs_list=kwargs_list)
        def parallel_kth_smallest(arr):
            algos.kth_smallest(arr, k)

        self.parallel_kth_smallest = parallel_kth_smallest

    def time_kth_smallest(self):
        self.parallel_kth_smallest()


class ParallelDatetimeFields:
    def setup(self):
        N = 10**6
        self.dti = date_range("1900-01-01", periods=N, freq="min")
        self.period = self.dti.to_period("D")

    def time_datetime_field_year(self):
        @run_parallel(num_threads=2)
        def run(dti):
            dti.year

        run(self.dti)

    def time_datetime_field_day(self):
        @run_parallel(num_threads=2)
        def run(dti):
            dti.day

        run(self.dti)

    def time_datetime_field_daysinmonth(self):
        @run_parallel(num_threads=2)
        def run(dti):
            dti.days_in_month

        run(self.dti)

    def time_datetime_field_normalize(self):
        @run_parallel(num_threads=2)
        def run(dti):
            dti.normalize()

        run(self.dti)

    def time_datetime_to_period(self):
        @run_parallel(num_threads=2)
        def run(dti):
            dti.to_period("s")

        run(self.dti)

    def time_period_to_datetime(self):
        @run_parallel(num_threads=2)
        def run(period):
            period.to_timestamp()

        run(self.period)


class ParallelRolling:
    params = ["median", "mean", "min", "max", "var", "skew", "kurt", "std"]
    param_names = ["method"]

    def setup(self, method):
        win = 100
        arr = np.random.rand(100000)
        if hasattr(DataFrame, "rolling"):
            df = DataFrame(arr).rolling(win)

            @run_parallel(num_threads=2)
            def parallel_rolling():
                getattr(df, method)()

            self.parallel_rolling = parallel_rolling
        elif have_rolling_methods:
            rolling = {
                "median": rolling_median,
                "mean": rolling_mean,
                "min": rolling_min,
                "max": rolling_max,
                "var": rolling_var,
                "skew": rolling_skew,
                "kurt": rolling_kurt,
                "std": rolling_std,
            }

            @run_parallel(num_threads=2)
            def parallel_rolling():
                rolling[method](arr, win)

            self.parallel_rolling = parallel_rolling
        else:
            raise NotImplementedError

    def time_rolling(self, method):
        self.parallel_rolling()


class ParallelReadCSV(BaseIO):
    number = 1
    repeat = 5
    params = ["float", "object", "datetime"]
    param_names = ["dtype"]

    def setup(self, dtype):
        rows = 10000
        cols = 50
        if dtype == "float":
            df = DataFrame(np.random.randn(rows, cols))
        elif dtype == "datetime":
            df = DataFrame(
                np.random.randn(rows, cols), index=date_range("1/1/2000", periods=rows)
            )
        elif dtype == "object":
            df = DataFrame(
                "foo", index=range(rows), columns=["object%03d" for _ in range(5)]
            )
        else:
            raise NotImplementedError

        self.fname = f"__test_{dtype}__.csv"
        df.to_csv(self.fname)

        @run_parallel(num_threads=2)
        def parallel_read_csv():
            read_csv(self.fname)

        self.parallel_read_csv = parallel_read_csv

    def time_read_csv(self, dtype):
        self.parallel_read_csv()


class ParallelFactorize:
    number = 1
    repeat = 5
    params = [2, 4, 8]
    param_names = ["threads"]

    def setup(self, threads):
        strings = Index([f"i-{i}" for i in range(100000)], dtype=object)

        @run_parallel(num_threads=threads)
        def parallel():
            factorize(strings)

        self.parallel = parallel

        def loop():
            factorize(strings)

        self.loop = loop

    def time_parallel(self, threads):
        self.parallel()

    def time_loop(self, threads):
        for i in range(threads):
            self.loop()


from .pandas_vb_common import setup  # noqa: F401 isort:skip
