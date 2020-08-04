from textwrap import dedent
from typing import Dict, Optional

from pandas.compat.numpy import function as nv
from pandas.util._decorators import Appender, Substitution, doc

from pandas.core.window.common import WindowGroupByMixin, _doc_template, _shared_docs
from pandas.core.window.rolling import _Rolling_and_Expanding


class Expanding(_Rolling_and_Expanding):
    """
    Provide expanding transformations.

    Parameters
    ----------
    min_periods : int, default 1
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    center : bool, default False
        Set the labels at the center of the window.
    axis : int or str, default 0

    Returns
    -------
    a Window sub-classed for the particular operation

    See Also
    --------
    rolling : Provides rolling window calculations.
    ewm : Provides exponential weighted functions.

    Notes
    -----
    By default, the result is set to the right edge of the window. This can be
    changed to the center of the window by setting ``center=True``.

    Examples
    --------
    >>> df = pd.DataFrame({"B": [0, 1, 2, np.nan, 4]})
    >>> df
         B
    0  0.0
    1  1.0
    2  2.0
    3  NaN
    4  4.0

    >>> df.expanding(2).sum()
         B
    0  NaN
    1  1.0
    2  3.0
    3  3.0
    4  7.0
    """

    _attributes = ["min_periods", "center", "axis"]

    def __init__(self, obj, min_periods=1, center=None, axis=0, **kwargs):
        super().__init__(obj=obj, min_periods=min_periods, center=center, axis=axis)

    @property
    def _constructor(self):
        return Expanding

    def _get_window(self, other=None, **kwargs):
        """
        Get the window length over which to perform some operation.

        Parameters
        ----------
        other : object, default None
            The other object that is involved in the operation.
            Such an object is involved for operations like covariance.

        Returns
        -------
        window : int
            The window length.
        """
        axis = self.obj._get_axis(self.axis)
        length = len(axis) + (other is not None) * len(axis)

        other = self.min_periods or -1
        return max(length, other)

    _agg_see_also_doc = dedent(
        """
    See Also
    --------
    pandas.DataFrame.aggregate : Similar DataFrame method.
    pandas.Series.aggregate : Similar Series method.
    """
    )

    _agg_examples_doc = dedent(
        """
    Examples
    --------
    >>> df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6], "C": [7, 8, 9]})
    >>> df
       A  B  C
    0  1  4  7
    1  2  5  8
    2  3  6  9

    >>> df.ewm(alpha=0.5).mean()
              A         B         C
    0  1.000000  4.000000  7.000000
    1  1.666667  4.666667  7.666667
    2  2.428571  5.428571  8.428571
    """
    )

    @doc(
        _shared_docs["aggregate"],
        see_also=_agg_see_also_doc,
        examples=_agg_examples_doc,
        versionadded="",
        klass="Series/Dataframe",
        axis="",
    )
    def aggregate(self, func, *args, **kwargs):
        return super().aggregate(func, *args, **kwargs)

    agg = aggregate

    @Substitution(name="expanding")
    @Appender(_shared_docs["count"])
    def count(self, **kwargs):
        return super().count(**kwargs)

    @Substitution(name="expanding")
    @Appender(_shared_docs["apply"])
    def apply(
        self,
        func,
        raw: bool = False,
        engine: Optional[str] = None,
        engine_kwargs: Optional[Dict[str, bool]] = None,
        args=None,
        kwargs=None,
    ):
        return super().apply(
            func,
            raw=raw,
            engine=engine,
            engine_kwargs=engine_kwargs,
            args=args,
            kwargs=kwargs,
        )

    @Substitution(name="expanding")
    @Appender(_shared_docs["sum"])
    def sum(self, *args, **kwargs):
        nv.validate_expanding_func("sum", args, kwargs)
        return super().sum(*args, **kwargs)

    @Substitution(name="expanding", func_name="max")
    @Appender(_doc_template)
    @Appender(_shared_docs["max"])
    def max(self, *args, **kwargs):
        nv.validate_expanding_func("max", args, kwargs)
        return super().max(*args, **kwargs)

    @Substitution(name="expanding")
    @Appender(_shared_docs["min"])
    def min(self, *args, **kwargs):
        nv.validate_expanding_func("min", args, kwargs)
        return super().min(*args, **kwargs)

    @Substitution(name="expanding")
    @Appender(_shared_docs["mean"])
    def mean(self, *args, **kwargs):
        nv.validate_expanding_func("mean", args, kwargs)
        return super().mean(*args, **kwargs)

    @Substitution(name="expanding")
    @Appender(_shared_docs["median"])
    def median(self, **kwargs):
        return super().median(**kwargs)

    @Substitution(name="expanding", versionadded="")
    @Appender(_shared_docs["std"])
    def std(self, ddof=1, *args, **kwargs):
        nv.validate_expanding_func("std", args, kwargs)
        return super().std(ddof=ddof, **kwargs)

    @Substitution(name="expanding", versionadded="")
    @Appender(_shared_docs["var"])
    def var(self, ddof=1, *args, **kwargs):
        nv.validate_expanding_func("var", args, kwargs)
        return super().var(ddof=ddof, **kwargs)

    @Substitution(name="expanding", func_name="skew")
    @Appender(_doc_template)
    @Appender(_shared_docs["skew"])
    def skew(self, **kwargs):
        return super().skew(**kwargs)

    _agg_doc = dedent(
        """
    Examples
    --------

    The example below will show an expanding calculation with a window size of
    four matching the equivalent function call using `scipy.stats`.

    >>> arr = [1, 2, 3, 4, 999]
    >>> import scipy.stats
    >>> print(f"{scipy.stats.kurtosis(arr[:-1], bias=False):.6f}")
    -1.200000
    >>> print(f"{scipy.stats.kurtosis(arr, bias=False):.6f}")
    4.999874
    >>> s = pd.Series(arr)
    >>> s.expanding(4).kurt()
    0         NaN
    1         NaN
    2         NaN
    3   -1.200000
    4    4.999874
    dtype: float64
    """
    )

    @Appender(_agg_doc)
    @Substitution(name="expanding")
    @Appender(_shared_docs["kurt"])
    def kurt(self, **kwargs):
        return super().kurt(**kwargs)

    @Substitution(name="expanding")
    @Appender(_shared_docs["quantile"])
    def quantile(self, quantile, interpolation="linear", **kwargs):
        return super().quantile(
            quantile=quantile, interpolation=interpolation, **kwargs
        )

    @Substitution(name="expanding", func_name="cov")
    @Appender(_doc_template)
    @Appender(_shared_docs["cov"])
    def cov(self, other=None, pairwise=None, ddof=1, **kwargs):
        return super().cov(other=other, pairwise=pairwise, ddof=ddof, **kwargs)

    @Substitution(name="expanding")
    @Appender(_shared_docs["corr"])
    def corr(self, other=None, pairwise=None, **kwargs):
        return super().corr(other=other, pairwise=pairwise, **kwargs)


class ExpandingGroupby(WindowGroupByMixin, Expanding):
    """
    Provide a expanding groupby implementation.
    """

    @property
    def _constructor(self):
        return Expanding
