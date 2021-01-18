from textwrap import dedent
from typing import Any, Callable, Dict, Optional, Tuple, Union

import numpy as np

from pandas._typing import FrameOrSeries
from pandas.compat.numpy import function as nv
from pandas.util._decorators import doc

from pandas.core.window.doc import (
    _shared_docs,
    doc_template,
    kwargs_compat,
    numba_notes,
    numpy_args_kwargs,
    window_agg_numba_args_kwargs_parameters,
    window_agg_numba_parameters,
    window_apply_parameters,
)
from pandas.core.window.indexers import BaseIndexer, ExpandingIndexer, GroupbyIndexer
from pandas.core.window.rolling import BaseWindowGroupby, RollingAndExpandingMixin


class Expanding(RollingAndExpandingMixin):
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
    method : str {'single', 'table'}, default 'single'
        Execute the rolling operation per single column or row (``'single'``)
        or over the entire object (``'table'``).

        This argument is only implemented when specifying ``engine='numba'``
        in the method call.

        .. versionadded:: 1.3.0

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

    _attributes = ["min_periods", "center", "axis", "method"]

    def __init__(
        self, obj, min_periods=1, center=None, axis=0, method="single", **kwargs
    ):
        super().__init__(
            obj=obj, min_periods=min_periods, center=center, axis=axis, method=method
        )

    def _get_window_indexer(self) -> BaseIndexer:
        """
        Return an indexer class that will compute the window start and end bounds
        """
        return ExpandingIndexer()

    def _get_cov_corr_window(
        self, other: Optional[Union[np.ndarray, FrameOrSeries]] = None, **kwargs
    ) -> int:
        """
        Get the window length over which to perform cov and corr operations.

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

    @doc(
        _shared_docs["aggregate"],
        see_also=dedent(
            """
        See Also
        --------
        pandas.DataFrame.aggregate : Similar DataFrame method.
        pandas.Series.aggregate : Similar Series method.
        """
        ),
        examples=dedent(
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
        ),
        klass="Series/Dataframe",
        axis="",
    )
    def aggregate(self, func, *args, **kwargs):
        return super().aggregate(func, *args, **kwargs)

    agg = aggregate

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="count of non NaN observations",
        parameters=kwargs_compat,
        numpy_args_kwargs="",
        agg_method="count",
        other_see_also="",
        notes="",
        examples="",
    )
    def count(self):
        return super().count()

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="custom aggregation function",
        parameters=window_apply_parameters,
        numpy_args_kwargs="",
        agg_method="apply",
        other_see_also="",
        notes=numba_notes,
        examples="",
    )
    def apply(
        self,
        func: Callable[..., Any],
        raw: bool = False,
        engine: Optional[str] = None,
        engine_kwargs: Optional[Dict[str, bool]] = None,
        args: Optional[Tuple[Any, ...]] = None,
        kwargs: Optional[Dict[str, Any]] = None,
    ):
        return super().apply(
            func,
            raw=raw,
            engine=engine,
            engine_kwargs=engine_kwargs,
            args=args,
            kwargs=kwargs,
        )

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="sum",
        parameters=window_agg_numba_parameters,
        numpy_args_kwargs=numpy_args_kwargs,
        agg_method="sum",
        other_see_also="",
        notes=numba_notes,
        examples="",
    )
    def sum(self, *args, engine=None, engine_kwargs=None, **kwargs):
        nv.validate_expanding_func("sum", args, kwargs)
        return super().sum(*args, engine=engine, engine_kwargs=engine_kwargs, **kwargs)

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="maximum",
        parameters=window_agg_numba_args_kwargs_parameters,
        numpy_args_kwargs="",
        agg_method="max",
        other_see_also="",
        notes=numba_notes,
        examples="",
    )
    def max(self, *args, engine=None, engine_kwargs=None, **kwargs):
        nv.validate_expanding_func("max", args, kwargs)
        return super().max(*args, engine=engine, engine_kwargs=engine_kwargs, **kwargs)

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="minimum",
        parameters=window_agg_numba_args_kwargs_parameters,
        numpy_args_kwargs="",
        agg_method="min",
        other_see_also="",
        notes=numba_notes,
        examples="",
    )
    def min(self, *args, engine=None, engine_kwargs=None, **kwargs):
        nv.validate_expanding_func("min", args, kwargs)
        return super().min(*args, engine=engine, engine_kwargs=engine_kwargs, **kwargs)

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="mean",
        parameters=window_agg_numba_args_kwargs_parameters,
        numpy_args_kwargs="",
        agg_method="mean",
        other_see_also="",
        notes=numba_notes,
        examples="",
    )
    def mean(self, *args, engine=None, engine_kwargs=None, **kwargs):
        nv.validate_expanding_func("mean", args, kwargs)
        return super().mean(*args, engine=engine, engine_kwargs=engine_kwargs, **kwargs)

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="median",
        parameters=window_agg_numba_parameters,
        numpy_args_kwargs="",
        agg_method="median",
        other_see_also="",
        notes=numba_notes,
        examples="",
    )
    def median(self, engine=None, engine_kwargs=None, **kwargs):
        return super().median(engine=engine, engine_kwargs=engine_kwargs, **kwargs)

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="standard deviation",
        parameters=dedent(
            """
        ddof : int, default 1
            Delta Degrees of Freedom.  The divisor used in calculations
            is ``N - ddof``, where ``N`` represents the number of elements.
        """
        ),
        numpy_args_kwargs=numpy_args_kwargs,
        agg_method="std",
        other_see_also="numpy.std : Equivalent method for Numpy array.\n",
        notes=dedent(
            """
        The default ``ddof`` of 1 used in :meth:`Series.std` is different
        than the default ``ddof`` of 0 in :func:`numpy.std`.

        A minimum of one period is required for the rolling calculation.
        """
        ),
        examples=dedent(
            """
        >>> s = pd.Series([5, 5, 6, 7, 5, 5, 5])

        >>> s.expanding(3).std()
        0         NaN
        1         NaN
        2    0.577350
        3    0.957427
        4    0.894427
        5    0.836660
        6    0.786796
        dtype: float64
        """
        ),
    )
    def std(self, ddof: int = 1, *args, **kwargs):
        nv.validate_expanding_func("std", args, kwargs)
        return super().std(ddof=ddof, **kwargs)

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="variance",
        parameters=dedent(
            """
        ddof : int, default 1
            Delta Degrees of Freedom.  The divisor used in calculations
            is ``N - ddof``, where ``N`` represents the number of elements.
        """
        ),
        numpy_args_kwargs=numpy_args_kwargs,
        agg_method="var",
        other_see_also="numpy.var : Equivalent method for Numpy array.\n",
        notes=dedent(
            """
        The default ``ddof`` of 1 used in :meth:`Series.var` is different
        than the default ``ddof`` of 0 in :func:`numpy.var`.

        A minimum of one period is required for the rolling calculation.
        """
        ),
        examples=dedent(
            """
        >>> s = pd.Series([5, 5, 6, 7, 5, 5, 5])

        >>> s.expanding(3).var()
        0         NaN
        1         NaN
        2    0.333333
        3    0.916667
        4    0.800000
        5    0.700000
        6    0.619048
        dtype: float64
        """
        ),
    )
    def var(self, ddof: int = 1, *args, **kwargs):
        nv.validate_expanding_func("var", args, kwargs)
        return super().var(ddof=ddof, **kwargs)

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="standard error of mean",
        parameters=dedent(
            """
        ddof : int, default 1
            Delta Degrees of Freedom.  The divisor used in calculations
            is ``N - ddof``, where ``N`` represents the number of elements.
        """
        ),
        numpy_args_kwargs="",
        agg_method="sem",
        other_see_also="",
        notes="A minimum of one period is required for the calculation.",
        examples=dedent(
            """
        >>> s = pd.Series([0, 1, 2, 3])

        >>> s.expanding().sem()
        0         NaN
        1    0.707107
        2    0.707107
        3    0.745356
        dtype: float64
        """
        ),
    )
    def sem(self, ddof: int = 1, *args, **kwargs):
        return super().sem(ddof=ddof, **kwargs)

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="unbiased skewness",
        parameters=kwargs_compat,
        numpy_args_kwargs=numpy_args_kwargs,
        agg_method="skew",
        other_see_also="scipy.stats.skew : Third moment of a probability density.\n",
        notes="A minimum of three periods is required for the rolling calculation.",
        examples="",
    )
    def skew(self, **kwargs):
        return super().skew(**kwargs)

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="Fisher's definition of kurtosis without bias",
        parameters=kwargs_compat,
        numpy_args_kwargs="",
        agg_method="kurt",
        other_see_also="scipy.stats.kurtosis : Reference SciPy method.\n",
        notes=dedent(
            """
        A minimum of four periods is required for the calculation.
        """
        ),
        examples=dedent(
            """
        The example below will show a rolling calculation with a window size of
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
        ),
    )
    def kurt(self, **kwargs):
        return super().kurt(**kwargs)

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="quantile",
        parameters=dedent(
            """
        quantile : float
            Quantile to compute. 0 <= quantile <= 1.
        interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
            This optional parameter specifies the interpolation method to use,
            when the desired quantile lies between two data points `i` and `j`:

                * linear: `i + (j - i) * fraction`, where `fraction` is the
                  fractional part of the index surrounded by `i` and `j`.
                * lower: `i`.
                * higher: `j`.
                * nearest: `i` or `j` whichever is nearest.
                * midpoint: (`i` + `j`) / 2.
        **kwargs
            For function compatibility and will not have an effect on the result.
        """
        ),
        numpy_args_kwargs="",
        agg_method="quantile",
        other_see_also="",
        notes="",
        examples="",
    )
    def quantile(
        self,
        quantile,
        interpolation="linear",
        **kwargs,
    ):
        return super().quantile(
            quantile=quantile,
            interpolation=interpolation,
            **kwargs,
        )

    @doc(
        doc_template,
        window_method="expanding",
        aggregation_description="sample covariance",
        parameters=dedent(
            """
        other : Series, DataFrame, or ndarray, optional
            If not supplied then will default to self and produce pairwise
            output.
        pairwise : bool, default None
            If False then only matching columns between self and other will be
            used and the output will be a DataFrame.
            If True then all pairwise combinations will be calculated and the
            output will be a MultiIndexed DataFrame in the case of DataFrame
            inputs. In the case of missing elements, only complete pairwise
            observations will be used.
        ddof : int, default 1
            Delta Degrees of Freedom.  The divisor used in calculations
            is ``N - ddof``, where ``N`` represents the number of elements.
        **kwargs
            For function compatibility and will not have an effect on the result.
        """
        ),
        numpy_args_kwargs="",
        agg_method="cov",
        other_see_also="",
        notes="",
        examples="",
    )
    def cov(
        self,
        other: Optional[Union[np.ndarray, FrameOrSeries]] = None,
        pairwise: Optional[bool] = None,
        ddof: int = 1,
        **kwargs,
    ):
        return super().cov(other=other, pairwise=pairwise, ddof=ddof, **kwargs)

    @doc(
        doc_template,
        window_method="rolling",
        aggregation_description="covariance",
        parameters=dedent(
            """
        other : Series, DataFrame, or ndarray, optional
            If not supplied then will default to self.
        pairwise : bool, default None
            Calculate pairwise combinations of columns within a
            DataFrame. If `other` is not specified, defaults to `True`,
            otherwise defaults to `False`.
            Not relevant for :class:`~pandas.Series`.
        **kwargs
            For function compatibility and will not have an effect on the result.
        """
        ),
        numpy_args_kwargs="",
        agg_method="corr",
        other_see_also=dedent(
            """
        cov : Similar method to calculate covariance.
        numpy.corrcoef : NumPy Pearson's correlation calculation.
        """
        ),
        notes=dedent(
            """
        This function uses Pearson's definition of correlation
        (https://en.wikipedia.org/wiki/Pearson_correlation_coefficient).

        When `other` is not specified, the output will be self correlation (e.g.
        all 1's), except for :class:`~pandas.DataFrame` inputs with `pairwise`
        set to `True`.

        Function will return ``NaN`` for correlations of equal valued sequences;
        this is the result of a 0/0 division error.

        When `pairwise` is set to `False`, only matching columns between `self` and
        `other` will be used.

        When `pairwise` is set to `True`, the output will be a MultiIndex DataFrame
        with the original index on the first level, and the `other` DataFrame
        columns on the second level.

        In the case of missing elements, only complete pairwise observations
        will be used.
        """
        ),
        examples="",
    )
    def corr(
        self,
        other: Optional[Union[np.ndarray, FrameOrSeries]] = None,
        pairwise: Optional[bool] = None,
        **kwargs,
    ):
        return super().corr(other=other, pairwise=pairwise, **kwargs)


class ExpandingGroupby(BaseWindowGroupby, Expanding):
    """
    Provide a expanding groupby implementation.
    """

    @property
    def _constructor(self):
        return Expanding

    def _get_window_indexer(self) -> GroupbyIndexer:
        """
        Return an indexer class that will compute the window start and end bounds

        Returns
        -------
        GroupbyIndexer
        """
        window_indexer = GroupbyIndexer(
            groupby_indicies=self._groupby.indices,
            window_indexer=ExpandingIndexer,
        )
        return window_indexer
