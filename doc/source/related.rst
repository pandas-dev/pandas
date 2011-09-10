************************
Related Python libraries
************************

la (larry)
----------

Keith Goodman's excellent `labeled array package
<http://pypi.python.org/pypi/la>`__ is very similar to pandas in many regards,
though with some key differences. The main philosophical design difference is
to be a wrapper around a single NumPy ``ndarray`` object while adding axis
labeling and label-based operations and indexing. Because of this, creating a
size-mutable object with heterogeneous columns (e.g. DataFrame) is not possible
with the ``la`` package.

  - Provide a single n-dimensional object with labeled axes with functionally
    analogous data alignment semantics to pandas objects
  - Advanced / label-based indexing similar to that provided in pandas but
    setting is not supported
  - Stays much closer to NumPy arrays than pandas-- ``larry`` objects must be
    homogeneously typed
  - GroupBy support is relatively limited, but a few functions are available:
    ``group_mean``, ``group_median``, and ``group_ranking``
  - It has a collection of analytical functions suited to quantitative
    portfolio construction for financial applications
  - It has a collection of moving window statistics implemented in
    `Bottleneck <http://pypi.python.org/pypi/Bottleneck>`__

scikits.statsmodels
-------------------

The main `statistics and econometrics library
<http://statsmodels.sourceforge.net>`__ for Python. pandas has become a
dependency of this library.

scikits.timeseries
------------------

`scikits.timeseries <http://pytseries.sourceforge.net/>`__ provides a data
structure for fixed frequency time series data based on the numpy.MaskedArray
class. For time series data, it provides some of the same functionality to the
pandas Series class. It has many more functions for time series-specific
manipulation. Also, it has support for many more frequencies, though less
customizable by the user (so 5-minutely data is easier to do with pandas for
example).

We are aiming to merge these libraries together in the near future.
