.. _performance_considerations:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import pandas as pd
   import numpy as np
   np.random.seed(123456)

******************
Performance Considerations
******************

Overview
--------
Pandas, for most use cases, will be fast; however, there are a few
anti-patterns to avoid. The following is a guide on achieving better
performance with pandas.

.. ipython:: python

   import sys

   import numpy as np
   import pandas as pd
   import matplotlib.pyplot as plt

   pd.options.display.max_rows = 10
   plt.style.use('default')

1) Using pandas when it is the wrong tool.

   * At least not for things it's not meant for.
   * Pandas is very fast at joins, reindex, factorization
   * Not as great at, say, matrix multiplications or problems that aren't vectorizable



2) Another optimization is to make sure to use native types instead of pandas types.

.. ipython:: python

   s1 = pd.Series(range(10000), dtype=object)
   s2 = pd.Series(range(10000), dtype=np.int64)

.. ipython:: python

   %timeit s1.sum()

.. ipython:: python

   %timeit s2.sum()


NumPy int64 datatype arrays are much faster than the python object version.

For other datatypes:

   * Strings - This is usually unavoidable. Pandas 2 will have a specialized
     string type, but for now you're stuck with python objects.
     If you have few distinct values (relative to the number of rows), you could
     use a Categorical.

   * Dates, Times - Pandas has implemented a specialized version of datetime.datetime,
     and datetime.timedelta, but not datetime.date or datetime.time. Depending on your
     application, you might be able to treat dates as datetimess, at midnight.

   * Decimal Types - Pandas uses floating-point arrays; there isn't a
     native arbitrary-precision Decimal type.


In certain cases, some things will automatically be object type:

   * Reading messy Excel files - read_excel will preserve the dtype
     of each cell in the spreadsheet. If you have a single column with
     an int, a float, and a datetime, pandas will have to store all of
     those as objects. This dataset probably isn't tidy though.

   * Integer NA - Unfortunately, pandas doesn't have real nullable types.
     To represent missingness, pandas uses NaN (not a number) which is a special
     floating point value. If you have to represent nullable integers, you can
     use object dtype.

.. ipython:: python

   s = pd.Series([1, 2, 3, np.nan, 5, 6, 7, 8, 9])
   s


.. ipython:: python

   type(s[0])

.. ipython:: python

   s = pd.Series([1, 2, 3, np.nan, 5, 6, 7, 8, 9], dtype=object)
   type(s[0])

