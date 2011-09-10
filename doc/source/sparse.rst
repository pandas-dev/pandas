.. currentmodule:: pandas
.. _sparse:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   from pandas import *
   import pandas.util.testing as tm
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   import matplotlib.pyplot as plt
   plt.close('all')

**********************
Sparse data structures
**********************

We have implemented "sparse" versions of Series, DataFrame, and Panel. These
are not sparse in the typical "mostly 0". You can view these objects as being
"compressed" where any data matching a specific value (NaN/missing by default,
though any value can be chosen) is omitted. A special ``SparseIndex`` object
tracks where data has been "sparsified". This will make much more sense in an
example. All of the standard pandas data structures have a ``to_sparse``
method:

.. ipython:: python

   ts = Series(randn(10))
   ts[2:-2] = np.nan
   sts = ts.to_sparse()
   sts

The ``to_sparse`` method takes a ``kind`` argument (for the sparse index, see
below) and a ``fill_value``. So if we had a mostly zero Series, we could
convert it to sparse with ``fill_value=0``:

.. ipython:: python

   ts.fillna(0).to_sparse(fill_value=0)

The sparse objects exist for memory efficiency reasons. Suppose you had a
large, mostly NA DataFrame:

.. ipython:: python

   df = DataFrame(randn(10000, 4))
   df.ix[:9998] = np.nan
   sdf = df.to_sparse()
   sdf
   sdf.density

As you can see, the density (% of values that have not been "compressed") is
extremely low. This sparse object takes up much less memory on disk (pickled)
and in the Python interpreter. Functionally, their behavior should be nearly
identical to their dense counterparts. If not, you should report any
inconsistencies as bugs on GitHub.

Kinds of ``SparseIndex`` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two kinds of ``SparseIndex`` are implemented, ``block`` and ``integer``. We
recommend using ``block`` as it's more memory efficient. The ``integer`` format
keeps an arrays of all of the locations where the data are not equal to the
fill value. The ``block`` format tracks only the locations and sizes of blocks
of data.
