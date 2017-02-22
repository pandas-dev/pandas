.. _remote_data:

.. currentmodule:: pandas

******************
Remote Data Access
******************

.. _remote_data.pandas_datareader:

DataReader
----------

The sub-package ``pandas.io.data`` is removed in favor of a separately
installable `pandas-datareader package
<https://github.com/pydata/pandas-datareader>`_. This will allow the data
modules to be independently updated to your pandas installation. The API for
``pandas-datareader v0.1.1`` is the same as in ``pandas v0.16.1``.
(:issue:`8961`)

   You should replace the imports of the following:

   .. code-block:: python

      from pandas.io import data, wb

   With:

   .. code-block:: python

      from pandas_datareader import data, wb
