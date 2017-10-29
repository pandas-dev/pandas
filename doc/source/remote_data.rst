.. _remote_data:

.. currentmodule:: pandas

******************
Remote Data Access
******************

.. _remote_data.pandas_datareader:

DataReader
----------

The sub-package ``pandas.io.data`` was deprecated in v.0.17 and removed in
`v.0.19 <http://pandas-docs.github.io/pandas-docs-travis/whatsnew.html#v0-19-0-october-2-2016>`__.
Instead there has been created a separately installable
`pandas-datareader package <https://github.com/pydata/pandas-datareader>`__.
This will allow the data modules to be independently updated on your pandas installation.

For code older than < 0.19 you should replace the imports of the following:

.. code-block:: python

   from pandas.io import data, wb

With:

.. code-block:: python

   from pandas_datareader import data, wb
