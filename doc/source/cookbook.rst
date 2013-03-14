.. _cookbook:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   import random
   import os
   np.random.seed(123456)
   from pandas import *
   import pandas as pd
   randn = np.random.randn
   randint = np.random.randint
   np.set_printoptions(precision=4, suppress=True)

********
Cookbook
********

This is a respository for *short and sweet* example and links for useful pandas recipes.
We encourage users to add to this documentation. This is a great *First Pull Request*.

Selection
---------
`Boolean Rows Indexing
<http://stackoverflow.com/questions/14725068/pandas-using-row-labels-in-boolean-indexing>`__

`Extending a panel along the minor axis
<http://stackoverflow.com/questions/15364050/extending-a-pandas-panel-frame-along-the-minor-axis>`__

`Using loc and iloc in selections
<https://github.com/pydata/pandas/issues/2904>`__

MultiIndexing
-------------

`Prepending a level to a multiindex
<http://stackoverflow.com/questions/14744068/prepend-a-level-to-a-pandas-multiindex>`__

`Slicing a multi-index with xs
<http://stackoverflow.com/questions/12590131/how-to-slice-multindex-columns-in-pandas-dataframes>`__

`Multi-index sorting
<http://stackoverflow.com/questions/14733871/mutli-index-sorting-in-pandas>`__

`Partial Selection, the need for sortedness
<https://github.com/pydata/pandas/issues/2995>`__

Grouping
--------

`Basic grouping with apply
<http://stackoverflow.com/questions/15322632/python-pandas-df-groupy-agg-column-reference-in-agg>`__

`Apply to different items in a group
<http://stackoverflow.com/questions/15262134/apply-different-functions-to-different-items-in-group-object-python-pandas>`__

`Replacing values with groupby means
<http://stackoverflow.com/questions/14760757/replacing-values-with-groupby-means>`__

`TimeGrouping of values grouped across time
<http://stackoverflow.com/questions/15297053/how-can-i-divide-single-values-of-a-dataframe-by-monthly-averages>`__

Merge
-----

Join
~~~~

`Joining a DataFrame to itself
<https://github.com/pydata/pandas/issues/2996>`__

Timeseries
----------

`Resample intraday frame without adding new days
<http://stackoverflow.com/questions/14898574/resample-intrday-pandas-dataframe-without-add-new-days>`__

Data In/Out
-----------

CSV
~~~

HDF5
~~~~

`Managing heteregenous data using a linked multiple table hierarchy
<https://github.com/pydata/pandas/issues/3032>`__

`Simple Queries with a Timestamp Index
<http://stackoverflow.com/questions/13926089/selecting-columns-from-pandas-hdfstore-table>`__

Miscellaneous
-------------

