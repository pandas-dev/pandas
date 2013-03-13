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

`Extending a panel along the minor axis
<http://stackoverflow.com/questions/15364050/extending-a-pandas-panel-frame-along-the-minor-axis>`__

Grouping
--------

`Basic grouping with apply
<http://stackoverflow.com/questions/15322632/python-pandas-df-groupy-agg-column-reference-in-agg>`__

`TimeGrouping of values grouped across time
<http://stackoverflow.com/questions/15297053/how-can-i-divide-single-values-of-a-dataframe-by-monthly-averages>`__

Merge
-----

Join
~~~~

`Joining a DataFrame to itself
<https://github.com/pydata/pandas/issues/2996>`__

Data In/Out
-----------

CSV
~~~

HDF5
~~~~

Managing heteregenous data using a linked multiple table hierarchy. 
See `here <https://github.com/pydata/pandas/issues/3032>`__.
