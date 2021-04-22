pandas DataFrames have a :meth:`~DataFrame.merge` method, which provides similar functionality. The
data does not have to be sorted ahead of time, and different join types are accomplished via the
``how`` keyword.

.. ipython:: python

   inner_join = df1.merge(df2, on=["key"], how="inner")
   inner_join

   left_join = df1.merge(df2, on=["key"], how="left")
   left_join

   right_join = df1.merge(df2, on=["key"], how="right")
   right_join

   outer_join = df1.merge(df2, on=["key"], how="outer")
   outer_join
