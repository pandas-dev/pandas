This doesn't work in pandas.  Instead, the :func:`pd.isna` or :func:`pd.notna` functions
should be used for comparisons.

.. ipython:: python

   outer_join[pd.isna(outer_join["value_x"])]
   outer_join[pd.notna(outer_join["value_x"])]

pandas also provides a variety of methods to work with missing data -- some of
which would be challenging to express in Stata. For example, there are methods to
drop all rows with any missing values, replacing missing values with a specified
value, like the mean, or forward filling from previous rows. See the
:ref:`missing data documentation<missing_data>` for more.

.. ipython:: python

   # Drop rows with any missing value
   outer_join.dropna()

   # Fill forwards
   outer_join.fillna(method="ffill")

   # Impute missing values with the mean
   outer_join["value_x"].fillna(outer_join["value_x"].mean())
