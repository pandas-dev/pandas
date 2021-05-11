In pandas, :meth:`Series.isna` and :meth:`Series.notna` can be used to filter the rows.

.. ipython:: python

   outer_join[outer_join["value_x"].isna()]
   outer_join[outer_join["value_x"].notna()]

pandas provides :ref:`a variety of methods to work with missing data <missing_data>`. Here are some examples:

Drop rows with missing values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python

   outer_join.dropna()

Forward fill from previous rows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python

   outer_join.fillna(method="ffill")

Replace missing values with a specified value
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the mean:

.. ipython:: python

   outer_join["value_x"].fillna(outer_join["value_x"].mean())
