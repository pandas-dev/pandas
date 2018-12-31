{{ header }}

.. _api.general_functions:

=================
General functions
=================
.. currentmodule:: pandas

Data manipulations
~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   melt
   pivot
   pivot_table
   crosstab
   cut
   qcut
   merge
   merge_ordered
   merge_asof
   concat
   get_dummies
   factorize
   unique
   wide_to_long

Top-level missing data
~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   isna
   isnull
   notna
   notnull

Top-level conversions
~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   to_numeric

Top-level dealing with datetimelike
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   to_datetime
   to_timedelta
   date_range
   bdate_range
   period_range
   timedelta_range
   infer_freq

Top-level dealing with intervals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   interval_range

Top-level evaluation
~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   eval

Hashing
~~~~~~~
.. autosummary::
   :toctree: generated/

   util.hash_array
   util.hash_pandas_object

Testing
~~~~~~~
.. autosummary::
   :toctree: generated/

   test
