{{ header }}

.. _api.general_functions:

=================
General functions
=================
.. currentmodule:: pandas

Data manipulations
~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

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
   from_dummies
   factorize
   unique
   lreshape
   wide_to_long

Top-level missing data
~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   isna
   isnull
   notna
   notnull

Top-level dealing with numeric data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   to_numeric

Top-level dealing with datetimelike data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   to_datetime
   to_timedelta
   date_range
   bdate_range
   period_range
   timedelta_range
   infer_freq

Top-level dealing with Interval data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   interval_range

Top-level evaluation
~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   eval

Datetime formats
~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   tseries.api.guess_datetime_format

Hashing
~~~~~~~
.. autosummary::
   :toctree: api/

   util.hash_array
   util.hash_pandas_object

Importing from other DataFrame libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   api.interchange.from_dataframe
