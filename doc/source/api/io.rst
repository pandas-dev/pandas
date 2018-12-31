{{ header }}

.. _api.io:

============
Input/Output
============
.. currentmodule:: pandas

Pickling
~~~~~~~~
.. autosummary::
   :toctree: generated/

   read_pickle

Flat File
~~~~~~~~~
.. autosummary::
   :toctree: generated/

   read_table
   read_csv
   read_fwf
   read_msgpack

Clipboard
~~~~~~~~~
.. autosummary::
   :toctree: generated/

   read_clipboard

Excel
~~~~~
.. autosummary::
   :toctree: generated/

   read_excel
   ExcelFile.parse

.. autosummary::
   :toctree: generated/
   :template: autosummary/class_without_autosummary.rst

   ExcelWriter

JSON
~~~~
.. autosummary::
   :toctree: generated/

   read_json

.. currentmodule:: pandas.io.json

.. autosummary::
   :toctree: generated/

   json_normalize
   build_table_schema

.. currentmodule:: pandas

HTML
~~~~
.. autosummary::
   :toctree: generated/

   read_html

HDFStore: PyTables (HDF5)
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   read_hdf
   HDFStore.put
   HDFStore.append
   HDFStore.get
   HDFStore.select
   HDFStore.info
   HDFStore.keys
   HDFStore.groups
   HDFStore.walk

Feather
~~~~~~~
.. autosummary::
   :toctree: generated/

   read_feather

Parquet
~~~~~~~
.. autosummary::
   :toctree: generated/

   read_parquet

SAS
~~~
.. autosummary::
   :toctree: generated/

   read_sas

SQL
~~~
.. autosummary::
   :toctree: generated/

   read_sql_table
   read_sql_query
   read_sql

Google BigQuery
~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   read_gbq

STATA
~~~~~
.. autosummary::
   :toctree: generated/

   read_stata

.. currentmodule:: pandas.io.stata

.. autosummary::
   :toctree: generated/

   StataReader.data
   StataReader.data_label
   StataReader.value_labels
   StataReader.variable_labels
   StataWriter.write_file
