{{ header }}

.. _api.io:

============
Input/output
============
.. currentmodule:: pandas

Pickling
~~~~~~~~
.. autosummary::
   :toctree: api/

   read_pickle

Flat file
~~~~~~~~~
.. autosummary::
   :toctree: api/

   read_table
   read_csv
   read_fwf
   read_msgpack

Clipboard
~~~~~~~~~
.. autosummary::
   :toctree: api/

   read_clipboard

Excel
~~~~~
.. autosummary::
   :toctree: api/

   read_excel
   ExcelFile.parse

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   ExcelWriter

JSON
~~~~
.. autosummary::
   :toctree: api/

   read_json

.. currentmodule:: pandas.io.json

.. autosummary::
   :toctree: api/

   json_normalize
   build_table_schema

.. currentmodule:: pandas

HTML
~~~~
.. autosummary::
   :toctree: api/

   read_html

HDFStore: PyTables (HDF5)
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

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
   :toctree: api/

   read_feather

Parquet
~~~~~~~
.. autosummary::
   :toctree: api/

   read_parquet

SAS
~~~
.. autosummary::
   :toctree: api/

   read_sas

SPSS
~~~~
.. autosummary::
   :toctree: api/

   read_spss

SQL
~~~
.. autosummary::
   :toctree: api/

   read_sql_table
   read_sql_query
   read_sql

Google BigQuery
~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   read_gbq

STATA
~~~~~
.. autosummary::
   :toctree: api/

   read_stata

.. currentmodule:: pandas.io.stata

.. autosummary::
   :toctree: api/

   StataReader.data_label
   StataReader.value_labels
   StataReader.variable_labels
   StataWriter.write_file
